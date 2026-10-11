"""Eye diagrams: the single-bit response against the time-domain engine, the statistical
eye against an enumeration of the intersymbol interference with noise, the worst case
against peak-distortion analysis and the pattern's eye, the equalisers, and jitter.
"""

from __future__ import annotations

import itertools
import math
import subprocess

import numpy as np
import pytest

import kairoq
import synthetic
from conftest import cli_path

NETWORK = {"f_max_hz": 40e9, "points": 800}  # period 20 ns
CHANNEL = [{"name": "C", "ports": ["R1.2", "R2.1"]}]


@pytest.fixture(scope="module")
def microstrip(tmp_path_factory):
    d = tmp_path_factory.mktemp("ms")
    synthetic.microstrip(d)
    return kairoq.open(d)


@pytest.fixture(scope="module")
def pair(tmp_path_factory):
    d = tmp_path_factory.mktemp("dp")
    synthetic.diff_pair(d)
    return kairoq.open(d)


def run_eye(board, channels=CHANNEL, **spec):
    e = {"name": "E", "channel": "C", "tx": [1], "rx": [2]}
    e.update(spec)
    # The eyes' arithmetic against closed forms: lines without microstrip dispersion, whose
    # response on this thick substrate would outlast the patterns these checks rely on.
    return kairoq.eye(board, {"channels": channels, "network": NETWORK, "eyes": [e], "dispersion": False})[0]


def superpose(r, bits, signs, t):
    """The waveform of a bit sequence from the single-bit response: sum_k s_k p(t - k UI)."""
    pt, pv = np.array(r["pulse_t"]), np.array(r["pulse_v"])
    return sum(s * np.interp(t - k * r["ui_s"], pt, pv, left=0.0, right=0.0) for k, s in zip(range(len(bits)), signs))


def test_pulse_against_waveforms(microstrip, validation):
    """Single-ended: 16 bits at 2.5 Gbit/s, 0 to 1 V behind 50 ohm into the line, 50 ohm at
    the far end. The superposition of the single-bit response against the convolution
    engine of the time-domain response (docs/solvers/timedomain.md), both at 3.125 ps."""
    bits = "0110100111010010"
    r = run_eye(microstrip, bit_rate=2.5e9, rise_ps=30, samples_per_ui=128)
    si = {"channels": CHANNEL, "network": NETWORK, "dispersion": False,
          "simulations": [{"name": "S", "channel": "C", "duration_ns": 7.0,
                           "drivers": [{"port": 1, "pattern": bits, "bit_rate": 2.5e9, "delay_ps": 0, "rise_ps": 30, "r_ohm": 50}],
                           "receivers": [{"port": 2, "r_ohm": 50}]}]}
    w = kairoq.simulate(microstrip, si)[0]
    # Until the pattern's end: then the driver holds its last level and the bits end.
    t = np.array(w["time_s"])
    keep = t <= len(bits) * r["ui_s"]
    expected = superpose(r, bits, [int(b) for b in bits], t[keep])
    worst = float(np.max(np.abs(expected - np.array(w["voltage"][1])[keep])))
    validation.check("microstrip 20 mm, 2.5 Gbit/s, 16 bits", "largest difference of the far-end waveform", worst, 0.0, abs_=2e-3,
                     unit="V", source="time-domain response")


def test_differential_pulse(pair, validation):
    """Differential: the pair driven by complementary legs (0 to 0.5 V each, 1 V of
    differential swing), terminated with 50 ohm per leg; the differential far-end
    waveform against the time-domain response."""
    ports = ["R1.1", "R1.2", "R2.1", "R2.2"]
    bits = "1101001110100101"
    comp = "".join("1" if b == "0" else "0" for b in bits)
    r = run_eye(pair, channels=[{"name": "C", "ports": ports}], tx=[1, 2], rx=[3, 4], bit_rate=2.5e9, rise_ps=30, samples_per_ui=128)
    drv = {"bit_rate": 2.5e9, "delay_ps": 0, "rise_ps": 30, "r_ohm": 50, "v_high": 0.5}
    si = {"channels": [{"name": "C", "ports": ports}], "network": NETWORK, "dispersion": False,
          "simulations": [{"name": "S", "channel": "C", "duration_ns": 7.0,
                           "drivers": [dict(port=1, pattern=bits, **drv), dict(port=2, pattern=comp, **drv)],
                           "receivers": [{"port": 3, "r_ohm": 50}, {"port": 4, "r_ohm": 50}]}]}
    w = kairoq.simulate(pair, si)[0]
    t = np.array(w["time_s"])
    keep = t <= len(bits) * r["ui_s"]
    # The legs' difference is +-0.5 V about 0: half the single-bit response per symbol.
    expected = superpose(r, bits, [0.5 if b == "1" else -0.5 for b in bits], t[keep])
    worst = float(np.max(np.abs(expected - (np.array(w["voltage"][2]) - np.array(w["voltage"][3]))[keep])))
    validation.check("differential pair 20 mm, 2.5 Gbit/s, 16 bits", "largest difference of the differential far-end waveform", worst, 0.0,
                     abs_=2e-3, unit="V", source="time-domain response")
    assert r["threshold_v"] == 0.0


def enumerated_height(cursors, pre, sigma, ber):
    """The eye height at a BER by enumerating the signs of the cursors other than the main
    one, the noise Gaussian: the level x where P(one below x) = BER, doubled."""
    main = 0.5 * cursors[pre]
    others = [0.5 * c for i, c in enumerate(cursors) if i != pre and abs(c) > 0]
    levels = np.array([main + sum(s * c for s, c in zip(signs, others)) for signs in itertools.product((-1, 1), repeat=len(others))])
    p = lambda x: float(np.mean(0.5 * np.vectorize(math.erfc)((levels - x) / (sigma * math.sqrt(2)))))
    a, b = levels.min() - 10 * sigma, levels.max()
    for _ in range(80):
        m = 0.5 * (a + b)
        a, b = (m, b) if p(m) <= ber else (a, m)
    return 2 * a


@pytest.mark.parametrize("sigma_mv", [5.0, 20.0])
def test_statistical_eye_with_noise(microstrip, validation, sigma_mv):
    """The statistical eye's height at the best phase against the enumeration of all signs
    of its cursors with Gaussian noise (10 Gbit/s: reflections at 2 UI)."""
    r = run_eye(microstrip, bit_rate=10e9, noise_mv=sigma_mv)
    ref = enumerated_height(r["cursors"], r["pre_cursors"], sigma_mv * 1e-3, 1e-12)
    validation.check(f"microstrip 20 mm, 10 Gbit/s, noise {sigma_mv:g} mV", "eye height at 1e-12", r["eye_height_v"], ref, abs_=1e-3,
                     unit="V", source="enumeration of the cursors' signs")


def test_noise_only_closed_form(microstrip, validation):
    """At 1 Gbit/s the line's reflections settle within the bit; with noise the eye height
    is 2 (A - Q sigma) for the main level A, Q = 7.034 at 1e-12."""
    sigma = 0.01
    r = run_eye(microstrip, bit_rate=1e9, noise_mv=10)
    a = 0.5 * r["cursors"][r["pre_cursors"]]
    isi = sum(0.5 * abs(c) for i, c in enumerate(r["cursors"]) if i != r["pre_cursors"])
    assert isi < 2e-3 * a
    validation.check("microstrip 20 mm, 1 Gbit/s, noise 10 mV", "eye height at 1e-12", r["eye_height_v"], 2 * (a - 7.034 * sigma), abs_=2 * isi + 1e-4,
                     unit="V", source="2 (A - Q sigma), the residual ISI as tolerance")


def test_worst_case_and_pattern(microstrip, validation):
    """Without noise or jitter, at 1e-12 the eye height is the peak-distortion opening
    main - sum |c_k| (fewer than 40 cursors), and the PRBS7 pattern reaches it: the
    response spans fewer than 7 UI."""
    r = run_eye(microstrip, bit_rate=10e9, pattern="prbs7", bits=2048)
    c, pre = r["cursors"], r["pre_cursors"]
    pda = c[pre] - sum(abs(x) for i, x in enumerate(c) if i != pre)
    validation.check("microstrip 20 mm, 10 Gbit/s, no noise", "eye height at 1e-12", r["eye_height_v"], pda, abs_=1e-3, unit="V",
                     source="peak-distortion analysis of the cursors")
    validation.check("microstrip 20 mm, 10 Gbit/s, no noise", "worst-case opening", r["peak_distortion_v"], pda, abs_=1e-9, unit="V",
                     source="peak-distortion analysis of the cursors")
    # The pattern's eye: the lowest one and the highest zero at the eye's centre (t = 0).
    centre = len(r["eye_t"]) // 2
    v, dens = np.array(r["eye_v"]), np.array(r["density"][centre])
    bin_v = v[1] - v[0]
    ones, zeros = v[(dens > 0) & (v > r["threshold_v"])], v[(dens > 0) & (v < r["threshold_v"])]
    validation.check("microstrip 20 mm, 10 Gbit/s, PRBS7", "inner opening of the pattern's eye", float(ones.min() - zeros.max()), pda,
                     abs_=1.5 * bin_v, unit="V", source="peak-distortion analysis (the density's bins)")


def test_dfe(microstrip, validation):
    """The DFE's taps are the post-cursors at the response's peak; they cancel them there."""
    r0 = run_eye(microstrip, bit_rate=10e9)
    r = run_eye(microstrip, bit_rate=10e9, dfe={"taps": 3})
    spu = 32
    pv, peak = r["pulse_v"], r["main_index"]
    for k, tap in enumerate(r["dfe"]):
        validation.check("microstrip 20 mm, 10 Gbit/s, DFE", f"tap {k + 1}", tap, pv[peak + (k + 1) * spu], abs_=1e-12, unit="V",
                         source="the single-bit response one UI apart")
    assert r["eye_height_v"] > r0["eye_height_v"]
    limited = run_eye(microstrip, bit_rate=10e9, dfe={"taps": 3, "limit": 0.01})
    assert max(abs(x) for x in limited["dfe"]) <= 0.01 * pv[peak] + 1e-12


def test_ffe_and_ctle(microstrip, validation):
    """A CTLE whose zero and pole coincide is its DC gain; the FFE shifts and weights the
    response by whole UI; the area of the single-bit response is its DC gain times the UI."""
    base = run_eye(microstrip, bit_rate=5e9)
    flat = run_eye(microstrip, bit_rate=5e9, ctle={"dc_gain_db": -6, "f_zero_hz": 3e9, "f_pole1_hz": 3e9})
    g = 10 ** (-6 / 20)
    worst = max(abs(a - g * b) for a, b in zip(flat["pulse_v"], base["pulse_v"]))
    validation.check("CTLE with coincident zero and pole, -6 dB", "largest difference from the scaled response", worst, 0.0, abs_=1e-12, unit="V",
                     source="the DC gain")
    peaking = run_eye(microstrip, bit_rate=5e9, ctle={"dc_gain_db": -6, "f_zero_hz": 1e9, "f_pole1_hz": 5e9, "f_pole2_hz": 20e9})
    dt = base["pulse_t"][1] - base["pulse_t"][0]
    ratio = sum(peaking["pulse_v"]) / sum(base["pulse_v"])
    validation.check("CTLE: -6 dB, zero 1 GHz, poles 5 and 20 GHz", "area ratio of the single-bit response", ratio, g, rel=1e-6,
                     source="the DC gain")
    assert sum(base["pulse_v"]) * dt == pytest.approx(base["threshold_v"] * 2 * base["ui_s"], rel=1e-6)
    ffe = run_eye(microstrip, bit_rate=5e9, ffe=[-0.1, 0.8, -0.1], ffe_main=2)
    pv = np.array(base["pulse_v"])
    expected = 0.8 * pv - 0.1 * np.roll(pv, -32) - 0.1 * np.roll(pv, 32)
    expected[:32] += 0.1 * pv[-32:]  # the shifts do not wrap
    expected[-32:] += 0.1 * pv[:32]
    assert np.max(np.abs(np.array(ffe["pulse_v"]) - expected)) < 1e-12


def test_jitter(microstrip, validation):
    """Dual-Dirac DJ of 4 samples narrows the eye by exactly 4 samples; Gaussian RJ by
    twice Q sigma (to a sample)."""
    spu, ui = 32, 1e-10
    r0 = run_eye(microstrip, bit_rate=10e9)
    dj = run_eye(microstrip, bit_rate=10e9, dj_ps=4 * ui / spu * 1e12)
    validation.check("microstrip 20 mm, 10 Gbit/s, DJ 12.5 ps", "eye width", dj["eye_width_ui"], r0["eye_width_ui"] - 4 / spu, abs_=1e-9, unit="UI",
                     source="the width without jitter less the DJ")
    sigma = 1.5e-12
    rj = run_eye(microstrip, bit_rate=10e9, rj_ps=sigma * 1e12)
    validation.check("microstrip 20 mm, 10 Gbit/s, RJ 1.5 ps", "eye width", rj["eye_width_ui"], r0["eye_width_ui"] - 2 * 7.034 * sigma / ui,
                     abs_=1.5 / spu, unit="UI", source="the width without jitter less 2 Q sigma")


def test_errors(microstrip):
    with pytest.raises(Exception, match="unknown key"):
        run_eye(microstrip, bitrate=1e9)
    with pytest.raises(Exception, match="port beyond"):
        run_eye(microstrip, rx=[3])
    with pytest.raises(Exception, match="PRBS"):
        run_eye(microstrip, pattern="edge")


@pytest.mark.skipif(cli_path() is None, reason="kairoq command not built")
def test_cli(tmp_path):
    synthetic.microstrip(tmp_path / "board")
    (tmp_path / "project.yaml").write_text(
        "board: [board]\nsi:\n  channels: [{name: C, ports: [R1.2, R2.1]}]\n  network: {f_max_hz: 40e9, points: 800}\n"
        "  eyes: [{name: E10, channel: C, tx: [1], rx: [2], bit_rate: 10e9, noise_mv: 5}]\n")
    r = subprocess.run([cli_path(), "eye", str(tmp_path / "project.yaml"), "-o", str(tmp_path / "out")], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert "eye E10" in r.stdout and "height" in r.stdout
    assert (tmp_path / "out" / "E10.svg").read_text().startswith("<svg")
    assert (tmp_path / "out" / "E10_pulse.csv").exists() and (tmp_path / "out" / "E10_contour.csv").exists()
