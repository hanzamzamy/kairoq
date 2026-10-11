"""Time-domain response: reflections on a line against the lattice diagram, the
convolution engine against an independent frequency-domain solution, IBIS buffers
against equivalent linear sources and their own waveforms, the IBIS reader and the bit
patterns.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import kairoq
import synthetic

VCC = 1.8


@pytest.fixture(scope="module")
def microstrip(tmp_path_factory):
    d = tmp_path_factory.mktemp("ms")
    synthetic.microstrip(d)
    b = kairoq.open(d)
    s = kairoq.analyse_si(b, {"nets": ["SIG"]}).samples[0]
    return b, s["z0"], 20.0 * s["delay_ps_per_mm"] * 1e-12


NETWORK = {"f_max_hz": 40e9, "points": 800}


def run(board, drivers, receivers=(), duration_ns=2.0, ports=("R1.2", "R2.1")):
    si = {"channels": [{"name": "C", "ports": list(ports)}], "network": NETWORK,
          "simulations": [{"name": "S", "channel": "C", "drivers": list(drivers), "receivers": list(receivers), "duration_ns": duration_ns}]}
    return kairoq.simulate(board, si)[0]


def value_at(r, port, t):
    k = min(range(len(r["time_s"])), key=lambda i: abs(r["time_s"][i] - t))
    return r["voltage"][port][k]


def test_lattice(microstrip, validation):
    """A 1 V step behind 50 ohm into the line, its far end open: the near end steps to
    Z / (Z + 50), the far end to twice that when the step arrives, both to 1 V at last."""
    b, z, td = microstrip
    r = run(b, [{"port": 1, "rise_ps": 50, "r_ohm": 50}])
    first = z / (z + 50)
    delay = 100e-12 + 1.25 * 25e-12  # the driver's delay and the middle of its edge
    validation.check("microstrip 20 mm, open end", "near-end step", value_at(r, 0, delay + td), first, rel=0.01, unit="V",
                     source="lattice diagram")
    # The line's loss (conductor and dielectric) lowers the step arriving at the far end.
    validation.check("microstrip 20 mm, open end", "far-end step", value_at(r, 1, delay + 1.5 * td), 2 * first, rel=0.02, unit="V",
                     source="lattice diagram (lossless)")
    validation.check("microstrip 20 mm, open end", "final value", r["voltage"][1][-1], 1.0, rel=0.002, unit="V", source="DC")


def test_against_frequency_domain(microstrip, validation):
    """The convolution with the drivers' Newton step against the exact linear solution in
    the frequency domain (numpy, from the channel's impedance matrix): a square wave of
    the sweep's period behind 50 ohm, the far end loaded with 100 ohm and 2 pF."""
    b, z, td = microstrip
    r = run(b, [{"port": 1, "rise_ps": 50, "r_ohm": 50}], [{"port": 2, "r_ohm": 100, "c_pf": 2.0}], duration_ns=4.0)
    ch = kairoq.analyse_channels(b, {"channels": [{"name": "C", "ports": ["R1.2", "R2.1"]}], "network": NETWORK})[0]
    f = np.array(ch["frequency_hz"])
    S = np.array(ch["s"]).reshape(-1, 2, 2)
    fmax, df = f[-1], f[1] - f[0]
    period = 1 / df
    w = np.where(f <= 0.75 * fmax, 1.0, 0.5 * (1 + np.cos(np.pi * (f - 0.75 * fmax) / (0.25 * fmax))))
    # Source: high from the delay for half the period, linear edges over 62.5 ps.
    t0, ramp = 100e-12, 62.5e-12
    k = np.arange(1, len(f) + 1)
    def edge(ts):  # Fourier coefficients of a unit ramp step at ts, periodic over the period
        x = 2j * np.pi * k * df
        return (np.exp(-x * ts) - np.exp(-x * (ts + ramp))) / (x * x * ramp) / period
    vs = edge(t0) - edge(t0 + period / 2)
    I = np.eye(2)
    v2 = np.zeros(len(f), complex)
    for i in range(len(f)):
        Z = 50 * (I + S[i]) @ np.linalg.inv(I - S[i])
        ZT = np.diag([50.0, 1 / (1 / 100 + 2j * np.pi * f[i] * 2e-12)])
        cur = np.linalg.solve(Z + ZT, np.array([vs[i], 0]))
        v2[i] = (Z @ cur)[1]
    t = np.array(r["time_s"])
    v_dc = 0.5 * 100 / (100 + 50)  # the square wave's mean through the resistive divider
    exact = v_dc + 2 * np.real((w * v2)[None, :] * np.exp(2j * np.pi * f[None, :] * t[:, None])).sum(1)
    worst = float(np.max(np.abs(exact - np.array(r["voltage"][1]))))
    validation.check("microstrip 20 mm, 100 ohm and 2 pF at the far end", "largest difference of the far-end waveform", worst, 0.0, abs_=0.01,
                     unit="V", source="frequency-domain solution")


def ibis_text(model="LIN", r_out=50.0, c_comp=0.0, ramp_ps=250.0, waveforms=None):
    """A buffer of output resistance r_out to VCC and to ground, clamps beyond the rails,
    switching linearly in ramp_ps; with waveforms, [Rising Waveform] tables."""
    v = np.linspace(-VCC, 2 * VCC, 55)
    clamp = np.where(v < -0.5, (v + 0.5) / 2.0, 0.0)  # into the pin, below -0.5 V
    rows = lambda cur: "\n".join(f"{a:.6f} {b:.9g}A NA NA" for a, b in zip(v, cur))
    text = f"""[IBIS Ver] 5.0
|comment line
[Component] SYN
[Manufacturer] none
[Package]
| typ min max
R_pkg 0 NA NA
L_pkg 0 NA NA
C_pkg 0 NA NA
[Pin] signal_name model_name R_pin L_pin C_pin
1 OUT {model}
[Model] {model}
Model_type Output
C_comp {c_comp * 1e12:.4f}pF NA NA
[Voltage Range] {VCC}V NA NA
[Pulldown]
{rows(v / r_out)}
[Pullup]
{rows(-v / r_out)}
[GND Clamp]
{rows(clamp)}
[POWER Clamp]
{rows(-clamp)}
[Ramp]
dV/dt_r {VCC:.3f}/{ramp_ps:.1f}p NA NA
dV/dt_f {VCC:.3f}/{ramp_ps:.1f}p NA NA
R_load = 50
"""
    for vf, table in waveforms or []:
        text += f"[Rising Waveform]\nR_fixture = 50\nV_fixture = {vf}\n" + "\n".join(f"{a * 1e12:.4f}p {b:.6f} NA NA" for a, b in table) + "\n"
    return text + "[End]\n"


def test_ibis_equals_linear_source(microstrip, tmp_path, validation):
    """A buffer with linear stages of 50 ohm switching in 250 ps, without package or die
    capacitance, is a source of 0 to 1.8 V behind 50 ohm with a 250 ps edge."""
    b, z, td = microstrip
    (tmp_path / "lin.ibs").write_text(ibis_text())
    receivers = [{"port": 2, "r_ohm": 50}]
    a = run(b, [{"port": 1, "ibis": {"file": str(tmp_path / "lin.ibs"), "model": "LIN", "package": False}}], receivers)
    c = run(b, [{"port": 1, "v_high": VCC, "r_ohm": 50, "rise_ps": 200}], receivers)
    worst = max(abs(x - y) for p in range(2) for x, y in zip(a["voltage"][p], c["voltage"][p]))
    validation.check("IBIS buffer of linear 50 ohm stages, 250 ps", "largest difference from the linear source", worst, 0.0, abs_=1e-3, unit="V",
                     source="equivalent source")


def test_ibis_switching_from_waveforms(tmp_path, validation):
    """Rising waveforms into 50 ohm to ground and to VCC, computed for stages of 40 ohm with
    1 pF at the die switching linearly in 200 ps: the switching coefficients recovered from
    them are that ramp."""
    r_out, cc, tr = 40.0, 1e-12, 200e-12
    ku = lambda t: min(1.0, max(0.0, (t - 50e-12) / tr))
    def waveform(vf):
        # The pin starts where the low state holds it against the fixture.
        h, t, v, out = 1e-14, 0.0, vf * r_out / (r_out + 50), []
        for n in range(80001):
            if n % 1000 == 0:
                out.append((t, v))
            k = ku(t)
            i_out = k * (VCC - v) / r_out - (1 - k) * v / r_out  # out of the pin
            v += h * (i_out - (v - vf) / 50) / cc
            t += h
        return out
    (tmp_path / "wf.ibs").write_text(ibis_text(model="WF", r_out=r_out, c_comp=cc, waveforms=[(0.0, waveform(0.0)), (VCC, waveform(VCC))]))
    sw = kairoq.ibis_switching(tmp_path / "wf.ibs", "WF", rising=True)
    worst = max(abs(k - ku(t)) for t, k in zip(sw["time_s"], sw["ku"]))
    worst_d = max(abs(k - (1 - ku(t))) for t, k in zip(sw["time_s"], sw["kd"]))
    validation.check("IBIS stages of 40 ohm, 1 pF, 200 ps ramp", "largest error of the recovered ku", worst, 0.0, abs_=0.02,
                     source="the ramp the waveforms were computed with")
    validation.check("IBIS stages of 40 ohm, 1 pF, 200 ps ramp", "largest error of the recovered kd", worst_d, 0.0, abs_=0.02,
                     source="the ramp the waveforms were computed with")


def test_ibis_reader(tmp_path):
    (tmp_path / "a.ibs").write_text(ibis_text(c_comp=1.5e-12).replace("[IBIS Ver] 5.0", "[IBIS Ver] 5.0\n[Comment Char] #_char")
                                    .replace("|comment line", "#comment line"))
    sw = kairoq.ibis_switching(tmp_path / "a.ibs", "LIN", rising=False)
    assert sw["type"] == "Output" and sw["c_comp"] == pytest.approx(1.5e-12) and sw["vpu"] == pytest.approx(VCC)
    # Ramp only: a linear fall over VCC / slope.
    assert sw["time_s"][-1] == pytest.approx(250e-12) and sw["ku"][0] == 1.0 and sw["kd"][-1] == 1.0
    assert any("[manufacturer]" in w for w in sw["warnings"])


def test_prbs():
    for name, n in (("prbs7", 7), ("prbs9", 9)):
        period = 2**n - 1
        bits = kairoq._core.prbs(name, 2 * period)
        assert bits[:period] == bits[period:]  # the period of a maximal sequence
        assert sum(bits[:period]) == 2 ** (n - 1)  # its ones
        assert all(bits[i : i + period] != bits[:period] for i in range(1, period))


def test_ibis_receiver(microstrip, tmp_path, validation):
    """An IBIS input with 2 pF at the die and its clamps off (within the rails) is a 2 pF
    load."""
    b, z, td = microstrip
    (tmp_path / "rx.ibs").write_text(ibis_text(model="RX", c_comp=2e-12))
    driver = [{"port": 1, "rise_ps": 50, "r_ohm": 50}]
    a = run(b, driver, [{"port": 2, "ibis": {"file": str(tmp_path / "rx.ibs"), "model": "RX", "package": False}}])
    c = run(b, driver, [{"port": 2, "c_pf": 2.0}])
    worst = max(abs(x - y) for x, y in zip(a["voltage"][1], c["voltage"][1]))
    validation.check("IBIS input of 2 pF, clamps off", "largest difference from a 2 pF load", worst, 0.0, abs_=1e-4, unit="V",
                     source="equivalent load")
