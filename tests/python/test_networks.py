"""Networks of a board's nets: S-parameters, mixed-mode parameters, TDR, Touchstone.

Against closed forms of uniform lines built from cross-sections of the same geometry,
against the analytic mode conversion of a skew, the quarter-wave resonance of a via stub,
and the reciprocity and passivity of a real channel.
"""

from __future__ import annotations

import cmath
import math

import numpy as np
import pytest

import kairoq
import synthetic
from conftest import TESTDATA

C0 = 299792458.0


def line_s(zc, gamma, length, z0):
    """S11, S21 of a uniform line between two ports of reference z0."""
    ch, sh = cmath.cosh(gamma * length), cmath.sinh(gamma * length)
    a, b, c, d = ch, zc * sh, sh / zc, ch
    den = a + b / z0 + c * z0 + d
    return (a + b / z0 - c * z0 - d) / den, 2 / den


def microstrip_cut(w=1.0, h=1.5, t=0.035, m=0.01, er=4.5, tand=0.02, lateral=12.0):
    """The cross-section of synthetic.microstrip as the board analysis cuts it: the strip
    with conformal mask over the laminate, the plane below, across the board's 24 mm, with
    the microstrip dispersion the analysis gives it."""
    zb, zt = 0.01 + 0.035 + h, 0.01 + 0.035 + h + t
    mask = [[-lateral, zb], [lateral, zb], [lateral, zb + m], [w / 2 + m, zb + m], [w / 2 + m, zt + m], [-w / 2 - m, zt + m],
            [-w / 2 - m, zb + m], [-lateral, zb + m]]
    return {"conductors": [{"name": "s", "rect_mm": [-w / 2, zb, w / 2, zt]},
                           {"name": "g", "rect_mm": [-lateral, 0.01, lateral, 0.045], "reference": True}],
            "dielectrics": [{"rect_mm": [-lateral, 0.045, lateral, zb], "epsilon_r": er, "loss_tangent": tand},
                            {"polygon_mm": mask, "epsilon_r": 3.8, "loss_tangent": 0.025}],
            "target_error": 0.02, "margin": 0.5,
            "microstrip": {"width_mm": w, "height_mm": h, "epsilon_r": er, "loss_tangent": tand}}


# The board's cuts at the hand-built cut's accuracy.
FINE = {"cut": {"target_error": 0.02}}


@pytest.fixture(scope="module")
def microstrip(tmp_path_factory):
    d = tmp_path_factory.mktemp("ms")
    synthetic.microstrip(d)
    return kairoq.open(d)


def channel(board, ports, **kw):
    spec = {"channels": [{"name": "C", "ports": ports, **kw.pop("channel", {})}], "network": kw.pop("network", {"f_max_hz": 10e9, "points": 100})}
    spec.update(kw)
    return kairoq.analyse_channels(board, spec)[0]


def test_uniform_line(microstrip, validation):
    """The track of synthetic.microstrip between its pads against the closed form of a
    uniform line with the per-length parameters of the same cross-section."""
    r = channel(microstrip, ["R1.2", "R2.1"], **FINE)
    spec = microstrip_cut()
    spec["frequencies_hz"] = r["frequency_hz"]
    x = kairoq.cross_section(spec)
    worst_mag = worst_phase = 0.0
    for f, s, fr in zip(r["frequency_hz"], r["s"], x["frequencies"]):
        w = 2 * math.pi * f
        z = complex(fr["r"][0], w * fr["l"][0])
        y = complex(fr["g"][0], w * fr["c"][0])
        s11, s21 = line_s(cmath.sqrt(z / y), cmath.sqrt(z * y), 0.020, 50.0)
        worst_mag = max(worst_mag, abs(20 * math.log10(abs(s[2])) - 20 * math.log10(abs(s21))))
        worst_phase = max(worst_phase, abs(math.degrees(cmath.phase(s[2] / s21))))
        assert abs(s[0] - s11) < 0.01
    validation.check("microstrip 20 mm, w 1, h 1.5 mm, to 10 GHz", "largest difference of |S21|", worst_mag, 0.0, abs_=0.03, unit="dB",
                     source="closed form, uniform line")
    validation.check("microstrip 20 mm, w 1, h 1.5 mm, to 10 GHz", "largest difference of the phase of S21", worst_phase, 0.0, abs_=0.5,
                     unit="degrees", source="closed form, uniform line")


def test_cascade_does_not_depend_on_the_step(microstrip):
    a = channel(microstrip, ["R1.2", "R2.1"], cut={"step_mm": 0.5})
    b = channel(microstrip, ["R1.2", "R2.1"], cut={"step_mm": 3.0})
    for sa, sb in zip(a["s"], b["s"]):
        assert max(abs(p - q) for p, q in zip(sa, sb)) < 1e-9


def test_matched_load(microstrip):
    """One port with a 50 ohm load at the far pad: the input reflection of a line into a
    resistor."""
    r = channel(microstrip, ["R1.2"], channel={"loads": [{"pad": "R2.1", "r_ohm": 50}]}, **FINE)
    spec = microstrip_cut()
    spec["frequencies_hz"] = r["frequency_hz"]
    x = kairoq.cross_section(spec)
    for f, s, fr in zip(r["frequency_hz"], r["s"], x["frequencies"]):
        w = 2 * math.pi * f
        z = complex(fr["r"][0], w * fr["l"][0])
        y = complex(fr["g"][0], w * fr["c"][0])
        zc, th = cmath.sqrt(z / y), cmath.tanh(cmath.sqrt(z * y) * 0.020)
        zin = zc * (50 + zc * th) / (zc + 50 * th)
        assert abs(s[0] - (zin - 50) / (zin + 50)) < 0.01


def test_tdr(microstrip, validation):
    """The TDR shows the line's impedance until the reflection of its end returns."""
    r = channel(microstrip, ["R1.2", "R2.1"], network={"f_max_hz": 20e9, "points": 400})
    si = kairoq.analyse_si(microstrip, {"nets": ["SIG"]})
    s = si.samples[0]
    t = r["tdr"][0]
    delay = 20.0 * s["delay_ps_per_mm"] * 1e-12
    plateau = [z for tt, z in zip(t["time_s"], t["impedance_ohm"]) if 0.4 * delay < tt < 1.6 * delay]
    validation.check("microstrip 20 mm", "TDR impedance on the line", sum(plateau) / len(plateau), s["z0"], rel=0.01, unit="ohm",
                     source="cross-section")
    # The end returns after twice the delay: the trace crosses halfway at that time.
    zmid = 0.5 * (s["z0"] + 50.0)
    cross = next(tt for tt, z in zip(t["time_s"], t["impedance_ohm"]) if tt > delay and z < zmid)
    validation.check("microstrip 20 mm", "TDR time of the far end", cross * 1e12, 2 * delay * 1e12, rel=0.05, unit="ps",
                     source="twice the line's delay")


def test_pair(tmp_path, validation):
    """A symmetric pair: its differential transmission against the closed form of the
    odd mode, and no conversion to the common mode."""
    synthetic.diff_pair(tmp_path)
    b = kairoq.open(tmp_path)
    r = channel(b, ["R1.1", "R1.2", "R2.1", "R2.2"], channel={"differential": [[1, 2], [3, 4]]})
    w, s, h, t, m = 0.2, 0.2, 0.2, 0.035, 0.01
    xp = (s / 2, s / 2 + w)
    zb, zt = 0.01 + 0.035 + h, 0.01 + 0.035 + h + t
    mask = [[-2.5, zb], [2.5, zb], [2.5, zb + m], [xp[1] + m, zb + m], [xp[1] + m, zt + m], [xp[0] - m, zt + m], [xp[0] - m, zb + m],
            [-xp[0] + m, zb + m], [-xp[0] + m, zt + m], [-xp[1] - m, zt + m], [-xp[1] - m, zb + m], [-2.5, zb + m]]
    spec = {"conductors": [{"name": "p", "rect_mm": [xp[0], zb, xp[1], zt]}, {"name": "n", "rect_mm": [-xp[1], zb, -xp[0], zt]},
                           {"name": "g", "rect_mm": [-2.5, 0.01, 2.5, 0.045], "reference": True}],
            "dielectrics": [{"rect_mm": [-2.5, 0.045, 2.5, zb], "epsilon_r": 4.5, "loss_tangent": 0.02},
                            {"polygon_mm": mask, "epsilon_r": 3.8, "loss_tangent": 0.025}],
            "frequencies_hz": r["frequency_hz"], "target_error": 0.03, "margin": 0.5,
            "microstrip": {"width_mm": w, "height_mm": h, "epsilon_r": 4.5, "loss_tangent": 0.02}}
    x = kairoq.cross_section(spec)
    worst = 0.0
    for f, mm, fr in zip(r["frequency_hz"], r["mixed"], x["frequencies"]):
        wv = 2 * math.pi * f
        z = lambda i: complex(fr["r"][i], wv * fr["l"][i])
        y = lambda i: complex(fr["g"][i], wv * fr["c"][i])
        zd, yd = 2 * (z(0) - z(1)), (y(0) - y(1)) / 2  # differential mode per length
        _, sdd21 = line_s(cmath.sqrt(zd / yd), cmath.sqrt(zd * yd), 0.020, 100.0)
        worst = max(worst, abs(mm[1 * 4 + 0] - sdd21))
        assert abs(mm[3 * 4 + 0]) < 1e-3  # SCD21
    validation.check("pair w = s = h = 0.2 mm, 20 mm, to 10 GHz", "largest difference of SDD21", worst, 0.0, abs_=0.01,
                     source="closed form, odd mode")


def test_skew(tmp_path, validation):
    """A skew converts the differential mode to the common mode: |SCD21| = |sin(w dt / 2)|
    for a delay difference dt between matched lines."""
    synthetic.diff_pair(tmp_path, extra_n=1.0)
    b = kairoq.open(tmp_path)
    r = channel(b, ["R1.1", "R1.2", "R2.1", "R2.2"], channel={"differential": [[1, 2], [3, 4]]})
    si = kairoq.analyse_si(b, {"nets": ["DP_N"]})
    alone = [s["delay_ps_per_mm"] for s in si.samples if s["z_differential"] == 0 and s["z0"] > 0]
    dt = 1.0 * sum(alone) / len(alone) * 1e-12
    for f_target in (2.5e9, 5e9):
        k = min(range(len(r["frequency_hz"])), key=lambda i: abs(r["frequency_hz"][i] - f_target))
        f = r["frequency_hz"][k]
        scd21 = 20 * math.log10(abs(r["mixed"][k][3 * 4 + 0]))
        ref = 20 * math.log10(abs(math.sin(2 * math.pi * f * dt / 2)))
        validation.check(f"pair with 1 mm skew, {f / 1e9:g} GHz", "SCD21", scd21, ref, abs_=1.0, unit="dB", source="sin(w dt / 2)")


def test_via_stub(tmp_path, validation):
    """Through vias leave stubs below the inner layer they connect: with the coaxial model,
    a notch where a stub is a quarter wave long (a check of the via topology)."""
    g = synthetic.via_stub_board(tmp_path)
    b = kairoq.open(tmp_path)
    r = channel(b, ["J1.1", "J2.1"], network={"f_max_hz": 20e9, "points": 400, "via_model": "coaxial"})
    s21 = [abs(m[2]) for m in r["s"]]
    k = min(range(len(s21)), key=lambda i: s21[i])
    quarter = C0 / (4 * g["stub_mm"] * 1e-3 * math.sqrt(4.5))
    validation.check("via stub 2.835 mm in er 4.5", "frequency of the S21 notch", r["frequency_hz"][k] / 1e9, quarter / 1e9, rel=0.02,
                     unit="GHz", source="quarter wave")
    assert r["vias"] == 2


def test_real_channel():
    """The PCIe transmit pair of cm5_minima: reciprocal and passive."""
    b = kairoq.open(TESTDATA / "kicad" / "cm5_minima" / "cm5_minima.ipc2581.zip")
    r = channel(b, ["Module302.122", "Module302.124", "J701.47", "J701.49"], channel={"differential": [[1, 2], [3, 4]]},
                network={"f_max_hz": 20e9, "points": 200})
    assert r["coupled_sections"] > 50 and r["vias"] >= 2
    for m in r["s"]:
        for i in range(4):
            for j in range(4):
                assert abs(m[i * 4 + j] - m[j * 4 + i]) < 1e-9
        # Passive: no column carries more power out than in.
        for j in range(4):
            assert sum(abs(m[i * 4 + j]) ** 2 for i in range(4)) <= 1 + 1e-9
    k = min(range(len(r["frequency_hz"])), key=lambda i: abs(r["frequency_hz"][i] - 1e9))
    si = kairoq.analyse_si(b, {"nets": ["/CM5/PCIE_PI.TX_*"]})
    loss = next(n["loss_db"] for n in si.nets if n["name"].endswith("TX_P"))
    # Matched, the insertion loss at 1 GHz is about the line's loss.
    assert -20 * math.log10(abs(r["mixed"][k][1 * 4 + 0])) == pytest.approx(loss, rel=0.25)


def test_touchstone_round_trip(microstrip, tmp_path):
    r = channel(microstrip, ["R1.2", "R2.1"])
    (tmp_path / "c.s2p").write_text(r["touchstone"])
    t = kairoq.read_touchstone(tmp_path / "c.s2p")
    assert t["frequency_hz"] == pytest.approx(r["frequency_hz"], rel=1e-9)
    for a, b in zip(t["s"], r["s"]):
        assert max(abs(p - q) for p, q in zip(a, b)) < 1e-9


# The notch of via_stack_board's stubs in openEMS (docs/validation/networks.md: cells of 0.025 mm,
# the pads' and antipads' circles traced by the mesh).
STACK_NOTCH_GHZ = 11.5


def test_via_model_against_fullwave(tmp_path, validation):
    """An 8-layer board: a track on L1 to a via, on L3 to a second via and back on L1; each
    via leaves a stub of 2 mm through four thin pairs of planes, with non-functional pads
    and four return vias. The default (physical) via model against the full-wave notch."""
    synthetic.via_stack_board(tmp_path)
    b = kairoq.open(tmp_path)
    r = channel(b, ["J1.1", "J2.1"], network={"f_max_hz": 20e9, "points": 200})
    s21 = [abs(m[2]) for m in r["s"]]
    k = min(range(len(s21)), key=s21.__getitem__)
    validation.check("8-layer via, 2 mm stub through four pairs, return vias", "frequency of the S21 notch", r["frequency_hz"][k] / 1e9,
                     STACK_NOTCH_GHZ, rel=0.1, unit="GHz", source="openEMS (stored)")


GERBER2EMS = TESTDATA / "gerber2ems"


def vna(case):
    rows = [line.split(",") for line in (GERBER2EMS / case / "vna.csv").read_text().splitlines()[1:]]
    return np.array([float(r[0]) for r in rows]) * 1e6, np.array([complex(float(r[1]), float(r[2])) for r in rows])


def test_registration_of_drills_and_placement():
    """The stub board's drill and placement files are relative to its lower-left corner: the
    via moves onto the ground pad of R4, the simulation ports onto the trace's ends."""
    b = kairoq.open(GERBER2EMS / "stub_short" / "fab")
    assert any("drill coordinates are relative" in w for w in b.warnings)
    via = b.drills[0]
    assert via["net"] == b.nets.index("GND") and via["position_mm"] == pytest.approx([116.576, -99.11], abs=1e-6)
    sp = {c["refdes"]: c["position_mm"] for c in b.components}
    assert sp["SP1"] == pytest.approx([100.8, -96.0], abs=1e-6) and sp["SP2"] == pytest.approx([137.0, -96.0], abs=1e-6)


@pytest.mark.parametrize("case,rms", [("stub_short", 0.12), ("stub_long", 0.15)])
def test_stub_against_measurement(validation, case, rms):
    """A microstrip with an open stub against its VNA measurement, 0.2 to 6 GHz, port SP2
    loaded with 50 ohm. The measurement carries a ripple of about 0.25 GHz from the set-up,
    averaged out over one period; the 0 ohm part at the stub's end is taken as not fitted, as
    the measurement shows (an open stub's resonance)."""
    b = kairoq.open(GERBER2EMS / case / "fab")
    sp = sorted((c["refdes"], c["position_mm"]) for c in b.components if c["refdes"].startswith("SP"))
    ports = ["%s@%.4f,%.4f" % (r, p[0], p[1]) for r, p in sp]
    r = kairoq.analyse_channels(b, {"channels": [{"name": "STUB", "ports": ports}], "network": {"f_max_hz": 6e9, "points": 1200}})[0]
    f, g = np.array(r["frequency_hz"]), np.array(r["s"])[:, 0]
    fm, gm = vna(case)
    n = int(round(0.25e9 / (fm[1] - fm[0])))
    gs = np.convolve(gm, np.ones(n) / n, mode="same")[n:-n]
    fv = fm[n:-n]
    gi = np.interp(fv, f, g.real) + 1j * np.interp(fv, f, g.imag)
    d = abs(gi - gs)
    validation.check(f"gerber2ems {case}, 0.2 to 6 GHz", "rms |S11 model - S11 measured|", float(np.sqrt(np.mean(d ** 2))), 0.0, abs_=rms,
                     source="VNA measurement (antmicro/gerber2ems)")
    if case == "stub_short":
        z, zm = 50 * (1 + gi) / (1 - gi), 50 * (1 + gs) / (1 - gs)
        validation.check("gerber2ems stub_short", "frequency of the largest |Z|", fv[np.argmax(abs(z))] / 1e9, fv[np.argmax(abs(zm))] / 1e9, rel=0.03,
                         unit="GHz", source="VNA measurement (antmicro/gerber2ems)")
        validation.check("gerber2ems stub_short", "frequency of the smallest |Z|", fv[np.argmin(abs(z))] / 1e9, fv[np.argmin(abs(zm))] / 1e9, rel=0.08,
                         unit="GHz", source="VNA measurement (antmicro/gerber2ems)")


def tee_board(out, split=False, gap_nm=0):
    """A 0.2 mm microstrip over a plane, 30 mm long, with a 5 mm stub from its middle: the
    stub ends on the line's inside, or the line is drawn in two pieces meeting there (with
    a gap of gap_nm between their ends)."""
    from synthetic import Gerber, _write
    top, bot, prof = Gerber("Copper,L1,Top"), Gerber("Copper,L2,Bot"), Gerber("Profile,NP")
    if split:
        top.track([(5.0, 10.0), (20.0 - gap_nm * 1e-6 / 2, 10.0)], 0.2, "SIG")
        top.track([(20.0 + gap_nm * 1e-6 / 2, 10.0), (35.0, 10.0)], 0.2, "SIG")
    else:
        top.track([(5.0, 10.0), (35.0, 10.0)], 0.2, "SIG")
    top.track([(20.0, 10.0), (20.0, 15.0)], 0.2, "SIG")
    bot.region([(0, 0), (40, 0), (40, 20), (0, 20)], "GND")
    prof.outline([(0, 0), (40, 0), (40, 20), (0, 20)])
    _write(out, "tee", {"tee-L1.gtl": top, "tee-L2.gbl": bot, "tee-Edge.gm1": prof}, {}, [0.035, 0.035], [0.2])


def tee_s(path):
    r = kairoq.analyse_channels(kairoq.open(path), {"channels": [{"name": "T", "ports": ["A@5,10", "B@35,10"]}],
                                                        "network": {"f_max_hz": 10e9, "points": 100}})[0]
    return np.array(r["s"])


def test_stub_ending_on_a_line(tmp_path):
    """A track ending on the inside of another meets it at a junction: the result equals
    that of the line drawn in two pieces; ends 2 nm apart across a micrometre are one node."""
    tee_board(tmp_path / "t")
    tee_board(tmp_path / "s", split=True)
    tee_board(tmp_path / "g", split=True, gap_nm=2)
    t, s, g = tee_s(tmp_path / "t"), tee_s(tmp_path / "s"), tee_s(tmp_path / "g")
    assert np.max(abs(t - s)) < 1e-9
    assert np.max(abs(g - s)) < 1e-6
    # The open stub is a quarter wave near 8 GHz: the line is not a plain through line.
    assert np.min(abs(t[:, 2])) < 0.5
