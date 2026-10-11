"""Transmission-line cross-sections and the signal-integrity analysis of boards.

Cross-sections against exact and published closed forms (coaxial line, stripline,
microstrip, edge-coupled stripline), the loss models against their definitions and an
independent finite difference, and the board analysis against its own cross-sections, a
plane slot and the impedance targets of a real design.
"""

from __future__ import annotations

import math

import pytest

import kairoq
import synthetic
from conftest import TESTDATA

EPS0 = 8.8541878128e-12
MU0 = 1.25663706212e-6
C0 = 299792458.0
ETA0 = 376.730313668


def ellipk(k: float) -> float:
    """Complete elliptic integral of the first kind, modulus k (arithmetic-geometric mean)."""
    a, b = 1.0, math.sqrt(1.0 - k * k)
    while abs(a - b) > 1e-15 * a:
        a, b = (a + b) / 2, math.sqrt(a * b)
    return math.pi / (2 * a)


def stripline_exact(w, b, er):
    """Centred stripline of zero thickness (conformal mapping; Pozar, Microwave Engineering)."""
    k = 1 / math.cosh(math.pi * w / (2 * b))
    kp = math.tanh(math.pi * w / (2 * b))
    return 30 * math.pi / math.sqrt(er) * ellipk(k) / ellipk(kp)


def coupled_stripline_exact(w, s, b, er):
    """Edge-coupled stripline of zero thickness, even and odd mode (Cohn 1955)."""
    ke = math.tanh(math.pi * w / (2 * b)) * math.tanh(math.pi * (w + s) / (2 * b))
    ko = math.tanh(math.pi * w / (2 * b)) / math.tanh(math.pi * (w + s) / (2 * b))
    z = lambda k: 30 * math.pi / math.sqrt(er) * ellipk(math.sqrt(1 - k * k)) / ellipk(k)
    return z(ke), z(ko)


def microstrip_hj(w, h, er):
    """Microstrip of zero thickness (Hammerstad and Jensen 1980): impedance and effective
    permittivity, stated accurate to 0.01 % and 0.2 %."""
    u = w / h
    fu = 6 + (2 * math.pi - 6) * math.exp(-((30.666 / u) ** 0.7528))
    z_air = ETA0 / (2 * math.pi) * math.log(fu / u + math.sqrt(1 + 4 / u**2))
    a = 1 + math.log((u**4 + (u / 52) ** 2) / (u**4 + 0.432)) / 49 + math.log(1 + (u / 18.1) ** 3) / 18.7
    b = 0.564 * ((er - 0.9) / (er + 3)) ** 0.053
    ee = (er + 1) / 2 + (er - 1) / 2 * (1 + 10 / u) ** (-a * b)
    return z_air / math.sqrt(ee), ee


def z_lossless(r, i=0, n=1):
    return math.sqrt(r["l_ext"][i * n + i] / r["frequencies"][0]["c"][i * n + i])


# Thin strips: a thickness of a thousandth of the width. Wheeler's thickness correction
# makes such a strip about 0.3 % wider electrically, which lowers its impedance by about
# 0.15 to 0.2 % against the zero-thickness forms; the tolerances allow for it.
T = 0.0002


def test_coaxial_line(validation):
    """Exact: C = 2 pi eps / ln(b/a), L = mu0 / (2 pi) ln(b/a), R = Rs / (2 pi) (1/a + 1/b)."""
    a, b, er, f = 0.5, 1.75, 2.25, 1e9
    spec = {"conductors": [{"name": "in", "circle_mm": [0, 0, a]}, {"name": "out", "ring_mm": [0, 0, b, b + 0.25], "reference": True}],
            "dielectrics": [{"circle_mm": [0, 0, b + 0.01], "epsilon_r": er}], "frequencies_hz": [f], "target_error": 0.02,
            "arc_tolerance_um": 0.05}
    r = kairoq.cross_section(spec)
    fr = r["frequencies"][0]
    ln = math.log(b / a)
    rs = math.sqrt(math.pi * f * MU0 / 5.8e7)
    validation.check("coaxial line a 0.5 mm, b 1.75 mm, er 2.25", "capacitance", fr["c"][0] * 1e12, 2 * math.pi * EPS0 * er / ln * 1e12,
                     rel=1e-3, unit="pF/m", source="exact")
    validation.check("coaxial line", "external inductance", r["l_ext"][0] * 1e9, MU0 / (2 * math.pi) * ln * 1e9, rel=1e-3, unit="nH/m",
                     source="exact")
    validation.check("coaxial line", "resistance at 1 GHz", fr["r"][0], rs / (2 * math.pi) * (1 / (a * 1e-3) + 1 / (b * 1e-3)),
                     rel=2e-3, unit="ohm/m", source="exact, skin depth 2.1 um")
    # The DC limit of the surface impedance model: the resistance of the two cross-sections.
    spec["frequencies_hz"] = [0.0]
    r0 = kairoq.cross_section(spec)["frequencies"][0]["r"][0]
    area_in = math.pi * (a * 1e-3) ** 2
    area_out = math.pi * (((b + 0.25) * 1e-3) ** 2 - (b * 1e-3) ** 2)
    assert r0 == pytest.approx(1 / (5.8e7 * area_in) + 1 / (5.8e7 * area_out), rel=2e-3)


def test_body_of_revolution(validation):
    """The axisymmetric option: concentric spheres, C = 4 pi eps0 a b / (b - a), and the
    capacitance per length of a coaxial line from the difference of two lengths."""
    def half(r0, r1, n=256):
        pts = [[r1 * math.sin(math.pi * k / n), r1 * math.cos(math.pi * k / n)] for k in range(n + 1)]
        return pts + [[r0 * math.sin(math.pi * k / n), r0 * math.cos(math.pi * k / n)] for k in range(n, -1, -1)] if r0 > 0 else pts
    a, b = 1.0, 2.0
    r = kairoq.cross_section({"axisymmetric": True, "margin": 0.2, "target_error": 0.01,
                                "conductors": [{"name": "in", "polygon_mm": half(0, a)},
                                               {"name": "out", "reference": True, "polygon_mm": half(b, b + 0.2)}]})
    validation.check("concentric spheres a 1 mm, b 2 mm", "capacitance", r["c0"][0] * 1e12, 4 * math.pi * EPS0 * a * b / (b - a) * 1e9,
                     rel=2e-3, unit="pF", source="exact")

    def coax(length):
        return kairoq.cross_section({"axisymmetric": True, "margin": 0.3, "target_error": 0.01,
                                       "conductors": [{"name": "in", "rect_mm": [0, 0, 0.15, length]},
                                                      {"name": "out", "reference": True, "rect_mm": [0.6, -1, 0.7, length + 1]}]})["c0"][0]
    validation.check("coaxial line a 0.15 mm, b 0.6 mm, axisymmetric", "capacitance of 4 mm more", (coax(8.0) - coax(4.0)) * 1e12,
                     2 * math.pi * EPS0 * 4e-3 / math.log(4.0) * 1e12, rel=2e-3, unit="pF", source="exact")
    # Filled with a dielectric: the capacitances at a frequency carry the permittivity.
    r = kairoq.cross_section({"axisymmetric": True, "margin": 0.3, "target_error": 0.01,
                                "conductors": [{"name": "in", "rect_mm": [0, 0, 0.15, 4.0]},
                                               {"name": "out", "reference": True, "rect_mm": [0.6, -1, 0.7, 5.0]}],
                                "dielectrics": [{"rect_mm": [0, -2, 2, 6], "epsilon_r": 4.5, "wideband": False}]})
    validation.check("coaxial line, axisymmetric, filled with er 4.5", "capacitance against the vacuum's", r["frequencies"][0]["c"][0] / r["c0"][0],
                     4.5, rel=1e-6, source="uniform fill")


def test_microstrip_dispersion(validation):
    """Kirschning-Jansen dispersion of a 50 ohm microstrip against scikit-rf's implementation
    of the same formulas, fed with the solver's static values; a stripline is unchanged."""
    mline = pytest.importorskip("skrf.media.mline")
    w, h, er, fs = 0.35, 0.2, 4.5, [0.1e9, 10e9, 20e9, 40e9]
    base = {"conductors": [{"name": "s", "rect_mm": [-w / 2, h, w / 2, h + 0.035]}, {"name": "g", "reference": True, "rect_mm": [-3, -0.035, 3, 0]}],
            "dielectrics": [{"rect_mm": [-3, 0, 3, h], "epsilon_r": er, "wideband": False}], "frequencies_hz": fs, "target_error": 0.01}
    qs = kairoq.cross_section(base)
    ms = kairoq.cross_section({**base, "microstrip": {"width_mm": w, "height_mm": h, "epsilon_r": er, "wideband": False}})
    for k, f in enumerate(fs):
        a, b = qs["frequencies"][k]["lines"][0], ms["frequencies"][k]["lines"][0]
        fn = min(f * h * 1e-9, 20.0)
        ef = mline.kirsching_er(w / h, fn, er, a["epsilon_eff"])
        zf, _ = mline.kirsching_zl(w / h, fn, er, a["epsilon_eff"], ef, a["z0"][0])
        case = "microstrip w 0.35, h 0.2 mm, er 4.5, %g GHz" % (f / 1e9)
        validation.check(case, "effective permittivity", b["epsilon_eff"], ef, rel=1e-3, source="scikit-rf, Kirschning-Jansen")
        validation.check(case, "impedance", b["z0"][0], zf, rel=1e-3, unit="ohm", source="scikit-rf, Jansen-Kirschning")


def test_stripline(validation):
    w, b, er, W = 0.2, 0.5, 4.0, 10.0
    spec = {"conductors": [{"name": "s", "rect_mm": [-w / 2, b / 2 - T / 2, w / 2, b / 2 + T / 2]},
                           {"name": "g1", "rect_mm": [-W, -0.035, W, 0], "reference": True},
                           {"name": "g2", "rect_mm": [-W, b, W, b + 0.035], "reference": True}],
            "dielectrics": [{"rect_mm": [-W, 0, W, b], "epsilon_r": er, "loss_tangent": 0.02, "wideband": False}],
            "frequencies_hz": [1e9], "target_error": 0.02, "margin": 0.2}
    r = kairoq.cross_section(spec)
    validation.check("stripline w 0.2 mm, b 0.5 mm, er 4", "impedance", z_lossless(r), stripline_exact(w, b, er), rel=4e-3, unit="ohm",
                     source="exact for zero thickness")
    fr = r["frequencies"][0]
    # Homogeneous dielectric: TEM in the material. (The thin strip's resistance changes
    # the lossy propagation constant; the material quantities come from C and G.)
    validation.check("stripline", "effective permittivity", fr["c"][0] / r["c0"][0], er, rel=1e-4, source="homogeneous")
    alpha_d = 20 / math.log(10) * math.pi * 1e9 * math.sqrt(er) * 0.02 / C0
    validation.check("stripline", "dielectric loss at 1 GHz", 20 / math.log(10) * fr["g"][0] * z_lossless(r) / 2, alpha_d, rel=4e-3,
                     unit="dB/m", source="pi f sqrt(er) tan d / c")


@pytest.mark.parametrize("w", [0.1, 0.3, 1.0])
def test_microstrip(w, validation):
    h, er = 0.2, 4.3
    W = 40 * h
    spec = {"conductors": [{"name": "s", "rect_mm": [-w / 2, h, w / 2, h + T]}, {"name": "g", "rect_mm": [-W, -0.035, W, 0], "reference": True}],
            "dielectrics": [{"rect_mm": [-W, 0, W, h], "epsilon_r": er}], "frequencies_hz": [1e9], "target_error": 0.02}
    r = kairoq.cross_section(spec)
    z, ee = microstrip_hj(w, h, er)
    case = f"microstrip w/h {w / h:g}, er 4.3"
    validation.check(case, "impedance", z_lossless(r), z, rel=4e-3, unit="ohm", source="Hammerstad and Jensen")
    validation.check(case, "effective permittivity", r["frequencies"][0]["c"][0] / r["c0"][0], ee, rel=3e-3, source="Hammerstad and Jensen")


def test_coupled_stripline(validation):
    w, s, b, er, W = 0.15, 0.15, 0.5, 4.0, 10.0
    spec = {"conductors": [{"name": "p", "rect_mm": [-s / 2 - w, b / 2 - T / 2, -s / 2, b / 2 + T / 2]},
                           {"name": "n", "rect_mm": [s / 2, b / 2 - T / 2, s / 2 + w, b / 2 + T / 2]},
                           {"name": "g1", "rect_mm": [-W, -0.035, W, 0], "reference": True},
                           {"name": "g2", "rect_mm": [-W, b, W, b + 0.035], "reference": True}],
            "dielectrics": [{"rect_mm": [-W, 0, W, b], "epsilon_r": er}], "frequencies_hz": [1e9], "target_error": 0.02, "margin": 0.2}
    p = kairoq.cross_section(spec)["frequencies"][0]["pair"]
    ze, zo = coupled_stripline_exact(w, s, b, er)
    case = "edge-coupled stripline w 0.15, s 0.15, b 0.5 mm"
    validation.check(case, "odd-mode impedance", p["z_odd"], zo, rel=4e-3, unit="ohm", source="Cohn, exact")
    validation.check(case, "even-mode impedance", p["z_even"], ze, rel=4e-3, unit="ohm", source="Cohn, exact")
    assert p["z_differential"] == pytest.approx(2 * p["z_odd"])


def test_wheeler_rule(validation):
    """The high-frequency resistance against an independent finite difference: Wheeler's
    rule R = Rs / mu0 dL/dn, with the inductance of the conductors grown and shrunk by 2 um
    on every face."""
    w, h, t, W, d = 0.2, 0.2, 0.035, 5.0, 0.002
    f = 1e11  # skin depth 0.2 um: the high-frequency limit
    def spec(g):
        return {"conductors": [{"name": "s", "rect_mm": [-w / 2 - g, h - g, w / 2 + g, h + t + g]},
                               {"name": "g", "rect_mm": [-W, -0.035, W, g], "reference": True}],
                "frequencies_hz": [f], "target_error": 0.006, "max_triangles": 600000}
    r = kairoq.cross_section(spec(0.0))
    lp, lm = (kairoq.cross_section(spec(s * d))["l_ext"][0] for s in (+1, -1))
    rs = math.sqrt(math.pi * f * MU0 / 5.8e7)
    wheeler = rs / MU0 * (lm - lp) / (2 * d * 1e-3)
    validation.check("microstrip w 0.2, h 0.2, t 0.035 mm in air", "resistance at 100 GHz", r["frequencies"][0]["r"][0], wheeler, rel=0.02,
                     unit="ohm/m", source="Wheeler's rule by finite difference")


def test_dielectric_model():
    from math import isclose
    spec = {"conductors": [{"name": "s", "rect_mm": [-0.1, 0.2, 0.1, 0.235]}, {"name": "g", "rect_mm": [-2, -0.035, 2, 0], "reference": True}],
            "dielectrics": [{"rect_mm": [-2, 0, 2, 0.5], "epsilon_r": 4.0, "loss_tangent": 0.02, "frequency_hz": 1e9}],
            "frequencies_hz": [1e6, 1e9, 1e10], "target_error": 0.05}
    # Fully filled: the effective permittivity is the material's, at each frequency.
    spec["dielectrics"][0]["rect_mm"] = [-20, -20, 20, 20]
    spec["margin"] = 0.1
    r = kairoq.cross_section(spec)
    ee = [f["c"][0] / r["c0"][0] for f in r["frequencies"]]
    tan = [f["g"][0] / (2 * math.pi * f["frequency_hz"] * f["c"][0]) for f in r["frequencies"]]
    assert ee[1] == pytest.approx(4.0, rel=1e-3) and tan[1] == pytest.approx(0.02, rel=1e-3)
    assert ee[0] > ee[1] > ee[2]  # the permittivity falls with frequency (causality with loss)
    assert all(isclose(t, 0.02, rel_tol=0.15) for t in tan)  # nearly constant loss tangent


def test_roughness_models():
    """The factors at their defining points: Hammerstad 1 + 2/pi atan(1.4) where the skin
    depth equals the RMS height; Huray 1 + ratio / (1 + d/a + d^2/2a^2) at skin depth d."""
    sigma = 5.8e7
    f_eq = 1 / (math.pi * MU0 * sigma * (1e-6) ** 2)  # skin depth 1 um
    base = {"conductors": [{"name": "s", "rect_mm": [-0.1, 0.2, 0.1, 0.235]}, {"name": "g", "rect_mm": [-2, -0.035, 2, 0], "reference": True}],
            "target_error": 0.05}
    def r_at(f, roughness):
        spec = dict(base, frequencies_hz=[f])
        spec["conductors"] = [dict(c, roughness=roughness) for c in base["conductors"]]
        return kairoq.cross_section(spec)["frequencies"][0]["r"][0]
    smooth = r_at(f_eq, {"model": "none"})
    rough = r_at(f_eq, {"model": "hammerstad", "rms_um": 1.0})
    assert rough / smooth == pytest.approx(1 + 2 / math.pi * math.atan(1.4), rel=2e-3)
    delta = 1 / math.sqrt(math.pi * 1e13 * MU0 * sigma)
    huray = 1 + 2.0 / (1 + delta / 0.5e-6 + delta**2 / (2 * 0.5e-6**2))
    assert r_at(1e13, {"model": "huray", "radius_um": 0.5, "ratio": 2.0}) / r_at(1e13, {"model": "none"}) == pytest.approx(huray, rel=2e-3)


# ------------------------------------------------------------------------ boards


def test_board_microstrip(tmp_path):
    """Every cut along a straight microstrip is the same cross-section; its delay sums
    over the length; the distances run from the first pad."""
    synthetic.microstrip(tmp_path)
    b = kairoq.open(tmp_path)
    r = kairoq.analyse_si(b, {"nets": ["SIG"]})
    (net,) = r["nets"]
    assert len(r["cross_sections"]) == 1 and r["cuts"] >= 18
    assert net["length_mm"] == pytest.approx(20.0, rel=1e-9)
    s = r["samples"][0]
    assert net["delay_ps"] == pytest.approx(20.0 * s["delay_ps_per_mm"], rel=1e-9)
    assert net["z_min"] == pytest.approx(net["z_max"])
    d = sorted(x["distance_mm"] for x in r["samples"])
    assert d[0] > 0 and d[-1] < 20.0 and all(b > a for a, b in zip(d, d[1:]))
    assert any("loss tangent" in w for w in r["warnings"])


def test_board_pair(tmp_path, validation):
    """A pair found by its names, its differential impedance against a cross-section of
    the same geometry (with mask and copper thickness as on the board)."""
    g = synthetic.diff_pair(tmp_path)
    b = kairoq.open(tmp_path)
    r = kairoq.analyse_si(b, {"targets": [{"ohm": 100, "differential": True}]})
    nets = {n["name"]: n for n in r["nets"]}
    assert nets["DP_P"]["partner"] == "DP_N" and nets["DP_N"]["partner"] == "DP_P"
    assert nets["DP_P"]["coupled_mm"] == pytest.approx(nets["DP_P"]["length_mm"])
    # The same cut by hand: tracks 35 um on 0.2 mm of er 4.5 over the plane, mask 10 um
    # of er 3.8 over all (the Gerber job file gives no mask permittivity).
    w, s, h, t, m = 0.2, 0.2, 0.2, 0.035, 0.01
    xp = (s / 2, s / 2 + w)
    zb, zt = 0.01 + 0.035 + h, 0.01 + 0.035 + h + t
    mask = [[-2.5, zb], [2.5, zb], [2.5, zb + m], [xp[1] + m, zb + m], [xp[1] + m, zt + m], [xp[0] - m, zt + m], [xp[0] - m, zb + m],
            [-xp[0] + m, zb + m], [-xp[0] + m, zt + m], [-xp[1] - m, zt + m], [-xp[1] - m, zb + m], [-2.5, zb + m]]
    spec = {"conductors": [{"name": "p", "rect_mm": [xp[0], zb, xp[1], zt]}, {"name": "n", "rect_mm": [-xp[1], zb, -xp[0], zt]},
                           {"name": "g", "rect_mm": [-2.5, 0.01, 2.5, 0.045], "reference": True}],
            "dielectrics": [{"rect_mm": [-2.5, 0.045, 2.5, zb], "epsilon_r": 4.5, "loss_tangent": 0.02},
                            {"polygon_mm": mask, "epsilon_r": 3.8, "loss_tangent": 0.025}],
            "frequencies_hz": [1e9], "target_error": 0.03, "margin": 0.5}
    ref = kairoq.cross_section(spec)["frequencies"][0]["pair"]["z_differential"]
    validation.check("synthetic pair w 0.2, s 0.2, h 0.2 mm", "differential impedance, board against cross-section",
                     nets["DP_P"]["zd_mean"], ref, rel=5e-3, unit="ohm", source="cross-section of the same geometry")


def test_board_pair_over_slot(tmp_path):
    synthetic.diff_pair(tmp_path, slot=1.0)
    b = kairoq.open(tmp_path)
    r = kairoq.analyse_si(b, {"nets": ["DP_*"], "cut": {"step_mm": 0.25}, "min_finding_mm": 0.2})
    kinds = {(f["net"], f["kind"]) for f in r["findings"]}
    assert ("DP_P", "no reference") in kinds and ("DP_N", "no reference") in kinds
    gap = [f for f in r["findings"] if f["kind"] == "no reference" and f["net"] == "DP_P"]
    assert abs(gap[0]["position_mm"][0] - 15.0) < 1.0 and gap[0]["length_mm"] == pytest.approx(1.0, abs=0.3)
    # The plane on either side of the slot is the same layer: no change of reference.
    assert not [f for f in r["findings"] if f["kind"] == "reference change"]


def class_pair(w, s, h=0.0994, t=0.035, m=0.01, er=4.45, em=3.3):
    """A coupled microstrip of the cm5_minima stackup: top copper on 0.0994 mm of 2116
    prepreg over the plane, with mask."""
    xp = (s / 2, s / 2 + w)
    zb, zt = h, h + t
    mask = [[-3, zb], [3, zb], [3, zb + m], [xp[1] + m, zb + m], [xp[1] + m, zt + m], [xp[0] - m, zt + m], [xp[0] - m, zb + m],
            [-xp[0] + m, zb + m], [-xp[0] + m, zt + m], [-xp[1] - m, zt + m], [-xp[1] - m, zb + m], [-3, zb + m]]
    spec = {"conductors": [{"name": "p", "rect_mm": [xp[0], zb, xp[1], zt]}, {"name": "n", "rect_mm": [-xp[1], zb, -xp[0], zt]},
                           {"name": "g", "rect_mm": [-3, -0.0152, 3, 0], "reference": True}],
            "dielectrics": [{"rect_mm": [-3, 0, 3, zb], "epsilon_r": er, "loss_tangent": 0.02}, {"polygon_mm": mask, "epsilon_r": em}],
            "frequencies_hz": [1e9], "target_error": 0.03}
    return kairoq.cross_section(spec)["frequencies"][0]["pair"]["z_differential"]


def test_design_classes(validation):
    """The net classes of cm5_minima: 100 ohm with tracks of 0.13 mm 0.25 mm apart, 90 ohm
    with 0.14 mm 0.154 mm apart, on its stackup."""
    validation.check("cm5_minima class 100ohm (w 0.13, s 0.25 mm)", "differential impedance", class_pair(0.13, 0.25), 100.0, rel=0.04,
                     unit="ohm", source="net class of the design")
    validation.check("cm5_minima class 90ohm (w 0.14, s 0.154 mm)", "differential impedance", class_pair(0.14, 0.154), 90.0, rel=0.04,
                     unit="ohm", source="net class of the design")


def test_real_board_pairs(validation):
    """cm5_minima: the differential impedance of its routed pairs, each group's mean over
    the coupled length, against the net classes of the design (100 ohm: Ethernet, HDMI,
    MIPI D-PHY; 90 ohm: USB, PCIe). The routes depart from the class geometry (ground
    between the members, other spacings, fan-outs); the tolerance is wide."""
    b = kairoq.open(TESTDATA / "kicad" / "cm5_minima" / "cm5_minima.ipc2581.zip")
    r = kairoq.analyse_si(b, {"nets": ["*_P", "*_N"]})
    groups = {"ETH": 100, "HDMI": 100, "DPHY": 100, "USB": 90, "PCIE": 90}
    for key, target in groups.items():
        nets = [n for n in r["nets"] if key in n["name"].upper() and n["zd_mean"] > 0]
        assert nets, key
        coupled = sum(n["coupled_mm"] for n in nets)
        mean = sum(n["zd_mean"] * n["coupled_mm"] for n in nets) / coupled
        validation.check(f"cm5_minima {key} pairs", "differential impedance, mean over the coupled length", mean, target, rel=0.10,
                         unit="ohm", source="net class of the design")


def test_etch_face(tmp_path, validation):
    """A symmetric stripline etched by 10 um per side has one impedance whichever face is
    the narrow one; the etch raises it."""
    synthetic.coupled_lines(tmp_path, stripline=True)
    b = kairoq.open(tmp_path)
    z = lambda **e: kairoq.analyse_si(b, {"nets": ["A"], "pairs": [], **e}).samples[0]["z0"]
    plain, top, bottom = z(), z(etch_um=10, etch_face={"L2": "top"}), z(etch_um=10, etch_face={"L2": "bottom"})
    validation.check("symmetric stripline, etch 10 um", "impedance, bottom face against top face", bottom, top, rel=1e-3, unit="ohm",
                     source="symmetry")
    assert top > plain * 1.01
    with pytest.raises(Exception, match="top or bottom"):
        z(etch_um=10, etch_face={"L2": "upper"})
