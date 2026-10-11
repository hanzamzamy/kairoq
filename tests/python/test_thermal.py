"""Thermal solver: analytical cases, heat balance and electrothermal coupling."""

from __future__ import annotations

import math

import numpy as np
import pytest

import kairoq
import synthetic
from conftest import ROOT, TESTDATA, cli_path

K_CU, K_XY, K_Z, K_MASK = 385.0, 0.8, 0.3, 0.25
SIGMA_SB = 5.670374419e-8


def surfaces(h=10.0, eps=0.0):
    return {"top": {"convection_w_m2k": h, "emissivity": eps}, "bottom": {"convection_w_m2k": h, "emissivity": eps}}


@pytest.fixture(scope="module")
def slab(tmp_path_factory):
    d = tmp_path_factory.mktemp("slab")
    return synthetic.slab(d), kairoq.open(d)


def two_node(g, q, h1, h2, gz, ta):
    a = np.array([[h1 + gz, -gz], [-gz, gz + h2]])
    return np.linalg.solve(a, [q, 0.0]) + ta


def test_slab_conduction_is_exact(slab, validation):
    """Uniform heat over the whole top: one-dimensional through the board."""
    g, b = slab
    p, ta = 2.0, 25.0
    r = kairoq.solve_thermal(b, {"ambient_c": ta, **surfaces(), "sources": [{"name": "H", "pad": "H1.1", "power_w": p}]})
    h1 = 1 / (g["mask_m"] / K_MASK + 1 / 10.0)
    t1, t2 = two_node(g, p / g["area_m2"], h1, h1, K_Z / g["dielectric_m"], ta)
    validation.check("slab 20 x 20 mm, 2 W", "top temperature", float(r.layer(1)["temperature_c"].max()), t1, rel=1e-10,
                     unit="C", source="series resistances")
    validation.check("slab 20 x 20 mm, 2 W", "bottom temperature", float(r.layer(2)["temperature_c"].min()), t2, rel=1e-10,
                     unit="C", source="series resistances")
    validation.check("slab 20 x 20 mm, 2 W", "heat out - heat in", r.heat_top_w + r.heat_bottom_w - p, 0.0, abs_=1e-9, unit="W",
                     source="energy balance")


def test_junction_resistance(slab):
    g, b = slab
    r = kairoq.solve_thermal(b, {**surfaces(), "sources": [{"name": "H", "pad": "H1.1", "power_w": 2.0, "junction_k_per_w": 5.0}]})
    s = r.sources["H"]
    assert s["junction_c"] == pytest.approx(s["mean_c"] + 10.0, rel=1e-9)


def test_radiation(slab, validation):
    """Radiation only (no convection): the slab's two-node problem with T^4 terms."""
    g, b = slab
    p, ta, eps = 2.0, 25.0, 0.9
    r = kairoq.solve_thermal(b, {"ambient_c": ta, **surfaces(h=0.0, eps=eps), "coupling_tolerance_k": 1e-4,
                                   "sources": [{"name": "H", "pad": "H1.1", "power_w": p}]})
    q, gz, rm = p / g["area_m2"], K_Z / g["dielectric_m"], g["mask_m"] / K_MASK
    tak = ta + 273.15

    def flux(t):  # sheet at t through the mask to a radiating surface; surface by Newton
        ts = t
        for _ in range(100):
            f = (t - ts) / rm - eps * SIGMA_SB * ((ts + 273.15) ** 4 - tak ** 4)
            df = -1 / rm - 4 * eps * SIGMA_SB * (ts + 273.15) ** 3
            ts -= f / df
        return (t - ts) / rm

    t = np.array([ta + 50.0, ta + 40.0])
    for _ in range(200):  # Newton on the two sheet temperatures
        f = np.array([q - flux(t[0]) - gz * (t[0] - t[1]), gz * (t[0] - t[1]) - flux(t[1])])
        e = 1e-6
        jac = np.zeros((2, 2))
        for k in range(2):
            dt = np.zeros(2)
            dt[k] = e
            tp = t + dt
            fp = np.array([q - flux(tp[0]) - gz * (tp[0] - tp[1]), gz * (tp[0] - tp[1]) - flux(tp[1])])
            jac[:, k] = (fp - f) / e
        t -= np.linalg.solve(jac, f)
    validation.check("slab, radiation only, emissivity 0.9", "top temperature", float(r.layer(1)["temperature_c"].max()), t[0],
                     rel=1e-5, unit="C", source="Stefan-Boltzmann, two nodes")
    assert r.heat_top_w + r.heat_bottom_w == pytest.approx(p, rel=1e-5)


@pytest.mark.parametrize("target", [0.1, 0.03])
def test_fin(tmp_path, validation, target):
    """A strip held at a fixed temperature at one end, cooled on both faces. The in-plane
    conduction is in the copper alone (the dielectric's would lie at a lower temperature
    than the closed form's lumped cross-section assumes). Elements of at most 0.5 mm keep
    the error of the lumped surface term, about (h m)^2 / 12, below 1e-4."""
    g = synthetic.fin(tmp_path)
    b = kairoq.open(tmp_path)
    core = next(l["name"] for l in b.stackup if l["type"] == "dielectric")
    ta, t0, h = 25.0, 75.0, 10.0
    r = kairoq.solve_thermal(b, {"ambient_c": ta, "target_error": target, **surfaces(h), "mesh": {"max_edge_mm": 0.5},
                                   "layers": {core: {"k_xy": 1e-9, "k_z": K_Z}},
                                   "sources": [{"name": "HOT", "pad": "J1.1", "temperature_c": t0}]})
    G = K_CU * 35e-6
    H = 1 / (g["mask_m"] / K_MASK + 1 / h) + 1 / (g["dielectric_m"] / K_Z + g["mask_m"] / K_MASK + 1 / h)
    m = math.sqrt(H / G)
    ell = g["length_m"] - g["pad_m"]
    theta = t0 - ta
    # Into the fin beyond the pad, plus the pad's own faces.
    q = G * g["width_m"] * theta * m * math.tanh(m * ell) + H * theta * g["width_m"] * g["pad_m"]
    validation.check(f"fin 5 x 60 mm, target {target}", "heat from the fixed end", r.sources["HOT"]["power_w"], q,
                     rel=max(r.estimated_error ** 2, 1e-4), unit="W", source="fin equation, insulated tip")
    validation.check(f"fin 5 x 60 mm, target {target}", "tip temperature", float(r.layer(1)["temperature_c"].min()),
                     ta + theta / math.cosh(m * ell), rel=max(r.estimated_error ** 2, 1e-4), unit="C", source="fin equation")
    assert r.heat_top_w + r.heat_bottom_w == pytest.approx(r.sources["HOT"]["power_w"], rel=1e-8)


def test_via_spreading(tmp_path, validation):
    """Both rims held at fixed temperatures; heat flows through copper and one via only.
    (Isothermal rims, as in the closed form: heat put in as a flux over the rim would add
    the conduction within the rim itself.)"""
    g = synthetic.via_spreading(tmp_path)
    b = kairoq.open(tmp_path)
    core = next(l["name"] for l in b.stackup if l["type"] == "dielectric")
    # The core conducts almost nothing (exactly nothing would leave the dielectric
    # outside the discs without any conductance); its share is below 1e-6 of the path.
    r = kairoq.solve_thermal(b, {**surfaces(h=0.0), "target_error": 0.03, "layers": {core: {"k_xy": 1e-3, "k_z": 1e-6}},
                                   "sources": [{"name": "IN", "pad": "P1.1", "temperature_c": 125.0},
                                               {"name": "OUT", "pad": "P2.1", "temperature_c": 25.0}]})
    t, tp, d = 35e-6, 18e-6, 2 * g["r_via"] * 1e-3
    resistance = 2 * math.log(g["r_inner_port"] / g["r_via"]) / (2 * math.pi * K_CU * t) + \
        g["barrel_length"] * 1e-3 / (K_CU * math.pi * tp * (d + tp))
    q = r.sources["IN"]["power_w"]
    validation.check("via spreading, copper only", "thermal resistance", 100.0 / q, resistance,
                     rel=max(r.estimated_error ** 2, 1e-3), unit="K/W", source="ln(r2/r1)/(pi k t) + barrel")
    assert r.sources["OUT"]["power_w"] == pytest.approx(-q, rel=1e-6)


# Three-dimensional reference (tools/ref3d.py, docs/validation/thermal.md): mean rise of the
# source copper in K, with the gaps of the outer copper filled with air and with solder
# mask, and the relative uncertainty of each reference.
# Heat sink cases are compared in docs/validation/thermal-transients.md: at the default target error their
# mesh error is of the order of the band, and a smaller target costs minutes.
REF3D = [
    ("trace 1 x 20 mm", (5.0, 14.5, 25.0, 15.5), False, 1.5, False, (83.09, 0.003), (81.61, 0.002)),
    ("square 2 mm over a plane", (14.0, 14.0, 16.0, 16.0), True, 1.5, False, (204.26, 0.006), (197.56, 0.004)),
    ("square 1 mm over a plane", (14.5, 14.5, 15.5, 15.5), True, 1.5, False, (423.03, 0.013), (398.64, 0.008)),
]


@pytest.mark.parametrize("case", REF3D, ids=[c[0] for c in REF3D])
def test_against_3d_reference(tmp_path, validation, case):
    """0.5 W in a pad or trace on a 30 x 30 mm board: the layered model lies between the two
    three-dimensional references."""
    name, pad, plane, dielectric, sink, (air, u_air), (mask, u_mask) = case
    synthetic.source_board(tmp_path, pad=pad, plane=plane, dielectric=dielectric)
    b = kairoq.open(tmp_path)
    sinks = [{"name": "S", "tim": {"thickness_mm": 0.5, "k": 3.0}, "temperature_c": 25.0}] if sink else []
    r = kairoq.solve_thermal(b, {**surfaces(), "sources": [{"name": "H", "pad": "H1.1", "power_w": 0.5}], "heatsinks": sinks})
    lo, hi = mask * (1 - u_mask), air * (1 + u_air)
    mid = 0.5 * (lo + hi)
    validation.check(f"{name}, 0.5 W", "mean rise of the source", r.sources["H"]["mean_c"] - 25.0, mid, rel=(hi - lo) / (2 * mid),
                     unit="K", source="3D reference, between air and mask in the gaps")


def test_electrothermal_strip(tmp_path, validation):
    """Long strip carrying 10 A: in the middle, Joule heat with rho(T) balances the faces."""
    g = synthetic.fin(tmp_path, width=5.0, length=400.0)
    b = kairoq.open(tmp_path)
    ta, h, i = 25.0, 10.0, 10.0
    r = kairoq.solve_thermal(b, {"ambient_c": ta, **surfaces(h)},
                               dc={"ports": [{"name": "A", "pad": "J1.1", "current_a": i}, {"name": "B", "pad": "J2.1", "current_a": -i}]})
    H = 1 / (g["mask_m"] / K_MASK + 1 / h) + 1 / (g["dielectric_m"] / K_Z + g["mask_m"] / K_MASK + 1 / h)
    w, t, alpha = g["width_m"], 35e-6, 0.00393
    a = i * i / 5.8e7 / (w * w * t * H)
    theta = a * (1 + alpha * (ta - 20)) / (1 - a * alpha)
    f = r.layer(1)
    mid = np.abs(f["points_mm"][:, 0] - 200.0) < 5.0
    validation.check("strip 5 mm x 400 mm, 10 A", "temperature at the middle", float(f["temperature_c"][mid].mean()), ta + theta,
                     rel=2e-3, unit="C", source="local balance with rho(T)")
    assert r.coupling_iterations > 1
    assert r.heat_top_w + r.heat_bottom_w == pytest.approx(r.heat_joule_w, rel=1e-4)


def test_coupled_board_balance(validation):
    r = kairoq.solve_thermal(kairoq.open(TESTDATA / "kicad" / "cm5_minima" / "gerber"),
                               str(ROOT / "examples" / "cm5_minima_thermal.yaml"))
    heat_in = r.heat_sources_w + r.heat_joule_w
    validation.check("cm5_minima, 3 A and 2 W", "(heat out - heat in) / heat in",
                     (r.heat_top_w + r.heat_bottom_w + r.heat_fixed_w - heat_in) / heat_in, 0.0, abs_=1e-3,
                     source="energy balance")
    assert r.estimated_error <= 0.1
    assert r.dc is not None and r.coupling_iterations >= 2


def test_settings_are_checked(slab):
    _, b = slab
    with pytest.raises(RuntimeError, match="either a power .power_w or power_schedule. or temperature_c"):
        kairoq.solve_thermal(b, {"sources": [{"name": "H", "pad": "H1.1"}]})
    with pytest.raises(RuntimeError, match="no path for the heat"):
        kairoq.solve_thermal(b, {**surfaces(h=0.0), "sources": [{"name": "H", "pad": "H1.1", "power_w": 1.0}]})
    with pytest.raises(RuntimeError, match="unknown key"):
        kairoq.solve_thermal(b, {"ambient": 25})


def test_coupling_acceleration_and_runaway(tmp_path, validation):
    """The strip of test_electrothermal_strip at a feedback gain a alpha of 0.8: the
    accelerated coupling settles in a few iterations on the plain iteration's temperature;
    above a gain of 1 there is no steady state, which is reported."""
    g = synthetic.fin(tmp_path, width=5.0, length=400.0)
    b = kairoq.open(tmp_path)
    ta, h, alpha = 25.0, 10.0, 0.00393
    H = 1 / (g["mask_m"] / K_MASK + 1 / h) + 1 / (g["dielectric_m"] / K_Z + g["mask_m"] / K_MASK + 1 / h)
    current = lambda gain: (gain / alpha * 5.8e7 * g["width_m"] ** 2 * 35e-6 * H) ** 0.5
    def run(gain, acc):
        i = current(gain)
        return kairoq.solve_thermal(b, {"ambient_c": ta, **surfaces(h), "coupling_acceleration": acc, "max_coupling_iterations": 60},
                                      dc={"ports": [{"name": "A", "pad": "J1.1", "current_a": i}, {"name": "B", "pad": "J2.1", "current_a": -i}]})
    plain, fast = run(0.8, False), run(0.8, True)
    mid = lambda r: float(r.layer(1)["temperature_c"][np.abs(r.layer(1)["points_mm"][:, 0] - 200.0) < 5.0].mean())
    validation.check("strip at a feedback gain of 0.8", "temperature at the middle, accelerated", mid(fast), mid(plain), abs_=0.5, unit="C",
                     source="plain fixed-point iteration")
    assert fast.coupling_iterations <= 12 < plain.coupling_iterations
    for acc in (False, True):
        assert any("thermal runaway" in w for w in run(1.05, acc).warnings)


@pytest.mark.parametrize("fill,ref", [("air", (200.3, 0.007)), ("mask", (165.8, 0.002))])
def test_thick_outer_copper_gap_fill(tmp_path, validation, fill, ref):
    """A 2 mm pad of 10 oz copper over a plane through 1.5 mm, 0.5 W: with the gaps of the
    outer copper filled with air or with mask, against the three-dimensional reference
    with the same fill (docs/validation/thermal.md), whose two values differ by 21 %. The model's
    own difference to the reference is about 2 %; a target error of 0.05 keeps the
    discretisation's share below 0.1 % (0.6 % at the default)."""
    synthetic.source_board(tmp_path, pad=(14.0, 14.0, 16.0, 16.0), plane=True, dielectric=1.5, copper=0.347)
    r = kairoq.solve_thermal(kairoq.open(tmp_path), {**surfaces(), "outer_gaps": fill, "target_error": 0.05,
                                                        "sources": [{"name": "H", "pad": "H1.1", "power_w": 0.5}]})
    validation.check(f"2 mm pad, 10 oz, gaps of {fill}", "mean rise of the source", r.sources["H"]["mean_c"] - 25.0, ref[0], rel=0.02 + ref[1],
                     unit="K", source=f"3D reference, {fill} in the gaps")
