"""Observed orders of convergence under uniform refinement (docs/validation/dc.md,
docs/validation/thermal.md).

A solution that agrees with a reference on one mesh may do so by chance. Here a quantity is
computed on three meshes whose element size halves, without adaptive refinement, and the
order of convergence is measured from the differences alone,

    p = log2 |(Q_h - Q_h/2) / (Q_h/2 - Q_h/4)|,

which needs no reference value. Linear elements converge in the energy, and so in a
resistance between fixed potentials, as h^2 where the solution is smooth. A corner of an
equipotential boundary with the interior angle alpha > pi (seen from the copper) makes the
solution singular as r^(pi / alpha) and the order 2 pi / alpha. Circles are polygons: a
circular port of n sides has corners of pi + 2 pi / n, so the order is 2n / (n + 2); a plated
hole is a 16-sided polygon (order 1.78), a circle of the import, at the arc tolerance, has
more sides.
"""

from __future__ import annotations

import math

import kairoq
import synthetic
from test_dc import spreading_exact
from test_thermal import K_CU, surfaces

ARC_TOLERANCE = 1e-3  # mm, the import's default


def observed_order(q):
    return math.log2(abs((q[0] - q[1]) / (q[1] - q[2])))


def polygon_order(radius_mm, tol_mm):
    """2n / (n + 2) for the polygon of a circle whose chords lie within tol of it."""
    n = math.pi / math.acos(1 - tol_mm / radius_mm)
    return 2 * n / (n + 2)


PORTS = [{"name": "A", "pad": "P1.1", "voltage_v": 1e-3}, {"name": "B", "pad": "P2.1", "voltage_v": 0.0}]


def test_dc_annulus_order(tmp_path, validation):
    """Radial current between two circular ports on one layer."""
    g = synthetic.annulus(tmp_path)
    b = kairoq.open(tmp_path)
    r = []
    for h in (0.1, 0.05, 0.025):
        s = kairoq.solve_dc(b, {"ports": PORTS, "max_passes": 1, "mesh": {"max_edge_mm": h}})
        r.append(1e-3 / s.ports["A"]["current_a"])
    assert r[0] < r[1] < r[2]  # a lower bound that rises with refinement
    case = "annulus, meshes of 0.1, 0.05, 0.025 mm"
    validation.check(case, "observed order of the resistance", observed_order(r), polygon_order(g["r_inner_port"], ARC_TOLERANCE),
                     abs_=0.1, source="2n / (n + 2), the inner port a polygon of n sides")
    exact = math.log(g["r_outer_port"] / g["r_inner_port"]) / (2 * math.pi * 5.8e7 * 35e-6)
    p = observed_order(r)
    validation.check(case, "extrapolated resistance", r[2] + (r[2] - r[1]) / (2 ** p - 1), exact, rel=2e-5, unit="ohm",
                     source="ln(r2/r1) / (2 pi sigma t)")


def test_dc_plated_hole_order(tmp_path, validation):
    """Radial spreading through a via of 1 mm: the hole's 16 corners set the order."""
    g = synthetic.via_spreading(tmp_path, drill=1.0)
    b = kairoq.open(tmp_path)
    r = []
    for h in (0.1, 0.05, 0.025):
        s = kairoq.solve_dc(b, {"ports": PORTS, "max_passes": 1, "mesh": {"max_edge_mm": h}})
        r.append(1e-3 / s.ports["A"]["current_a"])
    assert r[0] < r[1] < r[2]
    case = "via spreading, drill 1 mm, meshes of 0.1, 0.05, 0.025 mm"
    validation.check(case, "observed order of the resistance", observed_order(r), 2 * 16 / 18, abs_=0.1,
                     source="2n / (n + 2), the hole a polygon of 16 sides")
    spread, barrel = spreading_exact(g)
    validation.check(case, "resistance, finest mesh", r[2], spread + barrel, rel=2e-3, unit="ohm",
                     source="ln(r2/r1) / (pi sigma t) + barrel, circular hole")


def test_thermal_annulus_order(tmp_path, validation):
    """Radial heat flow between two circular rims at fixed temperatures, copper only."""
    g = synthetic.annulus(tmp_path)
    b = kairoq.open(tmp_path)
    core = next(l["name"] for l in b.stackup if l["type"] == "dielectric")
    q = []
    for h in (0.2, 0.1, 0.05):
        r = kairoq.solve_thermal(b, {**surfaces(h=0.0), "max_passes": 1, "mesh": {"max_edge_mm": h},
                                     "layers": {core: {"k_xy": 1e-3, "k_z": 1e-6}},
                                     "sources": [{"name": "IN", "pad": "P2.1", "temperature_c": 125.0},
                                                 {"name": "OUT", "pad": "P1.1", "temperature_c": 25.0}]})
        q.append(100.0 / r.sources["IN"]["power_w"])
    assert q[0] < q[1] < q[2]
    case = "annulus, meshes of 0.2, 0.1, 0.05 mm"
    validation.check(case, "observed order of the thermal resistance", observed_order(q),
                     polygon_order(g["r_inner_port"], ARC_TOLERANCE), abs_=0.1, source="2n / (n + 2), the inner rim a polygon of n sides")
    exact = math.log(g["r_outer_port"] / g["r_inner_port"]) / (2 * math.pi * K_CU * 35e-6)
    validation.check(case, "thermal resistance, finest mesh", q[2], exact, rel=2e-4, unit="K/W", source="ln(r2/r1) / (2 pi k t)")
