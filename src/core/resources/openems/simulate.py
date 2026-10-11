#!/usr/bin/env python3
"""Full-wave simulation of a loop exported by KairoQ, with openEMS (FDTD).

Reads model.json from the directory of this script, builds the structure with CSXCAD and
runs openEMS. Writes result.json: the impedance at the port and the loop inductance over
frequency, and the largest far field per ampere of port current at the model's distance,
with the loop estimate's values where the model carries them.

Requires openEMS with its Python interface (modules CSXCAD and openEMS), numpy and h5py.

usage: python simulate.py [--threads N] [--setup-only] [--sim-dir DIR]
"""
import argparse
import json
import os
import sys

import numpy as np

UNIT = 1e-3  # drawing unit of the simulation: millimetre
ETA0 = 376.730313668
C0 = 299792458.0


def mesh_lines(fixed, exact, lo, hi, min_cell, board_cell, max_cell, margin, smooth):
    """Lines through the exact coordinates, and through the other fixed coordinates at
    least a minimum cell from those already taken; smoothed to board_cell within
    [lo, hi] and to max_cell in the air margin beyond it."""
    taken = sorted(set(exact))
    for p in sorted(set([lo, hi] + [p for p in fixed if lo <= p <= hi])):
        i = np.searchsorted(taken, p)
        near = [taken[j] for j in (i - 1, i) if 0 <= j < len(taken)]
        if all(abs(p - q) >= 0.999 * min_cell for q in near):
            taken.insert(i, p)
    inner = smooth(taken, board_cell, 1.4)
    lines = np.unique(np.asarray(smooth([lo - margin] + list(inner) + [hi + margin], max_cell, 1.4)))
    # Smoothing moves lines by rounding; a sheet or port off its line by 1e-16 is dropped
    # from the grid, so lines near an exact coordinate take it exactly.
    ex = np.unique(np.asarray(exact, dtype=float))
    if len(ex):
        k = np.clip(np.searchsorted(ex, lines), 1, len(ex) - 1)
        near = np.where(np.abs(ex[k - 1] - lines) < np.abs(ex[k] - lines), ex[k - 1], ex[k])
        lines = np.where(np.abs(near - lines) < 1e-6 * min_cell, near, lines)
    return np.unique(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--threads", type=int, default=0, help="threads of openEMS (0: all)")
    ap.add_argument("--setup-only", action="store_true", help="build the structure and mesh, do not run")
    ap.add_argument("--sim-dir", default=None, help="directory of the simulation files (default: sim/ next to the model)")
    args = ap.parse_args()

    from CSXCAD import ContinuousStructure
    from CSXCAD.SmoothMeshLines import SmoothMeshLines
    from openEMS import openEMS

    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "model.json")) as f:
        model = json.load(f)
    s = model["simulation"]
    sim_dir = os.path.abspath(args.sim_dir or os.path.join(here, "sim"))
    d = lambda v: np.asarray(v, dtype=float) / UNIT  # metres to drawing units
    hmin, hboard, hmax, margin = (s[k] / UNIT for k in ("min_cell", "board_cell", "max_cell", "air_margin"))
    f_max = s["f_max_hz"]

    fdtd = openEMS(NrTS=int(s["max_timesteps"]), EndCriteria=10 ** (s["end_criterion_db"] / 10))
    fdtd.SetGaussExcite(f_max / 2, f_max / 2)
    fdtd.SetBoundaryCond(["PML_8"] * 6)
    csx = ContinuousStructure()
    fdtd.SetCSX(csx)
    grid = csx.GetGrid()
    grid.SetDeltaUnit(UNIT)

    fx, fy, fz = [], [], []  # fixed mesh coordinates
    exact = ([], [], [])     # coordinates that must be lines: the sheets' planes, the port's gap

    # Dielectric under the window, lowest priority.
    for k, slab in enumerate(model["dielectric"]):
        kappa = 2 * np.pi * (f_max / 2) * 8.8541878128e-12 * slab["epsilon_r"] * slab["loss_tangent"]
        mat = csx.AddMaterial("dielectric_%d" % k, epsilon=slab["epsilon_r"], kappa=kappa)
        z0, z1 = d(slab["z0"]), d(slab["z1"])
        for poly in model["window"]:
            p = d(poly).T
            mat.AddLinPoly(p, "z", z0, z1 - z0, priority=1)
        fz += [z0, z1]
    for poly in model["window"]:
        p = d(poly)
        fx += list(p[:, 0])
        fy += list(p[:, 1])

    metal = csx.AddMetal("copper")
    # Copper sheets with lines at the edges and one minimum cell to either side, and at
    # the sheets' planes with one cell above and below.
    for sheet in model["sheets"]:
        z = float(d(sheet["z"]))
        fz += [z - hmin, z + hmin]
        exact[2].append(z)
        for poly in sheet["polygons"]:
            p = d(poly)
            metal.AddPolygon(p.T, "z", z, priority=10)
            for c, lst in ((p[:, 0], fx), (p[:, 1], fy)):
                lst += list(c) + list(c - hmin) + list(c + hmin)
    for v in model["vias"]:
        x, y, r = d(v["x"]), d(v["y"]), d(v["radius"])
        metal.AddCylinder([x, y, d(v["z0"])], [x, y, d(v["z1"])], r, priority=10)
        fx += [x - r, x, x + r]
        fy += [y - r, y, y + r]

    # Elements: a bar of square section, side twice the radius, standing on the copper from
    # pad centre to pad centre; the port is a gap of one minimum cell in the middle of the
    # first element.
    port = None
    for e in model["elements"]:
        a, c = d(e["from"]), d(e["to"])
        zs, za = float(d(e["z_sheet"])), float(d(e["z_axis"]))
        ztop = zs + 2 * (za - zs)
        r = float(d(e["radius"]))
        along = 0 if abs(c[1] - a[1]) < 1e-6 else 1 if abs(c[0] - a[0]) < 1e-6 else None
        fz += [zs, ztop]
        if along is None:
            metal.AddCurve([[a[0], a[0], c[0], c[0]], [a[1], a[1], c[1], c[1]], [zs, za, za, zs]], priority=10)
            fx += [a[0], c[0]]
            fy += [a[1], c[1]]
            continue
        across = 1 - along
        lo, hi = sorted((a[along], c[along]))
        mid = a[across]
        def box(t0, t1):
            p0, p1 = [0.0, 0.0, min(zs, ztop)], [0.0, 0.0, max(zs, ztop)]
            p0[along], p1[along] = t0, t1
            p0[across], p1[across] = mid - r, mid + r
            return p0, p1
        (fx if along == 0 else fy).extend([lo, hi])
        (fy if along == 0 else fx).extend([mid - r, mid + r])
        if e.get("capacitance", 0) > 0:
            # A capacitor: the bar with a gap of one minimum cell in the middle, a lumped
            # capacitance across it on the element's axis.
            centre, g = 0.5 * (lo + hi), hmin
            metal.AddBox(*box(lo, centre - g / 2), priority=10)
            metal.AddBox(*box(centre + g / 2, hi), priority=10)
            q0, q1 = [0.0, 0.0, za], [0.0, 0.0, za]
            q0[along], q1[along] = centre - g / 2, centre + g / 2
            q0[across] = q1[across] = mid
            exact[along].extend([centre - g / 2, centre + g / 2])
            exact[across].append(mid)
            exact[2].append(za)
            cap = csx.AddLumpedElement("C_" + str(len(fz)), ny="xyz"[along], caps=False, C=float(e["capacitance"]))
            cap.AddBox(q0, q1, priority=20)
        elif e["port"]:
            centre, g = 0.5 * (lo + hi), hmin
            metal.AddBox(*box(lo, centre - g / 2), priority=10)
            metal.AddBox(*box(centre + g / 2, hi), priority=10)
            # The port is one line of mesh edges across the gap, on the element's axis.
            p0, p1 = [0.0, 0.0, za], [0.0, 0.0, za]
            p0[along], p1[along] = centre - g / 2, centre + g / 2
            p0[across] = p1[across] = mid
            exact[along].extend([centre - g / 2, centre + g / 2])
            exact[across].append(mid)
            exact[2].append(za)
            port = (p0, p1, along)
        else:
            metal.AddBox(*box(lo, hi), priority=10)
    if port is None:
        sys.exit("model.json has no port element")

    lo = [min(fx), min(fy), min(fz)]
    hi = [max(fx), max(fy), max(fz)]
    for ax, fixed in enumerate((fx, fy, fz)):
        grid.SetLines("xyz"[ax], mesh_lines(fixed, exact[ax], lo[ax], hi[ax], hmin, hboard, hmax, margin, SmoothMeshLines))
    p0, p1, along = port
    lumped = fdtd.AddLumpedPort(1, s["port_resistance_ohm"], p0, p1, "xyz"[along], excite=1.0, priority=20)

    # Near-field box ten lines inside the boundary (the absorbing layers take eight).
    lines = [grid.GetLines(ax) for ax in range(3)]
    nf_start = [lines[ax][10] for ax in range(3)]
    nf_stop = [lines[ax][-11] for ax in range(3)]
    nf2ff = fdtd.CreateNF2FFBox(start=nf_start, stop=nf_stop)
    ncells = np.prod([len(l) - 1 for l in lines])
    print("mesh: %d x %d x %d lines, %.3g million cells, smallest cell %.3g mm" % (
        len(lines[0]), len(lines[1]), len(lines[2]), ncells / 1e6, min(np.diff(l).min() for l in lines)), flush=True)

    fdtd.Run(sim_dir, cleanup=True, setup_only=args.setup_only, numThreads=args.threads)
    if args.setup_only:
        return

    # Impedance and inductance from 1 % of f_max up.
    freq = np.linspace(f_max / 100, f_max, 100)
    lumped.CalcPort(sim_dir, freq)
    z = lumped.uf_tot / lumped.if_tot
    out = {"loop": model["loop"], "frequency_hz": freq.tolist(), "impedance_ohm": [[v.real, v.imag] for v in z],
           "inductance_nh": (z.imag / (2 * np.pi * freq) * 1e9).tolist()}
    # Far field: the largest over the sphere, per ampere of port current; and the moment
    # a small loop would need to radiate it.
    ff = np.asarray(s["far_field_hz"], dtype=float)
    lumped.CalcPort(sim_dir, ff)
    theta, phi = np.arange(0.0, 181.0, 3.0), np.arange(0.0, 360.0, 3.0)
    centre = [0.5 * (lo[ax] + hi[ax]) * UNIT for ax in range(3)]
    res = nf2ff.CalcNF2FF(sim_dir, ff, theta, phi, radius=s["distance_m"], center=centre)
    out["far_field"] = []
    est = model.get("estimate")
    for i, f in enumerate(ff):
        e = float(np.max(res.E_norm[i])) / abs(lumped.if_tot[i])
        k = 2 * np.pi * f / C0
        row = {"frequency_hz": float(f), "field_v_per_m_per_a": e,
               "equivalent_moment_m2": e * 4 * np.pi * s["distance_m"] / (ETA0 * k * k)}
        if est and "far_field" in est:
            row["estimate_v_per_m_per_a"] = est["far_field"][i]["v_per_m_per_a"]
        elif est:
            m = float(np.linalg.norm(est["moment_m2"]))
            row["estimate_v_per_m_per_a"] = ETA0 * k * k * m / (4 * np.pi * s["distance_m"])
        out["far_field"].append(row)
    if est:
        out["estimate"] = est
    with open(os.path.join(here, "result.json"), "w") as f:
        json.dump(out, f, indent=1)
    i0 = int(np.argmin(np.abs(freq - f_max / 10)))
    print("inductance %.4g nH at %.4g MHz" % (out["inductance_nh"][i0], freq[i0] * 1e-6)
          + (" (estimate %.4g nH)" % est["inductance_nh"] if est else ""))
    for row in out["far_field"]:
        print("far field  %.4g MHz: %.4g V/m per ampere at %g m" % (row["frequency_hz"] * 1e-6, row["field_v_per_m_per_a"], s["distance_m"])
              + (" (estimate %.4g)" % row["estimate_v_per_m_per_a"] if est else ""))


if __name__ == "__main__":
    main()
