"""Three-dimensional reference for the layered (2.5D) thermal model.

A finite-volume solution of steady conduction on a graded tensor-product grid, written
independently of the library: every copper and dielectric layer is resolved through its
thickness, materials are anisotropic, the top and bottom faces lose heat by convection
through the solder mask, and the sides are adiabatic. The same boards are solved with the
library's 2.5D model and the two are compared (docs/validation/thermal.md).

Usage: python tools/ref3d.py [--quick] [--outer-gap air|mask] [--match TEXT] [--out logs/ref3d.json]
Requires numpy, scipy and pyamg (development tools, not package dependencies).
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys
import tempfile
import time

import numpy as np
import pyamg
import scipy.sparse as sp
import scipy.sparse.linalg as spla

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests" / "python"))
try:
    import kairoq
except ImportError:
    sys.path.insert(0, str(ROOT / "build" / "src" / "python"))
    import kairoq
import synthetic  # noqa: E402

K_CU = 385.0
K_DIEL = (0.8, 0.3)  # in-plane, through
K_MASK = 0.25
K_AIR = 0.0  # gaps of outer copper layers: air, not conducting
OUTER_GAP = "air"  # or "mask": the gaps filled with solder mask to the copper's height
T_MASK = 10e-6


# ---------------------------------------------------------------------------- grid


def graded(lo: float, hi: float, breaks: list[float], h_fine: float, h_max: float, ratio: float = 1.15) -> np.ndarray:
    """Cell edges on [lo, hi] including every break point, fine (h_fine) at the breaks and
    growing geometrically towards h_max between them."""
    pts = sorted({lo, hi, *[b for b in breaks if lo < b < hi]})
    edges = [pts[0]]
    for a, b in zip(pts[:-1], pts[1:]):
        # Grow from both ends towards the middle of [a, b].
        steps, h, s = [], h_fine, 0.0
        while s + h < (b - a) / 2:
            steps.append(h)
            s += h
            h = min(h * ratio, h_max)
        rest = (b - a) - 2 * s
        n_mid = max(1, int(np.ceil(rest / h_max)))
        seg = steps + [rest / n_mid] * n_mid + steps[::-1]
        for d in seg:
            edges.append(edges[-1] + d)
        edges[-1] = b
    return np.array(edges)


# ---------------------------------------------------------------------------- solver


def solve3d(x, y, layers, h, ta, sources, sink=None):
    """Steady temperature on a tensor grid.

    layers: list from the top: {"t": thickness, "n": cells, "copper": [rects]} (copper
    rectangles in a dielectric matrix; a dielectric layer has no "copper"). Between the copper
    of an outer layer is air, which is left out (its conduction is negligible), and the face
    is cooled where the board is, on the copper or on the dielectric below the gap; or, with
    OUTER_GAP = "mask", solder mask up to the height of the copper.
    sources: [{"layer": i, "rect": (x0, y0, x1, y1), "power": W}], spread over the
    copper volume of that layer inside the rectangle.
    sink: (resistance per area of the interface, temperature): the whole bottom face on a
    heat sink instead of convection.
    Returns the temperature (nx, ny, nz) and the z cell layer index.
    """
    nx, ny = len(x) - 1, len(y) - 1
    dx, dy = np.diff(x), np.diff(y)
    xc, yc = 0.5 * (x[1:] + x[:-1]), 0.5 * (y[1:] + y[:-1])
    dz, layer_of = [], []
    for i, l in enumerate(layers):
        dz += [l["t"] / l["n"]] * l["n"]
        layer_of += [i] * l["n"]
    dz = np.array(dz)
    nz = len(dz)
    kxy = np.empty((nx, ny, nz))
    kz = np.empty((nx, ny, nz))
    X, Y = np.meshgrid(xc, yc, indexing="ij")
    for k in range(nz):
        l = layers[layer_of[k]]
        cu = np.zeros((nx, ny), bool)
        for (x0, y0, x1, y1) in l.get("copper", []):
            cu |= (X > x0) & (X < x1) & (Y > y0) & (Y < y1)
        outer = "copper" in l and layer_of[k] in (0, len(layers) - 1)
        fill = ((K_AIR, K_AIR) if OUTER_GAP == "air" else (K_MASK, K_MASK)) if outer else K_DIEL
        kxy[:, :, k] = np.where(cu, K_CU, fill[0])
        kz[:, :, k] = np.where(cu, K_CU, fill[1])
    active = kz > K_AIR
    kxy = np.where(active, kxy, 1.0)  # placeholders; no conductance reaches air cells
    kz = np.where(active, kz, 1.0)
    idx = np.arange(nx * ny * nz).reshape(nx, ny, nz)
    rows, cols, vals = [], [], []

    def couple(a, b, g):
        rows.extend([a, b, a, b])
        cols.extend([a, b, b, a])
        vals.extend([g, g, -g, -g])

    # x faces
    DY, DZ = np.meshgrid(dy, dz, indexing="ij")
    for i in range(nx - 1):
        g = DY * DZ / (dx[i] / (2 * kxy[i]) + dx[i + 1] / (2 * kxy[i + 1])) * (active[i] & active[i + 1])
        couple(idx[i].ravel(), idx[i + 1].ravel(), g.ravel())
    DX, DZ2 = np.meshgrid(dx, dz, indexing="ij")
    for j in range(ny - 1):
        g = DX * DZ2 / (dy[j] / (2 * kxy[:, j]) + dy[j + 1] / (2 * kxy[:, j + 1])) * (active[:, j] & active[:, j + 1])
        couple(idx[:, j].ravel(), idx[:, j + 1].ravel(), g.ravel())
    A = np.outer(dx, dy)
    for k in range(nz - 1):
        g = A / (dz[k] / (2 * kz[:, :, k]) + dz[k + 1] / (2 * kz[:, :, k + 1])) * (active[:, :, k] & active[:, :, k + 1])
        couple(idx[:, :, k].ravel(), idx[:, :, k + 1].ravel(), g.ravel())
    rows, cols, vals = np.concatenate(rows), np.concatenate(cols), np.concatenate(vals)
    b = np.zeros(nx * ny * nz)
    # Convection through the mask on the top and bottom faces (from the cell centre), at
    # the first conducting cell of each column.
    diag_r, diag_v = [], []
    diag_t = []  # temperature each face exchanges with
    for bottom, order in ((False, range(nz)), (True, range(nz - 1, -1, -1))):
        done = np.zeros((nx, ny), bool)
        r_out, t_out = (sink[0], sink[1]) if bottom and sink else (1 / h, ta)
        for k in order:
            face = ~done & active[:, :, k]
            g = np.where(face, A / (dz[k] / (2 * kz[:, :, k]) + T_MASK / K_MASK + r_out), 0.0)
            diag_r.append(idx[:, :, k].ravel())
            diag_v.append(g.ravel())
            diag_t.append(t_out)
            np.add.at(b, idx[:, :, k].ravel(), (g * t_out).ravel())
            done |= face
            if done.all():
                break
    rows = np.concatenate([rows, *diag_r])
    cols = np.concatenate([cols, *diag_r])
    vals = np.concatenate([vals, *diag_v])
    K = sp.coo_matrix((vals, (rows, cols)), shape=(b.size, b.size)).tocsr()
    vol = np.einsum("i,j,k->ijk", dx, dy, dz)
    for s in sources:
        x0, y0, x1, y1 = s["rect"]
        sel = np.zeros((nx, ny, nz), bool)
        ks = [k for k in range(nz) if layer_of[k] == s["layer"]]
        inside = (X > x0) & (X < x1) & (Y > y0) & (Y < y1)
        for k in ks:
            sel[:, :, k] = inside & (kxy[:, :, k] == K_CU)
        w = np.where(sel, vol, 0.0)
        b += (s["power"] * w / w.sum()).ravel()
    a = active.ravel()
    K, b = K[a][:, a], b[a]
    # Symmetric diagonal scaling evens out the conductance jumps between copper,
    # dielectric and the faces before the multigrid hierarchy is built.
    d = 1.0 / np.sqrt(K.diagonal())
    D = sp.diags(d)
    Ks = (D @ K @ D).tocsr()
    ml = pyamg.ruge_stuben_solver(Ks)
    y, info = spla.cg(Ks, d * b, rtol=1e-12, maxiter=2000, M=ml.aspreconditioner(cycle="V"))
    Ta = d * y
    T = np.full(nx * ny * nz, np.nan)
    T[a] = Ta
    # Thin copper cells next to dielectric make the system stiff (condition number about
    # 1e8), so the residual stalls near 1e-7; the solution is accepted on the residual and
    # on the heat balance of the faces.
    residual = np.linalg.norm(b - K @ Ta) / np.linalg.norm(b)
    out = sum(float(np.sum(np.where(g > 0, g * (np.nan_to_num(T[r]) - tr), 0.0))) for r, g, tr in zip(diag_r, diag_v, diag_t))
    power = sum(s_["power"] for s_ in sources)
    if not (residual < 1e-6 and abs(out - power) < 1e-6 * power):
        raise RuntimeError(f"3D reference did not converge: relative residual {residual:.2e}, heat out {out:.9g} W of {power:.9g} W")
    return T.reshape(nx, ny, nz), np.array(layer_of)


# ---------------------------------------------------------------------------- cases


def board_25d(tmp: pathlib.Path, size, pad, plane, dielectric, copper, inner=False):
    """The same board for the library (tests/python/synthetic.py)."""
    synthetic.source_board(tmp, size, pad, plane, dielectric, copper, inner)
    return kairoq.open(tmp)


def ref_solution(size, pad, plane, dielectric, copper, inner, power, h, ta, fine, sink=None):
    """3D solution with cells of `fine` mm next to the source: in the plane around its
    edges and through the 0.1 mm of dielectric next to it. Returns the mean and maximum
    rise of the source copper, the cell count and the time."""
    mm = 1e-3
    x0, y0, x1, y1 = pad
    xs = graded(0, size * mm, [x0 * mm, x1 * mm], fine * mm, 1.0 * mm)
    ys = graded(0, size * mm, [y0 * mm, y1 * mm], fine * mm, 1.0 * mm)
    cu_t = copper * mm
    n_cu = max(2, int(np.ceil(copper / 0.0175)))
    src = [(x0 * mm, y0 * mm, x1 * mm, y1 * mm)]
    bottom = {"t": cu_t, "n": n_cu, "copper": [(0, 0, size * mm, size * mm)] if plane else []}
    d = dielectric / 2 if inner else dielectric
    near_t = min(0.1, d / 2)
    near = {"t": near_t * mm, "n": max(2, int(round(near_t / fine)))}
    rest = {"t": (d - near_t) * mm, "n": max(4, int(np.ceil((d - near_t) / 0.05)))}
    n_src = max(2, int(round(n_cu * 0.025 / fine)))
    source = {"t": cu_t, "n": n_src, "copper": src}
    if inner:
        layers = [{"t": cu_t, "n": n_cu, "copper": []}, rest, near, source, dict(near), dict(rest), bottom]
        src_layer = 3
    else:
        layers = [source, near, rest, bottom]
        src_layer = 0
    t0 = time.perf_counter()
    T, layer_of = solve3d(xs, ys, layers, h, ta, [{"layer": src_layer, "rect": src[0], "power": power}],
                          None if sink is None else (sink, ta))
    X, Y = np.meshgrid(0.5 * (xs[1:] + xs[:-1]), 0.5 * (ys[1:] + ys[:-1]), indexing="ij")
    inpad = (X > x0 * mm) & (X < x1 * mm) & (Y > y0 * mm) & (Y < y1 * mm)
    cells = T[:, :, layer_of == src_layer]
    vol = np.einsum("i,j->ij", np.diff(xs), np.diff(ys))
    mean = float((np.nanmean(cells[inpad], axis=1) * vol[inpad]).sum() / vol[inpad].sum())
    return mean - ta, float(np.nanmax(cells)) - ta, int(T.size), time.perf_counter() - t0


def run_case(name, size, pad, plane, dielectric, copper, inner=False, sink=None, power=0.5, h=10.0, ta=25.0):
    """sink: the whole bottom face on a heat sink at ambient through an interface of this
    resistance per area (m^2 K/W)."""
    x0, y0, x1, y1 = pad
    # Two resolutions; the edges of the source make the convergence first order, so the
    # reference is extrapolated as 2 T(h/2) - T(h), and the difference between the two
    # resolutions is reported as its uncertainty.
    m1, x1_, c1, s1 = ref_solution(size, pad, plane, dielectric, copper, inner, power, h, ta, 0.025, sink)
    m2, x2_, c2, s2 = ref_solution(size, pad, plane, dielectric, copper, inner, power, h, ta, 0.0125, sink)
    ref_mean, ref_max = 2 * m2 - m1, 2 * x2_ - x1_
    out = {
        "case": name, "size_mm": size, "pad_mm": [x1 - x0, y1 - y0], "plane": plane, "inner": inner, "sink": sink,
        "dielectric_mm": dielectric,
        "copper_mm": copper, "power_w": power, "rise_3d_max": ref_max, "rise_3d_mean": ref_mean,
        "rise_3d_mean_25um": m1, "rise_3d_mean_12um": m2, "uncertainty_3d": abs(m2 - m1) / ref_mean,
        "cells_3d": c2, "seconds_3d": round(s1 + s2, 1),
    }
    ref_mean += ta
    ref_max += ta
    # The library at its default dielectric slicing, and with one sheet per copper layer
    # only (each dielectric a single slice).
    with tempfile.TemporaryDirectory() as d:
        b = board_25d(pathlib.Path(d), size, pad, plane, dielectric, copper, inner)
        for key, extra in (("25d", {}), ("25d_single", {"dielectric_slice_mm": 100.0, "dielectric_first_slice_mm": 100.0})):
            t0 = time.perf_counter()
            r = kairoq.solve_thermal(b, {"ambient_c": ta, "target_error": 0.03, **extra,
                                           "top": {"convection_w_m2k": h, "emissivity": 0},
                                           "bottom": {"convection_w_m2k": h, "emissivity": 0},
                                           "sources": [{"name": "H", "pad": "H1.1", "power_w": power}],
                                           "heatsinks": [] if sink is None else
                                           [{"name": "S", "tim": {"resistance_k_mm2_per_w": sink * 1e6}, "temperature_c": ta}]})
            lay_max = float(r.layer(2 if inner else 1)["temperature_c"].max())
            lay_mean = r.sources["H"]["mean_c"]
            out.update({
                f"rise_{key}_max": lay_max - ta, f"rise_{key}_mean": lay_mean - ta,
                f"error_{key}_max": (lay_max - ref_max) / (ref_max - ta), f"error_{key}_mean": (lay_mean - ref_mean) / (ref_mean - ta),
                f"triangles_{key}": r.triangles, f"seconds_{key}": round(time.perf_counter() - t0, 1),
            })
    return out


def cases(quick: bool):
    c = 15.0  # board centre, 30 mm board
    out = []
    for s in ([1.0, 5.0] if quick else [0.5, 1.0, 2.0, 5.0]):
        out.append((f"square source {s} mm over a plane, 1.5 mm", 30.0, (c - s / 2, c - s / 2, c + s / 2, c + s / 2), True, 1.5, 0.035))
    for s in ([1.0] if quick else [0.5, 1.0, 2.0]):
        out.append((f"square source {s} mm over a plane, 0.2 mm", 30.0, (c - s / 2, c - s / 2, c + s / 2, c + s / 2), True, 0.2, 0.035))
    for w in ([1.0] if quick else [0.25, 1.0, 5.0]):
        out.append((f"trace {w} x 20 mm, no plane", 30.0, (5.0, c - w / 2, 25.0, c + w / 2), False, 1.5, 0.035))
    if not quick:
        out.append(("square source 2 mm, 10 oz copper, over a plane", 30.0, (c - 1, c - 1, c + 1, c + 1), True, 1.5, 0.35))
        out.append(("trace 1.0 x 20 mm, 3 oz, inner layer, no plane", 30.0, (5.0, c - 0.5, 25.0, c + 0.5), False, 1.5, 0.105, True))
        out.append(("trace 1.0 x 20 mm, 10 oz, inner layer, no plane", 30.0, (5.0, c - 0.5, 25.0, c + 0.5), False, 1.5, 0.35, True))
    # The bottom face on a heat sink through 0.5 mm of a 3 W/(m K) interface.
    tim = 0.5e-3 / 3.0
    out.append(("square source 2 mm on a sink, 1.5 mm, no plane", 30.0, (c - 1, c - 1, c + 1, c + 1), False, 1.5, 0.035, False, tim))
    out.append(("square source 2 mm on a sink, 0.2 mm, plane", 30.0, (c - 1, c - 1, c + 1, c + 1), True, 0.2, 0.035, False, tim))
    if not quick:
        out.append(("square source 0.5 mm on a sink, 0.2 mm, plane", 30.0, (c - .25, c - .25, c + .25, c + .25), True, 0.2, 0.035,
                    False, tim))
    return out


def slab_check():
    """The 3D solver itself on the one-dimensional slab (uniform heat, full copper)."""
    mm = 1e-3
    xs = np.linspace(0, 10 * mm, 11)
    layers = [{"t": 35e-6, "n": 2, "copper": [(0, 0, 10 * mm, 10 * mm)]}, {"t": 1.5 * mm, "n": 10},
              {"t": 35e-6, "n": 2, "copper": [(0, 0, 10 * mm, 10 * mm)]}]
    T, layer_of = solve3d(xs, xs, layers, 10.0, 25.0, [{"layer": 0, "rect": (0, 0, 10 * mm, 10 * mm), "power": 0.5}])
    q = 0.5 / (10 * mm) ** 2
    # Series resistances per area from the top copper centre (approximately: copper is thin).
    r_top = 35e-6 / 2 / K_CU + T_MASK / K_MASK + 1 / 10.0
    r_mid = 35e-6 / 2 / K_CU + 1.5 * mm / K_DIEL[1] + 35e-6 / 2 / K_CU
    r_bot = 35e-6 / 2 / K_CU + T_MASK / K_MASK + 1 / 10.0
    # q splits between the top path and the path through the board.
    t_top = 25.0 + q / (1 / r_top + 1 / (r_mid + r_bot))
    return float(T[:, :, layer_of == 0].mean()), t_top


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--outer-gap", choices=["air", "mask"], default="air")
    ap.add_argument("--match", default="", help="run only the cases whose name contains this text")
    ap.add_argument("--out", default=str(ROOT / "logs" / "ref3d.json"))
    args = ap.parse_args()
    global OUTER_GAP
    OUTER_GAP = args.outer_gap
    t3, t_exact = slab_check()
    print(json.dumps({"slab_check_3d": t3, "slab_exact": t_exact, "relative": (t3 - t_exact) / (t_exact - 25.0)}), flush=True)
    results = []
    for c in cases(args.quick):
        if args.match not in c[0]:
            continue
        res = run_case(*c)
        results.append(res)
        print(json.dumps(res), flush=True)
    pathlib.Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    pathlib.Path(args.out).write_text(json.dumps(results, indent=1))


if __name__ == "__main__":
    main()
