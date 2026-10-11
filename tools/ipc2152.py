"""Trace temperature rise against IPC-2152 data.

The IPC-2152 charts were measured on the test fixture of IPC-TM-650 2.5.4.1A: a 12 in
(304.8 mm) trace on a polyimide board without copper planes, in still air; the temperature
rise is that of the 6 in section in the middle, from its change in resistance. The charts
are used here through the curve fits of Brooks and Adam ("Trace Currents and Temperatures
Revisited", 2015, Table 3-1 and Appendix 6), with W and Th in mil and I in A:

    external, all weights:  dT = 215.3 I^2 W^-1.15 Th^-1.0
    internal, 1 oz:         dT = 200 I^1.9 W^-1.10 Th^-1.52
    internal, 2 oz:         dT = 300.3 I^2 W^-1.15 Th^-1.52
    internal, 3 oz:         dT = 300 I^1.9 W^-1.15 Th^-1.52 (225 at 10 mil)

The same paper gives simulations (TRM) of a 6 in, 1 oz, 200 mil trace at 15 A with planes
added, which are compared here as ratios to the case without a plane.

The baseline board is 0.070 in thick polyimide (Jouppi, "The Value of IPC-2152", IPC APEX
2012). Its conductivity is not published outside the standard; the traces are compared at
0.2 and 0.6 W/(m K). The lateral temperature profile around a 10 mil, 1 oz external trace
at a 10 C rise (Jouppi, Figure 10: 35, 29, 26.4 and 25.5 C at 0, 0.33, 0.6 and 0.9 in from
the trace, 25 C ambient) is compared at 0.2, 0.4 and 0.6 W/(m K); its decay length
constrains the board's in-plane conductance. The surface heat transfer coefficient is
10 W/(m^2 K), radiation included, the value the paper calibrated at small temperature
rises (it reports 11 to 14 at large ones), and the resistivity 1.7 uOhm cm at 20 C.

Usage: python tools/ipc2152.py [--quick] [--out logs/ipc2152.jsonl]
Results are appended per case; cases already in the output file are skipped.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys
import tempfile
import time

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests" / "python"))
try:
    import kairoq
except ImportError:
    sys.path.insert(0, str(ROOT / "build" / "src" / "python"))
    import kairoq
import synthetic  # noqa: E402

MIL = 0.0254  # mm
OZ = 0.035  # mm of copper per ounce
BOARD = (365.76, 142.24)  # mm, the fixture board
THICKNESS = 0.070 * 25.4  # mm
TRACE = 304.8  # mm
SECTION = 152.4  # mm, sensed in the middle
H = 10.0
RHO = 1.7e-8
K = (0.2, 0.4, 0.6)  # profile
K_SWEEP = (0.2, 0.6)  # trace and plane cases: the ends of the range


def fit(location: str, oz: int, width_mil: float, current: float) -> float:
    th = oz * OZ / MIL
    if location == "external":
        return 215.3 * current ** 2 * width_mil ** -1.15 * th ** -1.0
    if oz == 1:
        return 200.0 * current ** 1.9 * width_mil ** -1.10 * th ** -1.52
    if oz == 2:
        return 300.3 * current ** 2 * width_mil ** -1.15 * th ** -1.52
    c = 225.0 if width_mil <= 10 else 300.0
    return c * current ** 1.9 * width_mil ** -1.15 * th ** -1.52


def current_for(location, oz, width_mil, rise):
    lo, hi = 1e-3, 1e3
    for _ in range(200):
        mid = (lo * hi) ** 0.5
        lo, hi = (mid, hi) if fit(location, oz, width_mil, mid) < rise else (lo, mid)
    return (lo * hi) ** 0.5


def board(out: pathlib.Path, location: str, oz: int, width: float, length: float, pad: float, plane: str | None = None):
    """Trace along the middle of the fixture board. location: external (layer 1 of two) or
    internal (layer 2 of three, in the middle of the thickness). plane: None, "bottom"
    (layer 2 of two) or "near" (layer 2 of three, 10 mil below the trace)."""
    w, h = BOARD
    y = h / 2
    x0, x1 = (w - length) / 2, (w + length) / 2
    t = oz * OZ
    if location == "external" and plane != "near":
        names, coppers, diel = ["Top", "Bot"], [t, t], [THICKNESS - 2 * t]
        trace_layer = 1
    elif location == "external":
        names, coppers, diel = ["Top", "Inr", "Bot"], [t, t, t], [10 * MIL, THICKNESS - 3 * t - 10 * MIL]
        trace_layer = 1
    else:
        names, coppers, diel = ["Top", "Inr", "Bot"], [t, t, t], [(THICKNESS - 3 * t) / 2] * 2
        trace_layer = 2
    gs = [synthetic.Gerber(f"Copper,L{i + 1},{n}") for i, n in enumerate(names)]
    g = gs[trace_layer - 1]
    g.track([(x0, y), (x1, y)], width, "SIG")
    g.flash(f"R,{pad:.6f}X{pad:.6f}", x0, y, "SIG", ("J1", "1"), "SMDPad")
    g.flash(f"R,{pad:.6f}X{pad:.6f}", x1, y, "SIG", ("J2", "1"), "SMDPad")
    if plane:
        m = 2.0
        gs[1].region([(m, m), (w - m, m), (w - m, h - m), (m, h - m)], "GND")
    prof = synthetic.Gerber("Profile,NP")
    prof.outline([(0, 0), (w, 0), (w, h), (0, h)])
    files = {f"t-L{i + 1}.gbr": gg for i, gg in enumerate(gs)}
    files["t-Edge.gm1"] = prof
    synthetic._write(out, "t", files, {}, coppers, diel)
    return kairoq.open(out), trace_layer


def solve(b, k, current, h=H):
    names = [l["name"] for l in b.stackup if l["type"] == "dielectric"]
    th = {"ambient_c": 25.0, "target_error": 0.1, "coupling_tolerance_k": 0.1, "layers": {n: {"k_xy": k, "k_z": k} for n in names},
          "top": {"convection_w_m2k": h, "emissivity": 0.0}, "bottom": {"convection_w_m2k": h, "emissivity": 0.0}}
    dc = {"copper_conductivity_s_per_m": 1.0 / RHO,
          "ports": [{"name": "A", "pad": "J1.1", "current_a": current}, {"name": "B", "pad": "J2.1", "current_a": -current}]}
    return kairoq.solve_thermal(b, th, dc=dc)


def section_rise(r, layer, x_mid, half, ta=25.0):
    """Mean temperature rise of the copper over |x - x_mid| < half (mm), area weighted."""
    f = r.layer(layer)
    p, tri, t = f["points_mm"], f["triangles"], f["temperature_c"]
    cu = f["copper"].astype(bool)
    c = p[tri].mean(axis=1)
    sel = cu & (np.abs(c[:, 0] - x_mid) < half)
    u, v = p[tri[:, 1]] - p[tri[:, 0]], p[tri[:, 2]] - p[tri[:, 0]]
    a = 0.5 * np.abs(u[:, 0] * v[:, 1] - u[:, 1] * v[:, 0])
    tm = t[tri].mean(axis=1)
    return float((tm * a)[sel].sum() / a[sel].sum()) - ta


def ipc_cases(quick):
    widths = [10, 50, 200]
    out = []
    for location in ("external", "internal"):
        for oz in ([2] if quick else [2, 1, 3]):
            for wm in widths:
                for rise in ([30] if quick else [20, 60]):
                    for k in K_SWEEP:
                        out.append({"kind": "ipc", "location": location, "oz": oz, "width_mil": wm, "fit_rise": rise, "k": k})
    return out


def plane_cases():
    # Brooks and Adam, section 5: 6 in, 1 oz, 200 mil external trace, 15 A, HTC 11,
    # polyimide; midpoint rise 94.7 C without a plane (large pads at the ends), 47.9 C with
    # a plane on the bottom, 33.7 C with a plane 10 mil below the trace.
    return [{"kind": "plane", "plane": p, "trace_rise_trm": v, "k": k}
            for p, v in ((None, 94.7), ("bottom", 47.9), ("near", 33.7)) for k in K_SWEEP]


def profile_cases():
    return [{"kind": "profile", "k": k} for k in K]


def key(c):
    return json.dumps(c, sort_keys=True)


def run(c):
    t0 = time.perf_counter()
    with tempfile.TemporaryDirectory() as d:
        if c["kind"] == "ipc":
            width = c["width_mil"] * MIL
            current = current_for(c["location"], c["oz"], c["width_mil"], c["fit_rise"])
            b, layer = board(pathlib.Path(d), c["location"], c["oz"], width, TRACE, max(2.54, 2 * width))
            r = solve(b, c["k"], current)
            rise = section_rise(r, layer, BOARD[0] / 2, SECTION / 2)
            out = {"current_a": current, "rise": rise, "error": rise / c["fit_rise"] - 1}
        elif c["kind"] == "profile":
            current = current_for("external", 1, 10, 10.0)
            b, layer = board(pathlib.Path(d), "external", 1, 10 * MIL, TRACE, 2.54)
            r = solve(b, c["k"], current)
            f = r.layer(layer)
            p, t = f["points_mm"], f["temperature_c"]
            near = np.abs(p[:, 0] - BOARD[0] / 2) < 5.0
            out = {"current_a": current, "rise": section_rise(r, layer, BOARD[0] / 2, SECTION / 2), "profile": []}
            for dist_in, ref in ((0.0, 10.0), (0.33, 4.0), (0.6, 1.4), (0.9, 0.5)):
                dy = np.abs(np.abs(p[:, 1] - BOARD[1] / 2) - dist_in * 25.4)
                sel = near & (dy < 0.5)
                out["profile"].append({"distance_in": dist_in, "rise": float(t[sel].mean()) - 25.0, "rise_ref": ref})
        else:
            # The pads that stand for the connection to a plane: 1 in squares (not given).
            b, layer = board(pathlib.Path(d), "external", 1, 200 * MIL, 152.4, 25.4, c["plane"])
            r = solve(b, c["k"], 15.0, h=11.0)
            f = r.layer(layer)
            mid = np.abs(f["points_mm"][:, 0] - BOARD[0] / 2) < 1.0
            on = mid & (np.abs(f["points_mm"][:, 1] - BOARD[1] / 2) < 100 * MIL)
            out = {"current_a": 15.0, "rise": float(f["temperature_c"][on].max()) - 25.0}
    out.update({"triangles": r.triangles, "estimated_error": r.estimated_error, "coupling_iterations": r.coupling_iterations,
                "seconds": round(time.perf_counter() - t0, 1), "warnings": list(r.warnings)})
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--out", default=str(ROOT / "logs" / "ipc2152.jsonl"))
    args = ap.parse_args()
    path = pathlib.Path(args.out)
    path.parent.mkdir(parents=True, exist_ok=True)
    done = set()
    if path.exists():
        for line in path.read_text().splitlines():
            if line.strip():
                done.add(key(json.loads(line)["case"]))
    for c in profile_cases() + ipc_cases(args.quick) + plane_cases():
        if key(c) in done:
            continue
        res = {"case": c, **run(c)}
        with path.open("a") as f:
            f.write(json.dumps(res) + "\n")
        print(json.dumps(res), flush=True)


if __name__ == "__main__":
    main()
