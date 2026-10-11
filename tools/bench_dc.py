"""DC solver timing on the performance reference boards.

Usage: python tools/bench_dc.py [--target E] [--max-triangles N] [--case NAME]...

Cases: power100 (synthetic 100 x 100 mm, 4 layers, 2304 vias, 0.1 mm tracks and
clearances), cm5 (KiCad demo cm5_minima, 6 layers, 3 A supply and return). One JSON line
per case is printed; run under a memory cap (docs/solvers/dc.md, "Resources").
"""

from __future__ import annotations

import argparse
import json
import pathlib
import resource
import sys
import tempfile
import time

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests" / "python"))
try:
    import kairoq
except ImportError:
    sys.path.insert(0, str(ROOT / "build" / "src" / "python"))
    import kairoq
import synthetic  # noqa: E402

PORTS = {
    "power100": [
        {"name": "VIN", "pad": "J1.1", "voltage_v": 5.0}, {"name": "GND", "pad": "J1.2", "voltage_v": 0.0},
        {"name": "LOAD", "pad": "U1.1", "current_a": -10.0}, {"name": "RET", "pad": "U1.2", "current_a": 10.0},
    ],
    "cm5": [
        {"name": "VIN", "component": "J101", "net": "+5V", "voltage_v": 5.0},
        {"name": "LOAD", "component": "Module301", "net": "+5V", "current_a": -3.0},
        {"name": "GND_IN", "component": "J101", "net": "GND", "voltage_v": 0.0},
        {"name": "GND_RET", "component": "Module302", "net": "GND", "current_a": 3.0},
    ],
}


def board(case: str) -> kairoq.Board:
    if case == "power100":
        d = pathlib.Path(tempfile.mkdtemp())
        synthetic.power_board(d)
        return kairoq.open(d)
    return kairoq.open(ROOT / "testdata" / "kicad" / "cm5_minima" / "gerber")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", type=float, default=None)
    ap.add_argument("--max-triangles", type=int, default=None)
    ap.add_argument("--case", action="append", choices=sorted(PORTS), default=None)
    args = ap.parse_args()
    for case in args.case or sorted(PORTS):
        t0 = time.perf_counter()
        b = board(case)
        t1 = time.perf_counter()
        settings = {"ports": PORTS[case]}
        if args.target is not None:
            settings["target_error"] = args.target
        if args.max_triangles is not None:
            settings["max_triangles"] = args.max_triangles
        r = kairoq.solve_dc(b, settings)
        t2 = time.perf_counter()
        print(json.dumps({
            "case": case, "import_s": round(t1 - t0, 2), "solve_total_s": round(t2 - t1, 2),
            "mesh_s": round(r.seconds_mesh, 2), "assembly_s": round(r.seconds_assemble, 2), "solver_s": round(r.seconds_solve, 2),
            "passes": [(p["triangles"], round(p["estimated_error"], 4)) for p in r.passes],
            "unknowns": r.unknowns, "iterations": r.iterations, "estimated_error": round(r.estimated_error, 4),
            "ports": {k: round(v["voltage_v"], 7) for k, v in r.ports.items()}, "power_w": r.power_ports_w,
            "peak_rss_gb": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6, 2),
            "warnings": r.warnings,
        }), flush=True)


if __name__ == "__main__":
    main()
