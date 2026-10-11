"""The 3D references with Elmer. DC: a strip between equipotential pads against its closed
form, and a chain of vias against the 2.5D solver. Thermal: a source over a plane against
the finite-volume reference of tools/ref3d.py. Runs when an Elmer solver is found
(KAIROQ_ELMER_SOLVER, or ElmerSolver on the path).
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess

import pytest

import synthetic
from conftest import cli_path

SOLVER = os.environ.get("KAIROQ_ELMER_SOLVER") or shutil.which("ElmerSolver")
elmer = pytest.mark.skipif(not SOLVER or cli_path() is None, reason="no Elmer solver (KAIROQ_ELMER_SOLVER) or kairoq command")
GS = 5.8e7 * 35e-6  # sheet conductance of 35 um copper


def run(tmp_path, board_dir, ports, extra=""):
    (tmp_path / "project.yaml").write_text(
        f"board: {board_dir.name}\ndc:\n  elmer: {{solver: '{SOLVER}'{extra}}}\n  ports:\n" + "".join(f"    - {p}\n" for p in ports))
    r = subprocess.run([cli_path(), "elmer", str(tmp_path / "project.yaml"), "-o", str(tmp_path / "out"), "--run"], capture_output=True,
                       text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    rows = {}
    for line in r.stdout.splitlines():
        m = re.match(r"(\S+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)$", line)
        if m:
            rows[m.group(1)] = [float(m.group(k)) for k in range(2, 6)]  # V Elmer, I Elmer, V native, difference mV
    return rows


@elmer
def test_strip(tmp_path, validation):
    """48 x 1 mm of 35 um copper between equipotential pads: R = L / (sigma t w)."""
    synthetic.strip(tmp_path / "board", length=50.0, width=1.0)
    rows = run(tmp_path, tmp_path / "board", ["{name: A, pad: J1.1, current_a: 2.0}", "{name: B, pad: J2.1, voltage_v: 0.0}"])
    validation.check("strip 48 x 1 mm, Elmer 3D", "resistance", rows["A"][0] / 2.0, 48e-3 / (GS * 1e-3), rel=2e-3, unit="ohm",
                     source="R = L / (sigma t w)")


@elmer
def test_via_chain(tmp_path, validation):
    """A chain of vias joined by tracks on alternate layers: the 3D model against the 2.5D
    solver, which treats a barrel as a resistance between the layers' mid-planes and the
    copper as sheets."""
    synthetic.via_chain(tmp_path / "board")
    rows = run(tmp_path, tmp_path / "board", ["{name: A, at_mm: [5.0, 5.0], layer: 1, diameter_mm: 0.6, current_a: 1.0}",
                                              "{name: B, at_mm: [23.0, 5.0], layer: 1, diameter_mm: 0.6, voltage_v: 0.0}"])
    validation.check("via chain, Elmer 3D against the 2.5D solver", "resistance", rows["A"][2], rows["A"][0], rel=0.03, unit="ohm",
                     source="Elmer (3D)")


@elmer
def test_thermal_source_over_plane(tmp_path, validation):
    """0.5 W in a 2 mm square pad over a plane through 1.5 mm of FR-4, both faces at
    10 W/(m^2 K): the mean rise of the pad against the finite-volume 3D reference
    (docs/validation/thermal.md: 204.3 K, the layered model 200.4 K). At 0.1 mm Elmer is about
    3 % below its converged value (docs/validation/thermal.md, "Real board against Elmer")."""
    synthetic.source_board(tmp_path / "board", 30.0, (14.0, 14.0, 16.0, 16.0), True, 1.5, 0.035)
    (tmp_path / "project.yaml").write_text(
        "board: board\nthermal:\n  ambient_c: 25\n"
        "  top: {convection_w_m2k: 10, emissivity: 0}\n  bottom: {convection_w_m2k: 10, emissivity: 0}\n"
        "  sources:\n    - {name: H, pad: H1.1, power_w: 0.5}\n"
        f"  elmer: {{min_size_mm: 0.1, max_edge_mm: 2.0, solver: '{SOLVER}'}}\n")
    r = subprocess.run([cli_path(), "elmer", str(tmp_path / "project.yaml"), "-o", str(tmp_path / "out"), "--run", "--thermal"],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    m = re.search(r"^H\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)$", r.stdout, re.M)
    elmer_mean, native_mean = float(m.group(1)), float(m.group(3))
    validation.check("source 2 mm over a plane, Elmer 3D at 0.1 mm", "mean rise", elmer_mean - 25.0, 204.3, rel=0.05, unit="K",
                     source="finite-volume 3D reference")
    validation.check("source 2 mm over a plane, layered model against Elmer 3D", "mean rise", native_mean - 25.0, elmer_mean - 25.0,
                     rel=0.03, unit="K", source="Elmer (3D)")
