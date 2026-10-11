"""Command-line tool and display meshes."""

from __future__ import annotations

import json
import subprocess

import numpy as np
import pytest
import shapely.wkb

import kairoq
import synthetic
from conftest import TESTDATA, cli_path

CLI = cli_path()


@pytest.mark.skipif(CLI is None, reason="kairoq command not built")
def test_cli_import_and_info(tmp_path):
    out = tmp_path / "model"
    r = subprocess.run([CLI, "import", str(TESTDATA / "kicad" / "ecc83" / "gerber"), "-o", str(out)],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert (out / "board.json").exists() and (out / "copper" / "L1.wkb").exists()
    r = subprocess.run([CLI, "info", str(out), "--json"], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    data = json.loads(r.stdout)
    assert data["format"] == "kairoq-board" and data["version"] == kairoq.MODEL_VERSION
    assert len(data["copper_layers"]) == 2


@pytest.mark.skipif(CLI is None, reason="kairoq command not built")
def test_cli_errors():
    r = subprocess.run([CLI, "info", "/nonexistent/board"], capture_output=True, text=True)
    assert r.returncode == 1 and "kairoq:" in r.stderr
    r = subprocess.run([CLI, "frobnicate", "x"], capture_output=True, text=True)
    assert r.returncode != 0


def top_area(mesh, z):
    p, t = mesh["positions"], mesh["triangles"]
    tri = p[t]
    on_top = np.all(np.isclose(tri[:, :, 2], z, atol=1e-6), axis=1)
    a, b, c = tri[on_top, 0, :2], tri[on_top, 1, :2], tri[on_top, 2, :2]
    u, v = b - a, c - a
    return 0.5 * np.abs(u[:, 0] * v[:, 1] - u[:, 1] * v[:, 0]).sum()


def test_meshes_match_copper(tmp_path):
    synthetic.plane_with_clearance(tmp_path)
    b = kairoq.open(tmp_path)
    keys = b.mesh_keys()
    assert keys[:4] == ["copper/1", "copper/2", "copper/3", "copper/4"]
    assert "drills" in keys
    z = 0.0
    tops = []
    for layer in reversed(b.stackup):
        tops.append((layer, z + layer["thickness_mm"]))
        z += layer["thickness_mm"]
    ztop = {l.get("copper_index"): zt for l, zt in tops if l["type"] == "copper"}
    for L in (2, 3):
        m = b.mesh(f"copper/{L}")
        assert m["positions"].dtype == np.float32 and m["triangles"].shape[1] == 3
        assert top_area(m, ztop[L]) == pytest.approx(b.copper_area_mm2(L), rel=1e-5)
        # Every vertex carries the net id of its copper.
        assert set(np.unique(m["ids"])) <= {b.nets.index(n) for n in b.nets if n}
    # The dielectric mesh has the outline area less the drilled hole.
    for i, l in enumerate(b.stackup):
        if l["type"] == "dielectric":
            m = b.mesh(f"dielectric/{i}")
            outline = shapely.wkb.loads(b.outline_wkb()).area
            zt = [zt for ll, zt in tops if ll is l][0]
            assert top_area(m, zt) == pytest.approx(outline - np.pi * 0.15 ** 2, rel=1e-4)
            break
