"""Export of a loop to openEMS, and the full-wave comparison of the loop estimates.

The export tests check the written model against the board. The full-wave tests run the
exported simulation with openEMS and compare its inductance and far field with the
estimates; they need openEMS and run only when KAIROQ_OPENEMS_PYTHON names a Python
interpreter with the openEMS modules (minutes each).
"""

from __future__ import annotations

import json
import math
import os
import subprocess

import pytest
import shapely
import shapely.wkb

import kairoq
import synthetic

LOOP = [{"from": "R1.1", "to": "R1.2", "radius_mm": 0.25}]


def emi(**openems):
    return {"loops": [{"name": "RING", "elements": LOOP}], "loop_settings": {"max_triangles": 3000},
            "openems": {"f_max_hz": 1e9, **openems}}


def polys(model, layer):
    sheet = next(s for s in model["sheets"] if s["layer"] == layer)
    return [shapely.Polygon([(x * 1e3, y * 1e3) for x, y in p]) for p in sheet["polygons"]]


def test_export_ring(tmp_path):
    g = synthetic.ring_loop(tmp_path / "board")
    b = kairoq.open(tmp_path / "board")
    m = kairoq.export_openems(b, emi(), tmp_path / "export")
    assert (tmp_path / "export" / "model.json").exists()
    compile((tmp_path / "export" / "simulate.py").read_text(), "simulate.py", "exec")
    assert json.loads((tmp_path / "export" / "model.json").read_text())["loop"] == "RING"
    # The ring lies entirely inside the window: its sheet is the board's copper, with the
    # outline simplified as in the loop model (by 0.15 of the element size, 0.3 mm here).
    sheet = shapely.union_all(polys(m, 1))
    copper = shapely.union_all([shapely.wkb.loads(v) for v in b.copper_wkb(1).values()])
    assert sheet.hausdorff_distance(copper) < 0.15 * 0.3 + 1e-6
    assert sheet.area == pytest.approx(copper.area, rel=5e-3)
    # The element: from pad centre to pad centre, its axis one radius above the copper.
    (e,) = m["elements"]
    assert e["port"] and e["name"] == "R1.1 - R1.2"
    pads = {p["pin"]: p["position_mm"] for p in b.pads}
    assert [v * 1e3 for v in e["from"]] == pytest.approx(pads["1"], abs=1e-6)
    assert [v * 1e3 for v in e["to"]] == pytest.approx(pads["2"], abs=1e-6)
    assert e["z_axis"] - e["z_sheet"] == pytest.approx(0.25e-3, rel=1e-9)
    thickness = sum(s["z1"] - s["z0"] for s in m["dielectric"])
    copper_and_dielectric = sum(l["thickness_mm"] for l in b.stackup if l["type"] in ("copper", "dielectric")) * 1e-3
    assert thickness == pytest.approx(copper_and_dielectric, rel=1e-9)
    # The estimate is the loop solution of the same settings.
    r = kairoq.analyse_emi(b, emi()).loops["RING"]
    assert m["estimate"]["inductance_nh"] == pytest.approx(r["inductance_nh"], rel=1e-9)
    assert m["simulation"]["far_field_hz"] == [1e8, 2.5e8, 5e8]
    assert m["simulation"]["air_margin"] == pytest.approx(299792458.0 / 1e9 / 4)


def test_holes_are_cut_open(tmp_path):
    """CSXCAD polygons have no holes: a plane with openings is written as pieces that do
    not overlap and together leave the openings free."""
    openings = [(5.0, 5.0, 2.0), (12.0, 4.0, 1.5), (12.5, 6.0, 1.5)]
    synthetic.ring_loop(tmp_path / "board", plane=True, margin=25.0, openings=openings)
    b = kairoq.open(tmp_path / "board")
    m = kairoq.export_openems(b, emi(estimate=False), None)
    pieces = polys(m, 2)
    assert len(pieces) > 1 and all(p.is_valid for p in pieces)
    union = shapely.union_all(pieces)
    assert sum(p.area for p in pieces) == pytest.approx(union.area, rel=1e-9)
    window = shapely.union_all([shapely.Polygon([(x * 1e3, y * 1e3) for x, y in p]) for p in m["window"]])
    assert union.area == pytest.approx(window.area - sum(s * s for _, _, s in openings), rel=1e-6)
    x0 = y0 = 25.0
    for ox, oy, _ in openings:
        assert not union.contains(shapely.Point(x0 + ox, y0 + oy))
    assert "estimate" not in m



OPENEMS_PYTHON = os.environ.get("KAIROQ_OPENEMS_PYTHON")
fullwave = pytest.mark.skipif(not OPENEMS_PYTHON, reason="KAIROQ_OPENEMS_PYTHON names no Python with openEMS")


def simulate(tmp_path, **ring):
    synthetic.ring_loop(tmp_path / "board", **ring)
    b = kairoq.open(tmp_path / "board")
    model = kairoq.export_openems(b, emi(), tmp_path / "export")
    with open(tmp_path / "openems.log", "w") as log:
        subprocess.run([OPENEMS_PYTHON, str(tmp_path / "export" / "simulate.py")], stdout=log, stderr=subprocess.STDOUT,
                       check=True)
    return model, json.loads((tmp_path / "export" / "result.json").read_text())


def compare(case, model, result, validation, l_rel, field_rel):
    """The inductance where it is flat (at a tenth of f_max), and the far field where the
    loop is small against the wavelength (100 and 250 MHz)."""
    f = result["frequency_hz"]
    i = min(range(len(f)), key=lambda k: abs(f[k] - 1e8))
    low = [L for fk, L in zip(f, result["inductance_nh"]) if fk <= 1.1e8]
    assert max(low) - min(low) < 0.02 * low[-1]  # flat below a tenth of f_max, to the truncation ripple
    validation.check(case, "inductance (estimate against full-wave at 100 MHz)", model["estimate"]["inductance_nh"],
                     result["inductance_nh"][i], rel=l_rel, unit="nH", source="openEMS")
    m = math.hypot(*model["estimate"]["moment_m2"])
    for row in result["far_field"][:2]:
        validation.check(case, f"far field per ampere at 3 m, {row['frequency_hz'] / 1e6:g} MHz",
                         row["estimate_v_per_m_per_a"], row["field_v_per_m_per_a"], rel=field_rel, unit="V/m",
                         source="openEMS")
        assert row["equivalent_moment_m2"] == pytest.approx(m, rel=field_rel)


@fullwave
def test_fullwave_free_loop(tmp_path, validation):
    model, result = simulate(tmp_path)
    compare("rectangular loop 20 x 10 mm", model, result, validation, l_rel=0.06, field_rel=0.04)


@fullwave
def test_fullwave_loop_over_plane(tmp_path, validation):
    """The plane cancels most of the moment; what is left, from the plane's finite extent,
    is what radiates."""
    model, result = simulate(tmp_path, plane=True, margin=25.0)
    compare("rectangular loop 20 x 10 mm over a plane", model, result, validation, l_rel=0.06, field_rel=0.04)


SLOT_LOOP = [{"from": "J1.2", "to": "J1.1"}, {"from": "R1.1", "to": "R1.2"}]


@fullwave
def test_fullwave_stitch_across_a_slot(tmp_path, validation):
    """A trace over a slotted plane, with the slot open and with a strap across it beside the
    trace: the strap's effect on the inductance against the full-wave solutions. The
    estimates themselves lie 8 to 11 % above openEMS here: on a long trace close over a
    plane the loop solver converges slowly with the triangle budget (docs/validation/emi.md)."""
    res = {}
    for name, stitches in (("open", []), ("strap", [{"from": "C1.1", "to": "C1.2"}])):
        synthetic.slot_board(tmp_path / name / "board", stitch=bool(stitches))
        b = kairoq.open(tmp_path / name / "board")
        spec = {"loops": [{"name": "X", "elements": SLOT_LOOP, "stitches": stitches}], "openems": {"f_max_hz": 1e9}}
        model = kairoq.export_openems(b, spec, tmp_path / name / "export")
        with open(tmp_path / name / "openems.log", "w") as log:
            subprocess.run([OPENEMS_PYTHON, str(tmp_path / name / "export" / "simulate.py")], stdout=log, stderr=subprocess.STDOUT, check=True)
        result = json.loads((tmp_path / name / "export" / "result.json").read_text())
        f = result["frequency_hz"]
        i = min(range(len(f)), key=lambda k: abs(f[k] - 1e8))
        res[name] = (model["estimate"]["inductance_nh"], result["inductance_nh"][i])
        validation.check(f"slot crossing, {name}", "inductance (estimate against full-wave at 100 MHz)", res[name][0], res[name][1], rel=0.15,
                         unit="nH", source="openEMS")
    validation.check("slot crossing", "inductance removed by the strap", res["open"][0] - res["strap"][0], res["open"][1] - res["strap"][1],
                     rel=0.05, unit="nH", source="openEMS")
