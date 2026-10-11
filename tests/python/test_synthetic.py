"""Imports of synthetic boards against their analytical geometry."""

from __future__ import annotations

import shutil

import pytest
import shapely
import shapely.wkb

import kairoq
import synthetic


def area(wkb: bytes) -> float:
    return shapely.wkb.loads(wkb).area


def test_strip(tmp_path, validation):
    ref = synthetic.strip(tmp_path)
    b = kairoq.open(tmp_path)
    assert b.net_source == "x2"
    assert b.nets[1:] == ref["nets"]
    validation.check("strip", "L1 copper area", b.copper_area_mm2(1), ref["copper_area_l1"], rel=1e-6, unit="mm2")
    validation.check("strip", "outline area", area(b.outline_wkb()), ref["outline_area"], rel=1e-9, unit="mm2")
    # The core's area and an independent WKB decoding agree.
    net = b.nets.index("SIG")
    assert area(b.copper_wkb(1)[net]) == pytest.approx(b.copper_area_mm2(1), rel=1e-12)
    assert {(p["component"], p["pin"]) for p in b.pads} == {("J1", "1"), ("J2", "1")}
    assert all(b.nets[p["net"]] == "SIG" for p in b.pads)
    assert b.copper_area_mm2(2) == 0.0


@pytest.mark.parametrize("width", [0.1, 0.25, 3.0])
def test_strip_widths(tmp_path, validation, width):
    ref = synthetic.strip(tmp_path, length=20.0, width=width)
    b = kairoq.open(tmp_path)
    # Pads 2 x 2 mm at both ends; track wider than a pad extends past it.
    pad = 4.0
    if width <= 2.0:
        expected = 2 * pad + (20.0 - 2.0) * width
    else:
        r = width / 2
        import math
        # Stadium of the track, plus the pad corners outside it (none: pad lies within the track).
        expected = 20.0 * width + math.pi * r * r
    validation.check(f"strip w={width}", "L1 copper area", b.copper_area_mm2(1), expected, rel=2e-6, unit="mm2")


def test_via_chain(tmp_path, validation):
    ref = synthetic.via_chain(tmp_path)
    b = kairoq.open(tmp_path)
    assert b.nets[1:] == ref["nets"]
    vias = [d for d in b.drills if d["function"] == "via"]
    validation.check("via chain", "via count", len(vias), ref["vias"], abs_=0)
    assert all(d["plated"] and d["from_layer"] == 1 and d["to_layer"] == 2 for d in vias)
    assert all(b.nets[d["net"]] == "CHAIN" for d in vias)
    assert all(d["diameter_mm"] == pytest.approx(ref["via_drill"]) for d in vias)
    # The netlist agrees with the X2 nets: no warnings about it.
    assert not [w for w in b.warnings if "netlist" in w]


def test_via_chain_without_x2_uses_connectivity(tmp_path, validation):
    synthetic.via_chain(tmp_path)
    # Strip all X2 attributes and the job file: roles from file extensions, nets from
    # geometric connectivity named by the IPC-D-356 netlist.
    for f in tmp_path.iterdir():
        if f.suffix == ".gbrjob":
            f.unlink()
        elif f.suffix in (".gtl", ".gbl", ".gm1"):
            f.write_text("\n".join(line for line in f.read_text().splitlines()
                                   if not line.startswith(("%TO", "%TA", "%TD", "%TF"))))
    b = kairoq.open(tmp_path)
    assert b.net_source == "ipc356"
    assert b.data["stackup"]["source"] == "default"
    named = [n for n in b.nets if n]
    assert named == ["CHAIN"], named
    assert all(b.nets[d["net"]] == "CHAIN" for d in b.drills)


def test_plane_with_clear_polarity(tmp_path, validation):
    ref = synthetic.plane_with_clearance(tmp_path)
    b = kairoq.open(tmp_path)
    assert len(b.copper_layers) == 4
    assert [c["side"] for c in b.copper_layers] == ["top", "inner", "inner", "bottom"]
    validation.check("plane", "L2 area (plane minus clear square)", b.copper_area_mm2(2), ref["l2_area"], rel=1e-9,
                     unit="mm2")
    validation.check("plane", "L3 area (6 mm disc)", b.copper_area_mm2(3), ref["l3_area"], rel=1e-6, unit="mm2")
    thickness = sum(l["thickness_mm"] for l in b.stackup)
    validation.check("plane", "stackup thickness", thickness, ref["thickness"], rel=1e-9, unit="mm")
    # The through via joins TP1/TP2 to the plane: its net is GND.
    assert b.nets[b.drills[0]["net"]] == "GND"


def test_zip_input(tmp_path):
    synthetic.strip(tmp_path / "files")
    archive = shutil.make_archive(str(tmp_path / "fab"), "zip", tmp_path / "files")
    a = kairoq.open(tmp_path / "files")
    b = kairoq.open(archive)
    assert b.copper_area_mm2(1) == pytest.approx(a.copper_area_mm2(1), rel=1e-12)
    assert b.nets == a.nets


def test_arc_tolerance_option(tmp_path):
    synthetic.plane_with_clearance(tmp_path)
    coarse = kairoq.open(tmp_path, arc_tolerance_um=20.0)
    fine = kairoq.open(tmp_path, arc_tolerance_um=0.2)
    # The disc keeps its area at any tolerance (equal-area polygons) ...
    assert coarse.copper_area_mm2(3) == pytest.approx(fine.copper_area_mm2(3), rel=1e-5)
    # ... with fewer vertices when coarse.
    n = lambda b: len(shapely.get_coordinates(shapely.wkb.loads(b.copper_wkb(3)[b.nets.index("PWR")])))
    assert n(coarse) < n(fine)


def test_stackup_yaml_overrides_job_file(tmp_path):
    synthetic.strip(tmp_path)
    y = tmp_path.parent / "stack.yaml"
    y.write_text("layers:\n  - {name: top, type: copper, thickness_mm: 0.07}\n"
                 "  - {name: core, type: dielectric, thickness_mm: 0.8, material: FR-4}\n"
                 "  - {name: bot, type: copper, thickness_mm: 0.07}\n")
    b = kairoq.open(tmp_path, stackup=y)
    assert b.data["stackup"]["source"] == "yaml"
    assert [l["thickness_mm"] for l in b.stackup] == pytest.approx([0.07, 0.8, 0.07])


def test_crop(tmp_path):
    """The left half of a board with a plane and a 2 mm pad on its middle line."""
    synthetic.source_board(tmp_path, 30.0, (14.0, 14.0, 16.0, 16.0), True, 1.5, 0.035)
    whole = kairoq.open(tmp_path)
    half = kairoq.open(tmp_path, crop_mm=(0.0, 0.0, 15.0, 30.0))
    assert shapely.wkb.loads(half.outline_wkb()).area == pytest.approx(450.0, rel=1e-9)
    assert half.copper_area_mm2(2) == pytest.approx(whole.copper_area_mm2(2) / 2, rel=1e-9)
    assert half.copper_area_mm2(1) == pytest.approx(2.0, rel=1e-6)
    assert [p["component"] for p in half.pads] == ["H1"]
    assert half.warnings[-1].startswith("cropped to 0.000, 0.000 to 15.000, 30.000 mm")
    with pytest.raises(Exception, match="does not overlap"):
        kairoq.open(tmp_path, crop_mm=(40.0, 0.0, 50.0, 30.0))
