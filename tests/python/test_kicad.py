"""Cross-check of the two importers on KiCad demonstration boards.

The same board exported by KiCad as Gerber X2 (with drill files, job file, IPC-D-356
and placement) and as IPC-2581 must give the same canonical model, apart from known
differences of geometric discretisation. Both exports use the zone fill stored in the
board (tools/export_kicad_demos.sh).
"""

from __future__ import annotations

import collections
import math
import shutil
import time

import pytest
import shapely
import shapely.wkb

import kairoq
from conftest import ROOT, TESTDATA

BOARDS = ["ecc83", "pic_programmer", "cm5_minima"]
_cache: dict[str, tuple] = {}


def ipc2581_file(name):
    d = TESTDATA / "kicad" / name
    return d / f"{name}.xml" if (d / f"{name}.xml").exists() else d / f"{name}.ipc2581.zip"


def boards(name):
    if name not in _cache:
        t0 = time.perf_counter()
        g = kairoq.open(TESTDATA / "kicad" / name / "gerber")
        t1 = time.perf_counter()
        x = kairoq.open(ipc2581_file(name))
        t2 = time.perf_counter()
        _cache[name] = (g, x, t1 - t0, t2 - t1)
    return _cache[name][:2]


def layer_geometry(b, layer):
    return {b.nets[k]: shapely.wkb.loads(v) for k, v in b.copper_wkb(layer).items()}


@pytest.mark.parametrize("name", BOARDS)
def test_board_structure(name, validation):
    g, x = boards(name)
    assert g.net_source == "x2" and x.net_source == "ipc2581"
    assert len(g.copper_layers) == len(x.copper_layers)
    validation.check(name, "stackup thickness", sum(l["thickness_mm"] for l in x.stackup),
                     sum(l["thickness_mm"] for l in g.stackup), rel=1e-9, unit="mm", source="gbrjob")
    validation.check(name, "outline area", shapely.wkb.loads(x.outline_wkb()).area,
                     shapely.wkb.loads(g.outline_wkb()).area, rel=1e-4, unit="mm2", source="Gerber profile")
    assert set(g.nets) == set(x.nets)
    # The X2 nets agree with the IPC-D-356 netlist exported with them.
    assert not [w for w in g.warnings if "disagree" in w], g.warnings


@pytest.mark.parametrize("name", BOARDS)
def test_copper(name, validation):
    g, x = boards(name)
    for c in g.copper_layers:
        L = c["index"]
        gg, xx = layer_geometry(g, L), layer_geometry(x, L)
        ga, xa = shapely.union_all(list(gg.values())), shapely.union_all(list(xx.values()))
        diff = ga.symmetric_difference(xa)
        largest = max((p.area for p in getattr(diff, "geoms", [diff])), default=0.0)
        assert largest < 0.05, (L, largest)
        # The two exports place circle vertices at different angles, which leaves a
        # small symmetric difference along every round edge.
        validation.check(f"{name} L{L}", "copper area differing (XOR / area)", diff.area / ga.area, 0.0, abs_=5e-4,
                         source="Gerber X2 import")
        for net in set(gg) | set(xx):
            a = gg[net].area if net in gg else 0.0
            b = xx[net].area if net in xx else 0.0
            assert b == pytest.approx(a, rel=2e-3, abs=2e-3), (L, net)


@pytest.mark.parametrize("name", BOARDS)
def test_drills(name, validation):
    g, x = boards(name)
    validation.check(name, "drill count", len(x.drills), len(g.drills), abs_=0, source="Excellon")
    def key(d):
        # Slots are compared by centre and length: the formats order their ends differently.
        a = d["position_mm"]
        b = d.get("slot_end_mm", a)
        return (round((a[0] + b[0]) / 2, 3), round((a[1] + b[1]) / 2, 3), round(math.dist(a, b), 3),
                round(d["diameter_mm"], 3), d["plated"])

    gk = collections.Counter(map(key, g.drills))
    xk = collections.Counter(map(key, x.drills))
    assert gk == xk
    # Nets of plated holes: derived from copper geometry (Gerber) against stated (IPC-2581).
    xnet = {key(d): x.nets[d["net"]] for d in x.drills}
    plated = [d for d in g.drills if d["plated"]]
    agree = sum(g.nets[d["net"]] == xnet[key(d)] for d in plated)
    validation.check(name, "plated hole nets agreeing", agree, len(plated), abs_=0, source="IPC-2581 Set net")


@pytest.mark.parametrize("name", BOARDS)
def test_pads(name, validation):
    g, x = boards(name)

    def index(b):
        out = collections.defaultdict(list)
        for i, p in enumerate(b.pads):
            out[(p["component"], p["pin"], p["layer"])].append((i, p))
        return out

    gi, xi = index(g), index(x)
    named = {k for k in gi if k[1]}
    assert named <= set(xi), sorted(named - set(xi))[:5]
    worst = 0.0
    for k in named:
        a = shapely.union_all([shapely.wkb.loads(g.pad_wkb(i)) for i, _ in gi[k]])
        b = shapely.union_all([shapely.wkb.loads(x.pad_wkb(i)) for i, _ in xi[k]])
        worst = max(worst, a.symmetric_difference(b).area / a.area)
        assert {g.nets[p["net"]] for _, p in gi[k]} == {x.nets[p["net"]] for _, p in xi[k]}, k
    validation.check(name, "largest pad shape difference (XOR / area)", worst, 0.0, abs_=5e-3, source="IPC-2581 import")


@pytest.mark.parametrize("name", BOARDS)
def test_component_positions(name):
    g, x = boards(name)
    xc = {c["refdes"]: c for c in x.components}
    common = [c for c in g.components if c["refdes"] in xc]
    assert len(common) == len(g.components)
    for c in common:
        assert math.dist(c["position_mm"], xc[c["refdes"]]["position_mm"]) < 1e-3, c["refdes"]


@pytest.mark.parametrize("name", BOARDS)
def test_model_round_trip(name, tmp_path):
    g, _ = boards(name)
    g.save(tmp_path / "model.zip")
    h = kairoq.open(tmp_path / "model.zip")
    assert h.data == g.data
    for c in g.copper_layers:
        assert h.copper_wkb(c["index"]) == g.copper_wkb(c["index"])
    assert h.outline_wkb() == g.outline_wkb()


def test_import_times(capsys):
    for name in BOARDS:
        boards(name)
    with capsys.disabled():
        for name in BOARDS:
            _, _, tg, tx = _cache[name]
            print(f"\n  import {name}: Gerber set {tg:.2f} s, IPC-2581 {tx:.2f} s", end="")


@pytest.mark.parametrize("name", BOARDS)
def test_gerber_drill_files(name, validation):
    """Gerber X2 drill files give the same holes as the Excellon files."""
    g, _ = boards(name)
    src = TESTDATA / "kicad" / name
    files = [p for p in (src / "gerber").iterdir() if p.suffix != ".drl"] + list((src / "gerber-drill").iterdir())
    d = kairoq.open(files)

    def key(h):
        a = h["position_mm"]
        b = h.get("slot_end_mm", a)
        return (round((a[0] + b[0]) / 2, 3), round((a[1] + b[1]) / 2, 3), round(math.dist(a, b), 3),
                round(h["diameter_mm"], 3), h["plated"], h["function"], g.nets[h["net"]])

    validation.check(name, "holes from Gerber drill files", len(d.drills), len(g.drills), abs_=0, source="Excellon")
    assert collections.Counter(map(key, d.drills)) == collections.Counter(map(key, g.drills))


@pytest.mark.parametrize("name", BOARDS)
def test_tracks(name, validation):
    """The drawn tracks (centre line, width, net) agree between Gerber, IPC-2581 and ODB++,
    to the 1 um of the Gerber coordinates."""
    g, x = boards(name)
    o = kairoq.open(TESTDATA / "kicad" / name / f"{name}.odb.zip")

    def groups(b):
        out = collections.defaultdict(list)
        for t in b.tracks:
            out[(t["layer"], b.nets[t["net"]], round(t["width_mm"], 4), "centre_mm" in t)].append((t["start_mm"], t["end_mm"]))
        return out

    def same(a, b):
        # Each track of one import has a track of the other with both ends within 2 um.
        if len(a) != len(b):
            return False
        left = list(b)
        for s, e in a:
            for i, (s2, e2) in enumerate(left):
                if (math.dist(s, s2) < 2e-3 and math.dist(e, e2) < 2e-3) or (math.dist(s, e2) < 2e-3 and math.dist(e, s2) < 2e-3):
                    del left[i]
                    break
            else:
                return False
        return True

    validation.check(name, "tracks (Gerber)", len(g.tracks), len(x.tracks), abs_=0, source="IPC-2581")
    validation.check(name, "tracks (ODB++)", len(o.tracks), len(x.tracks), abs_=0, source="IPC-2581")
    gx = groups(x)
    for other in (g, o):
        go = groups(other)
        assert set(go) == set(gx)
        assert all(same(go[k], gx[k]) for k in gx)


def test_connectivity_nets_on_a_real_board(tmp_path):
    """cm5_minima's Gerber set without its X2 object attributes: nets from the copper's
    connectivity, named by the IPC-D-356 netlist, which KiCad writes in another origin and
    the importer places on the drills. Every plated hole's net is the reference net (the
    netlist keeps the last 14 characters of a name, in upper case)."""
    src = TESTDATA / "kicad" / "cm5_minima"
    ref = kairoq.open([src / "gerber", src / "gerber-drill"])
    for sub in ("gerber", "gerber-drill"):
        shutil.copytree(src / sub, tmp_path / sub)
        for f in (tmp_path / sub).iterdir():
            if f.suffix not in (".d356", ".csv", ".gbrjob"):
                f.write_text("\n".join(l for l in f.read_text().splitlines() if not l.startswith(("%TO", "; #@! TO"))) + "\n")
    b = kairoq.open([tmp_path / "gerber", tmp_path / "gerber-drill"])
    assert b.net_source == "ipc356"
    assert any("netlist is placed by an offset" in w for w in b.warnings)
    pairs = set()
    for dr, db in zip(ref.drills, b.drills):
        assert dr["position_mm"] == db["position_mm"]
        if dr["plated"]:
            pairs.add((ref.nets[dr["net"]], b.nets[db["net"]]))
    assert len(pairs) > 70
    for r, d in pairs:
        assert (d.startswith("N$") and r == "") or r.upper().endswith(d.split("#")[0].upper()), (r, d)
    # No net of the copper joins two reference nets.
    joined = {}
    for r, d in pairs:
        joined.setdefault(d, set()).add(r)
    assert all(len(v) == 1 for v in joined.values())


def test_de10_interface_board(validation):
    """The author's four-layer board: the import, and the DC drop of the 3.3 V video rail of
    examples/de10_interface.yaml (0.45 A from U305 to two LDOs and two converter ICs)."""
    b = kairoq.open(TESTDATA / "de10_interface" / "de10_interface.ipc2581.zip")
    assert (len(b.copper_layers), len(b.nets), len(b.drills), len(b.components)) == (4, 234, 619, 216)
    r = kairoq.solve_dc(b, ROOT / "examples" / "de10_interface.yaml")
    drops = {p: 3.3 - v["voltage_v"] for p, v in r.ports.items()}
    assert max(drops.values()) < 0.01
    validation.check("de10_interface +3V3_DOM2", "largest drop", max(drops.values()), 0.0049, rel=0.05, unit="V",
                     source="the same solution (regression)")


@pytest.mark.parametrize("file", ["ecc83.xml", "ecc83.odb.zip"])
def test_stackup_yaml_overrides_ipc2581_and_odbpp(file, tmp_path):
    """A stackup file given with the board replaces the stackup of an IPC-2581 or ODB++ file."""
    y = tmp_path / "stack.yaml"
    y.write_text("layers:\n  - {name: F.Mask, type: soldermask, thickness_mm: 0.02}\n"
                 "  - {name: top, type: copper, thickness_mm: 0.07}\n"
                 "  - {name: core, type: dielectric, thickness_mm: 0.8, material: FR-4}\n"
                 "  - {name: bot, type: copper, thickness_mm: 0.07}\n")
    b = kairoq.open(TESTDATA / "kicad" / "ecc83" / file, stackup=y)
    assert b.data["stackup"]["source"] == "yaml"
    assert [l["thickness_mm"] for l in b.stackup] == pytest.approx([0.02, 0.07, 0.8, 0.07])
    y.write_text("layers:\n  - {name: top, type: copper, thickness_mm: 0.035}\n")
    with pytest.raises(Exception, match="copper layers"):
        kairoq.open(TESTDATA / "kicad" / "ecc83" / file, stackup=y)
