"""Cross-check of the ODB++ importer on KiCad demonstration boards.

The ODB++ export of each board must give the same canonical model as its IPC-2581
export: copper, nets, pads, holes, stackup, outline and components. The archive forms an
ODB++ job comes in (directory, zip, tgz, members compressed with Unix compress) must give
the same board.
"""

from __future__ import annotations

import collections
import math
import shutil
import subprocess
import tarfile
import zipfile

import pytest
import shapely
import shapely.wkb

import kairoq
from conftest import TESTDATA
from test_kicad import BOARDS, ipc2581_file, layer_geometry

_cache: dict[str, tuple] = {}


def boards(name):
    if name not in _cache:
        d = TESTDATA / "kicad" / name
        _cache[name] = (kairoq.open(d / f"{name}.odb.zip"), kairoq.open(ipc2581_file(name)), kairoq.open(d / "gerber"))
    return _cache[name]


@pytest.mark.parametrize("name", BOARDS)
def test_structure(name, validation):
    o, x, _ = boards(name)
    assert o.net_source == "odbpp"
    assert not o.warnings or o.warnings == x.warnings, o.warnings
    assert len(o.copper_layers) == len(x.copper_layers)
    assert set(o.nets) == set(x.nets)
    validation.check(name, "stackup thickness (ODB++)", sum(l["thickness_mm"] for l in o.stackup),
                     sum(l["thickness_mm"] for l in x.stackup), rel=1e-6, unit="mm", source="IPC-2581")
    eps = lambda b: [(l["type"], l.get("epsilon_r")) for l in b.stackup if l["type"] in ("copper", "dielectric")]
    assert eps(o) == eps(x)
    validation.check(name, "outline area (ODB++)", shapely.wkb.loads(o.outline_wkb()).area,
                     shapely.wkb.loads(x.outline_wkb()).area, rel=1e-5, unit="mm2", source="IPC-2581")


@pytest.mark.parametrize("name", BOARDS)
def test_copper(name, validation):
    o, x, _ = boards(name)
    for c in x.copper_layers:
        L = c["index"]
        oo, xx = layer_geometry(o, L), layer_geometry(x, L)
        oa, xa = shapely.union_all(list(oo.values())), shapely.union_all(list(xx.values()))
        diff = oa.symmetric_difference(xa)
        validation.check(f"{name} L{L}", "copper area differing, ODB++ (XOR / area)", diff.area / xa.area, 0.0, abs_=5e-4,
                         source="IPC-2581 import")
        for net in set(oo) | set(xx):
            a = oo[net].area if net in oo else 0.0
            b = xx[net].area if net in xx else 0.0
            assert a == pytest.approx(b, rel=2e-3, abs=2e-3), (L, net)


@pytest.mark.parametrize("name", BOARDS)
def test_drills(name, validation):
    o, x, _ = boards(name)

    def key(b, d):
        a = d["position_mm"]
        e = d.get("slot_end_mm", a)
        return (round((a[0] + e[0]) / 2, 3), round((a[1] + e[1]) / 2, 3), round(math.dist(a, e), 3),
                round(d["diameter_mm"], 3), d["plated"], b.nets[d["net"]])

    validation.check(name, "drill count (ODB++)", len(o.drills), len(x.drills), abs_=0, source="IPC-2581")
    assert collections.Counter(key(o, d) for d in o.drills) == collections.Counter(key(x, d) for d in x.drills)


@pytest.mark.parametrize("name", BOARDS)
def test_pads(name, validation):
    o, x, _ = boards(name)

    def index(b):
        out = collections.defaultdict(list)
        for i, p in enumerate(b.pads):
            out[(p["component"], p["pin"], p["layer"])].append((i, p))
        return out

    oi, xi = index(o), index(x)
    named = {k for k in xi if k[1]}
    assert named == {k for k in oi if k[1]}
    worst = 0.0
    for k in named:
        a = shapely.union_all([shapely.wkb.loads(o.pad_wkb(i)) for i, _ in oi[k]])
        b = shapely.union_all([shapely.wkb.loads(x.pad_wkb(i)) for i, _ in xi[k]])
        worst = max(worst, a.symmetric_difference(b).area / b.area)
        assert {o.nets[p["net"]] for _, p in oi[k]} == {x.nets[p["net"]] for _, p in xi[k]}, k
    validation.check(name, "largest pad shape difference, ODB++ (XOR / area)", worst, 0.0, abs_=5e-3, source="IPC-2581 import")


@pytest.mark.parametrize("name", BOARDS)
def test_components(name):
    """ODB++ lists the parts with pads; IPC-2581 also lists footprints without any (logos),
    the placement file only the parts marked for assembly."""
    o, x, g = boards(name)
    oc, xc = {c["refdes"]: c for c in o.components}, {c["refdes"]: c for c in x.components}
    assert {c["refdes"] for c in g.components} <= set(oc) <= set(xc)
    for r, c in oc.items():
        assert math.dist(c["position_mm"], xc[r]["position_mm"]) < 1e-3, r
        assert c["side"] == xc[r]["side"], r
        assert (c["rotation_deg"] - xc[r]["rotation_deg"] + 180) % 360 - 180 == pytest.approx(0, abs=1e-6), r
    for c in g.components:
        assert (c["rotation_deg"] - oc[c["refdes"]]["rotation_deg"] + 180) % 360 - 180 == pytest.approx(0, abs=1e-6), c["refdes"]


def compress_lzw(data: bytes, maxbits: int = 16, clear: bool = False, stats: dict | None = None) -> bytes:
    """Unix compress (.Z), block mode. With clear, the table is cleared whenever it is
    full (ncompress clears when the ratio drops), so the stream holds CLEAR codes."""
    out = bytearray([0x1F, 0x9D, 0x80 | maxbits])
    acc = nacc = 0      # bit accumulator
    group = 0           # bits written in the current group of eight codes
    n_bits, maxcode, free = 9, 511, 257
    table: dict[tuple[int, int], int] = {}

    def emit(code, clearing=False):
        nonlocal acc, nacc, group, n_bits, maxcode
        acc |= code << nacc
        nacc += n_bits
        group = (group + n_bits) % (8 * n_bits)
        while nacc >= 8:
            out.append(acc & 0xFF)
            acc >>= 8
            nacc -= 8
        if free > maxcode or clearing:
            # The group is padded to its full length before the code width changes.
            if group:
                pad = 8 * n_bits - group
                nacc += pad
                while nacc >= 8:
                    out.append(acc & 0xFF)
                    acc >>= 8
                    nacc -= 8
                group = 0
            if clearing:
                n_bits, maxcode = 9, 511
            else:
                n_bits += 1
                maxcode = (1 << maxbits) if n_bits == maxbits else (1 << n_bits) - 1

    if not data:
        return bytes(out)
    ent = data[0]
    for c in data[1:]:
        k = (ent, c)
        if k in table:
            ent = table[k]
            continue
        emit(ent)
        ent = c
        if free < (1 << maxbits):
            table[k] = free
            free += 1
        elif clear:
            table.clear()
            free = 257
            emit(256, clearing=True)
            if stats is not None:
                stats["clears"] = stats.get("clears", 0) + 1
    emit(ent)
    if nacc:
        out.append(acc & 0xFF)
    return bytes(out)


def extract(name, tmp_path):
    with zipfile.ZipFile(TESTDATA / "kicad" / name / f"{name}.odb.zip") as z:
        z.extractall(tmp_path / "dir")
    return tmp_path / "dir"


def same_board(a, b):
    strip = lambda d: {k: v for k, v in d.items() if k not in ("name", "source")}
    assert strip(a.data) == strip(b.data)
    for c in a.copper_layers:
        assert a.copper_wkb(c["index"]) == b.copper_wkb(c["index"])


def test_archive_forms(tmp_path):
    ref = boards("ecc83")[0]
    job = extract("ecc83", tmp_path)
    same_board(kairoq.open(job), ref)
    for fmt, suffix in ((tarfile.PAX_FORMAT, "pax.tgz"), (tarfile.GNU_FORMAT, "gnu.tgz")):
        with tarfile.open(tmp_path / suffix, "w:gz", format=fmt) as t:
            # A root longer than the 100 characters of a tar header's name field.
            t.add(job, arcname="job_" + "x" * 120)
        same_board(kairoq.open(tmp_path / suffix), ref)


def test_compressed_members(tmp_path):
    dst = extract("cm5_minima", tmp_path)
    # Files the importer reads: the matrix, the netlist, the top copper and its components.
    chosen = [dst / "matrix" / "matrix", dst / "steps" / "pcb" / "eda" / "data",
              dst / "steps" / "pcb" / "layers" / "f.cu" / "features", dst / "steps" / "pcb" / "layers" / "comp_+_top" / "components"]
    largest = 0
    for p in chosen:
        raw = p.read_bytes()
        z = compress_lzw(raw)
        if shutil.which("gzip"):  # gzip reads .Z: an independent check of the compressor
            assert subprocess.run(["gzip", "-dc"], input=z, capture_output=True, check=True).stdout == raw
        p.with_name(p.name + ".Z").write_bytes(z)
        p.unlink()
        largest = max(largest, len(raw))
    # Large enough for every code width up to 16 bits.
    assert largest > 200_000
    same_board(kairoq.open(dst), boards("cm5_minima")[0])


def test_compressed_members_with_clears(tmp_path):
    """12-bit codes with the table cleared whenever it fills: the CLEAR code, the padding
    of its group and the return to 9-bit codes."""
    dst = extract("cm5_minima", tmp_path)
    p = dst / "steps" / "pcb" / "layers" / "f.cu" / "features"
    raw = p.read_bytes()
    stats = {}
    z = compress_lzw(raw, maxbits=12, clear=True, stats=stats)
    assert stats["clears"] >= 10
    if shutil.which("gzip"):
        assert subprocess.run(["gzip", "-dc"], input=z, capture_output=True, check=True).stdout == raw
    p.with_name(p.name + ".Z").write_bytes(z)
    p.unlink()
    same_board(kairoq.open(dst), boards("cm5_minima")[0])


def test_drill_plating_without_attributes(tmp_path):
    """Without the .drill attribute, plating comes from the layer's tools file (the tool of
    the finished size); without tools, from the layer's name, with a warning."""
    ref = boards("cm5_minima")[0]
    dst = extract("cm5_minima", tmp_path)
    layers = list((dst / "steps" / "pcb" / "layers").glob("drill*"))
    assert layers
    for d in layers:
        f = d / "features"
        f.write_text("\n".join(line.split(";")[0].rstrip() if line[:1] in "PL" else line for line in f.read_text().splitlines()) + "\n")
    b = kairoq.open(dst)
    same_board(b, ref)
    assert not any("taken from the layer's name" in w for w in b.data["warnings"])
    # The name tells plated from non-plated, not a via from a component hole.
    for d in layers:
        (d / "tools").unlink()
    b = kairoq.open(dst)
    key = lambda d: (d["position_mm"], d["diameter_mm"], d["plated"])
    assert [key(d) for d in b.data["drills"]] == [key(d) for d in ref.data["drills"]]
    assert any("taken from the layer's name" in w for w in b.data["warnings"])
