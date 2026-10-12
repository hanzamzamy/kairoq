"""Synthetic fabrication sets with known geometry, written independently of the core.

Each generator writes Gerber X2 copper and profile files, an Excellon drill file, a
Gerber job file and, optionally, an IPC-D-356 netlist into a directory, and returns the
analytical quantities the import must reproduce.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path


def _c(v_mm: float) -> str:
    """Coordinate in the 4.6 format (nm integer)."""
    return str(round(v_mm * 1e6))


@dataclass
class Gerber:
    function: str
    lines: list[str] = field(default_factory=list)
    _apertures: dict[str, int] = field(default_factory=dict)
    _defs: list[str] = field(default_factory=list)

    def aperture(self, template: str, function: str | None = None) -> int:
        key = f"{template}|{function}"
        if key not in self._apertures:
            code = 10 + len(self._apertures)
            self._apertures[key] = code
            if function:
                self._defs.append(f"%TA.AperFunction,{function}*%")
            self._defs.append(f"%ADD{code}{template}*%")
            if function:
                self._defs.append("%TD*%")
        return self._apertures[key]

    def net(self, name: str | None):
        self.lines.append(f"%TO.N,{name}*%" if name is not None else "%TD.N*%")

    def track(self, points, width: float, net: str | None = None):
        d = self.aperture(f"C,{width:.6f}", "Conductor")
        self.net(net)
        self.lines.append(f"D{d}*")
        x0, y0 = points[0]
        self.lines.append(f"X{_c(x0)}Y{_c(y0)}D02*")
        for x, y in points[1:]:
            self.lines.append(f"X{_c(x)}Y{_c(y)}D01*")
        self.lines.append("%TD*%")

    def flash(self, template: str, x: float, y: float, net: str | None = None, pad: tuple[str, str] | None = None,
              function: str = "ComponentPad", clear: bool = False):
        d = self.aperture(template, function)
        if clear:
            self.lines.append("%LPC*%")
        self.net(net)
        if pad:
            self.lines.append(f"%TO.P,{pad[0]},{pad[1]}*%")
        self.lines.append(f"D{d}*")
        self.lines.append(f"X{_c(x)}Y{_c(y)}D03*")
        self.lines.append("%TD*%")
        if clear:
            self.lines.append("%LPD*%")

    def region(self, points, net: str | None = None, clear: bool = False):
        if clear:
            self.lines.append("%LPC*%")
        self.net(net)
        self.lines.append("G36*")
        x0, y0 = points[0]
        self.lines.append(f"X{_c(x0)}Y{_c(y0)}D02*")
        for x, y in points[1:] + [points[0]]:
            self.lines.append(f"X{_c(x)}Y{_c(y)}D01*")
        self.lines.append("G37*")
        self.lines.append("%TD*%")
        if clear:
            self.lines.append("%LPD*%")

    def outline(self, points):
        d = self.aperture("C,0.100000", "Profile")
        self.lines.append(f"D{d}*")
        x0, y0 = points[0]
        self.lines.append(f"X{_c(x0)}Y{_c(y0)}D02*")
        for x, y in points[1:] + [points[0]]:
            self.lines.append(f"X{_c(x)}Y{_c(y)}D01*")

    def text(self) -> str:
        head = [
            "%TF.GenerationSoftware,kairoq,synthetic,1*%",
            f"%TF.FileFunction,{self.function}*%",
            "%TF.FilePolarity,Positive*%",
            "%FSLAX46Y46*%",
            "%MOMM*%",
            "%LPD*%",
            "G01*",
            "G75*",
        ]
        return "\n".join(head + self._defs + self.lines + ["M02*", ""])


def _excellon(holes, plated: bool, span=(1, 2), function="ViaDrill") -> str:
    lines = ["M48", f"; #@! TF.FileFunction,{'Plated' if plated else 'NonPlated'},{span[0]},{span[1]},{'PTH' if plated else 'NPTH'}",
             "METRIC"]
    tools = sorted({d for _, _, d in holes})
    for i, d in enumerate(tools, 1):
        lines.append(f"; #@! TA.AperFunction,{'Plated' if plated else 'NonPlated'},{'PTH' if plated else 'NPTH'},{function}")
        lines.append(f"T{i}C{d:.3f}")
    lines += ["%", "G90", "G05"]
    for i, d in enumerate(tools, 1):
        lines.append(f"T{i}")
        lines += [f"X{x:.4f}Y{y:.4f}" for x, y, dd in holes if dd == d]
    lines.append("M30")
    return "\n".join(lines) + "\n"


def _job(name: str, files: dict[str, str], copper_thickness_mm: list[float], dielectric_mm: list[float], epsilon_r="4.5",
         mask_mm=0.01) -> str:
    stack = [{"Type": "SolderMask", "Thickness": mask_mm, "Name": "top mask"}] if mask_mm > 0 else []
    for i, t in enumerate(copper_thickness_mm):
        stack.append({"Type": "Copper", "Thickness": t, "Name": f"L{i + 1}"})
        if i < len(dielectric_mm):
            stack.append({"Type": "Dielectric", "Thickness": dielectric_mm[i], "Material": "FR4",
                          "DielectricConstant": epsilon_r, "Name": f"D{i + 1}"})
    if mask_mm > 0:
        stack.append({"Type": "SolderMask", "Thickness": mask_mm, "Name": "bottom mask"})
    job = {
        "Header": {"GenerationSoftware": {"Vendor": "kairoq", "Application": "synthetic", "Version": "1"}},
        "GeneralSpecs": {"ProjectId": {"Name": name}, "LayerNumber": len(copper_thickness_mm),
                         "BoardThickness": sum(copper_thickness_mm) + sum(dielectric_mm) + 2 * mask_mm},
        "FilesAttributes": [{"Path": p, "FileFunction": f, "FilePolarity": "Positive"} for p, f in files.items()],
        "MaterialStackup": stack,
    }
    return json.dumps(job, indent=1)


def _write(out: Path, name: str, gerbers: dict[str, Gerber], drills: dict[str, str], copper_t, dielectric, extra=None, **job):
    out.mkdir(parents=True, exist_ok=True)
    for fname, g in gerbers.items():
        (out / fname).write_text(g.text())
    for fname, text in drills.items():
        (out / fname).write_text(text)
    for fname, text in (extra or {}).items():
        (out / fname).write_text(text)
    (out / f"{name}-job.gbrjob").write_text(_job(name, {f: g.function for f, g in gerbers.items()}, copper_t, dielectric, **job))


def strip(out: Path, length=50.0, width=1.0) -> dict:
    """One straight track between two square pads on a 2-layer board."""
    top, bot, prof = Gerber("Copper,L1,Top"), Gerber("Copper,L2,Bot"), Gerber("Profile,NP")
    top.track([(5.0, 10.0), (5.0 + length, 10.0)], width, "SIG")
    top.flash("R,2.000000X2.000000", 5.0, 10.0, "SIG", ("J1", "1"), "SMDPad")
    top.flash("R,2.000000X2.000000", 5.0 + length, 10.0, "SIG", ("J2", "1"), "SMDPad")
    prof.outline([(0, 0), (length + 10, 0), (length + 10, 20), (0, 20)])
    _write(out, "strip", {"strip-L1.gtl": top, "strip-L2.gbl": bot, "strip-Edge.gm1": prof}, {}, [0.035, 0.035], [1.5])
    # Track: rectangle between the pad centres plus two half discs, of which the
    # parts inside the pads do not add area.
    pad = 4.0
    track_outside = (length - 2.0) * width
    return {"copper_area_l1": 2 * pad + track_outside, "outline_area": (length + 10) * 20, "nets": ["SIG"]}


def via_chain(out: Path, n=10, pitch=2.0, via_drill=0.3, via_pad=0.6, width=0.25) -> dict:
    """Two-layer daisy chain: segments alternate between the layers through vias."""
    top, bot, prof = Gerber("Copper,L1,Top"), Gerber("Copper,L2,Bot"), Gerber("Profile,NP")
    holes = []
    for i in range(n):
        x = 5.0 + i * pitch
        holes.append((x, 5.0, via_drill))
        for g in (top, bot):
            g.flash(f"C,{via_pad:.6f}", x, 5.0, "CHAIN", None, "ViaPad")
        if i + 1 < n:
            (top if i % 2 == 0 else bot).track([(x, 5.0), (x + pitch, 5.0)], width, "CHAIN")
    w = 10.0 + (n - 1) * pitch
    prof.outline([(0, 0), (w, 0), (w, 10), (0, 10)])
    netlist = ["P  UNITS CUST 1"] + [
        f"317{'CHAIN':<14}   {'VIA':<6}      D{round(via_drill * 1000):04d}PA00X+{round((5.0 + i * pitch) * 1000):06d}Y+{5000:06d}"
        for i in range(n)] + ["999"]
    _write(out, "chain", {"chain-L1.gtl": top, "chain-L2.gbl": bot, "chain-Edge.gm1": prof},
           {"chain-PTH.drl": _excellon(holes, True)}, [0.035, 0.035], [1.5], {"chain.d356": "\n".join(netlist) + "\n"})
    return {"vias": n, "nets": ["CHAIN"], "via_drill": via_drill}


def plane_with_clearance(out: Path, size=30.0, hole=4.0) -> dict:
    """Four-layer board: a full plane on L2 with a cleared square, a disc on L3, slot cut-out."""
    l1, l2, l3, l4, prof = (Gerber("Copper,L1,Top"), Gerber("Copper,L2,Inr"), Gerber("Copper,L3,Inr"),
                            Gerber("Copper,L4,Bot"), Gerber("Profile,NP"))
    m = 1.0
    l2.region([(m, m), (size - m, m), (size - m, size - m), (m, size - m)], "GND")
    c = size / 2
    l2.region([(c - hole / 2, c - hole / 2), (c + hole / 2, c - hole / 2), (c + hole / 2, c + hole / 2),
               (c - hole / 2, c + hole / 2)], None, clear=True)
    l3.flash("C,6.000000", c, c, "PWR", None, "Conductor")
    l1.flash("R,1.000000X1.000000", 3, 3, "GND", ("TP1", "1"), "SMDPad")
    l4.flash("R,1.000000X1.000000", 3, 3, "GND", ("TP2", "1"), "SMDPad")
    prof.outline([(0, 0), (size, 0), (size, size), (0, size)])
    _write(out, "plane", {"plane-L1.gtl": l1, "plane-L2.g1": l2, "plane-L3.g2": l3, "plane-L4.gbl": l4,
                          "plane-Edge.gm1": prof},
           {"plane-PTH.drl": _excellon([(3.0, 3.0, 0.3)], True)}, [0.035, 0.0175, 0.0175, 0.035], [0.2, 1.0, 0.2])
    return {
        "l2_area": (size - 2 * m) ** 2 - hole ** 2,
        "l3_area": math.pi * 9.0,
        "thickness": 0.035 * 2 + 0.0175 * 2 + 1.4 + 0.02,
    }


def via_spreading(out: Path, r_outer=5.0, ring=0.5, drill=0.3, dielectric=1.51, copper=0.035) -> dict:
    """Disc of copper on both layers joined by one central via; annular ports at the rims.

    Current enters the top rim, spreads radially to the via, crosses the barrel, and leaves
    radially through the bottom rim.
    """
    top, bot, prof = Gerber("Copper,L1,Top"), Gerber("Copper,L2,Bot"), Gerber("Profile,NP")
    c = r_outer + 2.0
    for g, ref in ((top, "P1"), (bot, "P2")):
        g.flash(f"C,{2 * r_outer:.6f}", c, c, "N", None, "Conductor")
        # Annular port: a circular aperture with a hole.
        g.flash(f"C,{2 * r_outer:.6f}X{2 * (r_outer - ring):.6f}", c, c, "N", (ref, "1"), "ComponentPad")
    prof.outline([(0, 0), (2 * c, 0), (2 * c, 2 * c), (0, 2 * c)])
    _write(out, "spread", {"spread-L1.gtl": top, "spread-L2.gbl": bot, "spread-Edge.gm1": prof},
           {"spread-PTH.drl": _excellon([(c, c, drill)], True)}, [copper, copper], [dielectric])
    return {"r_inner_port": r_outer - ring, "r_via": drill / 2, "barrel_length": dielectric + copper, "copper": copper}


def annulus(out: Path, r_outer=5.0, ring=0.5, r_inner=1.0, copper=0.035, dielectric=1.5) -> dict:
    """Single copper layer: a disc with an annular port at its rim (P1.1) and a round port at
    its centre (P2.1); the field between them is radial."""
    top, prof = Gerber("Copper,L1,Top"), Gerber("Profile,NP")
    c = r_outer + 2.0
    top.flash(f"C,{2 * r_outer:.6f}", c, c, "N", None, "Conductor")
    top.flash(f"C,{2 * r_outer:.6f}X{2 * (r_outer - ring):.6f}", c, c, "N", ("P1", "1"), "ComponentPad")
    top.flash(f"C,{2 * r_inner:.6f}", c, c, "N", ("P2", "1"), "ComponentPad")
    prof.outline([(0, 0), (2 * c, 0), (2 * c, 2 * c), (0, 2 * c)])
    _write(out, "annulus", {"annulus-L1.gtl": top, "annulus-Edge.gm1": prof}, {}, [copper], [dielectric])
    return {"r_outer_port": r_outer - ring, "r_inner_port": r_inner}


def neck(out: Path, width=4.0, length=40.0, aperture=0.5, wall=0.02, pad=2.0) -> dict:
    """Strip with a thin wall across its middle, open over `aperture`; pads across both ends."""
    top, bot, prof = Gerber("Copper,L1,Top"), Gerber("Copper,L2,Bot"), Gerber("Profile,NP")
    x0, y0 = 5.0, 5.0
    xm, ym = x0 + length / 2, y0 + width / 2
    h = wall / 2
    a = aperture / 2
    # Outline of the strip with two notches leaving the aperture, counter-clockwise.
    pts = [(x0, y0), (xm - h, y0), (xm - h, ym - a), (xm + h, ym - a), (xm + h, y0), (x0 + length, y0),
           (x0 + length, y0 + width), (xm + h, y0 + width), (xm + h, ym + a), (xm - h, ym + a), (xm - h, y0 + width),
           (x0, y0 + width)]
    top.region(pts, "N")
    top.flash(f"R,{pad:.6f}X{width:.6f}", x0 + pad / 2, ym, "N", ("J1", "1"), "SMDPad")
    top.flash(f"R,{pad:.6f}X{width:.6f}", x0 + length - pad / 2, ym, "N", ("J2", "1"), "SMDPad")
    prof.outline([(0, 0), (length + 10, 0), (length + 10, width + 10), (0, width + 10)])
    _write(out, "neck", {"neck-L1.gtl": top, "neck-L2.gbl": bot, "neck-Edge.gm1": prof}, {}, [0.035, 0.035], [1.5])
    return {"free_length": length - 2 * pad, "width": width, "aperture": aperture, "wall": wall}


def power_board(out: Path, size=100.0, pitch=2.0, clearance=0.1, via_drill=0.3, via_pad=0.6) -> dict:
    """4-layer 100 x 100 mm board for timing: a 5 V plane (L3) and a ground plane (L2) with
    a stitching via every `pitch` mm, a ground pour on L4 around them with `clearance`
    gaps, and 0.1 mm wide power tracks on L1. A connector feeds the planes at one edge and
    a load draws the current at the other.
    """
    l1, l2, l3, l4, prof = (Gerber("Copper,L1,Top"), Gerber("Copper,L2,Inr"), Gerber("Copper,L3,Inr"),
                            Gerber("Copper,L4,Bot"), Gerber("Profile,NP"))
    m = 1.0
    plane = [(m, m), (size - m, m), (size - m, size - m), (m, size - m)]
    l2.region(plane, "GND")
    l3.region(plane, "+5V")
    l4.region(plane, "GND")
    gnd_vias, pwr_vias = [], []
    n = int((size - 4) / pitch)
    for i in range(n):
        for j in range(n):
            x, y = 3 + i * pitch, 3 + j * pitch
            if (i + j) % 2 == 0:
                gnd_vias.append((x, y))
            else:
                pwr_vias.append((x, y))
    anti = via_pad + 2 * clearance
    for x, y in pwr_vias:  # clearances in the ground layers around the power vias
        l2.flash(f"C,{anti:.6f}", x, y, None, None, "Other", clear=True)
        l4.flash(f"C,{anti:.6f}", x, y, None, None, "Other", clear=True)
    for x, y in gnd_vias:  # and in the power plane around the ground vias
        l3.flash(f"C,{anti:.6f}", x, y, None, None, "Other", clear=True)
    for x, y in pwr_vias:
        l1.flash(f"C,{via_pad:.6f}", x, y, "+5V", None, "ViaPad")
        l3.flash(f"C,{via_pad:.6f}", x, y, "+5V", None, "ViaPad")
    for x, y in gnd_vias:
        for g in (l1, l2, l4):
            g.flash(f"C,{via_pad:.6f}", x, y, "GND", None, "ViaPad")
    # Fine power tracks on L1 joining power vias along diagonals, between the ground vias.
    for start in range(-n, n):
        diag = [(3 + i * pitch, 3 + (i + start) * pitch) for i in range(n) if 0 <= i + start < n and (2 * i + start) % 2 == 1]
        if len(diag) > 1:
            l1.track(diag, 0.1, "+5V")
    # Connector: a 5 V pad on L1 over a column of its own vias (with clearances in the
    # ground layers), a ground pad on L4. Load: 1 mm pads on a power and a ground via.
    feed = [(1.5, size / 2 - 9 + 2 * k) for k in range(10)]
    for x, y in feed:
        l2.flash(f"C,{anti:.6f}", x, y, None, None, "Other", clear=True)
        l4.flash(f"C,{anti:.6f}", x, y, None, None, "Other", clear=True)
        l3.flash(f"C,{via_pad:.6f}", x, y, "+5V", None, "ViaPad")
    l1.flash("R,1.000000X20.000000", 1.5, size / 2, "+5V", ("J1", "1"), "SMDPad")
    l4.flash("R,1.000000X20.000000", 1.5, size / 2 + 25, "GND", ("J1", "2"), "SMDPad")
    px, py = 3 + 46 * pitch, 3 + 23 * pitch  # (i + j) odd: a power via
    l1.flash("R,1.000000X1.000000", px, py, "+5V", ("U1", "1"), "SMDPad")
    l4.flash("R,1.000000X1.000000", px - pitch, py, "GND", ("U1", "2"), "SMDPad")
    pwr_vias = pwr_vias + feed
    prof.outline([(0, 0), (size, 0), (size, size), (0, size)])
    _write(out, "power", {"power-L1.gtl": l1, "power-L2.g1": l2, "power-L3.g2": l3, "power-L4.gbl": l4, "power-Edge.gm1": prof},
           {"power-PTH.drl": _excellon([(x, y, via_drill) for x, y in gnd_vias + pwr_vias], True, (1, 4))},
           [0.035, 0.035, 0.035, 0.035], [0.2, 1.065, 0.2])
    return {"vias": len(gnd_vias) + len(pwr_vias)}


def slab(out: Path, size=20.0, copper=0.035, dielectric=1.5) -> dict:
    """Two full copper planes; a pad covering the whole top plane receives the heat."""
    top, bot, prof = Gerber("Copper,L1,Top"), Gerber("Copper,L2,Bot"), Gerber("Profile,NP")
    sq = [(0, 0), (size, 0), (size, size), (0, size)]
    top.region(sq, "A")
    bot.region(sq, "B")
    top.flash(f"R,{size:.6f}X{size:.6f}", size / 2, size / 2, "A", ("H1", "1"), "SMDPad")
    prof.outline(sq)
    _write(out, "slab", {"slab-L1.gtl": top, "slab-L2.gbl": bot, "slab-Edge.gm1": prof}, {}, [copper, copper], [dielectric])
    return {"area_m2": (size * 1e-3) ** 2, "dielectric_m": dielectric * 1e-3, "mask_m": 10e-6}


def fin(out: Path, width=5.0, length=60.0, pad=2.0, copper=0.035, dielectric=1.5) -> dict:
    """Single copper layer: a strip whose board outline is the strip itself, with pads
    (J1.1, J2.1) across both ends; one-dimensional along x."""
    top, prof = Gerber("Copper,L1,Top"), Gerber("Profile,NP")
    strip = [(0, 0), (length, 0), (length, width), (0, width)]
    top.region(strip, "N")
    top.flash(f"R,{pad:.6f}X{width:.6f}", pad / 2, width / 2, "N", ("J1", "1"), "SMDPad")
    top.flash(f"R,{pad:.6f}X{width:.6f}", length - pad / 2, width / 2, "N", ("J2", "1"), "SMDPad")
    prof.outline(strip)
    _write(out, "fin", {"fin-L1.gtl": top, "fin-Edge.gm1": prof}, {}, [copper], [dielectric])
    return {"width_m": width * 1e-3, "length_m": length * 1e-3, "pad_m": pad * 1e-3, "copper_m": copper * 1e-6 * 1e3,
            "dielectric_m": dielectric * 1e-3, "mask_m": 10e-6}


def source_board(out: Path, size=30.0, pad=(14.0, 14.0, 16.0, 16.0), plane=True, dielectric=1.5, copper=0.035,
                 inner=False) -> dict:
    """Square board with one rectangular pad H1.1 (a heat source) on the top layer, or on
    the middle layer of three when inner, and an optional plane on the bottom layer
    (tools/ref3d.py solves the same board in three dimensions)."""
    names = ["Top", "Inr", "Bot"] if inner else ["Top", "Bot"]
    gs = [Gerber(f"Copper,L{i + 1},{n}") for i, n in enumerate(names)]
    prof = Gerber("Profile,NP")
    x0, y0, x1, y1 = pad
    gs[1 if inner else 0].flash(f"R,{x1 - x0:.6f}X{y1 - y0:.6f}", (x0 + x1) / 2, (y0 + y1) / 2, "H", ("H1", "1"), "SMDPad")
    sq = [(0, 0), (size, 0), (size, size), (0, size)]
    if plane:
        gs[-1].region(sq, "GND")
    prof.outline(sq)
    files = {f"ref-L{i + 1}.gbr": g for i, g in enumerate(gs)}
    files["ref-Edge.gm1"] = prof
    _write(out, "ref", files, {}, [copper] * len(gs), [dielectric / 2] * 2 if inner else [dielectric])
    return {"pad": pad}


def ring_loop(out: Path, a=20.0, b=10.0, width=0.5, gap=1.0, plane=False, dielectric=1.5, margin=10.0,
              openings=()) -> dict:
    """A rectangular copper ring (centre line a x b mm) on the top layer, open at the middle
    of its left side, closed by a part R1 whose pads (1 below the gap, 2 above it) end the
    ring; optionally a full plane of another net on the bottom layer, with square openings
    (x, y, side in mm, relative to the ring's lower left corner)."""
    top, bot, prof = Gerber("Copper,L1,Top"), Gerber("Copper,L2,Bot"), Gerber("Profile,NP")
    m = margin  # board margin around the ring
    x0, y0 = m, m
    yc = y0 + b / 2
    pad = max(width, 0.6)
    path = [(x0, yc - gap / 2), (x0, y0), (x0 + a, y0), (x0 + a, y0 + b), (x0, y0 + b), (x0, yc + gap / 2)]
    top.track(path, width, "RING")
    top.flash(f"R,{pad:.6f}X{pad:.6f}", x0, yc - gap / 2, "RING", ("R1", "1"), "SMDPad")
    top.flash(f"R,{pad:.6f}X{pad:.6f}", x0, yc + gap / 2, "RING", ("R1", "2"), "SMDPad")
    w, h = a + 2 * m, b + 2 * m
    if plane:
        bot.region([(0, 0), (w, 0), (w, h), (0, h)], "GND")
        for ox, oy, side in openings:
            cx, cy, r = x0 + ox, y0 + oy, side / 2
            bot.region([(cx - r, cy - r), (cx + r, cy - r), (cx + r, cy + r), (cx - r, cy + r)], clear=True)
    prof.outline([(0, 0), (w, 0), (w, h), (0, h)])
    _write(out, "ring", {"ring-L1.gtl": top, "ring-L2.gbl": bot, "ring-Edge.gm1": prof}, {}, [0.035, 0.035], [dielectric])
    return {"a_m": a * 1e-3, "b_m": b * 1e-3, "width_m": width * 1e-3, "gap_m": gap * 1e-3, "pad_m": pad * 1e-3,
            "height_m": (dielectric + 0.035) * 1e-3}


def hairpin(out: Path, length=20.0, width=1.0, gap=1.0, margin=6.0) -> dict:
    """Two parallel strips on the top layer (width, gap between their edges), joined by
    copper at the far end; at the near end a part R1 bridges them (R1.1 on the lower strip,
    R1.2 on the upper)."""
    top, prof = Gerber("Copper,L1,Top"), Gerber("Profile,NP")
    x0, y0 = margin, margin
    pitch = width + gap
    ya, yb = y0 + width / 2, y0 + width / 2 + pitch
    top.region([(x0, y0), (x0 + length + width, y0), (x0 + length + width, y0 + 2 * width + gap), (x0, y0 + 2 * width + gap),
                (x0, y0 + width + gap), (x0 + length, y0 + width + gap), (x0 + length, y0 + width), (x0, y0 + width)], "LOOP")
    p = width
    top.flash(f"R,{p:.6f}X{p:.6f}", x0 + p / 2, ya, "LOOP", ("R1", "1"), "SMDPad")
    top.flash(f"R,{p:.6f}X{p:.6f}", x0 + p / 2, yb, "LOOP", ("R1", "2"), "SMDPad")
    w, h = x0 + length + width + margin, y0 + 2 * width + gap + margin
    prof.outline([(0, 0), (w, 0), (w, h), (0, h)])
    _write(out, "hairpin", {"hairpin-L1.gtl": top, "hairpin-Edge.gm1": prof}, {}, [0.035], [1.5])
    return {"length_m": length * 1e-3, "width_m": width * 1e-3, "gap_m": gap * 1e-3}


def microstrip(out: Path, length=20.0, width=1.0, dielectric=1.5, margin=10.0, pad=1.0) -> dict:
    """A strip on the top layer over a full plane (net GND) on the bottom layer. At each end
    a part joins the strip to a pad with a via to the plane: R1 (R1.1 on the plane pad,
    R1.2 on the strip) at the near end, R2 (R2.1 on the strip, R2.2 on the plane pad) at the
    far end. The pads are squares of side pad."""
    top, bot, prof = Gerber("Copper,L1,Top"), Gerber("Copper,L2,Bot"), Gerber("Profile,NP")
    x0 = margin
    w_b, h_b = length + 2 * margin, 2 * margin + 4.0
    yc = h_b / 2
    top.track([(x0, yc), (x0 + length, yc)], width, "SIG")
    top.flash(f"R,{pad:.6f}X{pad:.6f}", x0, yc, "SIG", ("R1", "2"), "SMDPad")
    top.flash(f"R,{pad:.6f}X{pad:.6f}", x0 + length, yc, "SIG", ("R2", "1"), "SMDPad")
    holes = []
    for xg, ref in ((x0, ("R1", "1")), (x0 + length, ("R2", "2"))):
        yg = yc + 2.0
        top.flash(f"R,{pad:.6f}X{pad:.6f}", xg, yg, "GND", ref, "SMDPad")
        holes.append((xg, yg, 0.4))
    bot.region([(0, 0), (w_b, 0), (w_b, h_b), (0, h_b)], "GND")
    prof.outline([(0, 0), (w_b, 0), (w_b, h_b), (0, h_b)])
    _write(out, "ms", {"ms-L1.gtl": top, "ms-L2.gbl": bot, "ms-Edge.gm1": prof}, {"ms-PTH.drl": _excellon(holes, True)},
           [0.035, 0.035], [dielectric])
    return {"length_m": length * 1e-3, "width_m": width * 1e-3, "height_m": (dielectric + 0.035) * 1e-3}


def diff_pair(out: Path, length=20.0, width=0.2, gap=0.2, dielectric=0.2, slot=0.0, margin=5.0, extra_n=0.0) -> dict:
    """Two coupled tracks (nets DP_P and DP_N, width and gap between their edges) on the
    top layer over a GND plane on the bottom; parts R1 (R1.1 on DP_P, R1.2 on DP_N) and
    R2 at the ends. With slot > 0 the plane has a slot of that width across the pair at its
    middle; with extra_n > 0 the DP_N track runs on alone by that length past the end of
    DP_P (a skew)."""
    top, bot, prof = Gerber("Copper,L1,Top"), Gerber("Copper,L2,Bot"), Gerber("Profile,NP")
    x0, yc = margin, margin + 2.0
    pitch = width + gap
    yp, yn = yc + pitch / 2, yc - pitch / 2
    for y, net, pin, extra in ((yp, "DP_P", "1", 0.0), (yn, "DP_N", "2", extra_n)):
        top.track([(x0, y), (x0 + length + extra, y)], width, net)
        top.flash(f"R,{width:.6f}X{width:.6f}", x0, y, net, ("R1", pin), "SMDPad")
        top.flash(f"R,{width:.6f}X{width:.6f}", x0 + length + extra, y, net, ("R2", pin), "SMDPad")
    w_b, h_b = length + extra_n + 2 * margin, 2 * margin + 4.0
    bot.region([(0, 0), (w_b, 0), (w_b, h_b), (0, h_b)], "GND")
    if slot > 0:
        xs = x0 + length / 2
        bot.region([(xs - slot / 2, 0.5), (xs + slot / 2, 0.5), (xs + slot / 2, h_b - 0.5), (xs - slot / 2, h_b - 0.5)], clear=True)
    prof.outline([(0, 0), (w_b, 0), (w_b, h_b), (0, h_b)])
    _write(out, "dp", {"dp-L1.gtl": top, "dp-L2.gbl": bot, "dp-Edge.gm1": prof}, {}, [0.035, 0.035], [dielectric])
    return {"length_mm": length, "width_mm": width, "gap_mm": gap, "dielectric_mm": dielectric, "x_slot_mm": x0 + length / 2}


def coupled_lines(out: Path, length=20.0, width=0.2, gap=0.2, dielectric=0.2, stripline=False, margin=5.0) -> dict:
    """Two parallel tracks, nets A and B (width, gap between their edges), with parts at the
    ends: R1 (R1.1 on A, R1.2 on B) and R2 (R2.1 on A, R2.2 on B). Microstrip on the top
    layer over a GND plane, or stripline on the middle of three layers between two GND
    planes, each the given dielectric away."""
    if stripline:
        layers = [Gerber("Copper,L1,Top"), Gerber("Copper,L2,Inr"), Gerber("Copper,L3,Bot")]
        sig, planes, diel = layers[1], (layers[0], layers[2]), [dielectric, dielectric]
    else:
        layers = [Gerber("Copper,L1,Top"), Gerber("Copper,L2,Bot")]
        sig, planes, diel = layers[0], (layers[1],), [dielectric]
    prof = Gerber("Profile,NP")
    x0, yc = margin, margin + 2.0
    pitch = width + gap
    for y, net, pin in ((yc + pitch / 2, "A", "1"), (yc - pitch / 2, "B", "2")):
        sig.track([(x0, y), (x0 + length, y)], width, net)
        sig.flash(f"R,{width:.6f}X{width:.6f}", x0, y, net, ("R1", pin), "SMDPad")
        sig.flash(f"R,{width:.6f}X{width:.6f}", x0 + length, y, net, ("R2", pin), "SMDPad")
    w_b, h_b = length + 2 * margin, 2 * margin + 4.0
    for pl in planes:
        pl.region([(0, 0), (w_b, 0), (w_b, h_b), (0, h_b)], "GND")
    prof.outline([(0, 0), (w_b, 0), (w_b, h_b), (0, h_b)])
    files = {f"cl-L{k + 1}.g{k + 1}": g for k, g in enumerate(layers)}
    files["cl-Edge.gm1"] = prof
    _write(out, "cl", files, {}, [0.035] * len(layers), diel)
    return {"length_mm": length, "width_mm": width, "gap_mm": gap, "dielectric_mm": dielectric}


def patch(out: Path, width=37.0, length=28.5, feed_width=2.8, feed_length=10.0, dielectric=1.5, margin=20.0) -> dict:
    """A rectangular patch (net ANT; length along x, the resonant dimension) on the top
    layer over a full plane (net GND) on the bottom, fed at the middle of an edge by a
    track of feed_length from the pad R1.1."""
    top, bot, prof = Gerber("Copper,L1,Top"), Gerber("Copper,L2,Bot"), Gerber("Profile,NP")
    x0 = margin
    w_b, h_b = feed_length + length + 2 * margin, width + 2 * margin
    yc = h_b / 2
    xp = x0 + feed_length
    top.region([(xp, yc - width / 2), (xp + length, yc - width / 2), (xp + length, yc + width / 2), (xp, yc + width / 2)], "ANT")
    top.track([(x0, yc), (xp, yc)], feed_width, "ANT")
    top.flash(f"R,{feed_width:.6f}X{feed_width:.6f}", x0, yc, "ANT", ("R1", "1"), "SMDPad")
    bot.region([(0, 0), (w_b, 0), (w_b, h_b), (0, h_b)], "GND")
    prof.outline([(0, 0), (w_b, 0), (w_b, h_b), (0, h_b)])
    _write(out, "patch", {"patch-L1.gtl": top, "patch-L2.gbl": bot, "patch-Edge.gm1": prof}, {}, [0.035, 0.035], [dielectric])
    return {"width_mm": width, "length_mm": length, "feed_length_mm": feed_length, "dielectric_mm": dielectric}


def via_stub_board(out: Path, width=0.3, drill=0.3, pad=0.6, antipad=1.2, dielectric=(0.2, 0.2, 2.8), segment=10.0,
                   margin=5.0) -> dict:
    """A four-layer board: a track (net SIG) on L1 from pad J1.1 to a via, on L3 to a second
    via, and on L1 again to pad J2.1; GND planes on L2 and L4 with antipads. The vias are
    drilled through all layers, so each leaves a stub from L3 to L4."""
    names = ("Copper,L1,Top", "Copper,L2,Inr", "Copper,L3,Inr", "Copper,L4,Bot")
    cu = [Gerber(n) for n in names]
    prof = Gerber("Profile,NP")
    x0, yc = margin, margin
    xs = [x0, x0 + segment, x0 + 2 * segment, x0 + 3 * segment]
    cu[0].track([(xs[0], yc), (xs[1], yc)], width, "SIG")
    cu[2].track([(xs[1], yc), (xs[2], yc)], width, "SIG")
    cu[0].track([(xs[2], yc), (xs[3], yc)], width, "SIG")
    cu[0].flash(f"R,{width:.6f}X{width:.6f}", xs[0], yc, "SIG", ("J1", "1"), "SMDPad")
    cu[0].flash(f"R,{width:.6f}X{width:.6f}", xs[3], yc, "SIG", ("J2", "1"), "SMDPad")
    for xv in (xs[1], xs[2]):
        for k in range(4):
            cu[k].flash(f"C,{pad:.6f}", xv, yc, "SIG", None, "ViaPad")
    w_b, h_b = xs[3] + margin, 2 * margin
    for k in (1, 3):
        cu[k].region([(0, 0), (w_b, 0), (w_b, h_b), (0, h_b)], "GND")
        for xv in (xs[1], xs[2]):
            r = antipad / 2
            n = 32
            cu[k].region([(xv + r * math.cos(2 * math.pi * i / n), yc + r * math.sin(2 * math.pi * i / n)) for i in range(n)], clear=True)
            cu[k].flash(f"C,{pad:.6f}", xv, yc, "SIG", None, "ViaPad")
    prof.outline([(0, 0), (w_b, 0), (w_b, h_b), (0, h_b)])
    files = {f"vs-L{k + 1}.g{k + 1}": g for k, g in enumerate(cu)}
    files["vs-Edge.gm1"] = prof
    _write(out, "vs", files, {"vs-PTH.drl": _excellon([(xs[1], yc, drill), (xs[2], yc, drill)], True, span=(1, 4))},
           [0.035] * 4, list(dielectric))
    return {"stub_mm": dielectric[2] + 0.035, "dielectric": list(dielectric)}


def plane_pair(out: Path, width=40.0, height=30.0, dielectric=(0.1, 0.2, 0.1), parts=(("U1", 12.0, 9.0), ("C1", 30.0, 22.0)),
               drill=0.3, antipad=0.8, pad=0.5, pitch=1.0) -> dict:
    """A four-layer board: GND plane on L2, PWR plane on L3, 0.2 mm apart. Each part has two
    pads on the top layer, x apart by pitch: pin 1 on PWR at (x, y) with a via to L3 that
    passes an antipad in L2, pin 2 on GND at (x + pitch, y) with a via to L2 that passes an
    antipad in L3; both vias are drilled through."""
    cu = [Gerber(n) for n in ("Copper,L1,Top", "Copper,L2,Inr", "Copper,L3,Inr", "Copper,L4,Bot")]
    prof = Gerber("Profile,NP")
    box = [(0, 0), (width, 0), (width, height), (0, height)]
    cu[1].region(box, "GND")
    cu[2].region(box, "PWR")
    holes = []
    r = antipad / 2
    circle = lambda x, y: [(x + r * math.cos(2 * math.pi * i / 32), y + r * math.sin(2 * math.pi * i / 32)) for i in range(32)]
    for ref, x, y in parts:
        cu[0].flash(f"C,{pad:.6f}", x, y, "PWR", (ref, "1"), "SMDPad")
        cu[0].flash(f"C,{pad:.6f}", x + pitch, y, "GND", (ref, "2"), "SMDPad")
        cu[1].region(circle(x, y), clear=True)
        cu[2].region(circle(x + pitch, y), clear=True)
        holes += [(x, y, drill), (x + pitch, y, drill)]
    prof.outline(box)
    files = {f"pp-L{k + 1}.g{k + 1}": g for k, g in enumerate(cu)}
    files["pp-Edge.gm1"] = prof
    _write(out, "pp", files, {"pp-PTH.drl": _excellon(holes, True, span=(1, 4))}, [0.035] * 4, list(dielectric))
    return {"height_mm": dielectric[1], "width_mm": width, "depth_mm": height, "drill_mm": drill}


def via_stack_board(out: Path, dielectric=(0.2, 0.2, 0.4, 0.4, 0.4, 0.4, 0.2), signal_layers=(1, 3), drill=0.3, pad=0.55,
                    antipad=0.9, width=0.3, segment=10.0, margin=5.0, return_vias=4, return_distance=1.0, blind=False) -> dict:
    """A board of thin dielectrics: a track (net SIG) on the first signal layer from pad J1.1
    to a via, on the second to a second via, and back on the first to pad J2.1; every other
    layer a GND plane with antipads, joined by return vias around each signal via. The vias
    are drilled through all layers, so each leaves a stub below the second signal layer;
    with blind, they end at the second signal layer."""
    n = len(dielectric) + 1
    cu = [Gerber(f"Copper,L{k + 1},{'Top' if k == 0 else 'Bot' if k == n - 1 else 'Inr'}") for k in range(n)]
    prof = Gerber("Profile,NP")
    x0, yc = margin, margin
    xs = [x0, x0 + segment, x0 + 2 * segment, x0 + 3 * segment]
    a, b = signal_layers
    cu[a - 1].track([(xs[0], yc), (xs[1], yc)], width, "SIG")
    cu[b - 1].track([(xs[1], yc), (xs[2], yc)], width, "SIG")
    cu[a - 1].track([(xs[2], yc), (xs[3], yc)], width, "SIG")
    cu[a - 1].flash(f"R,{width:.6f}X{width:.6f}", xs[0], yc, "SIG", ("J1", "1"), "SMDPad")
    cu[a - 1].flash(f"R,{width:.6f}X{width:.6f}", xs[3], yc, "SIG", ("J2", "1"), "SMDPad")
    w_b, h_b = xs[3] + margin, 2 * margin
    returns = []
    for xv in (xs[1], xs[2]):
        for k in range(n):
            cu[k].flash(f"C,{pad:.6f}", xv, yc, "SIG", None, "ViaPad")
        for i in range(return_vias):
            t = 2 * math.pi * (i + 0.5) / return_vias
            returns.append((xv + return_distance * math.cos(t), yc + return_distance * math.sin(t), drill))
    for k in range(n):
        if k + 1 in signal_layers:
            continue
        cu[k].region([(0, 0), (w_b, 0), (w_b, h_b), (0, h_b)], "GND")
        for xv in (xs[1], xs[2]):
            r = antipad / 2
            cu[k].region([(xv + r * math.cos(2 * math.pi * i / 32), yc + r * math.sin(2 * math.pi * i / 32)) for i in range(32)], clear=True)
            cu[k].flash(f"C,{pad:.6f}", xv, yc, "SIG", None, "ViaPad")
    prof.outline([(0, 0), (w_b, 0), (w_b, h_b), (0, h_b)])
    files = {f"vk-L{k + 1}.g{k + 1}": g for k, g in enumerate(cu)}
    files["vk-Edge.gm1"] = prof
    drills = {"vk-PTH.drl": _excellon(returns + ([] if blind else [(xs[1], yc, drill), (xs[2], yc, drill)]), True, span=(1, n))}
    if blind:
        drills["vk-PTH-1-%d.drl" % b] = _excellon([(xs[1], yc, drill), (xs[2], yc, drill)], True, span=(1, b))
    _write(out, "vk", files, drills, [0.035] * n, list(dielectric))
    stub = sum(dielectric[b - 1:]) + 0.035 * (n - b)
    return {"stub_mm": stub, "layers": n}


def mounting_board(out: Path, height=0.5, spacing=1.0, drill=0.3, pad=0.6, trace=0.0, trace_width=0.3, margin=10.0) -> dict:
    """A part C1 on the top layer over a ground plane on the bottom layer: two pads (C1.1,
    C1.2) spacing apart in x, each joined by a track of length trace (outwards, along y) to
    a via down to the plane; all of net GND, so the loop of the part, the pads, the tracks
    and the vias closes through the plane."""
    top, bot, prof = Gerber("Copper,L1,Top"), Gerber("Copper,L2,Bot"), Gerber("Profile,NP")
    x0, y0 = margin, margin
    holes = []
    for k, x in enumerate((x0, x0 + spacing)):
        top.flash(f"R,{pad:.6f}X{pad:.6f}", x, y0, "GND", ("C1", str(k + 1)), "SMDPad")
        if trace > 0:
            top.track([(x, y0), (x, y0 - trace)], trace_width, "GND")
            top.flash(f"C,{max(pad, drill + 0.2):.6f}", x, y0 - trace, "GND", None, "ViaPad")
        holes.append((x, y0 - trace, drill))
    w, h = x0 + spacing + margin, y0 + margin
    bot.region([(0, 0), (w, 0), (w, h), (0, h)], "GND")
    prof.outline([(0, 0), (w, 0), (w, h), (0, h)])
    _write(out, "mt", {"mt-L1.gtl": top, "mt-L2.gbl": bot, "mt-Edge.gm1": prof}, {"mt-PTH.drl": _excellon(holes, True, span=(1, 2))},
           [0.035, 0.035], [height])
    return {"height_mm": height + 0.035, "spacing_mm": spacing}


def return_path_board(out: Path, kind="solid", length=30.0, width=0.3, dielectric=0.2) -> dict:
    """A trace SIG on the top layer from J1.1 to J2.1 over a plane on the bottom layer:
    solid; with a slot 1 mm wide across the trace's middle (slot); split at the middle
    into GND and PWR by a 0.5 mm gap (split); with a row of overlapping antipads across the
    trace (antipads); or with one antipad under the trace's edge (notch)."""
    top, bot, prof = Gerber("Copper,L1,Top"), Gerber("Copper,L2,Bot"), Gerber("Profile,NP")
    W, H = length + 10.0, 20.0
    x0, x1, yc, xm = 5.0, 5.0 + length, H / 2, 5.0 + length / 2
    top.track([(x0, yc), (x1, yc)], width, "SIG")
    top.flash("R,1.000000X1.000000", x0, yc, "SIG", ("J1", "1"), "SMDPad")
    top.flash("R,1.000000X1.000000", x1, yc, "SIG", ("J2", "1"), "SMDPad")
    m = 1.0
    if kind == "split":
        bot.region([(m, m), (xm - 0.25, m), (xm - 0.25, H - m), (m, H - m)], "GND")
        bot.region([(xm + 0.25, m), (W - m, m), (W - m, H - m), (xm + 0.25, H - m)], "PWR")
    else:
        bot.region([(m, m), (W - m, m), (W - m, H - m), (m, H - m)], "GND")
        if kind == "slot":
            bot.region([(xm - 0.5, yc - 5), (xm + 0.5, yc - 5), (xm + 0.5, yc + 5), (xm - 0.5, yc + 5)], None, clear=True)
        elif kind == "antipads":
            for k in range(-4, 5):
                bot.flash("C,1.000000", xm, yc + 0.8 * k, None, None, "AntiPad", clear=True)
        elif kind == "notch":
            bot.flash("C,0.300000", xm, yc + width / 2, None, None, "AntiPad", clear=True)
    prof.outline([(0, 0), (W, 0), (W, H), (0, H)])
    _write(out, "rp", {"rp-L1.gtl": top, "rp-L2.gbl": bot, "rp-Edge.gm1": prof}, {}, [0.035, 0.035], [dielectric])
    return {"width_mm": width, "x_mid_mm": xm, "y_mm": yc}


def switch_node_board(out: Path, patch=(10.0, 8.0), shield=False, dielectric=1.5) -> dict:
    """A copper patch of net SW on the top layer of a two-layer board; with shield, a GND
    plane on the bottom layer under the left half of the patch."""
    top, bot, prof = Gerber("Copper,L1,Top"), Gerber("Copper,L2,Bot"), Gerber("Profile,NP")
    W, H = 40.0, 30.0
    a, b = patch
    x0, y0 = (W - a) / 2, (H - b) / 2
    top.region([(x0, y0), (x0 + a, y0), (x0 + a, y0 + b), (x0, y0 + b)], "SW")
    top.flash("R,1.000000X1.000000", 3.0, 3.0, "GND", ("TP1", "1"), "SMDPad")
    if shield:
        bot.region([(1.0, 1.0), (W / 2, 1.0), (W / 2, H - 1.0), (1.0, H - 1.0)], "GND")
    prof.outline([(0, 0), (W, 0), (W, H), (0, H)])
    _write(out, "sw", {"sw-L1.gtl": top, "sw-L2.gbl": bot, "sw-Edge.gm1": prof}, {}, [0.035, 0.035], [dielectric])
    return {"area_m2": a * b * 1e-6, "half_area_m2": a * b * 1e-6 / 2, "dielectric_m": dielectric * 1e-3}


def node_board(out: Path, outline, size=60.0, dielectric=1.5, copper=0.035, epsilon_r="4.5", mask_mm=0.01) -> dict:
    """A two-layer board size x size mm with one copper region of net SW on the top layer
    (outline in mm, relative to the board's centre) and nothing else."""
    top, bot, prof = Gerber("Copper,L1,Top"), Gerber("Copper,L2,Bot"), Gerber("Profile,NP")
    c = size / 2
    top.region([(c + x, c + y) for x, y in outline], "SW")
    top.flash("R,0.500000X0.500000", 1.0, 1.0, "GND", ("TP1", "1"), "SMDPad")
    prof.outline([(0, 0), (size, 0), (size, size), (0, size)])
    _write(out, "node", {"node-L1.gtl": top, "node-L2.gbl": bot, "node-Edge.gm1": prof}, {}, [copper, copper], [dielectric],
           epsilon_r=epsilon_r, mask_mm=mask_mm)
    return {"height_m": (dielectric + copper * 1.5) * 1e-3}


def slot_board(out: Path, slot_end=8.0, slot_width=1.0, dielectric=0.5, trace_width=0.5, stitch=True, size=(40.0, 30.0)) -> dict:
    """A trace of net SIG on the top layer from J1 to R1 along y = 20 mm over a GND plane on
    the bottom layer. Each part's second pad, 2 mm below its first, is GND with a via to the
    plane. The plane has a slot across the trace's path, slot_width wide at x = 20 mm, from
    the top edge down to y = slot_end, so the return current detours below it. With stitch,
    a part C1 across the slot beside the trace, its pads 1 mm either side of the slot's
    centre with vias to the plane."""
    top, bot, prof = Gerber("Copper,L1,Top"), Gerber("Copper,L2,Bot"), Gerber("Profile,NP")
    W, H = size
    yc, xs = 20.0, W / 2
    holes = []
    pad, drill = 0.6, 0.3
    top.flash(f"R,{pad:.6f}X{pad:.6f}", 5.0, yc, "SIG", ("J1", "1"), "SMDPad")
    top.flash(f"R,{pad:.6f}X{pad:.6f}", W - 5.0, yc, "SIG", ("R1", "1"), "SMDPad")
    top.track([(5.0, yc), (W - 5.0, yc)], trace_width, "SIG")
    parts = [("J1", 5.0), ("R1", W - 5.0)] + ([("C1", xs - 1.0), ("C1", xs + 1.0)] if stitch else [])
    for k, (ref, x) in enumerate(parts):
        pin = "2" if ref != "C1" else ("1" if x < xs else "2")
        top.flash(f"R,{pad:.6f}X{pad:.6f}", x, yc - 2.0, "GND", (ref, pin), "SMDPad")
        holes.append((x, yc - 2.0, drill))
    bot.region([(0, 0), (W, 0), (W, H), (0, H)], "GND")
    if slot_end < H:
        bot.region([(xs - slot_width / 2, slot_end), (xs + slot_width / 2, slot_end), (xs + slot_width / 2, H + 1), (xs - slot_width / 2, H + 1)],
                   clear=True)
    prof.outline([(0, 0), (W, 0), (W, H), (0, H)])
    _write(out, "slot", {"slot-L1.gtl": top, "slot-L2.gbl": bot, "slot-Edge.gm1": prof}, {"slot-PTH.drl": _excellon(holes, True, span=(1, 2))},
           [0.035, 0.035], [dielectric])
    return {"height_mm": dielectric + 0.035}
