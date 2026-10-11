#!/usr/bin/env python3
"""Full-wave S-parameters of a channel's structure exported by KairoQ, with openEMS (FDTD).

Reads model.json from the directory of this script, builds the structure with CSXCAD and
runs openEMS once per excited port. Writes result.sNp (Touchstone, the ports' resistance as
reference) and result.json; with far-field frequencies in the model, the radiation of
port 1: directivity, radiation efficiency, gain and realised gain.

Requires openEMS with its Python interface (modules CSXCAD and openEMS), numpy and h5py.

usage: python simulate.py [--threads N] [--setup-only] [--sim-dir DIR] [--ports 1,2,...]
"""
import argparse
import json
import os

import numpy as np

UNIT = 1e-3  # drawing unit of the simulation: millimetre
EPS0 = 8.8541878128e-12


def djordjevic_sarkar(eps_r, tand, f_ref, f, w1=1e4, w2=1e12):
    """Permittivity eps' - j eps'' of the Djordjevic-Sarkar model through eps_r and tand at
    f_ref (as in kairoq's cross-sections)."""
    w0, w = 2 * np.pi * f_ref, 2 * np.pi * np.asarray(f, dtype=float)
    span = np.log(w2 / w1)
    d = eps_r * tand * span / (np.arctan(w0 / w1) - np.arctan(w0 / w2))
    einf = eps_r - d / span * 0.5 * np.log((w2 ** 2 + w0 ** 2) / (w1 ** 2 + w0 ** 2))
    return einf + d / span * np.log((w2 + 1j * w) / (w1 + 1j * w))


def nnls(a, b, iterations=200):
    """Least squares with non-negative coefficients (Lawson and Hanson)."""
    n = a.shape[1]
    x, passive = np.zeros(n), np.zeros(n, bool)
    for _ in range(iterations):
        w = a.T @ (b - a @ x)
        if passive.all() or w[~passive].max() <= 1e-12 * max(1.0, np.abs(w).max()):
            break
        passive[np.argmax(np.where(passive, -np.inf, w))] = True
        while True:
            z = np.zeros(n)
            z[passive] = np.linalg.lstsq(a[:, passive], b, rcond=None)[0]
            if (z[passive] > 0).all():
                x = z
                break
            neg = passive & (z <= 0)
            alpha = np.min(x[neg] / (x[neg] - z[neg]))
            x = x + alpha * (z - x)
            passive &= x > 1e-14
    return x


def debye_poles(eps_r, tand, f_ref, f_lo, f_hi, poles=6):
    """Debye poles (eps_inf, [(delta, tau)]) fitted to the Djordjevic-Sarkar model over
    [f_lo, f_hi], relaxation frequencies spread over a little more; returns the fit's
    largest relative error too."""
    f = np.logspace(np.log10(f_lo), np.log10(f_hi), 80)
    target = djordjevic_sarkar(eps_r, tand, f_ref, f)
    fk = np.logspace(np.log10(f_lo / 3), np.log10(f_hi * 3), poles)
    tau = 1 / (2 * np.pi * fk)
    cols = [np.ones_like(f, dtype=complex)] + [1 / (1 + 1j * 2 * np.pi * f * t) for t in tau]
    a = np.array(cols).T
    c = nnls(np.vstack([a.real, a.imag]), np.concatenate([target.real, target.imag]))
    fit = a @ c
    return c[0], [(c[k + 1], tau[k]) for k in range(poles) if c[k + 1] > 0], float(np.max(np.abs(fit - target) / np.abs(target)))


def mesh_lines(fixed, exact, lo, hi, min_cell, board_cell, max_cell, margin, smooth):
    """Lines through the exact coordinates, and through the other fixed coordinates at
    least a minimum cell from those already taken; smoothed to board_cell within
    [lo, hi] and to max_cell in the air margin beyond it."""
    taken = sorted(set(exact))
    for p in sorted(set([lo, hi] + [p for p in fixed if lo <= p <= hi])):
        i = np.searchsorted(taken, p)
        near = [taken[j] for j in (i - 1, i) if 0 <= j < len(taken)]
        if all(abs(p - q) >= 0.999 * min_cell for q in near):
            taken.insert(i, p)
    inner = smooth(taken, board_cell, 1.4)
    lines = np.unique(np.asarray(smooth([lo - margin] + list(inner) + [hi + margin], max_cell, 1.4)))
    # Smoothing moves lines by rounding; a sheet or port off its line by 1e-16 is dropped
    # from the grid, so lines near an exact coordinate take it exactly.
    ex = np.unique(np.asarray(exact, dtype=float))
    if len(ex):
        k = np.clip(np.searchsorted(ex, lines), 1, len(ex) - 1)
        near = np.where(np.abs(ex[k - 1] - lines) < np.abs(ex[k] - lines), ex[k - 1], ex[k])
        lines = np.where(np.abs(near - lines) < 1e-6 * min_cell, near, lines)
    return np.unique(lines)


def openems_executable():
    """The openEMS program beside the library the Python module loads."""
    import openEMS.openEMS  # noqa: F401
    with open("/proc/self/maps") as f:
        for line in f:
            if "libopenEMS" in line:
                exe = os.path.join(os.path.dirname(line.split()[-1]), "..", "bin", "openEMS")
                if os.path.exists(exe):
                    return os.path.abspath(exe)
    import shutil
    return shutil.which("openEMS")


def run_debye(csx, model, sim_dir, threads, poles):
    """Run openEMS from an XML file in which the dielectrics are Debye materials (the Python
    interface has no Debye material)."""
    import re
    import shutil
    import subprocess
    s = model["simulation"]
    shutil.rmtree(sim_dir, ignore_errors=True)
    os.makedirs(sim_dir)
    geometry = os.path.join(sim_dir, "structure.xml")
    csx.Write2XML(geometry)
    text = open(geometry).read()
    for k, (einf, terms) in enumerate(poles):
        start = text.index('<Material ID="%s"' % re.search(r'<Material ID="(\d+)" Name="dielectric_%d"' % k, text).group(1))
        end = text.index("</Material>", start) + len("</Material>")
        block = text[start:end].replace("<Material ", "<DebyeMaterial ", 1).replace("</Material>", "</DebyeMaterial>")
        extra = "".join(' EpsilonDelta_{0}="{1:.9g},{1:.9g},{1:.9g}" EpsilonRelaxTime_{0}="{2:.9g},{2:.9g},{2:.9g}"'.format(i + 1, float(d), float(t))
                        for i, (d, t) in enumerate(terms))
        block = re.sub(r'(<Property )Epsilon="[^"]*"', r'\1Epsilon="{0:.9g},{0:.9g},{0:.9g}"{1}'.format(float(einf), extra), block, count=1)
        block = re.sub(r'(<Property [^>]*)Kappa="[^"]*"', r'\1Kappa="0,0,0"', block, count=1)
        text = text[:start] + block + text[end:]
    body = text[text.index("<ContinuousStructure"):]
    f_max = s["f_max_hz"]
    bc = " ".join('%s="3" PML_%s="8"' % (n, n) for n in ("xmin", "xmax", "ymin", "ymax", "zmin", "zmax"))
    xml = os.path.join(sim_dir, "simulation.xml")
    with open(xml, "w") as f:
        f.write('<?xml version="1.0" encoding="UTF-8" standalone="yes" ?>\n<openEMS>\n')
        f.write('  <FDTD NumberOfTimesteps="%d" endCriteria="%.9g" f_max="%.9g">\n' % (int(s["max_timesteps"]), 10 ** (s["end_criterion_db"] / 10), f_max))
        f.write('    <Excitation Type="0" f0="%.9g" fc="%.9g" />\n    <BoundaryCond %s />\n  </FDTD>\n' % (f_max / 2, f_max / 2, bc))
        f.write(body + "\n</openEMS>\n")
    cmd = [openems_executable(), "simulation.xml"] + (["--numThreads=%d" % threads] if threads else [])
    subprocess.run(cmd, cwd=sim_dir, check=True)


def move_to_pad(port, sim_dir, freq):
    """A microstrip port's voltage and current at the pad: from its measurement plane back
    along the line by its measured propagation constant, attenuation included (openEMS's
    own shift turns the phase only)."""
    port.CalcPort(sim_dir, freq)
    gamma = 1j * port.beta  # beta: the measured phase constant, its imaginary part the loss
    d = port.measplane_shift * port.CSX.GetGrid().GetDeltaUnit()
    z = port.Z_ref
    v, i = port.uf_tot, port.if_tot
    port.uf_tot = v * np.cosh(gamma * d) + z * i * np.sinh(gamma * d)
    port.if_tot = i * np.cosh(gamma * d) + v / z * np.sinh(gamma * d)


def build(model, excite, far_field, copper="pec", port_kind="auto"):
    """The structure, mesh and ports with port `excite` (0-based) driven."""
    from CSXCAD import ContinuousStructure
    from CSXCAD.SmoothMeshLines import SmoothMeshLines
    from openEMS import openEMS

    s = model["simulation"]
    d = lambda v: np.asarray(v, dtype=float) / UNIT
    hmin, hboard, hmax, margin = (s[k] / UNIT for k in ("min_cell", "board_cell", "max_cell", "air_margin"))
    f_max = s["f_max_hz"]
    fdtd = openEMS(NrTS=int(s["max_timesteps"]), EndCriteria=10 ** (s["end_criterion_db"] / 10))
    fdtd.SetGaussExcite(f_max / 2, f_max / 2)
    fdtd.SetBoundaryCond(["PML_8"] * 6)
    csx = ContinuousStructure()
    fdtd.SetCSX(csx)
    grid = csx.GetGrid()
    grid.SetDeltaUnit(UNIT)
    fx, fy, fz = [], [], []
    exact = ([], [], [])

    for k, slab in enumerate(model["dielectric"]):
        kappa = 2 * np.pi * (f_max / 2) * EPS0 * slab["epsilon_r"] * slab["loss_tangent"]
        mat = csx.AddMaterial("dielectric_%d" % k, epsilon=slab["epsilon_r"], kappa=kappa)
        z0, z1 = float(d(slab["z0"])), float(d(slab["z1"]))
        for poly in model["window"]:
            mat.AddLinPoly(d(poly).T, "z", z0, z1 - z0, priority=1)
        fz += [z0, z1]
    for poly in model["window"]:
        p = d(poly)
        fx += list(p[:, 0])
        fy += list(p[:, 1])

    # Copper: perfect conductor, or sheets of each layer's conductivity and thickness.
    sheet_props = {}

    def copper_of(sheet):
        if copper != "sheet":
            return sheet_props.setdefault("pec", csx.AddMetal("copper"))
        key = (sheet.get("conductivity", 5.8e7), sheet.get("thickness_m", 35e-6))
        if key not in sheet_props:
            sheet_props[key] = csx.AddConductingSheet("copper_%d" % len(sheet_props), conductivity=key[0], thickness=key[1])
        return sheet_props[key]
    barrels = csx.AddMetal("barrels")
    for sheet in model["sheets"]:
        metal = copper_of(sheet)
        z = float(d(sheet["z"]))
        exact[2].append(z)
        fz += [z - hmin, z + hmin]
        signal = sheet.get("signal", [True] * len(sheet["polygons"]))
        for poly, sig in zip(sheet["polygons"], signal):
            p = d(poly)
            metal.AddPolygon(p.T, "z", z, priority=10)
            # Lines where they resolve something: on the edges of the ports' nets' copper that
            # run along x or y, with a minimum cell to either side, and at the ends of its other
            # edges a board cell long or more (arcs' short chords resolve nothing); on other copper only along its long axis-parallel edges. Lines at
            # every vertex of a real board's pours and arcs would multiply the cells.
            q = np.roll(p, -1, axis=0)
            for a, b in zip(p, q):
                length = float(np.hypot(*(b - a)))
                for ax, lst in ((0, fx), (1, fy)):
                    along = abs(b[ax] - a[ax]) < 1e-6  # the edge runs along the other axis
                    if sig and along and length >= hmin:
                        lst += [a[ax], a[ax] - hmin, a[ax] + hmin]
                    elif sig and not along and length >= hboard:
                        lst += [a[ax], b[ax]]
                    elif along and length >= hboard:
                        lst.append(a[ax])
    for v in model["vias"]:
        x, y, r = float(d(v["x"])), float(d(v["y"])), float(d(v["radius"]))
        barrels.AddCylinder([x, y, float(d(v["z0"]))], [x, y, float(d(v["z1"]))], r, priority=10)
        fx += [x - r, x + r]
        fy += [y - r, y + r]
        # The channel's vias: lines that trace their pads' and antipads' circles (at 0, 22.5,
        # 45 and 67.5 degrees, so that the staircase stays within a few percent of the radius)
        # and halfway between, so that each ring of dielectric spans at least two cells.
        radii = sorted({r} | {float(d(q)) for q in v.get("rings", [])})
        for k, q in enumerate(radii):
            for c in ((1.0, 0.92388, 0.70711, 0.38268) if v.get("rings") else (1.0,)):
                fx += [x - q * c, x + q * c]
                fy += [y - q * c, y + q * c]
            if k:
                m = 0.5 * (q + radii[k - 1])
                fx += [x - m, x + m]
                fy += [y - m, y + m]

    # Ports: where a straight track leaves the pad along x or y, a microstrip port over its
    # first part (fed at the pad, measured along it, the reference plane moved back to the
    # pad, terminated in the port's resistance); elsewhere a vertical box from the reference
    # to the pad, square, as wide as the track, at the pad's centre (its current probe is a
    # loop around it at mid-height: a port of no extent across has a degenerate loop).
    specs = []
    for k, p in enumerate(model["ports"]):
        x, y, w = float(d(p["x"])), float(d(p["y"])), float(d(p["width"]))
        zp, zr = float(d(p["z_pad"])), float(d(p["z_ref"]))
        line = p.get("line")
        length = min(float(d(line["length"])), max(2.0, 10 * w)) if line else 0.0
        if port_kind == "auto" and line and length >= 8 * hmin:
            ax, sgn = line["axis"], line["dir"]
            start, stop = [x, y, zp], [x, y, zr]
            start[1 - ax] -= w / 2
            stop[1 - ax] += w / 2
            stop[ax] += sgn * length
            exact[1 - ax].extend([start[1 - ax], stop[1 - ax]])
            exact[ax].extend([x, x + sgn * length] if ax == 0 else [y, y + sgn * length])
            exact[2].extend([zr, zp])
            specs.append(("msl", k, start, stop, "xy"[ax], length, p["resistance_ohm"]))
        else:
            start, stop = [x - w / 2, y - w / 2, zr], [x + w / 2, y + w / 2, zp]
            exact[0].extend([x - w / 2, x + w / 2])
            exact[1].extend([y - w / 2, y + w / 2])
            exact[2].extend([zr, zp])
            specs.append(("lumped", k, start, stop, None, 0.0, p["resistance_ohm"]))

    lo = [min(fx), min(fy), min(fz)]
    hi = [max(fx), max(fy), max(fz)]
    for ax, fixed in enumerate((fx, fy, fz)):
        grid.SetLines("xyz"[ax], mesh_lines(fixed, exact[ax], lo[ax], hi[ax], hmin, hboard, hmax, margin, SmoothMeshLines))
    ports = []
    port_metal = None
    for kind, k, start, stop, prop, length, r in specs:
        ex = 1.0 if k == excite else 0.0
        if kind == "msl":
            port_metal = port_metal or csx.AddMetal("port_strips")
            port = fdtd.AddMSLPort(k + 1, port_metal, start, stop, prop, "z", excite=ex, MeasPlaneShift=length / 2, Feed_R=r, priority=5)
            port.kind = "msl"
        else:
            port = fdtd.AddLumpedPort(k + 1, r, start, stop, "z", excite=ex, priority=5)
            port.kind = "lumped"
        ports.append(port)
    nf2ff = None
    if far_field:
        lines = [grid.GetLines(ax) for ax in range(3)]
        nf2ff = fdtd.CreateNF2FFBox(start=[lines[ax][10] for ax in range(3)], stop=[lines[ax][-11] for ax in range(3)])
    return fdtd, ports, nf2ff, grid, (lo, hi)


def touchstone(path, freq, s, z0, names):
    n = s.shape[1]
    with open(path, "w") as f:
        f.write("! kairoq openEMS result, ports: %s\n# Hz S RI R %g\n" % (" ".join(names), z0))
        for k, fr in enumerate(freq):
            m = s[k]
            vals = [m[0, 0], m[1, 0], m[0, 1], m[1, 1]] if n == 2 else [m[i, j] for i in range(n) for j in range(n)]
            line = "%.9g" % fr
            for c, v in enumerate(vals):
                line += " %.12g %.12g" % (v.real, v.imag)
                if n > 2 and (c + 1) % 4 == 0 and c + 1 < len(vals):
                    line += "\n"
            f.write(line + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--threads", type=int, default=0, help="threads of openEMS (0: all)")
    ap.add_argument("--setup-only", action="store_true", help="build the structure and mesh of port 1, do not run")
    ap.add_argument("--sim-dir", default=None, help="directory of the simulation files (default: sim/ next to the model)")
    ap.add_argument("--ports", default=None, help="ports to excite, from 1 (default: all); others' columns follow by reciprocity")
    ap.add_argument("--port-kind", choices=("auto", "lumped"), default="auto",
                    help="auto: microstrip ports where a straight track leaves the pad along x or y, lumped boxes elsewhere")
    ap.add_argument("--losses", choices=("model", "simple"), default="model",
                    help="model: Debye dielectrics fitted to the dielectric model and copper sheets of their conductivity; "
                         "simple: a conductivity for the loss tangent at the band's middle, perfect copper")
    args = ap.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "model.json")) as f:
        model = json.load(f)
    s = model["simulation"]
    sim_root = os.path.abspath(args.sim_dir or os.path.join(here, "sim"))
    n = len(model["ports"])
    z0 = model["ports"][0]["resistance_ohm"]
    if any(abs(p["resistance_ohm"] - z0) > 1e-9 for p in model["ports"]):
        raise SystemExit("the ports need one resistance for one Touchstone reference")
    excite = sorted({int(p) - 1 for p in args.ports.split(",")}) if args.ports else list(range(n))
    f_max = s["f_max_hz"]
    freq = np.linspace(f_max / s["points"], f_max, s["points"])
    far = np.asarray(s["far_field_hz"], dtype=float)
    S = np.full((len(freq), n, n), np.nan, dtype=complex)
    out = {"channel": model["channel"], "ports": [p["name"] for p in model["ports"]], "frequency_hz": freq.tolist()}
    for j in excite:
        fdtd, ports, nf2ff, grid, (lo, hi) = build(model, j, len(far) > 0 and j == 0, "sheet" if args.losses == "model" else "pec",
                                                   args.port_kind)
        lines = [grid.GetLines(ax) for ax in range(3)]
        print("port %d: mesh %d x %d x %d lines, %.3g million cells, smallest cell %.3g mm" % (
            j + 1, len(lines[0]), len(lines[1]), len(lines[2]), np.prod([len(l) - 1 for l in lines]) / 1e6,
            min(np.diff(l).min() for l in lines)), flush=True)
        sim_dir = os.path.join(sim_root, "port%d" % (j + 1))
        os.makedirs(sim_dir, exist_ok=True)
        if args.losses == "model" and not args.setup_only:
            poles = []
            for slab in model["dielectric"]:
                einf, terms, err = debye_poles(slab["epsilon_r"], slab["loss_tangent"], f_max / 2, f_max / 200, f_max)
                poles.append((einf, terms))
            run_debye(fdtd.GetCSX() if hasattr(fdtd, "GetCSX") else None, model, sim_dir, args.threads, poles)
        else:
            fdtd.Run(sim_dir, cleanup=True, setup_only=args.setup_only, numThreads=args.threads)
        if args.setup_only:
            return
        # Waves against the reference resistance from each port's voltage and current at its
        # reference plane (a microstrip port's moved back to the pad along its line).
        for p in ports:
            if p.kind == "msl":
                move_to_pad(p, sim_dir, freq)
            else:
                p.CalcPort(sim_dir, freq, ref_impedance=z0)
        a_j = 0.5 * (ports[j].uf_tot + z0 * ports[j].if_tot)
        for i, p in enumerate(ports):
            S[:, i, j] = 0.5 * (p.uf_tot - z0 * p.if_tot) / a_j
        if nf2ff is not None:
            p = ports[0]
            if p.kind == "msl":
                move_to_pad(p, sim_dir, far)
            else:
                p.CalcPort(sim_dir, far, ref_impedance=z0)
            theta, phi = np.arange(0.0, 181.0, 2.0), np.arange(0.0, 360.0, 2.0)
            centre = [0.5 * (lo[ax] + hi[ax]) * UNIT for ax in range(3)]
            res = nf2ff.CalcNF2FF(sim_dir, far, theta, phi, center=centre)
            rows = []
            for k, f in enumerate(far):
                p_acc = float(np.real(0.5 * p.uf_tot[k] * np.conj(p.if_tot[k])))
                p_inc = float(0.5 * abs(0.5 * (p.uf_tot[k] + z0 * p.if_tot[k])) ** 2 / z0)
                d = float(res.Dmax[k])
                eff = float(res.Prad[k]) / p_acc
                rows.append({"frequency_hz": float(f), "directivity_dbi": 10 * np.log10(d), "radiation_efficiency": eff,
                             "gain_dbi": 10 * np.log10(d * eff), "realised_gain_dbi": 10 * np.log10(d * float(res.Prad[k]) / p_inc),
                             "radiated_power_w": float(res.Prad[k]), "accepted_power_w": p_acc})
                print("far field %.4g GHz: directivity %.2f dBi, efficiency %.3f, gain %.2f dBi, realised gain %.2f dBi" % (
                    f * 1e-9, rows[-1]["directivity_dbi"], eff, rows[-1]["gain_dbi"], rows[-1]["realised_gain_dbi"]))
            out["far_field"] = rows
    # Columns not excited follow by reciprocity where their rows are known.
    for j in range(n):
        if j in excite:
            continue
        for i in range(n):
            if i in excite:
                S[:, i, j] = S[:, j, i]
    out["excited"] = [j + 1 for j in excite]
    out["s"] = [[[v.real, v.imag] for v in m.ravel()] for m in S]
    with open(os.path.join(here, "result.json"), "w") as f:
        json.dump(out, f)
    if not np.isnan(S).any():
        touchstone(os.path.join(here, "result.s%dp" % n), freq, S, z0, out["ports"])
        print("wrote result.s%dp" % n)
    else:
        print("not every column of S was computed; result.json holds those that were")


if __name__ == "__main__":
    main()
