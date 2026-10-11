"""RF structures: the patch estimate against Balanis's worked example and the directivity
of a short slot, the openEMS export of a channel, and, with openEMS, the full-wave
S-parameters of a microstrip against the native channel and a patch against its estimate.
"""

from __future__ import annotations

import json
import math
import os
import subprocess

import numpy as np
import pytest

import kairoq
import synthetic
from conftest import cli_path

OPENEMS_PYTHON = os.environ.get("KAIROQ_OPENEMS_PYTHON")
fullwave = pytest.mark.skipif(not OPENEMS_PYTHON, reason="KAIROQ_OPENEMS_PYTHON names no Python with openEMS")


def test_patch_balanis(validation):
    """Balanis, Antenna Theory, chapter 14, worked examples: epsilon_r 2.2, h 0.1588 cm,
    10 GHz: W 1.186 cm, epsilon_eff 1.972, delta L 0.081 cm, L 0.906 cm; input resistance
    at the edge 228.3 ohm with the mutual conductance; 50 ohm inset 0.3126 cm."""
    r = kairoq.patch_antenna(2.2, 1.588, 10e9)
    case = "patch, epsilon_r 2.2, h 1.588 mm, 10 GHz"
    validation.check(case, "width", r["width_m"] * 100, 1.186, rel=2e-3, unit="cm", source="Balanis")
    validation.check(case, "effective permittivity", r["epsilon_eff"], 1.972, rel=1e-3, source="Balanis")
    validation.check(case, "fringing extension", r["delta_l_m"] * 100, 0.081, rel=0.01, unit="cm", source="Balanis")
    validation.check(case, "length", r["length_m"] * 100, 0.906, rel=2e-3, unit="cm", source="Balanis")
    validation.check(case, "edge resistance", r["edge_resistance_ohm"], 228.3, rel=5e-3, unit="ohm", source="Balanis")
    validation.check(case, "inset for 50 ohm", r["inset_m"] * 100, 0.3126, rel=5e-3, unit="cm", source="Balanis")
    assert r["resonant_hz"] == pytest.approx(10e9, rel=1e-9)


def test_short_slot_directivity(validation):
    """A short slot on an infinite ground plane radiates as a magnetic dipole and its image:
    directivity 3."""
    validation.check("slot of 1e-4 wavelengths", "directivity", kairoq._core.two_slot_directivity(1e-4, 1e-4), 3.0, rel=1e-4,
                     source="short magnetic dipole over a ground plane")


def microstrip_model(tmp_path):
    synthetic.microstrip(tmp_path / "board", width=0.35, dielectric=0.2, margin=4.0, pad=0.35)
    b = kairoq.open(tmp_path / "board")
    si = {"channels": [{"name": "C", "ports": ["R1.2", "R2.1"]}], "network": {"f_max_hz": 20e9, "points": 100},
          "openems": {"f_max_hz": 20e9, "points": 100, "margin_mm": 2.0}}
    return b, si


def test_channel_export(tmp_path):
    b, si = microstrip_model(tmp_path)
    m = kairoq.export_openems_channel(b, si, tmp_path / "export")
    assert (tmp_path / "export" / "simulate.py").read_text().startswith("#!/usr/bin/env python3")
    # Sheets on the faces towards the board: the strip 0.2 mm above the plane.
    p = m["ports"]
    assert len(p) == 2 and all(q["z_pad"] - q["z_ref"] == pytest.approx(0.2e-3) for q in p)
    assert all(q["width"] == pytest.approx(0.35e-3) for q in p)
    # The outer copper layers are air: one slab, the dielectric between the sheets.
    assert len(m["dielectric"]) == 1 and m["dielectric"][0]["z1"] - m["dielectric"][0]["z0"] == pytest.approx(0.2e-3)
    # Polygons without collinear vertices: the strip with its pads is a rectangle.
    top = next(s for s in m["sheets"] if s["layer"] == 1)
    assert sorted(len(q) for q in top["polygons"]) == [4, 4, 4]
    with pytest.raises(Exception, match="no channel"):
        kairoq.export_openems_channel(b, si, None, channel="X")


@pytest.mark.skipif(cli_path() is None, reason="kairoq command not built")
def test_cli(tmp_path):
    b, si = microstrip_model(tmp_path)
    (tmp_path / "project.yaml").write_text("board: [board]\nsi:\n" + "".join(f"  {k}: {json.dumps(v)}\n" for k, v in si.items()))
    r = subprocess.run([cli_path(), "openems", str(tmp_path / "project.yaml"), "--channel", "C", "-o", str(tmp_path / "out")],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert "port         R1.2" in r.stdout and (tmp_path / "out" / "model.json").exists()


def run(export_dir, log, *args):
    with open(log, "w") as f:
        subprocess.run([OPENEMS_PYTHON, str(export_dir / "simulate.py"), *args], stdout=f, stderr=subprocess.STDOUT, check=True)


@fullwave
def test_fullwave_microstrip(tmp_path, validation):
    """A 50 ohm microstrip of 20 mm on 0.2 mm, ports at its ends: the full-wave S-parameters
    against the native channel. openEMS has the dielectric model as Debye poles and copper
    sheets of the layer's conductivity; its ports are boxes as wide as the strip."""
    b, si = microstrip_model(tmp_path)
    kairoq.export_openems_channel(b, si, tmp_path / "export")
    run(tmp_path / "export", tmp_path / "openems.log")
    o = kairoq.read_touchstone(tmp_path / "export" / "result.s2p")
    n = kairoq.analyse_channels(b, si)[0]
    f, so = np.array(o["frequency_hz"]), np.array(o["s"]).reshape(-1, 2, 2)
    fn, sn = np.array(n["frequency_hz"]), np.array(n["s"]).reshape(-1, 2, 2)
    assert np.allclose(f, fn)
    band = (f >= 1e9) & (f <= 10e9)
    delay = lambda s: -np.polyfit(f[band], np.unwrap(np.angle(s[:, 1, 0]))[band], 1)[0] / (2 * np.pi)
    case = "microstrip 50 ohm, 20 mm, 0.2 mm dielectric"
    validation.check(case, "group delay 1-10 GHz", delay(so) * 1e12, delay(sn) * 1e12, rel=0.03, unit="ps", source="native channel")
    validation.check(case, "largest |S21| difference to 10 GHz", float(np.max(np.abs(np.abs(so[band, 1, 0]) - np.abs(sn[band, 1, 0])))), 0.0,
                     abs_=0.02, source="native channel")
    validation.check(case, "largest |S11| to 10 GHz", float(np.max(np.abs(so[band, 0, 0]))), 0.0, abs_=0.08, source="matched line")
    assert np.max(np.abs(so[:, 0, 1] - so[:, 1, 0])) < 0.01  # reciprocal


@fullwave
def test_fullwave_patch(tmp_path, validation):
    """An edge-fed patch of 37 x 28.5 mm on 1.5 mm of epsilon_r 4.5 (lossless), its ground
    plane 20 mm beyond the copper: the full-wave resonance (largest input resistance, the
    10 mm feed de-embedded with the native line), edge resistance and directivity against
    the transmission-line model, whose error on such substrates is several per cent."""
    synthetic.patch(tmp_path / "board")
    b = kairoq.open(tmp_path / "board")
    mat = {"D1": {"epsilon_r": 4.5, "loss_tangent": 0.0, "wideband": False}}
    feed = kairoq.analyse_si(b, {"nets": ["ANT"], "materials": mat}).samples[0]
    est = kairoq.patch_antenna(4.5, 1.535, width_mm=37.0, length_mm=28.5)
    si = {"materials": mat, "channels": [{"name": "A", "ports": ["R1.1"]}],
          "openems": {"f_max_hz": 4e9, "points": 400, "margin_mm": 20.0, "min_cell_mm": 0.3, "far_field_hz": [round(est["resonant_hz"], -6)]}}
    kairoq.export_openems_channel(b, si, tmp_path / "export")
    run(tmp_path / "export", tmp_path / "openems.log", "--losses", "simple", "--port-kind", "lumped")  # lossless, the energy balance; the lumped port of docs/validation/rf.md
    r = kairoq.deembed({"input": "export/result.s1p", "operations": [
        {"line": {"port": 1, "length_mm": -10.0, "z_ohm": feed["z0"], "epsilon_eff": feed["epsilon_eff"]}}]}, base_dir=tmp_path)
    f, s = np.array(r["frequency_hz"]), np.array(r["s"]).ravel()
    z = 50 * (1 + s) / (1 - s)
    k = int(np.argmax(z.real))
    ff = json.loads((tmp_path / "export" / "result.json").read_text())["far_field"][0]
    case = "edge-fed patch 37 x 28.5 mm, 1.5 mm, epsilon_r 4.5"
    validation.check(case, "resonant frequency", est["resonant_hz"] / 1e9, f[k] / 1e9, rel=0.05, unit="GHz", source="openEMS")
    validation.check(case, "edge resistance", est["edge_resistance_ohm"], z.real[k], rel=0.15, unit="ohm", source="openEMS")
    validation.check(case, "directivity", est["directivity_dbi"], ff["directivity_dbi"], abs_=1.5, unit="dBi", source="openEMS")
    validation.check(case, "radiation efficiency (lossless)", ff["radiation_efficiency"], 1.0, abs_=0.03, source="energy balance")
