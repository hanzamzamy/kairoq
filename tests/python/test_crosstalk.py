"""Crosstalk: coupling coefficients of the board scan, coupled sections of channels, and
the waveforms at a victim's ends.

Against the exact theory of the TEM directional coupler (coupled stripline in one
dielectric), the weak-coupling formulas of near- and far-end crosstalk, and the
cancellation of far-end crosstalk in a homogeneous dielectric.
"""

from __future__ import annotations

import math

import pytest

import kairoq
import synthetic
from conftest import TESTDATA

C0 = 299792458.0
PORTS = ["R1.1", "R1.2", "R2.1", "R2.2"]  # A near, B near, A far, B far

# A lossless, non-dispersive dielectric: the coupler theory is for TEM lines.
LOSSLESS = {"materials": {"D1": {"epsilon_r": 4.5, "loss_tangent": 0.0, "wideband": False},
                          "D2": {"epsilon_r": 4.5, "loss_tangent": 0.0, "wideband": False}},
            "copper": {"conductivity_s_per_m": 1e12}}


def coupled_stripline_modes(w=0.2, s=0.2, h=0.2, t=0.035, er=4.5):
    """Even- and odd-mode impedance of the coupled stripline of synthetic.coupled_lines."""
    zb, zt = h, h + t
    spec = {"conductors": [{"name": "a", "rect_mm": [s / 2, zb, s / 2 + w, zt]}, {"name": "b", "rect_mm": [-s / 2 - w, zb, -s / 2, zt]},
                           {"name": "g1", "rect_mm": [-5, -0.035, 5, 0], "reference": True},
                           {"name": "g2", "rect_mm": [-5, zt + h, 5, zt + h + 0.035], "reference": True}],
            "dielectrics": [{"rect_mm": [-5, 0, 5, zt + h], "epsilon_r": er, "wideband": False}],
            "frequencies_hz": [1e9], "target_error": 0.02, "margin": 0.2}
    p = kairoq.cross_section(spec)["frequencies"][0]["pair"]
    return p["z_even"], p["z_odd"]


def test_directional_coupler(tmp_path, validation):
    """Coupled stripline in one dielectric, ports of sqrt(Ze Zo): the coupled port's
    transmission peaks where the line is a quarter wave long, at (Ze - Zo) / (Ze + Zo);
    the isolated port receives nothing and the ports are matched (Oliver 1954)."""
    g = synthetic.coupled_lines(tmp_path, length=7.0, stripline=True)
    b = kairoq.open(tmp_path)
    ze, zo = coupled_stripline_modes()
    z0 = math.sqrt(ze * zo)
    si = {"channels": [{"name": "C", "ports": PORTS, "impedance_ohm": z0}], "network": {"f_max_hz": 10e9, "points": 1000},
          "cut": {"target_error": 0.02}, **LOSSLESS}
    r = kairoq.analyse_channels(b, si)[0]
    coupled = [abs(m[1 * 4 + 0]) for m in r["s"]]
    k = max(range(len(coupled)), key=lambda i: coupled[i])
    case = "coupled stripline w = s = 0.2 mm, 7 mm, er 4.5"
    validation.check(case, "largest coupling |S21|", coupled[k], (ze - zo) / (ze + zo), rel=0.01, source="TEM coupler, (Ze - Zo)/(Ze + Zo)")
    validation.check(case, "frequency of the largest coupling", r["frequency_hz"][k] / 1e9, C0 / (4 * 7e-3 * math.sqrt(4.5)) / 1e9, rel=0.02,
                     unit="GHz", source="quarter wave")
    assert max(abs(m[3 * 4 + 0]) for m in r["s"]) < 0.01  # isolation
    assert max(abs(m[0]) for m in r["s"]) < 0.01  # match


def test_far_end_cancels_in_one_dielectric(tmp_path, validation):
    """In a homogeneous dielectric the inductive and capacitive coupling are equal and the
    far-end crosstalk vanishes."""
    synthetic.coupled_lines(tmp_path, stripline=True)
    b = kairoq.open(tmp_path)
    r = kairoq.analyse_si(b, {"nets": ["A", "B"], "cut": {"target_error": 0.02}, **LOSSLESS})
    (c, _) = r.coupling
    validation.check("coupled stripline w = s = 0.2 mm, 20 mm", "far-end crosstalk, 100 ps", c["fext"], 0.0, abs_=2e-3,
                     source="homogeneous dielectric")
    assert c["length_mm"] == pytest.approx(20.0)


def test_microstrip_crosstalk(tmp_path, validation):
    """Coupled microstrip, ports matched to the lines: the waveforms at the victim's ends
    against the weak-coupling estimates of the scan, a saturated near end of Kb and a far
    end of Kf l / t_r (for the Gaussian step 1.02 Kf l / t_r)."""
    synthetic.coupled_lines(tmp_path)
    b = kairoq.open(tmp_path)
    scan = kairoq.analyse_si(b, {"nets": ["A", "B"], "crosstalk": {"rise_ps": 100}, **LOSSLESS})
    c = next(x for x in scan.coupling if x["victim"] == "B")
    z = next(n for n in scan.nets if n["name"] == "B")["z_mean"]
    si = {"channels": [{"name": "X", "ports": PORTS, "impedance_ohm": z, "drive": [1]}],
          "network": {"f_max_hz": 40e9, "points": 800, "rise_ps": 100}, **LOSSLESS}
    r = kairoq.analyse_channels(b, si)[0]
    (step,) = r["steps"]
    near = max(step["values"][1], key=abs)
    far = max(step["values"][3], key=abs)
    case = "coupled microstrip w = s = h = 0.2 mm, 20 mm, 100 ps"
    validation.check(case, "near-end crosstalk", near, c["kb_mean"], rel=0.05, source="Kb (saturated)")
    kf_l = c["kf_ps_per_mm"] * 1e-12 * 20.0
    validation.check(case, "far-end crosstalk", far, 1.0225 * kf_l / 100e-12, rel=0.10, source="1.02 Kf l / t_r")
    assert c["next"] == pytest.approx(c["kb_mean"])  # the run's round trip outlasts the rise
    assert far < 0  # microstrip: inductive coupling dominates at the far end


def test_board_scan():
    """The HDMI pairs of cm5_minima: coupling between neighbouring pairs, not within one."""
    b = kairoq.open(TESTDATA / "kicad" / "cm5_minima" / "cm5_minima.ipc2581.zip")
    r = kairoq.analyse_si(b, {"nets": ["/CM5/HDMI_PI.D*"]})
    assert r.coupling
    for c in r.coupling:
        assert c["victim"][:-1] != c["aggressor"][:-1]  # partners are not crosstalk
        assert 0 <= c["next"] < 0.5 and 0 <= c["fext"] < 0.5 and c["length_mm"] > 0
