"""PDN impedance: a rectangular plane pair with a load and a decoupling capacitor against the
closed-form modal sum of the pair with the same ports, plate and dielectric losses and
capacitor branch."""

from __future__ import annotations

import math

import numpy as np
import pytest

import kairoq
import synthetic


def j0(x):
    """Bessel function J0 from its integral (1/pi) int_0^pi cos(x sin t) dt; the midpoint rule
    on the periodic integrand is exact to 1e-15 for |x| < 40."""
    t = (np.arange(64) + 0.5) * np.pi / 64
    return np.mean(np.cos(np.multiply.outer(np.asarray(x, float), np.sin(t))), axis=-1)

MU0, EPS0 = 1.25663706212e-6, 8.8541878128e-12
CAP = {"component": "C1", "c_uf": 0.1, "esr_mohm": 20, "esl_nh": 0.4, "mount_nh": 0.5}


def rectangle_z(f, a, b, h, er, tand, ports, radius, sigma=5.8e7, t=35e-6, modes=(2000, 1500)):
    """Impedance matrix of a rectangular pair with open edges between circular ports: the
    double sum over its modes, each averaged over the port circles (J0), with the plates'
    surface impedance in the series impedance."""
    w = 2 * math.pi * f
    delta = 1 / math.sqrt(math.pi * f * MU0 * sigma)
    x = (1 + 1j) * t / delta
    zs = (1 + 1j) / (sigma * delta) / np.tanh(x)
    zser = 1j * w * MU0 * h + 2 * zs
    ysh = 1j * w * EPS0 * er * (1 - 1j * tand) / h
    k2 = -zser * ysh
    m = np.arange(modes[0])[:, None]
    n = np.arange(modes[1])[None, :]
    km, kn = m * math.pi / a, n * math.pi / b
    kmn = np.sqrt(km ** 2 + kn ** 2)
    e = np.where(m == 0, 1, 2) * np.where(n == 0, 1, 2)
    z = np.zeros((len(ports), len(ports)), complex)
    for i, (xi, yi) in enumerate(ports):
        for j, (xj, yj) in enumerate(ports):
            num = e * np.cos(km * xi) * np.cos(kn * yi) * np.cos(km * xj) * np.cos(kn * yj) * j0(kmn * radius) ** 2
            z[i, j] = zser / (a * b) * np.sum(num / (kmn ** 2 - k2))
    return z


def test_rectangular_pair_with_a_capacitor(tmp_path, validation):
    g = synthetic.plane_pair(tmp_path)
    b = kairoq.open(tmp_path)
    fs = [1e6, 30e6, 150e6, 600e6, 1.2e9, 2.2e9, 3.3e9]
    r = kairoq.pdn(b, {"power_net": "PWR", "ground_net": "GND", "f_min_hz": fs[0], "f_max_hz": 4e9, "points": 41, "loads": ["U1"],
                         "capacitors": [CAP], "dielectric": {"epsilon_r": 4.5, "loss_tangent": 0.02, "wideband": False}})
    assert len(r["pairs"]) == 1 and r["pairs"][0]["height_m"] == pytest.approx(0.2e-3, rel=1e-6)
    f = np.array(r["frequency_hz"])
    z = np.array(r["z"][0])
    worst = 0.0
    for k in range(0, len(f), 8):
        zc = rectangle_z(f[k], 0.04, 0.03, 0.2e-3, 4.5, 0.02, [(0.012, 0.009), (0.030, 0.022)], 0.15e-3)
        w = 2 * math.pi * f[k]
        ycap = 1 / (CAP["esr_mohm"] * 1e-3 + 1j * w * (CAP["esl_nh"] + CAP["mount_nh"]) * 1e-9 + 1 / (1j * w * CAP["c_uf"] * 1e-6))
        y = np.linalg.inv(zc) + np.diag([0, ycap])
        ref = np.linalg.inv(y)[0, 0]
        worst = max(worst, abs(z[k] - ref) / (abs(ref) + 0.3))
    # Relative to |Z| + 0.3 ohm: at the minima, where |Z| falls to a few tens of milliohms,
    # the difference is a few milliohms.
    validation.check("plane pair 40 x 30 mm, 0.2 mm, load and 0.1 uF, 1 MHz to 4 GHz", "largest difference of Z / (|Z| + 0.3 ohm)", worst, 0.0,
                     abs_=0.03, source="closed-form modal sum of the rectangular pair")
    # Low frequency: the capacitor and the plates' capacitance in parallel.
    c_plane = EPS0 * 4.5 * 0.04 * 0.03 / 0.2e-3
    validation.check("plane pair 40 x 30 mm, 1 MHz", "|Z|", abs(z[0]), 1 / (2 * math.pi * f[0] * (0.1e-6 + c_plane)), rel=0.01, unit="ohm",
                     source="capacitances in parallel")


def grover_rectangle(a, b, r):
    """Inductance of a rectangular loop a x b of round wire of radius r (Grover)."""
    d = math.hypot(a, b)
    return MU0 / math.pi * (-2 * (a + b) + 2 * d - a * math.log((a + d) / b) - b * math.log((b + d) / a) + a * math.log(2 * a / r) + b * math.log(2 * b / r))


def test_mounting_inductance(tmp_path, validation):
    """A capacitor with via-in-pad, pads and vias of one radius: the loop of the vias from the
    nearer plane's face and the pads between them, over the plane, is half a rectangle of
    twice the height with its image."""
    for pitch, top in ((2.0, 0.5), (3.0, 1.0)):
        d = tmp_path / f"m{pitch:g}"
        synthetic.plane_pair(d, pad=0.2, drill=0.1, pitch=pitch, dielectric=(top, 0.2, top), antipad=0.4)
        r = kairoq.pdn(kairoq.open(d), {"power_net": "PWR", "ground_net": "GND", "f_min_hz": 1e6, "f_max_hz": 1e8, "points": 3,
                                            "loads": ["U1"], "capacitors": [{"component": "C1", "c_uf": 0.1}]})
        m = next(b["mount_h"] for b in r["branches"] if b["component"] == "C1")
        h = (0.0175 + top) * 1e-3  # from the pads' layer to the plane's face
        validation.check(f"mounting loop, pitch {pitch:g} mm, {h * 1e3:.3f} mm over the plane", "inductance", m * 1e9,
                         grover_rectangle(pitch * 1e-3, 2 * h, 0.05e-3) / 2 * 1e9, rel=0.03, unit="nH", source="Grover, rectangle and image")
