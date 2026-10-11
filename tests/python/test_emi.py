"""EMI estimates: loop inductance against closed forms, magnetic moment and radiated field,
harmonics, return paths and switch-node coupling."""

from __future__ import annotations

import math

import numpy as np
import pytest

import kairoq
import synthetic

MU0, C0, ETA0, EPS0 = 4e-7 * math.pi, 299792458.0, 376.730313668, 8.8541878128e-12


def ellip_k(k):
    """Complete elliptic integral of the first kind, by the arithmetic-geometric mean."""
    a, b = 1.0, math.sqrt(1 - k * k)
    while abs(a - b) > 1e-15:
        a, b = (a + b) / 2, math.sqrt(a * b)
    return math.pi / (2 * a)


def filaments(l, d):
    """Partial mutual inductance of two equal, aligned parallel filaments (Grover); with d
    the radius, the partial self-inductance of a tube."""
    return MU0 / (2 * math.pi) * (l * math.asinh(l / d) - math.sqrt(l * l + d * d) + d)


def loop(board, elements, **settings):
    return kairoq.analyse_emi(board, {"loops": [{"name": "X", "elements": elements}], "loop_settings": settings}).loops["X"]


def per_length(make, elements, lengths, per_mm):
    """Inductance per length from two lengths at the same mesh density (the ends cancel)."""
    out = []
    for length in lengths:
        b = make(length)
        out.append(loop(b, elements, max_triangles=int(per_mm * length)))
    return (out[1]["inductance_nh"] - out[0]["inductance_nh"]) * 1e-9 / ((lengths[1] - lengths[0]) * 1e-3), out


def test_coplanar_strips(tmp_path, validation):
    """Two parallel strips 1 mm wide, 1 mm apart: L' = mu0 K(k)/K(k'), k = s/(s + 2w), in the
    high-frequency limit (conformal mapping; Simons, Coplanar Waveguide Circuits)."""
    def make(length):
        d = tmp_path / f"h{length:g}"
        synthetic.hairpin(d, length=length)
        return kairoq.open(d)

    lp, res = per_length(make, [{"from": "R1.2", "to": "R1.1"}], (8.0, 16.0), 200)
    k = 1.0 / 3.0
    ref = MU0 * ellip_k(k) / ellip_k(math.sqrt(1 - k * k))
    validation.check("coplanar strips w = s = 1 mm", "inductance per length", lp * 1e6, ref * 1e6, rel=0.04, unit="nH/mm",
                     source="conformal mapping, coplanar strips")
    for r in res:  # the last two adaptive passes agree closely
        assert abs(r["inductance_nh"] - r["inductance_coarse_nh"]) < 0.03 * r["inductance_nh"]


def test_microstrip(tmp_path, validation):
    """A 1 mm strip 1.535 mm over a plane: L' = Z0 / c0 with Z0 of Hammerstad and Jensen in
    air (accurate to 0.01 %)."""
    def make(length):
        d = tmp_path / f"m{length:g}"
        synthetic.microstrip(d, length=length)
        return kairoq.open(d)

    lp, _ = per_length(make, [{"from": "R1.1", "to": "R1.2"}, {"from": "R2.1", "to": "R2.2"}], (8.0, 16.0), 200)
    u = 1.0 / 1.535
    f = 6 + (2 * math.pi - 6) * math.exp(-(30.666 / u) ** 0.7528)
    ref = ETA0 / (2 * math.pi) * math.log(f / u + math.sqrt(1 + 4 / u ** 2)) / C0
    validation.check("microstrip w = 1 mm, h = 1.535 mm", "inductance per length", lp * 1e6, ref * 1e6, rel=0.04, unit="nH/mm",
                     source="Hammerstad and Jensen, air")


@pytest.fixture(scope="module")
def ring(tmp_path_factory):
    d = tmp_path_factory.mktemp("ring")
    g = synthetic.ring_loop(d)
    b = kairoq.open(d)
    return g, b, loop(b, [{"from": "R1.2", "to": "R1.1", "radius_mm": 0.3}], max_triangles=2500)


def test_rectangular_loop(ring, validation):
    """Rectangular loop 20 x 10 mm of a 0.5 mm strip: the partial inductances of four
    straight conductors of equivalent radius w/4 (Grover; Ruehli)."""
    g, _, r = ring
    a, b, w = g["a_m"], g["b_m"], g["width_m"]
    ref = 2 * filaments(a, w / 4) + 2 * filaments(b, w / 4) - 2 * filaments(a, b) - 2 * filaments(b, a)
    validation.check("rectangular loop 20 x 10 mm, strip 0.5 mm", "inductance", r["inductance_nh"], ref * 1e9, rel=0.03, unit="nH",
                     source="partial inductances of the sides")
    validation.check("rectangular loop 20 x 10 mm, strip 0.5 mm", "magnetic moment per ampere", r["moment_m2"][2] * 1e6, a * b * 1e6,
                     rel=0.005, unit="mm^2", source="area of the centre line")


def test_moment_is_independent_of_position(tmp_path, ring):
    """The current conserves charge, so its moment does not depend on the origin."""
    synthetic.ring_loop(tmp_path, margin=40.0)
    r2 = loop(kairoq.open(tmp_path), [{"from": "R1.2", "to": "R1.1", "radius_mm": 0.3}], max_triangles=2500)
    assert r2["moment_m2"][2] == pytest.approx(ring[2]["moment_m2"][2], rel=2e-3)


def test_loop_over_plane(tmp_path, ring, validation):
    """The same loop 1.535 mm over a plane of another net reaching 20 mm beyond it: the
    plane's induced current acts as the loop's image, L = L_free - M(loop, image at 2h), and
    cancels most of the moment (all of it for an infinite plane)."""
    g = synthetic.ring_loop(tmp_path, plane=True, margin=25.0)
    r = loop(kairoq.open(tmp_path), [{"from": "R1.2", "to": "R1.1", "radius_mm": 0.3}], max_triangles=3000, window_margin_mm=20)
    a, b, D = g["a_m"], g["b_m"], 2 * g["height_m"]

    def m_par(l, d):
        return filaments(l, d)

    image = 2 * (m_par(a, D) - m_par(a, math.sqrt(b * b + D * D))) + 2 * (m_par(b, D) - m_par(b, math.sqrt(a * a + D * D)))
    free = ring[2]["inductance_nh"] * 1e-9
    validation.check("rectangular loop 20 x 10 mm over a plane, h = 1.535 mm", "inductance", r["inductance_nh"], (free - image) * 1e9,
                     rel=0.03, unit="nH", source="image in the plane")
    validation.check("rectangular loop 20 x 10 mm over a plane, h = 1.535 mm", "moment left / moment without the plane",
                     abs(r["moment_m2"][2]) / ring[2]["moment_m2"][2], 0.0, abs_=0.1, source="image cancels the moment")


def test_harmonics(validation):
    """Trapezoid with equal edges: c_n = 2 A d |sinc(n pi d)| |sinc(n pi tr/T)| (Paul); with
    unequal edges, against the discrete Fourier transform of the sampled waveform."""
    a, f0, duty, tr = 2.0, 1e6, 0.3, 20e-9
    h = kairoq._core.trapezoid_harmonics(a, f0, duty, tr, -1.0, 50e6)
    sinc = lambda x: 1.0 if x == 0 else math.sin(x) / x
    worst = max(abs(c - 2 * a * duty * abs(sinc(n * math.pi * duty)) * abs(sinc(n * math.pi * tr * f0))) for n, (f, c) in enumerate(h, 1))
    validation.check("trapezoid 2 A, 1 MHz, d = 0.3, 20 ns edges", "largest deviation of a harmonic", worst, 0.0, abs_=1e-12,
                     unit="A", source="sinc envelope")
    tf = 50e-9
    h2 = kairoq._core.trapezoid_harmonics(a, f0, duty, tr, tf, 20e6)
    N = 1 << 18
    t = np.arange(N) / N / f0
    t3 = duty / f0 + 0.5 * (tr - tf)
    x = np.clip(np.minimum(t / tr, (t3 + tf - t) / tf), 0, 1) * a
    spec = 2 * np.abs(np.fft.rfft(x)) / N
    err = max(abs(c - spec[n]) for n, (f, c) in enumerate(h2, 1))
    assert err < 1e-6


def test_small_loop_field(ring, validation):
    """Small-loop far field: E = 1.316e-14 f^2 A I / r (Ott), with the computed moment."""
    _, b, r = ring
    e = kairoq.analyse_emi(b, {"loops": [{"name": "X", "elements": [{"from": "R1.2", "to": "R1.1", "radius_mm": 0.3}],
                                            "waveform": {"amplitude_a": 3.0, "frequency_hz": 2e6, "duty": 0.4, "rise_ns": 5}}],
                                 "loop_settings": {"max_triangles": 2500, "distance_m": 3.0, "max_frequency_hz": 300e6,
                                                   "field": "small_loop"}}).loops["X"]
    m = math.sqrt(sum(v * v for v in e["moment_m2"]))
    f, i, field = e["harmonic_frequency_hz"], e["harmonic_current_a"], e["harmonic_field_v_per_m"]
    ott = 1.316e-14 * f ** 2 * m * i / 3.0
    validation.check("loop 20 x 10 mm, 3 A, 2 MHz, 5 ns", "field / Ott's formula, largest deviation", float(np.max(np.abs(field / ott - 1))), 0.0,
                     abs_=2e-3, source="small-loop far field")
    assert e["max_field_v_per_m"] == pytest.approx(float(field.max()))


@pytest.mark.parametrize("kind,expect", [("solid", None), ("slot", "crossing"), ("split", "split"), ("antipads", "crossing"), ("notch", "partial")])
def test_return_path(tmp_path, validation, kind, expect):
    g = synthetic.return_path_board(tmp_path, kind)
    e = kairoq.analyse_emi(kairoq.open(tmp_path), {"return_path": True})
    found = [f for f in e.return_path if f["net"] == "SIG"]
    if expect is None:
        assert found == []
        return
    assert [f["kind"] for f in found] == [expect]
    f = found[0]
    assert f["position_mm"][0] == pytest.approx(g["x_mid_mm"], abs=0.05)
    if kind == "split":
        assert f["reference_nets"] == ["GND", "PWR"]
    if kind in ("slot", "split"):
        gap = 1.0 if kind == "slot" else 0.5
        validation.check(f"trace 0.3 mm over a {gap:g} mm {kind}", "span over the gap", f["span_mm"], gap, rel=0.03, unit="mm",
                         source="drawn gap")


@pytest.mark.parametrize("shield", [False, True])
def test_switch_node_coupling(tmp_path, validation, shield):
    """Parallel-plate capacitance of a 10 x 8 mm patch to a plate under the board through
    the dielectric (4.5), the bottom copper layer and mask (3.8) and a 0.5 mm gap (4.0); a
    plane under half of the patch shields that half. With fringing the capacitance is larger,
    and the shield takes a share of it."""
    g = synthetic.switch_node_board(tmp_path, shield=shield)
    e = kairoq.analyse_emi(kairoq.open(tmp_path), {"plates": [{"name": "CASE", "face": "bottom", "gap_mm": 0.5, "epsilon_r": 4.0}],
                                                       "switch_nodes": [{"net": "SW", "plate": "CASE", "dv_dt_v_per_ns": 10}]})
    c = e.coupling[0]
    d = g["dielectric_m"] / 4.5 + 35e-6 / 3.8 + 10e-6 / 3.8 + 0.5e-3 / 4.0
    area = g["half_area_m2"] if shield else g["area_m2"]
    validation.check(f"patch 10 x 8 mm over a plate{', half shielded' if shield else ''}", "parallel-plate capacitance", c["parallel_plate_pf"],
                     EPS0 * area / d * 1e12, rel=1e-9, unit="pF", source="parallel plates, layered")
    assert c["common_mode_current_ma"] == pytest.approx(10 * c["capacitance_pf"])
    assert c["capacitance_pf"] > 1.2 * c["parallel_plate_pf"]
    assert c["total_pf"] == pytest.approx(c["capacitance_pf"] + sum(pf for _, pf in c["nets_pf"]), rel=1e-9)
    if shield:
        assert c["neighbours"][0][:2] == ("GND", 2)
        assert c["neighbours"][0][2] == pytest.approx(EPS0 * 4.5 * g["half_area_m2"] / g["dielectric_m"] * 1e12, rel=1e-9)
        assert c["nets_pf"][0][0] == "GND" and c["nets_pf"][0][1] > c["capacitance_pf"]


def disc(radius, n=720):
    return [(radius * math.cos(2 * math.pi * k / n), radius * math.sin(2 * math.pi * k / n)) for k in range(n)]


def axisymmetric_disc(radius, z_mm, layers, thickness=0.0005, plane=30.0):
    """Capacitance of a thin disc at height z over a plane, in horizontal dielectric layers
    ((z0, z1, epsilon_r) in mm), by the axisymmetric field solver."""
    r = kairoq.cross_section({"axisymmetric": True, "margin": 2.0, "target_error": 0.005,
                                "conductors": [{"name": "disc", "rect_mm": [0, z_mm - thickness / 2, radius, z_mm + thickness / 2]},
                                               {"name": "plane", "reference": True, "rect_mm": [0, -0.1, plane, 0]}],
                                "dielectrics": [{"rect_mm": [0, z0, plane, z1], "epsilon_r": er, "wideband": False} for z0, z1, er in layers]})
    return r["frequencies"][0]["c"][0]


@pytest.mark.parametrize("case", ["vacuum", "layered"])
def test_switch_node_with_fringing(tmp_path, validation, case):
    """A disc of 5 mm radius over a plate, against the axisymmetric field solution of the
    same stack: in vacuum 0.5 mm over the plate; on FR-4 (4.5) of 0.5 mm with masks of
    10 um (3.8) over an air gap of 0.5 mm."""
    t, h = 0.002, 0.5
    if case == "vacuum":
        synthetic.node_board(tmp_path, disc(5.0), dielectric=h, copper=t, epsilon_r="1.0", mask_mm=0)
        gap, mask, layers = 0.0, 0.0, []
    else:
        synthetic.node_board(tmp_path, disc(5.0), dielectric=h, copper=t, epsilon_r="4.5", mask_mm=0.01)
        gap, mask = 0.5, 0.01
        z, layers = gap, []
        for th, er in ((mask, 3.8), (t, 4.5), (h, 4.5), (t, 4.5), (mask, 3.8)):
            layers.append((z, z + th, er))
            z += th
    e = kairoq.analyse_emi(kairoq.open(tmp_path), {"plates": [{"name": "P", "face": "bottom", "gap_mm": gap, "epsilon_r": 1.0}],
                                                       "switch_nodes": [{"net": "SW", "plate": "P"}]})
    c = e.coupling[0]
    ref = axisymmetric_disc(5.0, gap + mask + 1.5 * t + h, layers) * 1e12
    validation.check(f"disc 5 mm radius over a plate, {case}", "capacitance", c["capacitance_pf"], ref, rel=0.005, unit="pF",
                     source="axisymmetric field solution")
    assert c["coarse_pf"] < c["mesh_pf"] < c["capacitance_pf"]  # the Galerkin charges converge from below


def test_switch_node_strip_per_length(tmp_path, validation):
    """A strip 1 mm wide on FR-4 (4.5) 0.5 mm over a plate: the capacitance per length from
    strips of 10 and 20 mm against the cross-section solver, whose conductors of 5 and 2 um
    extrapolate to zero thickness."""
    w, h, t = 1.0, 0.5, 0.002
    cs = []
    for length in (10.0, 20.0):
        d = tmp_path / f"s{length:g}"
        synthetic.node_board(d, [(-length / 2, -w / 2), (length / 2, -w / 2), (length / 2, w / 2), (-length / 2, w / 2)],
                             dielectric=h, copper=t, epsilon_r="4.5", mask_mm=0)
        e = kairoq.analyse_emi(kairoq.open(d), {"plates": [{"name": "P", "face": "bottom"}], "switch_nodes": [{"net": "SW", "plate": "P"}]})
        cs.append(e.coupling[0]["capacitance_pf"])
    zc = 1.5 * t + h

    def cut(th):
        return kairoq.cross_section({"conductors": [{"name": "s", "rect_mm": [-w / 2, zc - th / 2, w / 2, zc + th / 2]},
                                                      {"name": "g", "reference": True, "rect_mm": [-15, -0.1, 15, 0]}],
                                       "dielectrics": [{"rect_mm": [-15, 0, 15, h + 2 * t], "epsilon_r": 4.5, "wideband": False}],
                                       "target_error": 0.005})["frequencies"][0]["c"][0] * 1e9
    c5, c2 = cut(0.005), cut(0.002)
    validation.check("strip 1 mm on FR-4 0.5 mm over a plate", "capacitance per length", (cs[1] - cs[0]) / 10, c2 - (c5 - c2) * 2 / 3,
                     rel=0.01, unit="pF/mm", source="cross-section solver, zero thickness")


def test_settings_are_checked(ring):
    _, b, _ = ring
    with pytest.raises(RuntimeError, match="no pad"):
        loop(b, [{"from": "R9.1", "to": "R1.1"}])
    with pytest.raises(RuntimeError, match="unknown key"):
        kairoq.analyse_emi(b, {"loops": [], "colour": 1})
    with pytest.raises(RuntimeError, match="free, ground or small_loop"):
        kairoq.analyse_emi(b, {"loops": [], "loop_settings": {"field": "chamber"}})


def test_moment_independent_of_window(tmp_path, validation):
    """Over a finite plane the residual moment depends on how much of the plane the window
    holds (37 mm2 at a 5 mm margin, 14 mm2 at 20 mm): the window grows until the moment
    settles or the window holds the board, whatever margin it starts from."""
    synthetic.ring_loop(tmp_path, plane=True, margin=25.0)
    b = kairoq.open(tmp_path)
    el = [{"from": "R1.2", "to": "R1.1", "radius_mm": 0.3}]
    small = loop(b, el, max_triangles=1500, window_margin_mm=5)
    large = loop(b, el, max_triangles=1500, window_margin_mm=20)
    validation.check("rectangular loop over a plane 25 mm beyond it", "moment from margins of 5 and 20 mm",
                     abs(small["moment_m2"][2]), abs(large["moment_m2"][2]), rel=0.05, unit="m2", source="the same board")
    assert small["window_margin_mm"] > 5


def test_loop_resistance_and_inductance_over_frequency(tmp_path, validation):
    """Two coplanar strips 1 mm wide, 1 mm apart, 35 um thick: per length (difference of two
    lengths), the resistance at DC is 2 / (sigma t w); at 1 GHz the inductance approaches the
    high-frequency limit and the resistance that of zero-thickness strips by Wheeler's rule,
    Rs / mu0 dL/dn over both faces."""
    def make(length):
        d = tmp_path / f"h{length:g}"
        synthetic.hairpin(d, length=length)
        return kairoq.open(d)

    fs = [1e3, 1e9]
    out = []
    for length in (8.0, 16.0):
        r = kairoq.analyse_emi(make(length), {"loops": [{"name": "X", "elements": [{"from": "R1.2", "to": "R1.1"}]}],
                                                "loop_settings": {"max_triangles": int(400 * length), "frequencies_hz": fs,
                                                                  "moment_tolerance": 0}}).loops["X"]
        out.append(r["impedance"])
    per = lambda k, i: (out[1][k][i] - out[0][k][i]) / 8e-3
    validation.check("coplanar strips w = s = 1 mm, 35 um, 1 kHz", "resistance per length", per(0, 1), 2 / (5.8e7 * 35e-6 * 1e-3), rel=0.01,
                     unit="ohm/m", source="2 / (sigma t w)")
    k = 1.0 / 3.0
    l_hf = MU0 * ellip_k(k) / ellip_k(math.sqrt(1 - k * k))
    validation.check("coplanar strips w = s = 1 mm, 1 GHz", "inductance per length", per(1, 2) * 1e-3, l_hf * 1e6, rel=0.04,
                     unit="nH/mm", source="conformal mapping, coplanar strips")
    rs = math.sqrt(math.pi * 1e9 * MU0 / 5.8e7)
    lw = lambda w, s: MU0 * ellip_k(s / (s + 2 * w)) / ellip_k(math.sqrt(1 - (s / (s + 2 * w)) ** 2))
    dn = 1e-8
    wheeler = 2 * rs / MU0 * (lw(1e-3 - 2 * dn, 1e-3 + 2 * dn) - lw(1e-3, 1e-3)) / dn
    validation.check("coplanar strips w = s = 1 mm, 1 GHz", "resistance per length", per(1, 1), wheeler, rel=0.1, unit="ohm/m",
                     source="Wheeler's rule, zero-thickness strips, both faces")


def small_loop_exact(f, m, r):
    """Largest field of a small loop of moment m at distance r, near-field term included:
    E = eta0 k^2 m / (4 pi r) |1 + 1 / (j k r)| in the loop's plane."""
    k = 2 * math.pi * f / C0
    return ETA0 * k * k * m / (4 * math.pi * r) * math.sqrt(1 + 1 / (k * r) ** 2)


def test_field_of_the_currents(ring, validation):
    """The field of the solved currents at 3 m in free space, against the small loop's
    with its near-field term, from 1 MHz (k r = 0.06) to 300 MHz."""
    _, b, _ = ring
    e = kairoq.analyse_emi(b, {"loops": [{"name": "X", "elements": [{"from": "R1.2", "to": "R1.1", "radius_mm": 0.3}],
                                            "waveform": {"amplitude_a": 1.0, "frequency_hz": 1e6, "duty": 0.5, "rise_ns": 5}}],
                                 "loop_settings": {"max_triangles": 2500, "max_frequency_hz": 300e6}}).loops["X"]
    m = math.sqrt(sum(v * v for v in e["moment_m2"]))
    worst = max(abs(v / small_loop_exact(f, m, 3.0) - 1) for f, v in e["field_per_ampere"])
    validation.check("loop 20 x 10 mm at 3 m, 1-300 MHz", "field per ampere / small loop with near field, largest deviation", worst, 0.0,
                     abs_=1e-3, source="small loop, near and far field")
    f, i, field = e["harmonic_frequency_hz"], e["harmonic_current_a"], e["harmonic_field_v_per_m"]
    ratio = field / np.array([small_loop_exact(x, m, 3.0) * a for x, a in zip(f, i)])
    assert np.max(np.abs(ratio - 1)) < 2e-3  # interpolated to the harmonics


def test_field_of_a_large_loop(tmp_path, validation):
    """A loop 100 x 100 mm, one wavelength around at 750 MHz: the far field of the solved
    currents against the radiation integral of a uniform current on the centre line, the
    largest over the sphere (2-degree grid)."""
    synthetic.ring_loop(tmp_path, a=100.0, b=100.0)
    e = loop(kairoq.open(tmp_path), [{"from": "R1.2", "to": "R1.1", "radius_mm": 0.3}], max_triangles=4000, moment_tolerance=0)
    P, M = np.array(e["current_positions_m"]), np.array(e["current_moments_a_m"])
    a, x0 = 0.1, 0.01
    corners = np.array([(x0, x0), (x0 + a, x0), (x0 + a, x0 + a), (x0, x0 + a), (x0, x0)])
    u = (np.arange(500) + 0.5) / 500
    pts = np.concatenate([c0 + u[:, None] * (c1 - c0) for c0, c1 in zip(corners[:-1], corners[1:])])
    dls = np.concatenate([np.repeat([(c1 - c0) / 500], 500, axis=0) for c0, c1 in zip(corners[:-1], corners[1:])])
    pts, dls = np.c_[pts, np.zeros(len(pts))], np.c_[dls, np.zeros(len(dls))] * np.sign(e["moment_m2"][2])
    th, ph = np.meshgrid(np.radians(np.arange(0, 181, 2.0)), np.radians(np.arange(0, 360, 2.0)))
    dirs = np.stack([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)], -1).reshape(-1, 3)
    for f in (3e8, 7.5e8, 1e9):
        k = 2 * math.pi * f / C0
        n = np.exp(1j * k * (dirs @ (pts - pts.mean(axis=0)).T)) @ dls
        nt = n - np.sum(n * dirs, 1)[:, None] * dirs
        ref = 2 * math.pi * f * MU0 / (4 * math.pi * 3.0) * np.max(np.linalg.norm(nt, axis=1))
        got = kairoq._core.loop_field_per_ampere(P.tolist(), M.tolist(), f, "free", far_only=True)
        validation.check(f"loop 100 x 100 mm, {f / 1e6:g} MHz", "far field per ampere at 3 m", got, ref, rel=2e-3, unit="V/m",
                         source="radiation integral, uniform current on the centre line")
    small = ETA0 * (2 * math.pi * 1e9 / C0) ** 2 * abs(e["moment_m2"][2]) / (4 * math.pi * 3.0)
    assert small > 1.15 * got  # the small-loop estimate is 20 % high here


def test_field_over_a_reflecting_ground(validation):
    """Measurement site: a small loop 0.8 m over a reflecting ground, the antenna 3 m away at
    heights of 1 to 4 m; against the loop and its image (reversed moment), both with their
    near-field terms."""
    a, n = 0.01, 40
    pos, mom = [], []
    for (x0, y0), (dx, dy) in [((-a / 2, -a / 2), (1, 0)), ((a / 2, -a / 2), (0, 1)), ((a / 2, a / 2), (-1, 0)), ((-a / 2, a / 2), (0, -1))]:
        for i in range(n):
            t = (i + 0.5) / n * a
            pos.append((x0 + dx * t, y0 + dy * t, 0.0))
            mom.append((dx * a / n, dy * a / n, 0.0))
    for f in (1e6, 3e7, 1e8, 3e8):
        k = 2 * math.pi * f / C0
        ref = 0.0
        for h in np.linspace(1, 4, 13):
            e = 0
            for z, sign in ((0.8, 1), (-0.8, -1)):
                R = math.hypot(3.0, h - z)
                e += sign * ETA0 * k * k * a * a / (4 * math.pi * R) * (1 + 1 / (1j * k * R)) * (3.0 / R) * np.exp(-1j * k * R)
            ref = max(ref, abs(e))
        got = kairoq._core.loop_field_per_ampere(pos, mom, f, "ground")
        validation.check(f"loop 10 x 10 mm 0.8 m over ground, 3 m, {f / 1e6:g} MHz", "largest field per ampere over the height scan",
                         got, ref, rel=5e-4, unit="V/m", source="loop and image, near and far field")


SLOT_LOOP = [{"from": "J1.2", "to": "J1.1"}, {"from": "R1.1", "to": "R1.2"}]


def test_stitch_across_a_slot(tmp_path, validation):
    """A trace over a plane with a slot across its path: the return detours around the
    slot's end. A part across the slot beside the trace offers a shorter path: as a strap it
    carries most of the return; as a capacitor its share follows the two-mesh circuit of the
    loop and the stitch, whose inductances the solutions with the stitch open and strapped
    give: Z = j w L11 - (j w L12)^2 / (j w L22 + 1 / (j w C))."""
    synthetic.slot_board(tmp_path / "plain", slot_end=40.0, stitch=False)
    synthetic.slot_board(tmp_path / "slot", stitch=False)
    synthetic.slot_board(tmp_path / "stitch")
    plain = loop(kairoq.open(tmp_path / "plain"), SLOT_LOOP, moment_tolerance=0)
    slot = loop(kairoq.open(tmp_path / "slot"), SLOT_LOOP, moment_tolerance=0)
    b = kairoq.open(tmp_path / "stitch")

    def run(c_nf, freqs):
        st = {"from": "C1.1", "to": "C1.2", **({"c_nf": c_nf} if c_nf else {})}
        return kairoq.analyse_emi(b, {"loops": [{"name": "X", "elements": SLOT_LOOP, "stitches": [st]}],
                                        "loop_settings": {"frequencies_hz": freqs, "moment_tolerance": 0}}).loops["X"]
    strap = run(None, [1e9])
    assert slot["inductance_nh"] > 1.3 * plain["inductance_nh"]
    assert plain["inductance_nh"] < strap["inductance_nh"] < 0.8 * slot["inductance_nh"]
    assert strap["stitch_current"][0] == pytest.approx(strap["impedance"][0][3][0].real, rel=1e-3)
    assert abs(strap["stitch_current"][0]) > 0.8
    open_ = run(1e-9, [1e9])
    alpha = strap["impedance"][0][3][0].real
    L11, Ls = open_["impedance"][0][2] * 1e-9, strap["impedance"][0][2] * 1e-9
    L12 = (Ls - L11) / alpha
    L22 = -L12 / alpha
    C = 0.1e-9
    for f, R, L, sc in run(0.1, [1e8, 3e8, 1e9])["impedance"]:
        w = 2 * math.pi * f
        z2 = 1j * w * L22 + 1 / (1j * w * C)
        Z = 1j * w * L11 - (1j * w * L12) ** 2 / z2
        I2 = -1j * w * L12 / z2
        validation.check(f"slot crossing, 100 pF stitch, {f / 1e6:g} MHz", "inductance", L, Z.imag / w * 1e9, rel=0.01, unit="nH",
                         source="two-mesh circuit from the open and strapped stitch")
        validation.check(f"slot crossing, 100 pF stitch, {f / 1e6:g} MHz", "stitch current per ampere", abs(sc[0]), abs(I2), rel=0.01,
                         source="two-mesh circuit from the open and strapped stitch")


def test_streamed_partial_inductances(tmp_path):
    """Above dense_matrix_mb the partial inductances are formed by blocks of rows and
    multiplied out: the same solution as with the dense matrix, frequency sweep and
    stitches included."""
    synthetic.slot_board(tmp_path)
    b = kairoq.open(tmp_path)
    spec = {"name": "X", "elements": SLOT_LOOP, "stitches": [{"from": "C1.1", "to": "C1.2", "c_nf": 0.1}]}
    out = [kairoq.analyse_emi(b, {"loops": [spec], "loop_settings": {"max_triangles": 1500, "frequencies_hz": [1e8, 1e9], "dense_matrix_mb": mb}}).loops["X"]
           for mb in (2048, 0)]
    assert out[1]["inductance_nh"] == pytest.approx(out[0]["inductance_nh"], rel=1e-10)
    assert np.allclose(out[1]["moment_m2"], out[0]["moment_m2"], rtol=1e-9, atol=1e-18)
    for za, zb in zip(out[0]["impedance"], out[1]["impedance"]):
        assert zb[1] == pytest.approx(za[1], rel=1e-9) and zb[2] == pytest.approx(za[2], rel=1e-9)
        assert abs(zb[3][0] - za[3][0]) < 1e-9 * abs(za[3][0])


def test_budget_lowers_the_inductance(tmp_path):
    """The inductance converges from above as the triangle budget grows (docs/solvers/emi.md,
    "Accuracy and cost"); at the default budget the slotted board's first mesh has the
    coarsest elements allowed at the loop, which the result reports and warns of."""
    synthetic.slot_board(tmp_path, stitch=False)
    b = kairoq.open(tmp_path)
    runs = [loop(b, SLOT_LOOP, moment_tolerance=0, max_triangles=n) for n in (4000, 8000)]
    assert runs[1]["inductance_nh"] < runs[0]["inductance_nh"] < runs[0]["inductance_coarse_nh"]
    assert runs[1]["element_size_mm"] < runs[0]["element_size_mm"] == pytest.approx(1.2, rel=1e-6)
    assert any("coarsest elements" in w for w in runs[0]["warnings"])
    assert not any("coarsest elements" in w for w in runs[1]["warnings"])
