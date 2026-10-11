"""Transient thermal solutions: exact discrete cases, diffusion, the steady limit and the
heat balance."""

from __future__ import annotations

import math

import numpy as np
import pytest

import kairoq
import synthetic

K_CU, K_Z, K_MASK = 385.0, 0.3, 0.25
C_CU, C_D, C_MASK = 3.45e6, 2.2e6, 1.8e6


def faces(h=10.0, eps=0.0):
    return {"top": {"convection_w_m2k": h, "emissivity": eps}, "bottom": {"convection_w_m2k": h, "emissivity": eps}}


ONE_SLICE = {"dielectric_slice_mm": 100.0, "dielectric_first_slice_mm": 100.0}


def chain(g, c, load, x0, times):
    """Backward Euler on a lumped network: conductance matrix g (W/K), capacities c (J/K),
    load(t) (W), from x0 over the given times."""
    x = np.array(x0, float)
    out = [x.copy()]
    for t0, t1 in zip(times[:-1], times[1:]):
        dt = t1 - t0
        x = np.linalg.solve(np.diag(c / dt) + g, c / dt * x + load(t1))
        out.append(x.copy())
    return np.array(out)


def test_slab_network_is_exact(tmp_path, validation):
    """Uniform slab, one dielectric slice: junction, top sheet, bottom sheet and a heat sink
    body form a chain of four nodes, which the solution must follow step by step."""
    g_ = synthetic.slab(tmp_path)
    b = kairoq.open(tmp_path)
    ta, rjc, cj, rsa, cs, tim = 25.0, 4.0, 0.5, 3.0, 20.0, (0.5e-3, 3.0)
    sched = [[0.0, 0.0], [0.0, 2.0], [40.0, 2.0], [40.0, 0.5]]
    r = kairoq.solve_thermal(b, {
        "ambient_c": ta, **faces(), **ONE_SLICE,
        "sources": [{"name": "H", "pad": "H1.1", "power_schedule": sched, "junction_k_per_w": rjc, "heat_capacity_j_per_k": cj}],
        "heatsinks": [{"name": "S", "face": "bottom", "tim": {"thickness_mm": tim[0] * 1e3, "k": tim[1]},
                       "resistance_k_per_w": rsa, "heat_capacity_j_per_k": cs}],
        "transient": {"duration_s": 100.0, "step_s": 0.5, "max_step_s": 8.0, "steps_per_doubling": 3}})
    a, td, t_cu, tm = g_["area_m2"], g_["dielectric_m"], 35e-6, g_["mask_m"]
    g_top = a / (tm / K_MASK + 1 / 10.0)
    g_mid = a * K_Z / td
    g_sink = a / (tm / K_MASK + tim[0] / tim[1])
    G = np.zeros((4, 4))
    for i, j, gv in ((0, 1, 1 / rjc), (1, 2, g_mid), (2, 3, g_sink)):
        G[i, i] += gv
        G[j, j] += gv
        G[i, j] -= gv
        G[j, i] -= gv
    G[1, 1] += g_top
    G[3, 3] += 1 / rsa
    sheet = a * (C_CU * t_cu + C_D * td / 2 + C_MASK * tm)
    c = np.array([cj, sheet, sheet, cs])

    def power(t):
        return 2.0 if t <= 40.0 else 0.5

    def load(t):
        return np.array([power(t), g_top * ta, 0.0, ta / rsa])

    tr = r.transient
    times = tr["time_s"]
    ref = chain(G, c, load, [ta] * 4, times)
    worst = max(np.max(np.abs(tr["source_junction_c"][0] - ref[:, 0])), np.max(np.abs(tr["heatsink_c"][0] - ref[:, 3])),
                np.max(np.abs(tr["source_mean_c"][0] - ref[:, 1])))
    validation.check("slab with junction and sink body, 100 s", "largest deviation over all steps", worst, 0.0, abs_=1e-7,
                     unit="K", source="backward Euler on the four-node network")
    assert tr["steps"] == len(times) - 1
    # The schedule's step at 40 s is a breakpoint, with the step length reset.
    assert 40.0 in set(np.round(times, 12))
    assert np.min(np.diff(times)[times[1:] > 40.0][:1]) == pytest.approx(0.5)


def test_time_convergence(tmp_path, validation):
    """The same slab without junction and sink against the exact two-node solution: the
    error of backward Euler halves with the step."""
    g_ = synthetic.slab(tmp_path)
    b = kairoq.open(tmp_path)
    a, td, tm = g_["area_m2"], g_["dielectric_m"], g_["mask_m"]
    h_eff = a / (tm / K_MASK + 1 / 10.0)
    g_mid = a * K_Z / td
    G = np.array([[h_eff + g_mid, -g_mid], [-g_mid, g_mid + h_eff]])
    sheet = a * (C_CU * 35e-6 + C_D * td / 2 + C_MASK * tm)
    lam, vec = np.linalg.eigh(G / sheet)
    steady = np.linalg.solve(G, [2.0, 0.0])
    t_end = 60.0

    def exact(t):  # rise above ambient
        return steady - vec @ (np.exp(-lam * t) * (vec.T @ steady))

    errs = []
    for dt in (2.0, 1.0, 0.5):
        r = kairoq.solve_thermal(b, {**faces(), **ONE_SLICE, "sources": [{"name": "H", "pad": "H1.1", "power_w": 2.0}],
                                       "transient": {"duration_s": t_end, "step_s": dt, "max_step_s": dt}})
        errs.append(abs(r.transient["source_mean_c"][0][-1] - 25.0 - exact(t_end)[0]))
    validation.check("slab, 2 W step, 60 s", "error ratio for halved steps", errs[1] / errs[0], 0.5, rel=0.05,
                     source="first-order convergence of backward Euler")
    assert errs[2] / errs[1] == pytest.approx(0.5, rel=0.05)


def test_diffusion_into_a_strip(tmp_path, validation):
    """Copper strip with an end held at a fixed temperature, no losses: the temperature
    spreads as the semi-infinite solution erfc(x / 2 sqrt(alpha t))."""
    g_ = synthetic.fin(tmp_path)
    b = kairoq.open(tmp_path)
    core = next(l["name"] for l in b.stackup if l["type"] == "dielectric")
    t0, ta, t_end = 75.0, 25.0, 1.0
    r = kairoq.solve_thermal(b, {
        "ambient_c": ta, **faces(h=0.0), "mesh": {"max_edge_mm": 0.25}, "heat_capacity": {"mask": 0.0},
        "layers": {core: {"k_xy": 1e-9, "k_z": 1e-3, "heat_capacity": 1e-3}},
        "sources": [{"name": "HOT", "pad": "J1.1", "temperature_c": t0}],
        "transient": {"duration_s": t_end, "step_s": 1e-3, "max_step_s": 5e-3, "snapshots_s": [t_end]}})
    alpha = K_CU / C_CU
    f = r.layer(1)
    x = f["points_mm"][:, 0] - g_["pad_m"] * 1e3
    tr = r.transient
    temp = r.snapshot(len(tr["snapshot_time_s"]) - 1, 1)
    assert tr["snapshot_time_s"][-1] == pytest.approx(t_end)
    for xm in (2.0, 5.0, 10.0, 20.0):
        sel = np.abs(x - xm) < 0.13
        ref = ta + (t0 - ta) * math.erfc(xm * 1e-3 / (2 * math.sqrt(alpha * t_end)))
        validation.check("copper strip, end at 75 C, t = 1 s", f"rise {xm:g} mm from the end", float(temp[sel].mean()) - ta, ref - ta,
                         rel=0.015, unit="K", source="erfc, semi-infinite")


def test_steady_limit(tmp_path, validation):
    """Long after a step the transient reaches the steady solution."""
    synthetic.fin(tmp_path)
    b = kairoq.open(tmp_path)
    base = {**faces(), "sources": [{"name": "HOT", "pad": "J1.1", "temperature_c": 75.0}]}
    steady = kairoq.solve_thermal(b, base)
    r = kairoq.solve_thermal(b, {**base, "transient": {"duration_s": 4000.0, "step_s": 1.0, "max_step_s": 200.0}})
    validation.check("fin, step to 75 C, 4000 s", "heat from the fixed end against steady", r.sources["HOT"]["power_w"],
                     steady.sources["HOT"]["power_w"], rel=1e-4, unit="W", source="steady solution")
    assert np.max(np.abs(r.layer(1)["temperature_c"] - steady.layer(1)["temperature_c"])) < 2e-3


def test_balance_with_sink_junction_and_radiation(tmp_path, validation):
    """Heat in = heat out + rate of change of the stored heat, at every step, on a board with a
    pad source over a plane, a partial heat sink, a junction and radiation: exact with steps
    solved to the solver tolerance, within the step tolerance by default."""
    synthetic.source_board(tmp_path, pad=(13.0, 13.0, 17.0, 17.0))
    b = kairoq.open(tmp_path)

    def run(step_tolerance):
        return kairoq.solve_thermal(b, {
            **faces(eps=0.9),
            "sources": [{"name": "H", "pad": "H1.1", "power_schedule": [[0, 0], [1, 3.0], [30, 3.0], [31, 0.0]],
                         "junction_k_per_w": 5.0, "heat_capacity_j_per_k": 0.2}],
            "heatsinks": [{"name": "S", "face": "bottom", "area": [{"rect_mm": [5, 5, 25, 25]}],
                           "tim": {"thickness_mm": 0.3, "k": 2.0}, "resistance_k_per_w": 4.0, "heat_capacity_j_per_k": 15.0}],
            "transient": {"duration_s": 60.0, "step_s": 0.25, "max_step_s": 4.0, "step_tolerance": step_tolerance}})

    def imbalance(r):
        tr = r.transient
        return float(np.max(np.abs((tr["heat_in_w"][1:] - tr["heat_out_w"][1:] - tr["heat_stored_w"][1:]) / 3.0)))
    r = run(0.0)
    validation.check("30 x 30 mm board, 3 W pulse, sink and radiation", "largest (in - out - stored) / peak power",
                     imbalance(r), 0.0, abs_=1e-6, source="energy balance of each step, steps solved to the solver tolerance")
    validation.check("30 x 30 mm board, 3 W pulse, sink and radiation", "the same with steps to 1e-4", imbalance(run(1e-4)), 0.0,
                     abs_=1e-3, source="energy balance of each step, default step tolerance")
    tr = r.transient
    j = tr["source_junction_c"][0]
    assert np.argmax(j) > 0 and tr["time_s"][np.argmax(j)] <= 31.0 + 1e-9  # heats while powered, then cools
    assert j[-1] < j.max()


def test_coupled_strip_balance_and_limit(tmp_path, validation):
    """5 A through a strip from ambient: Joule heat with rho(T) renewed during the transient;
    the balance closes and the end state approaches the coupled steady solution."""
    synthetic.fin(tmp_path, width=5.0, length=100.0)
    b = kairoq.open(tmp_path)
    dc = {"ports": [{"name": "A", "pad": "J1.1", "current_a": 5.0}, {"name": "B", "pad": "J2.1", "current_a": -5.0}]}
    steady = kairoq.solve_thermal(b, {**faces()}, dc=dc)
    r = kairoq.solve_thermal(b, {**faces(), "transient": {"duration_s": 5000.0, "step_s": 1.0, "max_step_s": 250.0,
                                                            "dc_update_k": 0.2, "step_tolerance": 0.0}}, dc=dc)
    tr = r.transient
    assert tr["dc_solutions"] > 3
    imbalance = (tr["heat_in_w"][1:] - tr["heat_out_w"][1:] - tr["heat_stored_w"][1:]) / tr["heat_in_w"][1:]
    validation.check("strip 5 x 100 mm, 5 A", "largest relative (in - out - stored)", float(np.max(np.abs(imbalance))), 0.0, abs_=1e-6,
                     source="energy balance of each step")
    validation.check("strip 5 x 100 mm, 5 A, 5000 s", "maximum temperature against steady", r.max_c, steady.max_c, rel=2e-3,
                     unit="C", source="coupled steady solution")


def test_transient_settings_are_checked(tmp_path):
    synthetic.slab(tmp_path)
    b = kairoq.open(tmp_path)
    with pytest.raises(RuntimeError, match="duration_s"):
        kairoq.solve_thermal(b, {"sources": [{"name": "H", "pad": "H1.1", "power_w": 1.0}], "transient": {"step_s": 1.0}})
    with pytest.raises(RuntimeError, match="must not decrease"):
        kairoq.solve_thermal(b, {"sources": [{"name": "H", "pad": "H1.1", "power_schedule": [[1, 1.0], [0, 2.0]]}],
                                   "transient": {"duration_s": 1.0}})
