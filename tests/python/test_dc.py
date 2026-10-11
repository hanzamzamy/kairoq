"""DC solver: analytical cases, conservation laws and the cross-check of importers."""

from __future__ import annotations

import math
import subprocess

import pytest

import kairoq
import synthetic
from conftest import ROOT, TESTDATA, cli_path

SIGMA = 5.8e7  # S/m, default copper and plating conductivity
T_CU = 35e-6  # m
GS = SIGMA * T_CU  # sheet conductance, S


def solve(board, ports, **settings):
    return kairoq.solve_dc(board, {"ports": ports, **settings})


@pytest.fixture(scope="module")
def strip(tmp_path_factory):
    d = tmp_path_factory.mktemp("strip")
    synthetic.strip(d, length=50.0, width=1.0)
    return kairoq.open(d)


def test_strip_resistance_is_exact(strip, validation):
    # Between the pad edges the track is a rectangle with equipotential ends.
    r = solve(strip, [{"name": "A", "pad": "J1.1", "voltage_v": 1.0}, {"name": "B", "pad": "J2.1", "voltage_v": 0.0}])
    exact = 48e-3 / (GS * 1e-3)
    validation.check("strip 48 x 1 mm", "resistance", 1.0 / r.ports["A"]["current_a"], exact, rel=1e-9, unit="ohm",
                     source="R = L / (sigma t w)")
    assert r.ports["B"]["current_a"] == pytest.approx(-r.ports["A"]["current_a"], rel=1e-12)
    assert r.power_ports_w == pytest.approx(r.power_copper_w, rel=1e-9)


def test_strip_with_current_ports_only(strip, validation):
    r = solve(strip, [{"name": "A", "pad": "J1.1", "current_a": 2.0}, {"name": "B", "pad": "J2.1", "current_a": -2.0}])
    exact = 48e-3 / (GS * 1e-3)
    validation.check("strip, current ports", "resistance", (r.ports["A"]["voltage_v"] - r.ports["B"]["voltage_v"]) / 2.0,
                     exact, rel=1e-9, unit="ohm", source="R = L / (sigma t w)")
    assert any("relative to port A" in w for w in r.warnings)
    # Currents that do not balance on floating copper are an error.
    with pytest.raises(RuntimeError, match="do not sum to zero"):
        solve(strip, [{"name": "A", "pad": "J1.1", "current_a": 2.0}, {"name": "B", "pad": "J2.1", "current_a": -1.0}])


def test_copper_temperature(strip, validation):
    ports = [{"name": "A", "pad": "J1.1", "voltage_v": 1.0}, {"name": "B", "pad": "J2.1", "voltage_v": 0.0}]
    r20 = solve(strip, ports)
    r75 = solve(strip, ports, temperature_c=75.0)
    validation.check("strip at 75 C", "R(75 C) / R(20 C)", r20.ports["A"]["current_a"] / r75.ports["A"]["current_a"],
                     1 + 0.00393 * 55, rel=1e-9, source="rho(T) = rho0 (1 + 0.00393 (T - 20))")


def test_settings_are_checked(strip):
    with pytest.raises(RuntimeError, match="unknown key"):
        kairoq.solve_dc(strip, {"ports": [], "plating_um": 18})
    with pytest.raises(RuntimeError, match="either voltage_v or current_a"):
        solve(strip, [{"name": "A", "pad": "J1.1"}])
    with pytest.raises(RuntimeError, match="no pad"):
        solve(strip, [{"name": "A", "pad": "J9.1", "voltage_v": 0}])


def spreading_exact(g):
    spread = 2 * math.log(g["r_inner_port"] / g["r_via"]) / (2 * math.pi * GS)
    tp, d = 18e-6, 2 * g["r_via"] * 1e-3
    barrel = g["barrel_length"] * 1e-3 / (SIGMA * math.pi * tp * (d + tp))
    return spread, barrel


@pytest.mark.parametrize("target", [0.1, 0.03])
def test_via_spreading(tmp_path, validation, target):
    """Radial spreading on both layers through one plated via (analytical)."""
    g = synthetic.via_spreading(tmp_path)
    b = kairoq.open(tmp_path)
    r = solve(b, [{"name": "A", "pad": "P1.1", "voltage_v": 1e-3}, {"name": "B", "pad": "P2.1", "voltage_v": 0.0}],
              target_error=target)
    spread, barrel = spreading_exact(g)
    validation.check("via spreading", "barrel resistance", r.barrels[0]["resistance_ohm"], barrel, rel=1e-9, unit="ohm",
                     source="L / (sigma pi t (d + t))")
    # The energy-norm estimate e bounds the relative error of a two-port resistance by e^2.
    bound = r.estimated_error ** 2
    validation.check(f"via spreading, target {target}", "resistance", 1e-3 / r.ports["A"]["current_a"], spread + barrel,
                     rel=bound, unit="ohm", source="ln(r2/r1) / (pi sigma t) + barrel")
    assert r.estimated_error <= target


@pytest.mark.parametrize("target", [0.1, 0.03])
def test_neck_constriction(tmp_path, validation, target):
    """Strip with a thin wall and an aperture: conformal-mapping constriction resistance."""
    g = synthetic.neck(tmp_path, aperture=0.5, wall=0.005)
    b = kairoq.open(tmp_path)
    r = solve(b, [{"name": "A", "pad": "J1.1", "voltage_v": 1e-3}, {"name": "B", "pad": "J2.1", "voltage_v": 0.0}],
              target_error=target)
    w, a = g["width"], g["aperture"]
    base = g["free_length"] / w / GS
    constriction = 2 / (math.pi * GS) * math.log(1 / math.sin(math.pi * a / (2 * w)))
    wall = g["wall"] * (1 / a - 1 / w) / GS  # first-order length of the aperture
    validation.check(f"neck, target {target}", "resistance", 1e-3 / r.ports["A"]["current_a"], base + constriction + wall,
                     rel=max(r.estimated_error ** 2, 2e-3), unit="ohm", source="Schwarz-Christoffel, zero-thickness wall")


def test_via_chain_plating(tmp_path, validation):
    """The plating thickness changes only the barrels: R(t1) - R(t2) = n (Rb(t1) - Rb(t2))."""
    synthetic.via_chain(tmp_path)
    b = kairoq.open(tmp_path)
    ports = [{"name": "A", "at_mm": [5.0, 5.0], "layer": 1, "diameter_mm": 0.6, "voltage_v": 1e-3},
             {"name": "B", "at_mm": [23.0, 5.0], "layer": 1, "diameter_mm": 0.6, "voltage_v": 0.0}]
    res = {}
    for tp in (18.0, 25.0):
        r = solve(b, ports, plating_thickness_um=tp, max_passes=1)  # same mesh for both
        res[tp] = (1e-3 / r.ports["A"]["current_a"], r.barrels)
    length = (1.5 + 0.035) * 1e-3

    def rb(t_um):
        t = t_um * 1e-6
        return length / (SIGMA * math.pi * t * (0.3e-3 + t))

    carrying = [h for h in res[18.0][1] if abs(h["current_a"]) > 1e-9]
    assert len(carrying) == 8
    validation.check("via chain, 8 vias in series", "R(18 um) - R(25 um)", res[18.0][0] - res[25.0][0],
                     8 * (rb(18) - rb(25)), rel=1e-6, unit="ohm", source="barrel formula")


@pytest.fixture(scope="module")
def cm5():
    return kairoq.open(TESTDATA / "kicad" / "cm5_minima" / "gerber")


CM5_PORTS = [
    {"name": "VIN", "component": "J101", "net": "+5V", "voltage_v": 5.0},
    {"name": "LOAD", "component": "Module301", "net": "+5V", "current_a": -3.0},
    {"name": "GND_IN", "component": "J101", "net": "GND", "voltage_v": 0.0},
    {"name": "GND_RET", "component": "Module302", "net": "GND", "current_a": 3.0},
]


@pytest.fixture(scope="module")
def cm5_solution(cm5):
    return solve(cm5, CM5_PORTS)


def test_conservation_on_a_real_board(cm5_solution, validation):
    r = cm5_solution
    # Exact up to the linear solver's tolerance (relative residual 1e-10 of the right-hand side).
    validation.check("cm5_minima 5 V, 3 A", "(port power - Joule heat) / port power",
                     (r.power_ports_w - r.power_copper_w - r.power_barrels_w) / r.power_ports_w, 0.0, abs_=1e-4,
                     source="energy balance")
    total = sum(abs(p["current_a"]) for p in r.ports.values())
    validation.check("cm5_minima 5 V, 3 A", "sum of port currents / sum of magnitudes",
                     sum(p["current_a"] for p in r.ports.values()) / total, 0.0, abs_=1e-6, source="charge conservation")
    assert r.estimated_error <= 0.1
    assert 0.0 < 5.0 - r.ports["LOAD"]["voltage_v"] < 0.1


def test_importers_agree_on_dc(cm5_solution, validation):
    """The same board from IPC-2581 gives the same drops (conductivity set equal: the
    IPC-2581 stackup states one, the Gerber job file does not)."""
    x = kairoq.open(TESTDATA / "kicad" / "cm5_minima" / "cm5_minima.ipc2581.zip")
    g = kairoq.open(TESTDATA / "kicad" / "cm5_minima" / "gerber")
    rg = solve(g, CM5_PORTS, copper_conductivity_s_per_m=5.8e7)
    rx = solve(x, CM5_PORTS, copper_conductivity_s_per_m=5.8e7)
    drop_g = 5.0 - rg.ports["LOAD"]["voltage_v"]
    drop_x = 5.0 - rx.ports["LOAD"]["voltage_v"]
    validation.check("cm5_minima, IPC-2581 against Gerber", "+5V drop", drop_x, drop_g, rel=0.01, unit="V",
                     source="Gerber X2 import")


@pytest.mark.skipif(cli_path() is None, reason="kairoq command not built")
def test_cli_dc_example():
    r = subprocess.run([cli_path(), "dc", str(ROOT / "examples" / "cm5_minima_5v.yaml")], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert "LOAD" in r.stdout and "estimated error" in r.stdout
