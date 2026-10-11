"""Heat sinks through a thermal interface, and component power from a table."""

from __future__ import annotations

import subprocess

import numpy as np
import pytest

import kairoq
import synthetic
from conftest import cli_path

K_Z, K_MASK = 0.3, 0.25


def faces(h=10.0, eps=0.0):
    return {"top": {"convection_w_m2k": h, "emissivity": eps}, "bottom": {"convection_w_m2k": h, "emissivity": eps}}


@pytest.fixture(scope="module")
def slab(tmp_path_factory):
    d = tmp_path_factory.mktemp("slab")
    return synthetic.slab(d), kairoq.open(d)


def series(g_, sink_r, ta, p, t_sink=None, rsa=None):
    """Slab with 2 W on the top: convection on the top, the whole bottom on a sink through
    mask and interface; the sink fixed or with a resistance to ambient. Returns the top,
    bottom and sink temperatures."""
    a = g_["area_m2"]
    g_top = a / (g_["mask_m"] / K_MASK + 0.1)
    g_mid = a * K_Z / g_["dielectric_m"]
    g_sink = a / (g_["mask_m"] / K_MASK + sink_r)
    if t_sink is not None:
        G = np.array([[g_top + g_mid, -g_mid], [-g_mid, g_mid + g_sink]])
        t = np.linalg.solve(G, [p + g_top * ta, g_sink * t_sink])
        return t[0], t[1], t_sink
    G = np.array([[g_top + g_mid, -g_mid, 0], [-g_mid, g_mid + g_sink, -g_sink], [0, -g_sink, g_sink + 1 / rsa]])
    return tuple(np.linalg.solve(G, [p + g_top * ta, 0.0, ta / rsa]))


def test_sink_at_fixed_temperature(slab, validation):
    g_, b = slab
    tim = 0.5e-3 / 3.0
    r = kairoq.solve_thermal(b, {**faces(), "sources": [{"name": "H", "pad": "H1.1", "power_w": 2.0}],
                                   "heatsinks": [{"name": "PLATE", "face": "bottom", "tim": {"thickness_mm": 0.5, "k": 3.0},
                                                  "temperature_c": 40.0}]})
    t1, t2, _ = series(g_, tim, 25.0, 2.0, t_sink=40.0)
    validation.check("slab on a plate at 40 C through 0.5 mm of 3 W/(m K)", "top temperature",
                     float(r.layer(1)["temperature_c"].mean()), t1, rel=1e-9, unit="C", source="series resistances")
    validation.check("slab on a plate at 40 C through 0.5 mm of 3 W/(m K)", "bottom temperature",
                     float(r.layer(2)["temperature_c"].mean()), t2, rel=1e-9, unit="C", source="series resistances")
    h = r.heatsinks["PLATE"]
    assert h["temperature_set"] and h["temperature_c"] == 40.0
    assert r.heat_bottom_w == pytest.approx(0.0, abs=1e-12)  # covered by the sink
    assert r.heat_top_w + h["heat_w"] == pytest.approx(2.0, rel=1e-9)


def test_sink_body_with_resistance(slab, validation):
    g_, b = slab
    r = kairoq.solve_thermal(b, {**faces(), "sources": [{"name": "H", "pad": "H1.1", "power_w": 2.0}],
                                   "heatsinks": [{"name": "S", "tim": {"resistance_k_mm2_per_w": 50.0}, "resistance_k_per_w": 5.0}]})
    t1, t2, ts = series(g_, 50e-6, 25.0, 2.0, rsa=5.0)
    validation.check("slab on a heat sink of 5 K/W", "sink temperature", r.heatsinks["S"]["temperature_c"], ts, rel=1e-9, unit="C",
                     source="series resistances")
    validation.check("slab on a heat sink of 5 K/W", "top temperature", float(r.layer(1)["temperature_c"].mean()), t1, rel=1e-9,
                     unit="C", source="series resistances")
    assert r.heatsinks["S"]["heat_w"] == pytest.approx((ts - 25.0) / 5.0, rel=1e-8)


def test_partial_sink_balance_and_area(tmp_path, validation):
    """A rectangular sink under part of the bottom face: the heat balance closes, and with
    the board and interface far better conductors than the faces the sink takes heat
    through exactly its area."""
    synthetic.source_board(tmp_path, pad=(13.0, 13.0, 17.0, 17.0))
    b = kairoq.open(tmp_path)
    r = kairoq.solve_thermal(b, {**faces(eps=0.9), "sources": [{"name": "H", "pad": "H1.1", "power_w": 3.0}],
                                   "heatsinks": [{"name": "S", "area": [{"rect_mm": [5, 7, 25, 22]}],
                                                  "tim": {"thickness_mm": 0.3, "k": 2.0}, "resistance_k_per_w": 4.0}]})
    out = r.heat_top_w + r.heat_bottom_w + r.heat_sinks_w
    validation.check("30 x 30 mm board, sink under 20 x 15 mm", "(heat out - heat in) / heat in", (out - 3.0) / 3.0, 0.0, abs_=1e-6,
                     source="energy balance")
    # Isothermal board: k huge in the dielectric, no faces; the sink then takes all the heat
    # through its area, and the board's rise over the sink is P r / A.
    core = next(l["name"] for l in b.stackup if l["type"] == "dielectric")
    r2 = kairoq.solve_thermal(b, {**faces(h=0.0), "layers": {core: {"k_xy": 1e6, "k_z": 1e6}},
                                    "sources": [{"name": "H", "pad": "H1.1", "power_w": 3.0}],
                                    "heatsinks": [{"name": "S", "area": [{"rect_mm": [5, 7, 25, 22]}],
                                                   "tim": {"thickness_mm": 0.3, "k": 2.0}, "temperature_c": 25.0}]})
    rise = float(r2.layer(2)["temperature_c"].mean()) - 25.0
    area = 20e-3 * 15e-3
    validation.check("isothermal board on a 20 x 15 mm sink", "rise over the sink", rise, 3.0 * (10e-6 / K_MASK + 0.3e-3 / 2.0) / area,
                     rel=1e-4, unit="K", source="P (mask + interface) / area")


def test_sink_settings_are_checked(slab):
    _, b = slab
    src = [{"name": "H", "pad": "H1.1", "power_w": 1.0}]
    with pytest.raises(RuntimeError, match="either temperature_c or resistance_k_per_w"):
        kairoq.solve_thermal(b, {"sources": src, "heatsinks": [{"name": "S", "temperature_c": 30, "resistance_k_per_w": 2}]})
    with pytest.raises(RuntimeError, match="outside the board"):
        kairoq.solve_thermal(b, {"sources": src, "heatsinks": [{"name": "S", "area": [{"rect_mm": [50, 50, 60, 60]}],
                                                                  "temperature_c": 30}]})
    with pytest.raises(RuntimeError, match="defined twice"):
        kairoq.solve_thermal(b, {"sources": src, "heatsinks": [{"name": "S", "temperature_c": 30},
                                                                 {"name": "S", "face": "top", "temperature_c": 30}]})


# ---------------------------------------------------------------------- component power


@pytest.fixture(scope="module")
def strip(tmp_path_factory):
    d = tmp_path_factory.mktemp("strip")
    synthetic.strip(d)
    return kairoq.open(d), d


def test_component_power_table(strip, tmp_path):
    b, _ = strip
    (tmp_path / "power.csv").write_text('Reference,Value,Power (W),Rjc\n"J1, J2",conn,250 mW,8\nX9,unused,,\n')
    r = kairoq.solve_thermal(b, {**faces(), "component_power": {"file": "power.csv", "junction_column": "Rjc"}},
                               base_dir=tmp_path)
    assert set(r.sources) == {"J1", "J2"}
    assert r.heat_sources_w == pytest.approx(0.5)
    assert r.sources["J1"]["junction_c"] == pytest.approx(r.sources["J1"]["mean_c"] + 0.25 * 8, rel=1e-9)
    assert r.heat_top_w + r.heat_bottom_w == pytest.approx(0.5, rel=1e-6)


def test_component_power_table_checks(strip, tmp_path):
    b, _ = strip
    (tmp_path / "semi.csv").write_text("Designator;Power\nJ1;0.1\nX9;0.2\n")
    with pytest.raises(RuntimeError, match="1 components not on the board .X9."):
        kairoq.solve_thermal(b, {**faces(), "component_power": {"file": "semi.csv"}}, base_dir=tmp_path)
    r = kairoq.solve_thermal(b, {**faces(), "component_power": {"file": "semi.csv", "ignore_unknown": True}}, base_dir=tmp_path)
    assert set(r.sources) == {"J1"} and any("X9" in w for w in r.warnings)
    (tmp_path / "twice.csv").write_text("Ref,P\nJ1,0.1\nJ1,0.2\n")
    with pytest.raises(RuntimeError, match="listed twice"):
        kairoq.solve_thermal(b, {**faces(), "component_power": {"file": "twice.csv"}}, base_dir=tmp_path)
    (tmp_path / "bad.csv").write_text("Ref,P\nJ1,0.1 A\n")
    with pytest.raises(RuntimeError, match="unknown unit"):
        kairoq.solve_thermal(b, {**faces(), "component_power": {"file": "bad.csv"}}, base_dir=tmp_path)
    # A source in the project file takes precedence.
    (tmp_path / "both.csv").write_text("Ref,P\nJ1,0.1\nJ2,0.2\n")
    r = kairoq.solve_thermal(b, {**faces(), "sources": [{"name": "J1", "component": "J1", "power_w": 1.0}],
                                   "component_power": {"file": "both.csv"}}, base_dir=tmp_path)
    assert r.sources["J1"]["power_w"] == 1.0 and r.sources["J2"]["power_w"] == pytest.approx(0.2)


@pytest.mark.skipif(cli_path() is None, reason="kairoq command not built")
def test_cli_project_with_table_sink_and_transient(tmp_path):
    """A project file with a component power table, a heat sink and a transient, run from
    the command line."""
    synthetic.strip(tmp_path / "board")
    (tmp_path / "power.csv").write_text("Ref,Power\nJ1,0.3\nJ2,0.2\n")
    (tmp_path / "project.yaml").write_text(
        "board: board\n"
        "thermal:\n"
        "  component_power: {file: power.csv}\n"
        "  heatsinks: [{name: PLATE, face: bottom, tim: {thickness_mm: 0.5, k: 3}, resistance_k_per_w: 8, heat_capacity_j_per_k: 5}]\n"
        "  transient: {duration_s: 120, step_s: 1, max_step_s: 20}\n")
    r = subprocess.run([cli_path(), "thermal", str(tmp_path / "project.yaml")], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert "heat sinks" in r.stdout and "PLATE" in r.stdout
    assert "transient    120 s" in r.stdout
    assert "J1" in r.stdout and "J2" in r.stdout
