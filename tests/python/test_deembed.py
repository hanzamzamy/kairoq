"""De-embedding: lines removed from ports against the shorter line, the change of
reference against the matched line, and the half of a symmetric 2x-thru removed from both
ports of a fixture-DUT-fixture cascade, all against closed forms in numpy.
"""

from __future__ import annotations

import subprocess

import numpy as np
import pytest

import kairoq
from conftest import cli_path

C0 = 299792458.0
F = np.linspace(0.1e9, 40e9, 400)


def line(zc, eeff, length, loss_db_1ghz=0.0):
    """ABCD of a quasi-TEM line per frequency (attenuation scaling with sqrt f)."""
    alpha = loss_db_1ghz / (20 / np.log(10)) * np.sqrt(F / 1e9)
    g = (alpha + 2j * np.pi * F * np.sqrt(eeff) / C0) * length
    ch, sh = np.cosh(g), np.sinh(g)
    return np.stack([np.stack([ch, zc * sh], -1), np.stack([sh / zc, ch], -1)], -2)


def shunt_c(c):
    y = 2j * np.pi * F * c
    one, zero = np.ones_like(F), np.zeros_like(F)
    return np.stack([np.stack([one, zero], -1), np.stack([y + 0j, one + 0j], -1)], -2)


def s_of(abcd, z0=50.0):
    a, b, c, d = abcd[:, 0, 0], abcd[:, 0, 1], abcd[:, 1, 0], abcd[:, 1, 1]
    den = a + b / z0 + c * z0 + d
    return np.stack([np.stack([(a + b / z0 - c * z0 - d) / den, 2 * (a * d - b * c) / den], -1),
                     np.stack([2 / den, (-a + b / z0 - c * z0 + d) / den], -1)], -2)


def write(path, s, z0=50.0):
    with open(path, "w") as f:
        f.write(f"# Hz S RI R {z0}\n")
        for fi, m in zip(F, s):
            vals = [m[0, 0], m[1, 0], m[0, 1], m[1, 1]]
            f.write(f"{fi:.9g} " + " ".join(f"{v.real:.15g} {v.imag:.15g}" for v in vals) + "\n")
    return path


def worst(result, s):
    got = np.array(result["s"]).reshape(-1, 2, 2)
    return float(np.max(np.abs(got - s)))


def test_line_removed(tmp_path, validation):
    """A lossy 60 ohm line of 30 mm, 10 mm removed at port 1 and 5 mm at port 2: the line
    of 15 mm."""
    seg = dict(z_ohm=60, epsilon_eff=3.2, loss_db_per_m_1ghz=4.0)
    write(tmp_path / "l30.s2p", s_of(line(60, 3.2, 30e-3, 4.0)))
    r = kairoq.deembed({"input": "l30.s2p", "output": "l15.s2p",
                          "operations": [{"line": dict(port=1, length_mm=-10, **seg)}, {"line": dict(port=2, length_mm=-5, **seg)}]},
                         base_dir=tmp_path)
    validation.check("lossy line 30 mm, 15 mm removed", "largest S difference", worst(r, s_of(line(60, 3.2, 15e-3, 4.0))), 0.0, abs_=1e-9,
                     source="the 15 mm line")
    assert (tmp_path / "l15.s2p").exists() and len(r["steps"]) == 2


def test_renormalise(tmp_path, validation):
    """A lossless 60 ohm line referred to 60 ohm reflects nothing; back to 50 ohm it is the
    original."""
    s = s_of(line(60, 3.2, 30e-3))
    write(tmp_path / "l.s2p", s)
    r = kairoq.deembed({"input": "l.s2p", "operations": [{"renormalise_ohm": 60}]}, base_dir=tmp_path)
    s11 = max(abs(m[0]) for m in r["s"])
    validation.check("lossless 60 ohm line referred to 60 ohm", "largest |S11|", s11, 0.0, abs_=1e-12, source="matched line")
    back = kairoq.deembed({"input": "l.s2p", "operations": [{"renormalise_ohm": 60}, {"renormalise_ohm": 50}]}, base_dir=tmp_path)
    assert worst(back, s) < 1e-12 and back["reference_ohm"] == 50


def test_2x_thru(tmp_path, validation):
    """Fixture: 4 mm of 40 ohm, 0.2 pF, 4 mm of 40 ohm (symmetric, several wavelengths at
    40 GHz). The fixture-DUT-fixture cascade with the 2x-thru de-embedded at both ports is
    the DUT, a 70 ohm line of 12 mm with loss."""
    fix = line(40, 4.0, 4e-3) @ shunt_c(0.2e-12) @ line(40, 4.0, 4e-3)
    dut = line(70, 3.0, 12e-3, 6.0)
    write(tmp_path / "thru.s2p", s_of(fix @ fix))
    write(tmp_path / "meas.s2p", s_of(fix @ dut @ fix))
    r = kairoq.deembed({"input": "meas.s2p", "operations": [{"thru": {"port": 1, "file": "thru.s2p"}}, {"thru": {"port": 2, "file": "thru.s2p"}}]},
                         base_dir=tmp_path)
    validation.check("symmetric fixture with 0.2 pF, 2x-thru removed at both ports", "largest S difference", worst(r, s_of(dut)), 0.0, abs_=1e-6,
                     source="the DUT's closed form")


def test_errors(tmp_path):
    write(tmp_path / "l.s2p", s_of(line(60, 3.2, 30e-3)))
    with pytest.raises(Exception, match="port 3"):
        kairoq.deembed({"input": "l.s2p", "operations": [{"line": {"port": 3, "length_mm": 1, "epsilon_eff": 2}}]}, base_dir=tmp_path)
    with pytest.raises(Exception, match="unknown operation"):
        kairoq.deembed({"input": "l.s2p", "operations": [{"shift": 1}]}, base_dir=tmp_path)


@pytest.mark.skipif(cli_path() is None, reason="kairoq command not built")
def test_cli(tmp_path):
    write(tmp_path / "l.s2p", s_of(line(60, 3.2, 30e-3)))
    (tmp_path / "spec.yaml").write_text("input: l.s2p\noutput: out.s2p\noperations:\n  - renormalise_ohm: 60\n")
    r = subprocess.run([cli_path(), "deembed", str(tmp_path / "spec.yaml")], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert "reference changed to 60 ohm" in r.stdout
    assert kairoq.read_touchstone(tmp_path / "out.s2p")["reference_ohm"] == 60
