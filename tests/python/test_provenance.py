"""Provenance of results: digests of the inputs, the build, and their records (docs/provenance.md)."""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import tomllib

import pytest

import kairoq
from conftest import ROOT, TESTDATA, cli_path

CLI = cli_path()


def sha256(path) -> str:
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def directory_sha256(path) -> str:
    """The digest of a directory as docs/provenance.md defines it, written independently."""
    lines = []
    for base, _, files in os.walk(path):
        for f in files:
            full = os.path.join(base, f)
            lines.append(os.path.relpath(full, path).replace(os.sep, "/") + " " + sha256(full) + "\n")
    return hashlib.sha256("".join(sorted(lines)).encode()).hexdigest()


def test_build_info():
    b = kairoq.build_info()
    version = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]["version"]
    assert kairoq.__version__ == b["version"] == version
    assert b["commit"] == "unknown" or len(b["commit"].split("+")[0]) == 40


def test_board_sources_survive_the_model(tmp_path):
    gerber = TESTDATA / "kicad" / "ecc83" / "gerber"
    board = kairoq.open(gerber)
    d = board.data["source"]["digests"]
    assert d == [{"path": "gerber", "sha256": directory_sha256(gerber), "bytes": d[0]["bytes"], "directory": True}]
    board.save(tmp_path / "model.zip")
    assert kairoq.open(tmp_path / "model.zip").data["source"]["digests"] == d  # not the model's own digest


def test_python_provenance(tmp_path):
    board = kairoq.open(TESTDATA / "de10_interface" / "de10_interface.ipc2581.zip")
    project = ROOT / "examples" / "de10_interface.yaml"
    p = kairoq.provenance(board, project)
    names = {i["path"]: i["sha256"] for i in p["inputs"]}
    assert names[str(project)] == sha256(project)
    assert names["../testdata/de10_interface/de10_interface.ipc2581.zip"] == sha256(TESTDATA / "de10_interface" / "de10_interface.ipc2581.zip")
    assert p["board_sources"][0]["sha256"] == sha256(TESTDATA / "de10_interface" / "de10_interface.ipc2581.zip")
    assert p["kairoq"] == kairoq.build_info() and p["time_utc"].endswith("Z")


@pytest.mark.skipif(CLI is None, reason="kairoq command not built")
def test_cli_record(tmp_path):
    project = ROOT / "examples" / "de10_interface.yaml"
    record = tmp_path / "dc.json"
    r = subprocess.run([CLI, "dc", str(project), "--record", str(record)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert r.stdout.splitlines()[-1].startswith("kairoq " + kairoq.__version__ + ", commit ")
    j = json.loads(record.read_text())
    assert j["summary"] in r.stdout and j["command"][1:3] == ["dc", str(project)]
    assert {i["sha256"] for i in j["inputs"]} == {sha256(project), sha256(TESTDATA / "de10_interface" / "de10_interface.ipc2581.zip")}


@pytest.mark.skipif(CLI is None, reason="kairoq command not built")
def test_cli_touchstone_header(tmp_path):
    """A written Touchstone file carries its provenance as comments; the file the run writes
    is not one of its inputs, and the result still reads as Touchstone."""
    (tmp_path / "in.s2p").write_text("# GHz S RI R 50\n1 0.1 0 0.9 0 0.9 0 0.1 0\n2 0.2 0 0.8 0 0.8 0 0.2 0\n")
    (tmp_path / "spec.yaml").write_text("input: in.s2p\noutput: out.s2p\noperations:\n  - renormalise_ohm: 75\n")
    r = subprocess.run([CLI, "deembed", str(tmp_path / "spec.yaml"), "--record", str(tmp_path / "rec.json")], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    text = (tmp_path / "out.s2p").read_text()
    assert text.startswith("! kairoq " + kairoq.__version__)
    assert "! input in.s2p sha256 " + sha256(tmp_path / "in.s2p") in text
    assert [i["path"] for i in json.loads((tmp_path / "rec.json").read_text())["inputs"]] == [str(tmp_path / "spec.yaml"), "in.s2p"]
    t = kairoq.read_touchstone(str(tmp_path / "out.s2p"))
    assert t["reference_ohm"] == 75
