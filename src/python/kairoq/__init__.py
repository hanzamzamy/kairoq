"""PCB analysis: board import, canonical model and 3D viewer.

Lengths are in millimetres. Geometry is exchanged as OGC Well-Known Binary
MultiPolygons, readable with ``shapely.wkb.loads``.
"""

from __future__ import annotations

import builtins
import json
import os
from typing import Iterable

from . import _core

__all__ = ["Board", "DcResult", "EmiResult", "SiResult", "ThermalResult", "analyse_channels", "analyse_emi", "analyse_si", "cross_section", "deembed", "export_openems", "export_openems_channel", "eye", "open", "patch_antenna", "pdn", "serve", "solve_dc",
           "ibis_switching", "read_touchstone", "simulate", "solve_thermal", "MODEL_VERSION"]

MODEL_VERSION = _core.MODEL_VERSION
builtins_open = builtins.open


class Board:
    """An imported board: the canonical model in memory."""

    def __init__(self, core: _core.Board):
        self._core = core
        self._data = json.loads(core.json())

    # Board data ------------------------------------------------------------

    @property
    def data(self) -> dict:
        """board.json content without geometry references (see docs/canonical-model.md)."""
        return self._data

    @property
    def name(self) -> str:
        return self._data["name"]

    @property
    def nets(self) -> list[str]:
        """Net names by id; id 0 is "no net"."""
        return self._data["nets"]

    @property
    def net_source(self) -> str:
        return self._data["net_source"]

    @property
    def stackup(self) -> list[dict]:
        return self._data["stackup"]["layers"]

    @property
    def copper_layers(self) -> list[dict]:
        return self._data["copper_layers"]

    @property
    def pads(self) -> list[dict]:
        return self._data["pads"]

    @property
    def drills(self) -> list[dict]:
        return self._data["drills"]

    @property
    def tracks(self) -> list[dict]:
        """Drawn tracks: layer, net, width_mm, start_mm, end_mm (centre_mm and ccw for arcs)."""
        return self.data["tracks"]

    @property
    def components(self) -> list[dict]:
        return self._data["components"]

    @property
    def warnings(self) -> list[str]:
        return self._data["warnings"]

    # Geometry --------------------------------------------------------------

    def copper_wkb(self, layer: int) -> dict[int, bytes]:
        """Copper of a layer (1 = top) per net id, as WKB MultiPolygons in mm."""
        return dict(self._core.copper_wkb(layer))

    def copper_area_mm2(self, layer: int, net: int | None = None) -> float:
        return self._core.copper_area_mm2(layer, net)

    def outline_wkb(self) -> bytes:
        return self._core.outline_wkb()

    def mask_wkb(self, side: str) -> bytes | None:
        """Solder mask openings of "top" or "bottom"."""
        return self._core.mask_wkb(side)

    def pad_wkb(self, index: int) -> bytes:
        return self._core.pad_wkb(index)

    def mesh_keys(self) -> list[str]:
        return self._core.mesh_keys()

    def mesh(self, key: str) -> dict:
        """Display mesh: positions (n, 3) in mm relative to origin_mm, triangles (m, 3), ids (n,)."""
        return self._core.mesh(key)

    # Output ----------------------------------------------------------------

    def save(self, path: str | os.PathLike) -> None:
        """Write the canonical model to a directory, or to a .zip archive."""
        self._core.save(os.fspath(path))

    def summary(self) -> str:
        return self._core.summary()

    def __repr__(self) -> str:
        return f"<kairoq.Board {self.name!r}: {len(self.copper_layers)} copper layers, {len(self.nets) - 1} nets>"


def open(
    inputs: str | os.PathLike | Iterable[str | os.PathLike],
    *,
    arc_tolerance_um: float = 1.0,
    stackup: str | os.PathLike | None = None,
    crop_mm: tuple[float, float, float, float] | None = None,
) -> Board:
    """Open a canonical model, an IPC-2581 file, or a set of fabrication files.

    ``inputs`` is a path or a list of paths (files, directories or zip archives).
    ``arc_tolerance_um`` is the largest chord deviation from arcs; ``stackup`` a
    stackup YAML file overriding the job file; ``crop_mm`` (x0, y0, x1, y1) keeps the
    part of the board inside that rectangle, with open edges.
    """
    if isinstance(inputs, (str, os.PathLike)):
        inputs = [inputs]
    paths = [os.fspath(p) for p in inputs]
    return Board(_core.open(paths, arc_tolerance_um, None if stackup is None else os.fspath(stackup),
                            None if crop_mm is None else tuple(float(v) for v in crop_mm)))


class DcResult:
    """DC solution (see docs/solvers/dc.md)."""

    def __init__(self, board: Board, core: _core.DcResult):
        self.board = board
        self._core = core

    def __getattr__(self, name):
        return getattr(self._core, name)

    @property
    def ports(self) -> dict[str, dict]:
        """Port results by name: voltage_v, current_a (into the copper), voltage_set, nodes."""
        return {p["name"]: p for p in self._core.ports}

    def layer(self, index: int) -> dict:
        """Solution on copper layer ``index``: points_mm, triangles, voltage_v (per point),
        sheet_current_a_per_mm, power_w and net (per triangle)."""
        return self._core.layer(index)

    def summary(self) -> str:
        return _core.dc_summary(self.board._core, self._core)


def solve_dc(board: Board, settings: dict | str | os.PathLike) -> DcResult:
    """Solve the DC current distribution.

    ``settings`` is the ``dc`` section of a project file as a dict, or the path of a
    project file. Example::

        solve_dc(board, {"ports": [{"name": "IN", "pad": "J1.1", "voltage_v": 5},
                                    {"name": "LOAD", "pad": "U1.3", "current_a": -2}]})
    """
    if isinstance(settings, dict):
        text, base = json.dumps({"dc": settings}), "."
    else:
        path = os.fspath(settings)
        with builtins_open(path) as f:
            text = f.read()
        base = os.path.dirname(os.path.abspath(path))
    return DcResult(board, _core.solve_dc(board._core, text, base))


class ThermalResult:
    """Thermal solution, steady or transient (see docs/solvers/thermal.md). For a transient
    the fields and source results are those at the end; ``transient`` holds the time series
    and ``snapshot(i, layer)`` the temperature of a layer at ``transient["snapshot_time_s"][i]``."""

    def __init__(self, board: Board, core):
        self.board = board
        self._core = core

    def __getattr__(self, name):
        return getattr(self._core, name)

    @property
    def sources(self) -> dict[str, dict]:
        """Results by source name: power_w, mean_c, max_c, junction_c, temperature_set."""
        return {s["name"]: s for s in self._core.sources}

    @property
    def heatsinks(self) -> dict[str, dict]:
        """Results by heat sink name: heat_w (into the sink), temperature_c, temperature_set."""
        return {h["name"]: h for h in self._core.heatsinks}

    @property
    def dc(self) -> "DcResult | None":
        d = self._core.dc
        return None if d is None else DcResult(self.board, d)

    def layer(self, index: int) -> dict:
        """Temperature on copper layer ``index``: points_mm, triangles, temperature_c (per
        point), copper (per triangle)."""
        return self._core.layer(index)

    def summary(self) -> str:
        return _core.thermal_summary(self.board._core, self._core)


def solve_thermal(board: Board, thermal: dict | str | os.PathLike, dc: dict | None = None,
                  base_dir: str | os.PathLike = ".") -> ThermalResult:
    """Solve the temperature, steady or (with a ``transient`` section) over time, coupled
    with the DC current when ``dc`` has ports.

    ``thermal`` is the ``thermal`` section of a project file as a dict, or the path of a
    project file (whose ``dc`` section then gives the coupling). Files named in a dict
    (component power tables) are relative to ``base_dir``.
    """
    if isinstance(thermal, dict):
        doc = {"thermal": thermal}
        if dc:
            doc["dc"] = dc
        text, base = json.dumps(doc), os.fspath(base_dir)
    else:
        path = os.fspath(thermal)
        with builtins_open(path) as f:
            text = f.read()
        base = os.path.dirname(os.path.abspath(path))
    return ThermalResult(board, _core.solve_thermal(board._core, text, base))


class EmiResult:
    """EMI estimates (see docs/solvers/emi.md): ``loops`` by name, ``return_path`` findings,
    ``coupling`` of switch nodes to plates."""

    def __init__(self, board: Board, core):
        self.board = board
        self._core = core

    def __getattr__(self, name):
        return getattr(self._core, name)

    @property
    def loops(self) -> dict[str, dict]:
        """Loop results by name: inductance_nh, inductance_coarse_nh, moment_m2 (per ampere),
        harmonic_frequency_hz / harmonic_current_a / harmonic_field_v_per_m, sheets."""
        return {l["name"]: l for l in self._core.loops}

    def summary(self) -> str:
        return _core.emi_summary(self.board._core, self._core)


def analyse_emi(board: Board, emi: dict | str | os.PathLike, base_dir: str | os.PathLike = ".") -> EmiResult:
    """Loop inductance and radiated estimate, return-path checks and switch-node coupling.

    ``emi`` is the ``emi`` section of a project file as a dict, or the path of a project file.
    """
    if isinstance(emi, dict):
        text, base = json.dumps({"emi": emi}), os.fspath(base_dir)
    else:
        path = os.fspath(emi)
        with builtins_open(path) as f:
            text = f.read()
        base = os.path.dirname(os.path.abspath(path))
    return EmiResult(board, _core.analyse_emi(board._core, text, base))


def export_openems(board: Board, emi: dict | str | os.PathLike, out_dir: str | os.PathLike | None,
                   loop: str | None = None, base_dir: str | os.PathLike = ".") -> dict:
    """Export a loop of the ``emi`` settings (the first, or the one named) to openEMS: writes
    model.json and simulate.py into ``out_dir`` and returns the model. With ``out_dir`` None
    nothing is written. See docs/solvers/openems.md."""
    if isinstance(emi, dict):
        text, base = json.dumps({"emi": emi}), os.fspath(base_dir)
    else:
        path = os.fspath(emi)
        with builtins_open(path) as f:
            text = f.read()
        base = os.path.dirname(os.path.abspath(path))
    model, summary = _core.export_openems(board._core, text, base, "" if out_dir is None else os.fspath(out_dir), loop or "")
    model = json.loads(model)
    model["summary"] = summary
    return model


class SiResult:
    """Signal integrity of the board's nets (see docs/solvers/si.md): ``nets``, ``samples``
    along the tracks, ``tracks``, ``findings`` and ``cross_sections``; ``result["nets"]``
    and the other keys index the same data."""

    def __init__(self, board: Board, core):
        self.board = board
        self._core = core
        self.data = json.loads(core.json(board._core))

    def __getitem__(self, key):
        return self.data[key]

    def __getattr__(self, name):
        if name in ("nets", "samples", "tracks", "coupling", "findings", "cross_sections", "warnings", "cuts", "seconds"):
            return self.data[name]
        raise AttributeError(name)

    def summary(self) -> str:
        return self._core.summary(self.board._core)


def analyse_si(board: Board, si: dict | str | os.PathLike, base_dir: str | os.PathLike = ".") -> SiResult:
    """Impedance profile, delay and loss of the board's nets, differential impedance of
    pairs, and findings (docs/solvers/si.md). ``si`` is the ``si`` section of a project
    file as a dict, or the path of a project file."""
    if isinstance(si, dict):
        text, base = json.dumps({"si": si}), os.fspath(base_dir)
    else:
        path = os.fspath(si)
        with builtins_open(path) as f:
            text = f.read()
        base = os.path.dirname(os.path.abspath(path))
    return SiResult(board, _core.analyse_si(board._core, text, base))


def analyse_channels(board: Board, si: dict | str | os.PathLike, base_dir: str | os.PathLike = ".") -> list[dict]:
    """S-parameters, mixed-mode parameters and TDR of the channels in ``si["channels"]``
    (docs/solvers/networks.md). Each result has ``frequency_hz``, ``s`` and ``mixed`` (per
    frequency, the matrix row major as complex numbers), ``tdr``, and its Touchstone texts
    ``touchstone`` and ``touchstone_mixed``; the list's summary is in each ``summary``."""
    if isinstance(si, dict):
        text, base = json.dumps({"si": si}), os.fspath(base_dir)
    else:
        path = os.fspath(si)
        with builtins_open(path) as f:
            text = f.read()
        base = os.path.dirname(os.path.abspath(path))
    data, summary, ts, tm = _core.analyse_channels(board._core, text, base)
    out = json.loads(data)
    for r, a, b in zip(out, ts, tm):
        r["s"] = [[complex(*c) for c in m] for m in r["s"]]
        r["mixed"] = [[complex(*c) for c in m] for m in r["mixed"]]
        r["touchstone"], r["touchstone_mixed"], r["summary"] = a, b, summary
    return out


def simulate(board: Board, si: dict | str | os.PathLike, base_dir: str | os.PathLike = ".") -> list[dict]:
    """Waveforms of the simulations in ``si["simulations"]``, each on a channel of
    ``si["channels"]``, with linear or IBIS drivers and receivers
    (docs/solvers/timedomain.md): ``time_s``, ``voltage`` per port, ``metrics``; the
    printable summary in each ``summary``."""
    if isinstance(si, dict):
        text, base = json.dumps({"si": si}), os.fspath(base_dir)
    else:
        path = os.fspath(si)
        with builtins_open(path) as f:
            text = f.read()
        base = os.path.dirname(os.path.abspath(path))
    data, summary = _core.simulate(board._core, text, base)
    out = json.loads(data)
    for r in out:
        r["summary"] = summary
    return out


def export_openems_channel(board: Board, si: dict | str | os.PathLike, out_dir: str | os.PathLike | None,
                           channel: str | None = None, base_dir: str | os.PathLike = ".") -> dict:
    """Export a channel of ``si["channels"]`` (the first, or the one named) to openEMS for
    full-wave S-parameters, with ``si["openems"]`` settings: writes model.json and
    simulate.py into ``out_dir`` and returns the model (docs/solvers/rf.md)."""
    if isinstance(si, dict):
        text, base = json.dumps({"si": si}), os.fspath(base_dir)
    else:
        path = os.fspath(si)
        with builtins_open(path) as f:
            text = f.read()
        base = os.path.dirname(os.path.abspath(path))
    model, summary = _core.export_openems_channel(board._core, text, base, "" if out_dir is None else os.fspath(out_dir), channel or "")
    model = json.loads(model)
    model["summary"] = summary
    return model


def pdn(board: Board, spec: dict | str | os.PathLike) -> dict:
    """Impedance of a power distribution network (docs/solvers/pdn.md): ``spec`` is the
    project's ``pdn`` section or a project file. Returns ``frequency_hz``, ``loads``,
    ``z`` (complex lists, per load), ``pairs``, ``branches``, ``summary`` and ``warnings``."""
    if isinstance(spec, dict):
        text = json.dumps({"pdn": spec})
    else:
        with builtins_open(os.fspath(spec)) as f:
            text = f.read()
    r = json.loads(_core.pdn(board._core, text))
    r["z"] = [[complex(a, b) for a, b in zl] for zl in r.pop("z_ohm")]
    return r


def eye(board: Board, si: dict | str | os.PathLike, base_dir: str | os.PathLike = ".") -> list[dict]:
    """Eye diagrams of ``si["eyes"]``, each on a channel of ``si["channels"]``
    (docs/solvers/eye.md): the single-bit response (``pulse_t``, ``pulse_v``), the cursors,
    the pattern's eye (``eye_t``, ``eye_v``, ``density``), the statistical eye at the BER
    (``phase_t``, ``upper_v``, ``lower_v``, ``eye_height_v``, ``eye_width_ui``), the SVG
    in ``svg`` and the printable summary in ``summary``."""
    if isinstance(si, dict):
        text, base = json.dumps({"si": si}), os.fspath(base_dir)
    else:
        path = os.fspath(si)
        with builtins_open(path) as f:
            text = f.read()
        base = os.path.dirname(os.path.abspath(path))
    data, summary, svg = _core.eye(board._core, text, base)
    out = json.loads(data)
    for r, s in zip(out, svg):
        r["summary"], r["svg"] = summary, s
    return out


def ibis_switching(path: str | os.PathLike, model: str, rising: bool = True, corner: str = "typ") -> dict:
    """The switching coefficients of an IBIS model's output stage over an edge: time_s, ku
    (pullup), kd (pulldown), and the reader's warnings."""
    t, ku, kd, warnings, kind, c_comp, vpu = _core.ibis_switching(os.fspath(path), model, rising, {"typ": 0, "min": 1, "max": 2}[corner])
    return {"time_s": list(t), "ku": list(ku), "kd": list(kd), "warnings": list(warnings), "type": kind, "c_comp": c_comp, "vpu": vpu}


def deembed(spec: dict | str | os.PathLike, base_dir: str | os.PathLike = ".") -> dict:
    """De-embedding of a Touchstone file (docs/solvers/rf.md): ``input``, ``operations``
    in order (``{"renormalise_ohm": z}``, ``{"line": {...}}``, ``{"thru": {...}}``) and
    an optional ``output`` file. Returns frequency_hz, s, reference_ohm and the steps."""
    if isinstance(spec, dict):
        text, base = json.dumps(spec), os.fspath(base_dir)
    else:
        path = os.fspath(spec)
        with builtins_open(path) as f:
            text = f.read()
        base = os.path.dirname(os.path.abspath(path))
    freq, s, ref, steps, output = _core.deembed(text, base)
    return {"frequency_hz": list(freq), "s": [[complex(*c) for c in m] for m in s], "reference_ohm": ref, "steps": list(steps),
            "output": output}


def patch_antenna(epsilon_r: float, height_mm: float, frequency_hz: float = 0.0, width_mm: float = 0.0, length_mm: float = 0.0,
                  feed_ohm: float = 50.0) -> dict:
    """A rectangular microstrip patch by the transmission-line model (docs/solvers/rf.md):
    designed for frequency_hz where its width or length is not given. Returns the
    dimensions (m), effective permittivity, fringing extension, resonant frequency, slot
    conductances, edge resistance, the inset for feed_ohm, directivity and bandwidth."""
    return json.loads(_core.patch_antenna(epsilon_r, height_mm * 1e-3, frequency_hz, width_mm * 1e-3, length_mm * 1e-3, feed_ohm))


def read_touchstone(path: str | os.PathLike) -> dict:
    """A Touchstone (version 1) file of S-parameters: frequency_hz, s (per frequency, the
    matrix row major), reference_ohm."""
    path = os.fspath(path)
    ext = os.path.splitext(path)[1].lower()
    ports = int(ext[2:-1])
    with builtins_open(path) as f:
        freq, s, ref = _core.parse_touchstone(f.read(), ports)
    return {"frequency_hz": list(freq), "s": [[complex(*c) for c in m] for m in s], "reference_ohm": ref}


def cross_section(spec: dict) -> dict:
    """Per-unit-length parameters of a transmission-line cross-section (docs/solvers/si.md):
    ``conductors`` and ``dielectrics`` with shapes in millimetres, ``frequencies_hz``.
    Returns the matrices (r, l, g, c per frequency, row major) and the line quantities."""
    return json.loads(_core.cross_section(json.dumps(spec)))


def serve(board: Board, *, host: str = "127.0.0.1", port: int = 8150, open_browser: bool = True,
          dc: DcResult | None = None, thermal: ThermalResult | None = None, emi: EmiResult | None = None,
          si: SiResult | None = None) -> None:
    """Serve the 3D viewer for the board, with the fields of a DC and a thermal solution,
    the EMI estimates and the signal-integrity results, until interrupted."""
    if dc is None and thermal is not None and thermal.dc is not None:
        dc = thermal.dc
    _core.serve(board._core, host, port, open_browser, None if dc is None else dc._core,
                None if thermal is None else thermal._core, None if emi is None else emi._core,
                None if si is None else si._core)
