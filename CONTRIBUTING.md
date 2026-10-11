# Contributing

## Layout

| path | content |
| --- | --- |
| `src/core/` | the library: geometry, importers, canonical model, meshes, viewer server |
| `src/cli/` | the `kairoq` command |
| `src/python/` | Python bindings (`_core`, nanobind) and the `kairoq` package |
| `viewer/` | the browser viewer, compiled into the library |
| `tests/cpp/` | C++ unit tests (doctest) |
| `tests/python/` | Python tests: synthetic boards, importer cross-checks, validation report |
| `testdata/` | exported test boards (`testdata/README.md`) |
| `tools/` | scripts that produce test data, and benchmarks |
| `docs/` | model specification, importer and viewer documentation, validation notes |

## Development build

```
python3 -m venv .venv
.venv/bin/pip install nanobind scikit-build-core numpy shapely pytest
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release \
      -DKAIROQ_PYTHON=ON -DPython_EXECUTABLE=$PWD/.venv/bin/python
cmake --build build
```

## Tests

```
ctest --test-dir build              # C++ unit tests
.venv/bin/python -m pytest          # Python tests, using the module in build/
```

The Python tests use the package installed in the environment if there is one, and the
development build in `build/src/python` otherwise. They write the validation report to
`build/validation/`.

Solver runs can need several gigabytes; start benchmarks with a memory limit of their own
(for example `systemd-run --user --scope -p MemoryMax=10G ...` on Linux), and see
`tools/bench_dc.py` for the timing cases.

Comparisons that take longer than a test run are scripts in `tools/`, with their results
in the validation notes: `tools/ref3d.py` (an independent three-dimensional reference for
the thermal model; needs `scipy` and `pyamg`) and `tools/ipc2152.py` (trace temperature
rise against IPC-2152 data). These tools are development dependencies only.

The full-wave comparisons in `tests/python/test_openems.py` run openEMS, a few minutes
each, and only when `KAIROQ_OPENEMS_PYTHON` names a Python interpreter with the openEMS
modules (CSXCAD, openEMS) and h5py; otherwise they are skipped.

A result enters the documentation only with a test that checks it against a reference:
an analytical solution, published data, an independent implementation or another export of
the same board. Each test states its tolerance and, through the `validation` fixture, records
the compared quantity in the report.

## Conventions

- Internal lengths are integer nanometres (`kairoq::Point64`); exchange formats use
  millimetres. Name quantities with their unit where it could be in doubt (`thickness_nm`,
  `diameter_mm`).
- Code, comments, documentation and commit messages are in neutral technical English.
- A commit is one coherent change, with a subject line of at most 72 characters in the
  imperative and, if useful, up to three short lines of body. No attribution trailers
  (see "AI-assisted contributions" for how tool use is disclosed).
- Dependencies are pinned in `cmake/Dependencies.cmake` with their SHA-256 and listed in
  `THIRD_PARTY.md`.
- Test boards are produced by scripts in `tools/` and documented with their source and
  licence.

## AI-assisted contributions

Contributions made with AI coding assistants are welcome on the same terms as any other:

- Disclose the use in the pull request description: the tool, and how much of the change
  it produced (for example "tests written with ...; the solver change by hand"). Spelling
  and search need no disclosure.
- You must understand the change and be able to explain it, and how it fits the rest of
  the code, without the tool. You are responsible for it and for its licence.
- The validation rule above applies unchanged: a result needs a test against a reference
  outside the code under test.
- Issues and reports must be checked by a person before they are filed: a defect
  reproduced, a disagreement with a reference stated with the reference's source. Reports
  that are not are closed.

`AGENTS.md` gives the instructions for the assistants themselves.
