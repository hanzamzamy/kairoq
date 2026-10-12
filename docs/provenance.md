# Provenance of results

A result is traceable to the build that computed it and to the exact files it was computed
from. Two runs whose provenance agrees used the same program on the same inputs; where it
differs, it shows what changed.

## What is recorded

| field | content |
| --- | --- |
| `kairoq` | `version`; `commit`, the git commit of the source it was built from, followed by `+modified` when the source tree had uncommitted changes, or `unknown` when built outside git; `compiler`; `platform` |
| `time_utc` | when the result was computed |
| `command` | the command line |
| `inputs` | the files given, the project file, and every file the project's values name (boards, stackups, Touchstone and IBIS files), each with its SHA-256 and size; a project's files are named as the project names them |
| `board_sources` | the files the board was imported from, by file name, with their SHA-256 |
| `summary` | the printed summary of the result |

The digest of a file is the SHA-256 of its bytes. The digest of a directory is the SHA-256
of the lines `<relative path> <sha256>` and a line feed, one per file in the directory and its
subdirectories, sorted by relative path, with `/` as separator: it changes when any file is
added, removed, renamed or changed.

The settings of an analysis are those of the project file, whose digest is recorded, and
the defaults of the recorded version and commit. Values under an `output` key of a project
name files the run writes and are not inputs.

## Where it is written

- **Command line.** Every analysis prints the version and commit as the last line of its
  summary. With `-o dir`, `dir/provenance.json` holds the record; with `--record file.json`,
  the record is written there (for analyses without `-o` too). Touchstone files carry it as
  comment lines (`!`) before the data.
- **Canonical model.** A model records the digests of the files it was imported from
  (`source.digests`, `docs/canonical-model.md`); a board opened from a model keeps them, so
  a result computed from the model traces to the original fabrication files.
- **Python.** `kairoq.provenance(board, project, files)` returns the record for a result
  computed now, to be stored with it; `kairoq.build_info()` and `kairoq.__version__` give
  the build.

```python
board = kairoq.open("fab.zip")
r = kairoq.solve_dc(board, "project.yaml")
record = kairoq.provenance(board, "project.yaml")
```
