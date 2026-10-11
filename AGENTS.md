# Instructions for AI coding assistants

This file is addressed to AI coding assistants (agents) working on KairoQ, and to the
people who direct them. People contributing without such tools need only
`CONTRIBUTING.md`, which applies to agents as well. How the project itself uses these tools
is stated in `DISCLOSURE.md`.

## Process

- Follow `CONTRIBUTING.md`: the layout, the build, the tests and the conventions.
- Read the documentation of the part being changed before changing it: its description in
  `docs/solvers/` or `docs/importers.md`, and its validation note in `docs/validation/`.
  Read whole documents, not fragments found by searching for keywords.
- Code is licensed GPL-3.0-or-later. Do not copy code whose licence is unknown or
  incompatible; third-party code enters only as a pinned dependency
  (`cmake/Dependencies.cmake`, `THIRD_PARTY.md`).

## Responsibility

- An agent does not act outside the working copy on its own: it does not push, open or
  comment on pull requests or issues, publish packages or send reports. The person
  directing it does, after reviewing the work.
- The person submitting a change must have reviewed it, be able to explain it without the
  tool, and take responsibility for it, as for any other contribution.
- Commits carry no attribution trailers. Use of an assistant is disclosed in the pull
  request description: the tool and the extent of its use (`CONTRIBUTING.md`).
- Committed text (code, comments, documentation, commit messages) is written for any
  reader of the repository: it does not refer to the assistant, the conversation or the
  person who asked for the change.

## Results and tests

- A result enters the documentation only with a test that checks it against a reference
  outside the code under test (`CONTRIBUTING.md`).
- Never change a reference value, widen a tolerance or skip a test to make it pass. If a
  test fails, report the failure with its output; if a tolerance is wrong, say why, and
  change it as a change of its own.
- Report what was run and what was not: the commands, the results, any step skipped. State
  plainly when something could not be built, run or verified.

## Defects

To fix a defect, found by the agent or reported:

1. Read the relevant documentation and validation note first.
2. Reproduce it: a minimal board (as in `tests/python/synthetic.py`) or a test that shows
   it. If it cannot be reproduced, stop and say so; do not fix a defect that is not shown
   to exist.
3. For a result that disagrees with a reference, name the reference and where it comes
   from (formula, measurement, solver and its settings), and check the reference too.
4. Fix it with a test that fails before the fix and passes after it.
5. Build and run the C++ and Python tests (`CONTRIBUTING.md`, "Tests"). A fix must not add
   compiler warnings.
6. Write the commit message as `CONTRIBUTING.md` asks, with the cause in the body if it is
   not evident from the change.
7. State what is not covered: a defect reproduced only on synthetic boards, a fix not tested
   with the external solvers, an accuracy change not measured on a real board.

## Resources

Solver runs can need several gigabytes and long times. Run heavy jobs with a memory limit
of their own, in the background, with their output written to a file, and stop processes by
their process identifier rather than by name.
