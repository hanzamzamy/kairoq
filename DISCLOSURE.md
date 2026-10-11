# Use of automated tools

KairoQ is developed with AI coding assistants based on large language models. This
document states how they are used and what the project's results rest on. The rules for
contributions made with such tools are in `CONTRIBUTING.md` ("AI-assisted contributions")
and, for the tools themselves, in `AGENTS.md`.

## How they are used

Much of the code, the tests and the documentation was drafted with these tools, under the
direction of the author: implementation, literature search, technical discussion,
exploration of alternative designs, documentation, and repository operations such as
organising commits. The commits do not mark individual changes as tool-assisted; this
document covers the whole history.

## Who is responsible

The author decides the architecture, the physical models and numerical methods, the
validation method and its tolerances, and accepts or rejects every result. The author
reviews the changes, is responsible for them and for their licensing, and answers for the
project. The tools are not authors of the project and carry no responsibility for it.

## What the results rest on

Trust in KairoQ does not rest on how its code was written. Every result in the
documentation is checked by a test against a reference outside the code under test: an
analytical solution, published or measured data, an independent implementation, or an
external solver (Elmer, openEMS). The tests state their tolerances and record every
compared quantity in a report (`docs/validation/`). A result without such a check is
stated as unvalidated, and the known limits of each analysis are listed with it
(`docs/solvers/`).
