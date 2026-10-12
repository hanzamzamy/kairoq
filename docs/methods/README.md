# Methods

These documents derive each analysis from its physical model to the numbers it reports:
the governing equations and their assumptions, the weak or integral form, the
discretisation, the linear algebra and the estimates of error. The code implements them;
each document ends with a table of where each step is implemented, so that the code can be
read against the mathematics. Where the two disagree, the code is wrong or the document is
incomplete, and either is a defect.

The other documents answer other questions: `docs/solvers/` how to use an analysis and its
limits, `docs/validation/` how well it agrees with references. Methods taken from the
literature are credited where they are used and listed in `docs/references.bib`.

| analysis | method | user guide | validation |
| --- | --- | --- | --- |
| DC current | [dc.md](dc.md) | [solvers/dc.md](../solvers/dc.md) | [validation/dc.md](../validation/dc.md) |
| EMI loops: inductance, impedance, field | [emi.md](emi.md) | [solvers/emi.md](../solvers/emi.md) | [validation/emi.md](../validation/emi.md) |
| transmission lines | [si.md](si.md) | [solvers/si.md](../solvers/si.md) | [validation/si.md](../validation/si.md) |
| switch-node capacitance | [capacitance.md](capacitance.md) | [solvers/emi.md](../solvers/emi.md) | [validation/emi.md](../validation/emi.md) |
| thermal | [thermal.md](thermal.md) | [solvers/thermal.md](../solvers/thermal.md) | [validation/thermal.md](../validation/thermal.md), [thermal-transients.md](../validation/thermal-transients.md) |

Each document follows the same order: symbols; physical model and its assumptions;
variational or integral form; discretisation; solution of the discrete system; reported
quantities; error; implementation; references. Equations are numbered within a document.
