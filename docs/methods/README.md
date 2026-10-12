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
| thermal | [thermal.md](thermal.md) | [solvers/thermal.md](../solvers/thermal.md) | [validation/thermal.md](../validation/thermal.md), [thermal-transients.md](../validation/thermal-transients.md) |
| EMI loops: inductance, impedance, field | [emi.md](emi.md) | [solvers/emi.md](../solvers/emi.md) | [validation/emi.md](../validation/emi.md) |
| switch-node capacitance | [capacitance.md](capacitance.md) | [solvers/emi.md](../solvers/emi.md) | [validation/emi.md](../validation/emi.md) |
| transmission lines | [si.md](si.md) | [solvers/si.md](../solvers/si.md) | [validation/si.md](../validation/si.md) |
| crosstalk | [crosstalk.md](crosstalk.md) | [solvers/si.md](../solvers/si.md) | [validation/crosstalk.md](../validation/crosstalk.md) |
| networks: S-parameters, mixed mode, TDR | [networks.md](networks.md) | [solvers/networks.md](../solvers/networks.md) | [validation/networks.md](../validation/networks.md) |
| time-domain response | [timedomain.md](timedomain.md) | [solvers/timedomain.md](../solvers/timedomain.md) | [validation/timedomain.md](../validation/timedomain.md) |
| eye diagrams | [eye.md](eye.md) | [solvers/eye.md](../solvers/eye.md) | [validation/eye.md](../validation/eye.md) |
| PDN impedance, plane-pair cavities | [pdn.md](pdn.md) | [solvers/pdn.md](../solvers/pdn.md) | [validation/pdn.md](../validation/pdn.md) |
| RF: de-embedding, patch antennas | [rf.md](rf.md) | [solvers/rf.md](../solvers/rf.md) | [validation/rf.md](../validation/rf.md) |
| openEMS and Elmer references | [external.md](external.md) | [solvers/openems.md](../solvers/openems.md), [solvers/elmer.md](../solvers/elmer.md) | [validation/emi.md](../validation/emi.md), [validation/rf.md](../validation/rf.md), [validation/dc.md](../validation/dc.md), [validation/thermal.md](../validation/thermal.md) |

## Evidence

What each analysis has been checked against, from its validation note. *Closed forms*
are exact or published analytical results and invariants (conservation, reciprocity,
passivity), which verify that the code solves its model. *Independent solver* is another
implementation of the same or a more complete physical model, which also tests the
model's assumptions. *Measurement* is measured data. An empty cell means no such check.

| analysis | closed forms | independent solver | measurement | known gaps |
| --- | --- | --- | --- | --- |
| DC current | exact finite-element cases, constrictions, observed order of convergence, charge and energy balance, agreement of importers | Elmer in 3D | | |
| thermal | slabs, fins, series resistances, energy balance, observed order of convergence | 3D finite differences (`tools/ref3d.py`); Elmer on a region of a real board | IPC-2152 trace temperature rise | radiation not in the Elmer comparison |
| thermal transients | the discrete scheme, diffusion, steady limit, energy per step | 3D finite differences for heat sinks | | |
| EMI loops | strips, microstrip, rectangular loop, image theory, moment invariance | openEMS: inductance and far field | | |
| switch-node capacitance | parallel plates | axisymmetric field solution of the same stack | | |
| transmission lines | coaxial line, stripline, microstrip, coupled stripline; loss at defining points | finite differences for conductor loss; openEMS for delay and dispersion | | |
| crosstalk | directional coupler, far-end cancellation in one dielectric | | | only against the native channel |
| networks | uniform lines, skew mode conversion, stub resonance, reciprocity, passivity | openEMS: via stubs (7 % and 19 %) | two stub lines (published VNA data) | real channels not measured |
| time-domain response | lattice diagram, linear limits of IBIS buffers | frequency-domain solution (numpy) | | IBIS packages and supply tables not against a reference simulator |
| eye diagrams | sign enumeration, noise-limited eye, peak distortion | native time-domain response only | | no statistical simulator |
| PDN impedance | modal sum of a rectangle, eigenvalue bounds and order of convergence, static limit, mounting inductance by images | | | no full-wave or measured check |
| RF | de-embedding identities, textbook patch example | openEMS: microstrip, patch | | no measured S-parameters or patterns |

The comparisons with Elmer and openEMS run only where those solvers are installed
(`KAIROQ_OPENEMS_PYTHON`, an Elmer solver on the path); a test run without them skips them
and says so. Their recorded figures are those of the last run with the solvers.

Each document follows the same order: symbols; physical model and its assumptions;
variational or integral form; discretisation; solution of the discrete system; reported
quantities; error; implementation; references. Equations are numbered within a document.
