# Validation: heat sinks, component power and transients

## Scope

Heat sinks, component power and transients in the thermal solver (`docs/solvers/thermal.md`)
are checked for:

1. heat sinks: exactness against series resistances, the contact area of a partial sink,
   the energy balance, and the spreading through the board into a sink against the
   three-dimensional reference;
2. transients: exactness of the discrete scheme, convergence in time, diffusion against a
   closed form, the steady limit and the energy balance of every step, with and without
   electrothermal coupling;
3. component power tables: reading, grouping, units and the checks;
4. the example board, the viewer and the time taken.

The tests are in `tests/python/test_heatsink.py` and `tests/python/test_transient.py`; they
write their compared quantities to `build/validation/report.md`. The three-dimensional
comparison is `tools/ref3d.py` (logs in `logs/`, not tracked).

Material values unless stated otherwise: as in `docs/validation/thermal.md`; volumetric heat
capacities $3.45 \times 10^6$ (copper), $2.2 \times 10^6$ (dielectric) and
$1.8 \times 10^6$ J/(m³ K) (mask).

## Heat sinks

| case | quantity | reference | relative error |
| --- | --- | --- | --- |
| slab 20 × 20 mm, 2 W on the top, the bottom on a plate at 40 °C through 0.5 mm of 3 W/(m K) | top and bottom temperature | series resistances | $1.3 \times 10^{-12}$ |
| same, the plate a body of 5 K/W to ambient through 50 K mm²/W | sink and top temperature | series resistances | $1.6 \times 10^{-12}$ |
| 30 × 30 mm board, 3 W pad, sink of 4 K/W under 20 × 15 mm, radiation | (heat out − heat in) / heat in | 0 | $2.6 \times 10^{-11}$ |
| isothermal board on the same 20 × 15 mm region, sink at 25 °C | rise over the sink | $P\,(r_\mathrm{mask} + r_\mathrm{TIM})/A$ | $2.2 \times 10^{-5}$ |

The isothermal case gives the dielectric a conductivity of $10^6$ W/(m K) and switches the
faces off, so that the whole power crosses the interface over exactly the sink's area; the
remaining error, $2 \times 10^{-5}$, is the temperature variation left in the board (it was
$2 \times 10^{-3}$ at $10^4$ W/(m K)). It confirms that the region is a constraint of the
mesh and is counted by area.

**Against the three-dimensional reference.** The cases of `docs/validation/thermal.md`, with the
whole bottom face on a heat sink at ambient through 0.5 mm of a 3 W/(m K) interface instead
of convection, 0.5 W in the pad, references extrapolated from two resolutions, with air
and with mask in the gaps of the outer copper:

| case | 3D, air | 3D, mask | sheets | error, sheets | error, single sheet |
| --- | --- | --- | --- | --- | --- |
| 2 mm square, 1.5 mm, no plane | 175.6 ± 0.7 % | 169.9 ± 0.5 % | 171.9 | −2.1 % / +1.2 % | −42.4 % / −40.4 % |
| 2 mm square, 0.2 mm, plane | 68.5 ± 1.0 % | 65.6 ± 0.6 % | 67.6 | −1.2 % / +3.1 % | −12.6 % / −8.8 % |
| 0.5 mm square, 0.2 mm, plane | 592.9 ± 3.3 % | 523.7 ± 2.0 % | 555.0 | −6.4 % / +6.0 % | −31.9 % / −22.9 % |

As without a sink, the layered model lies between the two references, within about 3 % for
the 2 mm pads; one sheet per copper layer is far below. A sink does not change where the
layered approximation holds.

These comparisons were run at a target error of 0.03. At the default target, 0.1, the
0.2 mm case on a sink comes out 2 % higher (69.1 K): its heat flows mostly through the
thickness into the sink, where the error estimate, which covers the in-plane conduction
only, sees little, and a tighter target then refines from 0.24 to 4.1 million triangles at
once. This case is therefore not part of the test suite; the three cases without a sink
are.

## Transients

| case | quantity | reference | result |
| --- | --- | --- | --- |
| slab with a junction (4 K/W, 0.5 J/K) and a sink body (3 K/W, 20 J/K), power 2 W then 0.5 W at 40 s, 100 s | largest deviation of junction, sheet and sink temperatures over all steps | backward Euler on the four-node network | $3.0 \times 10^{-11}$ K |
| slab, 2 W step, 60 s, steps of 2, 1 and 0.5 s | ratio of the errors for halved steps | 0.5 (first order) | 0.502; the next halving within 5 % of 0.5 |
| copper strip, one end held at 75 °C, no losses, 1 s | rise 2, 5, 10, 20 mm from the end | $\theta_0\, \mathrm{erfc}\big(x / 2\sqrt{\alpha t}\big)$ | within 0.14 % |
| fin held at 75 °C at one end, 4000 s | heat from the fixed end against the steady solution | steady solution | $1.8 \times 10^{-8}$ |
| 30 × 30 mm board, 3 W pulse of 30 s, partial sink, junction, radiation | largest (in − out − stored) / peak power over all steps | 0 | $9.1 \times 10^{-8}$ |
| strip 5 × 100 mm carrying 5 A from ambient, coupled, 5000 s | largest relative (in − out − stored) | 0 | $6.9 \times 10^{-11}$ |
| same | maximum temperature at the end against the coupled steady solution | steady solution | $2.9 \times 10^{-4}$ |

- **The discrete scheme.** With one dielectric slice and a uniform field the model reduces
  to four lumped nodes: junction, top sheet, bottom sheet and sink body, with the capacities
  of copper, half the dielectric and mask on each sheet. The solver follows backward Euler
  on that network to rounding, over a schedule with a step of power and the step lengths
  the solver chose (reset to the first step at the schedule's breakpoint). This checks the
  capacities, the junction and sink nodes, the interface and the schedule at once.
- **Convergence in time.** Against the exact (exponential) solution of the same network the
  error halves with the step, as it must for backward Euler.
- **Diffusion.** In a copper strip with negligible dielectric the temperature spreads with
  $\alpha = k / (\rho c) = 1.12 \times 10^{-4}$ m²/s; at 1 s the computed profile agrees with the
  semi-infinite solution within 0.14 % (0.25 mm elements, steps up to 5 ms).
- **Steady limit and balance.** Long after a step the transient reaches the steady
  solution on the same mesh. Every step closes its heat balance: in = out + rate of change of
  the stored heat, with the faces as linearised in the step; with radiation the
  linearisation uses the slopes of the peak state, so the balance is exact for the
  discrete problem solved.
- **Electrothermal coupling.** With the DC solution renewed after every 0.2 K of change,
  the end state is within $3 \times 10^{-4}$ of the coupled steady solution, a difference set by
  that renewal threshold.

## Component power tables

`test_component_power_table` and `test_component_power_table_checks` read tables with
comma and semicolon separators, quoted groups of references ("J1, J2", each dissipating the
stated power), units (`250 mW`), an optional junction column, empty power cells (skipped),
and check the errors: references not on the board (or a warning when they are to be
ignored), a reference listed twice, an unknown unit; a source given in the project file
replaces the table's entry. `test_cli_project_with_table_sink_and_transient` runs a
project with a table, a heat sink and a transient from the command line.

## Example board

`examples/cm5_minima_transient.yaml`: the `cm5_minima` carrier switched on at $t = 0$:
0.6 W in the converter U702 (15 K/W to its junction), 0.25 W in its inductor and 80 mW in
each of two load switches, from `examples/cm5_minima_power.csv`; 3 A through the 5 V path
(coupled); a plate of 6 K/W and 90 J/K under 60 × 30 mm of the bottom face through a 1 mm
pad of 3 W/(m K); 600 s.

After 600 s the converter's junction is at 68.1 °C and still rising; the plate is at
28.4 °C, takes 0.83 W of the 1.07 W through the pad and passes 0.57 W on to ambient, the rest
heating it (its time constant is about 540 s).

The viewer shows the fields kept during the run on a time control, with the highest
temperature and the sources over time; checked headless on a smaller board (temperatures
change with the time selected, the probe reads them, no console errors).

## Performance

A 20-thread x86-64 computer, Release build.

| case | triangles | unknowns | steps | time | peak memory |
| --- | --- | --- | --- | --- | --- |
| `cm5_minima`, 600 s, coupled, sink | 1.37 M | 0.64 M | 52 | 4 min 49 s | 3.5 GB |

Of that time, 55 s are the steady solutions (the mesh at the peak power), 16 s the seven
multigrid set-ups (one per step length), about 125 s the 1961 conjugate-gradient
iterations of the steps (38 per step) and 32 s the nine DC solutions. Each step was then
solved to the same relative tolerance as a steady solution ($10^{-10}$). A step now solves
for the change of the temperature, to a tolerance relative to the residual of the previous
temperature (`step_tolerance`, $10^{-4}$): on the diffusion into a strip, 6.0 instead of
12.7 iterations per step, the temperatures unchanged to four digits. The energy
balance of each step then closes to the step tolerance rather than to $10^{-6}$: on the
board with a sink, a junction and radiation, $3 \times 10^{-4}$ of the peak power; the
balance tests solve their steps to the solver tolerance (`step_tolerance: 0`).

The test suite (`tests/python`) runs in 98 s.

## Findings

- Where the heat flows mostly through the thickness (a thin dielectric on a heat sink),
  the 0.2 mm case above is 0.9 % high at the default target error and converges with a
  tighter one (`docs/validation/thermal.md`, "Error estimate where the heat leaves through
  the thickness").
- The steady part of a transient (the mesh) can cost as much as the time steps; for
  repeated transients on one board a mesh kept from an earlier run would save it.
