# Validation: DC current

## Scope

The DC solver (`docs/solvers/dc.md`) is checked for:

1. exactness where the finite element solution is exact;
2. convergence to closed-form solutions where it is not;
3. the conservation of charge and energy on a real board;
4. agreement between the two importers of the same board;
5. the meaning of its error estimate;
6. time and memory.

The tests are in `tests/python/test_dc.py` and `tests/cpp/test_mesh2d.cpp`; the Python
tests write their compared quantities to `build/validation/report.md`. The timing and
convergence figures come from `tools/bench_dc.py` (logs in `logs/`, not tracked).

Conductivity in the analytical cases: $\sigma = 5.8 \times 10^7\ \mathrm{S/m}$, copper
$t = 35\ \mu\mathrm{m}$, plating $18\ \mu\mathrm{m}$.

## Exact cases

| case | quantity | reference | relative error |
| --- | --- | --- | --- |
| straight track 48 × 1 mm between equipotential pads | resistance | $L/(\sigma t w)$ | $7 \times 10^{-14}$ |
| same, with set currents only (floating copper) | resistance | $L/(\sigma t w)$ | $2 \times 10^{-14}$ |
| same at 75 °C | $R(75)/R(20)$ | $1 + 0.00393 \cdot 55$ | $2 \times 10^{-13}$ |
| one via, 0.3 mm, 1.545 mm between layer mid-planes | barrel resistance | $L/(\sigma \pi t_p (d + t_p))$ | 0 |
| chain of 8 vias, plating 18 and 25 µm, same mesh | $R(18) - R(25)$ | $8\,(R_b(18) - R_b(25))$ | $1 \times 10^{-11}$ |

The track is exact because the field between equipotential ends of a rectangle is uniform,
which linear elements represent exactly; the plating difference is exact because the copper
part of the chain does not depend on the plating.

## Convergence to closed-form solutions

**Radial spreading through a via.** A copper disc on each layer, a plated via
($d = 0.3$ mm) at the centre and an annular port at the rim (inner radius $r_2 = 4.5$ mm):

$$R = \frac{2 \ln(r_2 / r_1)}{2 \pi \sigma t} + R_\mathrm{barrel}.$$

**Constriction.** A strip of width $W = 4$ mm with a wall across it, open over $a = 0.5$ mm.
For a wall of zero thickness, conformal mapping gives the additional resistance

$$\Delta R = \frac{2}{\pi \sigma t} \ln \csc \frac{\pi a}{2W};$$

the wall of 5 µm adds about $w (1/a - 1/W)/(\sigma t)$. With finer meshes the computed
constriction approaches the formula as the wall gets thinner (+1.3 % of $\Delta R$ for a
20 µm wall, +0.3 % for 5 µm, at 1.6 million triangles), which supports both the formula and
the solver.

| case | target $e$ | estimated $e$ | relative error of $R$ | bound $e^2$ |
| --- | --- | --- | --- | --- |
| via spreading | 0.1 | 0.069 | $-6.7 \times 10^{-4}$ | $4.8 \times 10^{-3}$ |
| via spreading | 0.03 | 0.021 | $-9.8 \times 10^{-5}$ | $4.5 \times 10^{-4}$ |
| constriction | 0.1 | 0.087 | $-6.9 \times 10^{-3}$ | $7.5 \times 10^{-3}$ |
| constriction | 0.03 | 0.028 | $-1.4 \times 10^{-3}$ | $7.6 \times 10^{-4}$ |

The error decreases with the target. The squared estimate bounds the error of the
two-port resistance in three of the four cases; for the constriction at the tighter target
the error is about twice $e^2$. The re-entrant corners of the wall are field singularities,
where gradient recovery is known to underestimate; the test therefore allows
$\max(e^2, 2 \times 10^{-3})$ for this case. The resistance converges from below, as it must for a
conforming finite element solution with set potentials.

## Real board

KiCad demonstration board `cm5_minima` (6 layers, 444 vias), 3 A from the input connector
J101 to the module connectors, supply (+5V) and return (GND) at once
(`examples/cm5_minima_5v.yaml`).

| target | triangles | +5V drop | GND rise at the load | time |
| --- | --- | --- | --- | --- |
| 0.1 (default) | 0.45 M | 19.017 mV | 0.369 mV | 3.3 s |
| 0.05 | 4.0 M | 19.058 mV | 0.378 mV | 40 s |
| 0.03 | 4.7 M | 19.063 mV | 0.379 mV | 46 s |
| 0.02 (triangle budget) | 5.0 M | 19.065 mV | 0.380 mV | 129 s |

At the default target the supply drop is within 0.3 % of the finest solution; the ground
rise, whose path carries 1.9 % of the dissipation, within 3 %. This is the limitation of a
global estimate stated in `docs/solvers/dc.md`.

Conservation at the default target: port power and Joule heat in copper and plating agree
to $1.3 \times 10^{-5}$ of the dissipation, the port currents sum to $1.3 \times 10^{-8}$ of
their magnitudes; both are limited by the linear solver tolerance.

The same board imported from IPC-2581 gives a +5V drop within $2.2 \times 10^{-5}$ of the
Gerber import, with the copper conductivity set equal (the IPC-2581 stackup states
$5.959 \times 10^7\ \mathrm{S/m}$, the Gerber job file none).

## Mesh

Mesh tests: a plate with a hole and a port on its edge is covered exactly (area to
$10^{-12}$), the port is an exact union of triangles, no triangle is degenerate and the
smallest angle is at least 19.9°; refinement seeds that coincide or lie on the outline are
dropped.

Defects found and corrected during this validation, each now covered by a test or a
measurement above:

- outline simplification moved a straight copper edge off a coinciding pad edge by
  0.3 µm, which produced slivers and slowed the linear solver; constraint edges within the
  simplification tolerance of the outline are now merged into it;
- collinear outline points produced zero-area triangles outside the copper whose centroid
  lay on the outline; degenerate triangles are now removed;
- a hole wall inside a port pad was merged into the pad and left the barrel unconnected;
  constraints are now intersected with the copper one by one;
- copper inside equipotential areas set a floor to the error estimate; it is now excluded;
- holes joined the first copper at their centre whatever its net, which coupled nets on a
  board with shorts; a hole now joins copper of its own net only.

## Performance

A 20-thread x86-64 computer, Release build, default settings
(`tools/bench_dc.py`):

| board | size | triangles | unknowns | time | peak memory |
| --- | --- | --- | --- | --- | --- |
| `cm5_minima`, 2 nets | 61 × 64 mm, 6 layers | 0.45 M | 0.22 M | 1.5 to 1.9 s | 0.8 GB |
| synthetic power board, 2 nets | 100 × 100 mm, 4 layers, 2304 vias, 0.1 mm tracks and clearances | 2.3 M | 1.1 M | 8.8 to 8.9 s | 2.4 GB |

The synthetic board takes about 9 s: mesh 2.8 s, assembly 1.1 s, solver 4.1 to 4.3 s. The
set-current ports are free unknowns of one system, and refinement aims at 0.9 of the
target, which avoids a third pass. The mesher clips each port and hole wall against the
outline rings near it only, and its point-in-region test ends at the first grid cell that
no boundary crosses; clipping against the whole plane, whose outline has a ring per
clearance, took 3 s per pass on a plane with about 1150 holes, and the mesh 8.0 to 8.8 s
in all. The largest plane takes about 1.7 s per pass, most of it in the
triangulation's vertex insertion and quality refinement; meshing a plane in parallel tiles
would shorten that, at the cost of seams that must stay conforming under refinement, and
is not done.

The multigrid preconditioner uses classical (Ruge–Stüben) coarsening. With smoothed
aggregation instead, the same solutions need 31 instead of 24 conjugate-gradient iterations
on the synthetic power board and 30 in both cases on `cm5_minima`, with identical port
voltages.

## The error estimate on more geometries

The estimate is meant to bound the relative error of the dissipation by about its square.
Against fine solutions of the same problems (target 0.01; voltages set at both ports), the
ratio $\lvert \Delta P \rvert / P$ to the squared estimate over the adaptive passes:

| case | passes (triangles) | ratio |
| --- | --- | --- |
| strip with a 0.5 mm aperture in a 5 µm wall | 587 to 1.3 million | 1.58, 1.24, 0.84, 0.73, 0.57, 0.30 |
| chain of 8 vias | 1655 to 307 000 | 0.36, 0.30, 0.26, 0.20 |
| cm5_minima +5V, connector to module | 48 000 to 3.7 million | 0.57, 0.47, 0.17 |
| cm5_minima GND, connector to module | 276 000 to 4.7 million | 0.68, 0.61, 0.39 (then within the reference's own error) |

The bound holds within a factor of 1.6 where the singular aperture is still coarsely
resolved and is conservative elsewhere. On the aperture, the energy-norm error itself
(from $\sqrt{(P_h - P)/P}$) is estimated within factors of 0.8 to 1.3 down to 4 %.
The dissipation compared is the copper's; the barrels' resistance does not depend on the
mesh.

## 3D reference (Elmer)

Elmer solves the same DC problems in three dimensions (`docs/solvers/elmer.md`;
`tests/python/test_elmer.py`, run when an Elmer solver is found). Port voltages for the
currents given:

| case | quantity | reference | Elmer (3D) | KairoQ (2.5D) |
| --- | --- | --- | --- | --- |
| strip 48 × 1 mm, 2 A | voltage | $IL/(\sigma t w)$: 47.290 mV | 47.306 mV (+0.03 %) | 47.291 mV |
| chain of 8 vias with tracks on alternate layers, 1 A | voltage | | 41.45, 41.56, 41.66 mV (three refinements) | 40.73 mV; 40.88 mV refined |
| `cm5_minima`, 3 A from the 5 V connector to the module and back | drop of +5V | | 19.042 mV | 19.017 mV (−0.13 %) |
| the same | rise of GND at the module | | 0.380 mV | 0.369 mV (−2.9 %) |

- **Strip.** Elmer's coarse default mesh is within 0.03 % of the closed form.
- **Via chain.** The 3D solution rises by about 0.1 mV at each refinement (the
  re-entrant corners where the current turns from the copper into a hole converge slowly),
  and lies 1.9 % or more above the 2.5D solution: about 0.1 mΩ per via, the constriction
  of the current where it leaves the 35 µm copper for the hole, which sheets and a
  barrel between mid-planes do not see. A through hole carrying a large share of a path's
  resistance is underestimated by the 2.5D model by this amount.
- **Real board.** On `cm5_minima` (Elmer: triangles of 0.2 to 1 mm, 4.0 million nodes, four
  solutions of conjugate gradients) the supply's drop agrees within 0.13 %; the ground's
  small rise, which crosses the board through many vias, is 2.9 % higher in 3D, in the
  direction of the via chain. The 12 plated slots, left out of the 3D model, would lower
  its resistance.
- **Defects found.** The first exports gave wrong currents: Elmer's diffusive flux at a
  port is computed from element gradients and misses current where it turns sharply into
  the copper (12.6 instead of 42.3 S); the sum of the nodal loads is used instead. A fine
  mesh failed to converge: nearly coincident outlines of different layers gave slivers
  (angles below 0.03°), and BiCGStab broke down without an error. The outlines are now
  made identical where they nearly coincide, the solver is conjugate gradients, and a
  solution that does not converge stops the run.

## Limits of this validation

- The analytical cases are planar; the through-thickness approximation is checked against
  the 3D reference above on a strip and a via chain.
- On a real board the external reference is Elmer's 3D model above, at one mesh.
- The error estimate is checked against closed forms on two geometries, and against fine
  solutions on four (below); its reliability on every geometry is not proven.
