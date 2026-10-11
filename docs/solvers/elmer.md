# 3D references with Elmer

The DC solver (`docs/solvers/dc.md`) treats copper as sheets and plated holes as
resistances between the layers' mid-planes. Elmer (<https://www.elmerfem.org>), an
open-source finite element package, solves the same problem in three dimensions as a
reference: the current through the copper's thickness, into and along the holes. KairoQ
writes the model and Elmer's input, runs Elmer when asked, and solves the ports' conditions
from Elmer's results. Elmer is not part of KairoQ; it is found on the path or named in
the project. Validation: `docs/validation/dc.md`, section "3D reference".

## Model

- **One triangulation** of the footprint of the ports' nets (the union of their copper on
  all layers and of their plated holes) carries every layer's copper outline, the ports and
  the holes as constraints. Outlines that nearly coincide (a track's edge where it joins a
  pad on one layer and the pad's arc on the next) are made identical first: vertices
  within a micrometre become one, and every edge takes the vertices that lie on it.
- **Copper.** The triangulation is stacked at the faces of the copper layers, with
  `copper_divisions` element layers through each; a triangle inside a layer's copper gives
  a 6-node wedge through that layer, of the layer's conductivity at the DC temperature.
- **Plated holes.** A triangle inside a plated hole gives a column of wedges over the
  hole's span, with the dielectric gaps divided into steps of at most four minimum sizes.
  The column conducts along the hole as the plating spread over it,
  $\sigma_z = \sigma_p\, \pi t_p (d + t_p) / (\pi d^2/4)$, and across it as copper, as the
  thin wall does around the hole. Plated slots are left out and reported.
- **Ports** are the outer faces of the copper inside each port's area: the top face, or
  the bottom face on the bottom layer. Every layer shares the triangulation's nodes, so the
  volumes join without coupling conditions.

## Solution

For $n$ ports, Elmer solves the static current problem $n$ times (`port_<k>.sif`), each with
port $k$ at 1 V and the others at 0 V, by conjugate gradients with an incomplete LU
preconditioner, aborting if it does not converge. The sum of the nodal loads on each port
is the current into it, exact for the discrete solution; together they form the
conductance matrix $G$ of the equipotential ports. The project's conditions (set
voltages, set currents) give the ports' voltages from $I = G V$; a group of connected ports
without a set voltage takes its first port as 0 V, as in the DC solver.

```yaml
dc:
  ports: [...]
  elmer: {max_edge_mm: 1.0, min_size_mm: 0.2, grading: 0.5, copper_divisions: 1, solver: ElmerSolver}
```

| key of `dc.elmer` or `thermal.elmer` | default | meaning |
| --- | --- | --- |
| `max_edge_mm`, `min_size_mm`, `grading` | 0.5, 0.05, 0.3 | the triangulation: largest edge, size at outlines, growth away from them |
| `copper_divisions` | 1 | element layers through each copper layer |
| `solver` | `ElmerSolver` | the Elmer solver executable |

```
kairoq elmer project.yaml [<input>...] -o dir [--run]
```

writes `dir/mesh/` and the cases; with `--run` it runs them, solves the DC problem with
KairoQ too and prints both sets of port voltages.

## Thermal reference

The thermal solver (`docs/solvers/thermal.md`) treats every copper layer and every slice of
dielectric as a conducting sheet. The thermal export resolves the whole board through its
thickness instead:

- **One triangulation** of the board, without its non-plated holes, carries every layer's
  copper outline, the plated holes and the sources' areas as constraints, as for DC.
- **Layers.** Every copper layer gives `copper_divisions` element layers of copper or, between
  the copper, of the material that fills the gaps: the dielectric below on inner layers, air
  (no elements) or mask on outer layers, as `thermal.outer_gaps` says. Every dielectric gives
  element layers of at most two minimum sizes, with its in-plane and through-thickness
  conductivities.
- **Plated holes** are columns over the hole's span that conduct along the hole as the
  plating spread over it and across it as copper; plated slots are left out and reported.
- **Sources** heat the copper inside their areas uniformly by volume, at their set power.
- **Faces.** The top and bottom of every column lose heat by convection through the solder
  mask, which is a series resistance, $h_\mathrm{eff} = 1/(1/h + t_\mathrm{mask}/k_\mathrm{mask})$,
  as in the layered model. The sides are adiabatic.

Radiation is left out (set `emissivity: 0` in both models to compare them). Heat sinks, fixed
temperatures, junction resistances, Joule heat and transients are not exported. Elmer
solves the steady heat equation by IDR(s) with an ILU1 preconditioner, which on these stiff
systems (thin copper, conductivities a thousand times apart) is about fifteen times faster
than conjugate gradients, and writes the temperature in the gmsh format. KairoQ reads it
and reports, for each source, the mean and the largest temperature of its copper, and the
largest temperature of each copper layer, next to its own solution.

```yaml
thermal:
  top: {convection_w_m2k: 10, emissivity: 0}
  bottom: {convection_w_m2k: 10, emissivity: 0}
  sources: [...]
  elmer: {min_size_mm: 0.1, max_edge_mm: 1.0, grading: 0.3, solver: ElmerSolver}
```

```
kairoq elmer project.yaml [<input>...] -o dir --thermal [--run] [--crop x0,y0,x1,y1]
```

A whole board gives tens of millions of elements; `--crop` keeps the part of the board
inside a rectangle (mm), whose cut edges are adiabatic in both models, so that a region of
a real board can be compared. Validation: `docs/validation/thermal.md`, section "Real board
against Elmer".

## Limits

- The model has one triangulation for all layers: a board whose nets cover it on many
  layers gives millions of elements (cm5_minima's supply and ground: 4.5 million nodes at
  1 mm and 0.2 mm). The sizes are set by hand; there is no adaptive refinement.
- Ports are equipotential faces on the copper's outer side; contacts on inner layers take
  the top face.
- DC: the temperature is uniform (no electrothermal coupling); plated slots are left out.
- Thermal: steady state, linear (no radiation), set powers only.
