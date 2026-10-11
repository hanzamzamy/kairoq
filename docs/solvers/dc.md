# DC current solver

The DC solver computes the steady current distribution in the copper of a board: the
voltage on every copper layer, the current density, the Joule heating, the current in each
plated hole and the voltage or current at each port. Validation: `docs/validation/dc.md`.

## Model

### Copper layers

Each copper layer is a conducting sheet of thickness $t$ and conductivity $\sigma$. The
potential $V$ satisfies

$$\nabla \cdot (\sigma t\, \nabla V) = 0$$

in the copper of the layer, with no current across the copper edges. The sheet current
density is $\mathbf{K} = -\sigma t\, \nabla V$ (A/m, reported in A/mm of width), the
volume current density $\mathbf{J} = \mathbf{K}/t$, and the Joule heating per unit area
$\sigma t\, |\nabla V|^2$. The current is uniform through the thickness of a layer; this
holds while the copper is thin compared with the distances over which the current changes
direction.

The conductivity at temperature $T$ is

$$\sigma(T) = \frac{\sigma_{20}}{1 + \alpha\,(T - 20\ \mathrm{^\circ C})}, \qquad \alpha = 0.00393\ \mathrm{K^{-1}},$$

with $\sigma_{20}$ from the stackup when it states one (IPC-2581), and otherwise
$5.8 \times 10^7\ \mathrm{S/m}$ (annealed copper, IEC 60028). The temperature is uniform
(`temperature_c`), or, when the DC problem is solved with the thermal one, the computed
temperature field (`docs/solvers/thermal.md`).

### Plated holes

A plated hole joins the copper islands of its net that it passes through. Between two
consecutive connected layers $a$ and $b$, the barrel is a conductance

$$G = \frac{\sigma_p\, A}{L}, \qquad A = \pi\, t_p\,(d + t_p) \;(+\, 2\, t_p\, \ell \text{ for a slot of length } \ell),$$

where $d$ is the finished hole diameter, $t_p$ the plating thickness (default
$18\ \mu\mathrm{m}$), $L$ the distance between the mid-planes of the two layers and
$\sigma_p$ the plating conductivity (default that of copper). The hole wall on each layer is
an equipotential ring; layers the hole passes without copper of its net (antipads) are not
connected.

### Ports

A port is an equipotential area of copper with either a set voltage or a set current
(positive into the copper). Its area is selected by pad (`REF.PIN`: every pad of the pin,
on every layer), by component and net (every pad of the component on the net), or as a
contact disc at a point on a layer. All copper of a port is at one potential, as for a
soldered pin.

Copper reached only by set-current ports takes the potential of its first port as 0 V; the
currents must then sum to zero. Copper not connected to any port carries no current and is
not solved.

## Discretisation

Only the copper islands connected to a port are meshed. Each island is triangulated with a
constrained Delaunay triangulation with quality refinement (smallest angle about
$20^\circ$); hole walls and port outlines are mesh edges, so that every equipotential area
is an exact union of triangles. Copper outlines are simplified within $5\ \mu\mathrm{m}$
and hole walls are polygons of at least 16 sides of equal area.

The potential is approximated by linear finite elements. Hole walls and ports tie their
nodes to one unknown; set potentials are boundary values, and a set-current port is a free
unknown loaded with its current, so that one solution gives all ports. The linear system
is solved by conjugate gradients preconditioned with classical (Ruge-Stüben) algebraic
multigrid (AMGCL), to a relative residual of $10^{-10}$.

### Adaptive refinement and the error estimate

The solution is refined where it is inaccurate. After each solution the error is estimated
per triangle by gradient recovery (Zienkiewicz-Zhu): the recovered gradient at a node is
the area-weighted mean of the element gradients around it, and the estimate on a triangle
$T$ is

$$\eta_T^2 = \sigma t \int_T |\mathbf{G}^* - \nabla V_h|^2 \, dA .$$

Triangles inside equipotential areas are excluded. The relative error in the energy norm
is estimated as

$$e = \sqrt{\frac{\sum_T \eta_T^2}{\sum_T \sigma t \int_T |\nabla V_h|^2 \, dA}} .$$

While $e$ exceeds `target_error` (default 0.1), triangles are marked in descending order
of $\eta_T$ until the predicted total meets the target, and are refined by new vertices on
a lattice inside them. The mesh is rebuilt from the previous vertices and the new ones.

The relative error of the total dissipated power, and of the resistance between two ports,
is about $e^2$: 1 % at the default target. Gradient recovery underestimates the error at
sharp re-entrant corners, where it can be about twice $e^2$ (`docs/validation/dc.md`). The estimate is global: a quantity of a path that
carries a small share of the power (a ground return next to a supply plane, for example)
can be less accurate. A smaller target, or a solution with only that path's ports, gives
it more accuracy. Measured accuracy: `docs/validation/dc.md`.

## Project file

```yaml
board: ../fab/                 # inputs, relative to this file (a path or a list)
dc:
  temperature_c: 20
  plating_thickness_um: 18
  ports:
    - {name: VIN,  pad: J1.1, voltage_v: 5.0}
    - {name: LOAD, component: U3, net: +3V3, current_a: -2.0}
    - {name: TP,   at_mm: [41.2, 17.0], layer: 1, diameter_mm: 1.0, current_a: 0.5}
```

| key under `dc` | default | meaning |
| --- | --- | --- |
| `temperature_c` | 20 | copper temperature |
| `copper_conductivity_s_per_m` | from the stackup, else $5.8 \times 10^7$ | overrides the conductivity of every layer (at 20 °C) |
| `plating_thickness_um` | 18 | hole wall plating |
| `plating_conductivity_s_per_m` | $5.8 \times 10^7$ | at 20 °C |
| `target_error` | 0.1 | relative error in the energy norm at which refinement stops |
| `max_passes` | 6 | adaptive passes |
| `max_triangles` | 5 000 000 | refinement stops before the mesh exceeds this size |
| `mesh` | | initial mesh: `max_edge_mm` (1), `min_angle_deg` (20), `boundary_tolerance_um` (5), `min_edge_um` (0.5); `grading` and `min_size_um` add a size field graded from the copper edges (off by default) |
| `solver` | | `tolerance` ($10^{-10}$), `max_iterations` (2000) |
| `ports` | | list of ports, as above |

A port gives exactly one of `pad`, `pads` (a list), `component` with `net`, or `at_mm` with
`layer` and `diameter_mm`; and exactly one of `voltage_v` or `current_a`. Unknown keys are
errors.

## Use

```
kairoq dc project.yaml [<input>...] [--view]
```

prints the port voltages and currents, the dissipation, the holes with the largest
currents and the places with the largest current densities, with the mesh size, the passes
and the estimated error. `--view` opens the viewer with the solution: voltage (scaled per
net), current density and power density on each layer, a probe for values at a point, and
the port results.

```python
board = kairoq.open("fab/")
r = kairoq.solve_dc(board, {"ports": [...]})   # or the path of a project file
r.ports["LOAD"]["voltage_v"]
layer = r.layer(1)        # points_mm, triangles, voltage_v, sheet_current_a_per_mm, power_w, net
kairoq.serve(board, dc=r)
```

## Resources

Memory grows with the mesh: about 1 kB per triangle in the measurements of
`docs/validation/dc.md` (2.4 GB for 2.3 million triangles). `max_triangles` bounds it; when
the budget stops refinement, the result says so and reports the error estimate reached.
Large runs are best started with a memory limit of their own, so that the operating system
stops only that run if the budget is set too high.

## Limits

- The current is uniform through the thickness of each layer; vertical current inside thick
  copper and the current distribution within a pad (as opposed to through it) are not
  modelled. Ports are equipotential.
- Hole walls are equipotential on each layer. Filled or capped vias, and a component pin in
  a hole, are modelled as the plating alone.
- Alone, the solver assumes a uniform temperature; the temperature field comes from the
  coupled thermal solution.
- Copper without a net attribute is treated as copper of no net; plated holes join copper of
  their net only, so a short between nets in the copper data is not conducted through holes
  (the importer reports such shorts).
