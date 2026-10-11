# Thermal solver

The thermal solver computes the temperature of a board, steady or over time: on every
copper layer, at every heat source and at the junctions behind them, with the heat leaving
through each face and into heat sinks. Heat comes from components and from the Joule
heating of the DC current, whose copper conductivity depends on the temperature
(electrothermal coupling). Validation: `docs/validation/thermal.md` (steady state) and
`docs/validation/thermal-transients.md` (heat sinks and transients).

## Model

### Conduction sheets

The board is represented by a stack of conduction sheets, each a two-dimensional field
$T(x, y)$ over the whole board outline (non-plated holes and cut-outs removed):

- one sheet on every copper layer;
- in every dielectric, further sheets: the slices start at `dielectric_first_slice_mm`
  (default 0.025 mm) at both faces of the dielectric and double in thickness towards its
  middle, up to `dielectric_slice_mm` (default 0.2 mm), which fills the middle;
- a sheet on an outer face that is dielectric rather than copper.

A slice of thickness $\tau$ and conductivities $k_{xy}$ (in-plane) and $k_z$ (through the
thickness) gives half of its in-plane conductance $k_{xy} \tau$ to the sheet above it and
half to the sheet below it, and joins them by the conductance per area $k_z / \tau$. On a
copper sheet, the copper adds $k_\mathrm{Cu}\, t$ where there is copper; between the copper
of an inner layer, the resin adds $k_{xy}\, t$ (with the conductivity of the dielectric below
the layer). Between the copper of an outer layer is air, which adds nothing, or with
`outer_gaps: mask` solder mask filling it to the copper's thickness, which adds
$k_\mathrm{mask}\, t$. Each sheet $i$
therefore satisfies

$$-\nabla \cdot (G_i \nabla T_i) + \frac{T_i - T_{i-1}}{R_{i-1}} + \frac{T_i - T_{i+1}}{R_i} = q_i ,$$

with $G_i$ the in-plane conductance (W/K), $R_i$ the resistance per area to the next sheet
(m² K/W) and $q_i$ the heat input per area.

With one sheet per copper layer only, the model is the usual layered (2.5D) approximation.
It treats the dielectric between two copper layers as a single vertical resistance and
underestimates the temperature of sources that are small compared with the dielectric
thickness, by two thirds for a 0.5 mm source over 1.5 mm of dielectric. The intermediate
sheets resolve the spreading through the thickness. Heat from a small source spreads
mostly within a distance of the order of its size from the copper, hence the thin slices
there. At the default slicing the model is within 2.5 % of a three-dimensional solution
for traces and for sources of 2 mm and more (`docs/validation/thermal.md`). For a 0.5 mm pad
the three-dimensional temperature itself depends by about 11 % on what fills the gaps of the
outer copper (mask or air), and for a 2 mm pad of 10 oz copper by 21 %; with the matching
`outer_gaps` the model is within 2 % of either for sources of 2 mm and more
(`docs/validation/thermal.md`).

### Surfaces

The outermost sheets lose heat through the solder mask (series resistance
$t_\mathrm{mask} / k_\mathrm{mask}$) to the surface, which is cooled by convection and
radiation:

$$q = h\,(T_s - T_a) + \varepsilon \sigma_\mathrm{SB}\,(T_s^4 - T_a^4) ,$$

with $h$ (`convection_w_m2k`) and $\varepsilon$ (`emissivity`) set per face; the surface
temperature $T_s$ follows from the flux balance through the mask. The sides of the board are
adiabatic. Radiation is to surroundings at the ambient temperature; view factors and
reflections are not modelled.

### Plated holes

A plated hole joins, by its barrel, the copper of the layers it touches (copper at its
centre). Between consecutive touched layers $a$ and $b$ the barrel is the conductance

$$G = \frac{k_\mathrm{Cu}\, \pi\, t_p\,(d + t_p)}{L_{ab}}$$

($+\,2 t_p \ell$ in the area for a slot of length $\ell$), $L_{ab}$ being the distance between
the layer mid-planes. The hole wall is isothermal on each layer. Fill and pins in holes
are not modelled.

### Heat sources and fixed temperatures

A source is an area selected as a DC port is (pad, component with or without net, or a
disc at a point on a layer), with either

- a power: distributed over the area in proportion to the area, or, with
  `junction_k_per_w`, entering a junction node joined to the area by that resistance (the
  junction temperature is reported); or
- a temperature: the area is isothermal at that temperature (a clamp, a soldered
  connection to a large body); the heat it takes out is reported.

For transients, a source can follow a power schedule (below) and carry a heat capacity: of
its junction node when it has one, otherwise spread over its area. The heat capacity of a
package (its body and die) is what delays the junction temperature after a step of power.

Component power can also come from a table (a bill of materials exported with a power
column, or a list written for the purpose): one source per component, over all its pads,
with the power and optionally the junction-to-case resistance and the heat capacity read
from named columns. A cell may name several components ("R1, R2"), each of which then
dissipates the stated power; a value may carry a unit (`250 mW`). A component named in the
table but absent from the board is an error, unless such entries are to be ignored. A
source written in the project file replaces the table's entry for the same component.

### Heat sinks

A heat sink is a body pressed onto a face of the board, over the whole face or over
rectangles and polygons of it, through a thermal interface material (a pad, a paste, a gap
filler) of resistance per area

$$r_\mathrm{TIM} = \frac{t_\mathrm{TIM}}{k_\mathrm{TIM}} .$$

Under the sink the face exchanges heat with the sink through the solder mask and the
interface, $q = (T - T_\mathrm{sink}) / (r_\mathrm{mask} + r_\mathrm{TIM})$, instead of with the
air. The sink is either held at a temperature (a cold plate, a chassis much larger than the
board) or a body with one temperature, a resistance to ambient $R_\mathrm{sa}$ (as heat sink
datasheets state it, convection and radiation included) and, for transients, a heat
capacity. One sink may cover several rectangles and polygons of its face. The heat into
each sink and its temperature are reported. The sink body is isothermal: spreading within a
plate is not modelled.

### Electrothermal coupling

When the project has DC ports, the DC problem is solved with the conductivity at the local
copper temperature,

$$\sigma(T) = \frac{\sigma_{20}}{1 + \alpha\,(T - 20\ \mathrm{^\circ C})}, \qquad \alpha = 0.00393\ \mathrm{K^{-1}},$$

its Joule heat (copper and barrels) is added to the thermal problem, and the two are
iterated until the temperature changes by less than `coupling_tolerance_k` (default
0.02 K) between iterations. The DC mesh of the first iteration is kept for the following
ones. A current that heats the copper beyond the point where the losses rise faster than the
cooling has no steady state; the iteration then does not settle and says so.

## Discretisation

Each sheet is meshed separately over the board outline with linear triangles; copper
outlines, hole walls and source areas are mesh edges on the copper sheets. The sheets are
coupled on non-matching meshes: every node couples, over half its lumped area, to the point
directly below (or above) it on the neighbouring sheet, so that the coupling matrix is
symmetric. Hole walls and fixed-temperature areas tie their nodes to one unknown.

The radiation and the mask are linearised about the previous temperature (Newton's
method) and the linear system is solved again until the temperature settles, together with
the electrothermal iteration. The linear systems are solved by conjugate gradients with
algebraic multigrid, as for the DC solver.

### Transients

The heat capacity per area of each sheet is that of its copper (where there is copper), of
the resin between the copper of an inner layer, and half of each adjacent dielectric slice;
the outer sheets also carry the solder mask. Junctions and sink bodies carry their own
capacities. With $C$ the lumped (diagonal) capacity and $K$ the conductance matrix, each step
of length $\Delta t$ solves (backward Euler)

$$\Big(\frac{C}{\Delta t} + K\Big)\, T^{n+1} = \frac{C}{\Delta t}\, T^{n} + f(t^{n+1}),$$

which is stable for any step and first-order accurate. The faces' convection and radiation
are taken as $q(T^n) + g\,(T^{n+1} - T^n)$ with a fixed slope $g$ per node (that of the steady
solution at the peak power), so that the matrix depends on the step length only and its
multigrid hierarchy serves every step of that length; with convection alone this is exact.

The run starts from ambient or from the steady solution at $t = 0$. The mesh is that of the
steady solution at the peak power of every source. Steps start at `step_s` and double
every `steps_per_doubling` steps up to `max_step_s`; they land on every point of a power
schedule, after which they start again at `step_s`. With electrothermal coupling the DC
solution is renewed whenever the temperature has changed by more than `dc_update_k` since
the last one. Each step records the temperature of every source (mean, maximum, junction),
of every sink, the highest temperature and the heat balance (heat in, heat out, rate of
change of the stored heat); the fields of all copper layers are kept at chosen times.

### Adaptive refinement

The mesh is refined adaptively by gradient recovery, as in the DC solver
(`docs/solvers/dc.md`); the gradient is recovered per material (copper and dielectric
separately), since the in-plane conductance jumps at the copper edges. The vertical
coupling between the sheets is lumped at the nodes, and its error is estimated too: for a
linear temperature difference $w$ between two sheets over a triangle of area $A$, the lumped
energy exceeds the exact one by $G A \sum_{i<j} (w_i - w_j)^2 / 12$, with $G$ the
conductance per area between the sheets. This term dominates near small sources, where the
thin dielectric slices next to the copper, stiffly coupled and nearly without in-plane
conduction, change steeply under coarse elements: on a region of a real board, the mean
rise of a 2 × 2 mm package is 1.4 % low at the default target, and 6.4 % low on the mesh
that the in-plane estimate alone gives (`docs/validation/thermal.md`, "Real board against
Elmer"). The error
of the lumped surface terms, of order $(h_e/\ell)^2$ for elements of size $h_e$ and a
thermal length $\ell$, is not included. The mesh of the last adaptive pass is kept for the
coupling iterations.

## Project file

```yaml
board: ../fab/
dc:                                   # optional: electrothermal coupling
  ports: [...]
thermal:
  ambient_c: 25
  top:    {convection_w_m2k: 10, emissivity: 0.9}
  bottom: {convection_w_m2k: 10, emissivity: 0.9}
  dielectric: {k_xy: 0.8, k_z: 0.3}
  layers: {D1: {k_xy: 0.35, k_z: 0.35}}
  sources:
    - {name: U1, component: U1, power_w: 2.0, junction_k_per_w: 10}
    - {name: Q3, component: Q3, power_schedule: [[0, 0], [0, 3.0], [20, 3.0], [20, 0]],
       junction_k_per_w: 1.5, heat_capacity_j_per_k: 0.4}
    - {name: CLAMP, pad: J1.1, temperature_c: 40}
  component_power: {file: power.csv, junction_column: Rjc}
  heatsinks:
    - name: PLATE
      face: bottom
      area: [{rect_mm: [10, 10, 60, 40]}]      # omitted: the whole face
      tim: {thickness_mm: 0.5, k: 3.0}          # or {resistance_k_mm2_per_w: 170}
      resistance_k_per_w: 2.0                   # or temperature_c: 40
      heat_capacity_j_per_k: 120
  transient:
    duration_s: 600
    step_s: 0.5
    max_step_s: 20
```

| key under `thermal` | default | meaning |
| --- | --- | --- |
| `ambient_c` | 25 | ambient temperature (air and surroundings) |
| `top`, `bottom` | $h = 10$, $\varepsilon = 0.9$ | `convection_w_m2k`, `emissivity` of each face |
| `copper_k` | 385 | W/(m K) |
| `dielectric` | `k_xy` 0.8, `k_z` 0.3 | default dielectric (FR-4), W/(m K) |
| `layers` | | conductivities per stackup layer name |
| `mask_k` | 0.25 | solder mask, W/(m K) |
| `outer_gaps` | `air` | what fills the space between the copper of the outer layers: `air` or `mask` |
| `dielectric_first_slice_mm` | 0.025 | slices next to the faces of a dielectric |
| `dielectric_slice_mm` | 0.2 | largest slice in a dielectric |
| `heat_capacity` | `copper` $3.45 \times 10^6$, `dielectric` $2.2 \times 10^6$, `mask` $1.8 \times 10^6$ | J/(m³ K); per layer as `heat_capacity` under `layers` |
| `plating_thickness_um` | 18 | hole wall plating |
| `sources` | | heat sources and fixed temperatures, as above |
| `component_power` | | a table, or a list of them: `file`, `reference_column`, `power_column`, `junction_column`, `heat_capacity_column`, `ignore_unknown` |
| `heatsinks` | | heat sinks, as above; `face` is `bottom` unless given |
| `transient` | | `duration_s` (required), `step_s` (`max_step_s` / 64), `max_step_s` (`duration_s` / 50), `steps_per_doubling` (8), `start` (`ambient` or `steady`), `snapshots_s` (20 times, logarithmically spaced), `dc_update_k` (1), `step_tolerance` ($10^{-4}$: each step solves for the change of the temperature to this tolerance relative to its residual; 0: the whole state to the solver tolerance) |
| `coupling` | true | electrothermal coupling when the project has DC ports |
| `coupling_tolerance_k`, `max_coupling_iterations` | 0.02, 30 | iteration of the coupling and of radiation |
| `target_error`, `max_passes`, `max_triangles` | 0.1, 6, 5 000 000 | adaptive refinement, as for the DC solver |
| `mesh` | `max_edge_mm` 2 | initial mesh, keys as for the DC solver |
| `solver` | | `tolerance` ($10^{-10}$), `max_iterations` (2000) |

A source gives a power (`power_w`, `power_schedule` or both) or `temperature_c`;
`junction_k_per_w` and `heat_capacity_j_per_k` go with a power. A schedule lists
`[time_s, power_w]` points, linear in between and constant beyond its ends; a step is two
points at the same time. A steady solution uses `power_w`, or else the schedule's value at
$t = 0$. A component selected without a net covers all its pads.

The reference column of a table is found as `Reference`, `Designator`, `Ref`, `RefDes` or
`References`, the power column as `Power`, `Power_W`, `Power (W)`, `P` or `P_W`, unless the
project names them; fields are separated by commas, semicolons or tabs. Unknown keys are
errors.

## Use

```
kairoq thermal project.yaml [<input>...] [--view]
```

prints the mesh, the iterations, the hottest point, the heat balance, per source the mean
and maximum temperature of its area and its junction temperature, and per heat sink its
heat and temperature. For a transient these are the values at the end, followed by the
highest temperature reached by each source and when. `--view` shows the temperature on
each copper layer in the viewer, with the DC solution when coupled; for a transient, a time
control selects the kept fields and a plot shows the temperatures over time.

```python
board = kairoq.open("fab/")
r = kairoq.solve_thermal(board, {"sources": [...]}, dc={"ports": [...]})  # or a project file
r.sources["U1"]["junction_c"]
r.heatsinks["PLATE"]["heat_w"]
layer = r.layer(1)        # points_mm, triangles, temperature_c, copper
tr = r.transient          # time_s, max_c, source_mean_c / max_c / junction_c [source, step],
                          # heatsink_c / heatsink_w, heat_in_w / out_w / stored_w, snapshot_time_s
r.snapshot(3, 1)          # temperature of layer 1 at tr["snapshot_time_s"][3]
```

## Resources

A transient costs one multigrid set-up per distinct step length (a few per run) and one
solution per step, plus a steady solution for the mesh. Every sheet is a mesh of the whole
board, so the size grows with the number of sheets:
a two-layer board of 1.6 mm has thirteen. Sheets away from the sources stay coarse, since
refinement follows the error. Memory is about 2 kB per triangle in a coupled solution
(2.5 GB for 1.3 million triangles, `docs/validation/thermal.md`); `max_triangles` bounds it.

## Limits

- Transients use the mesh of the steady solution at the peak power: sharp transient
  gradients shortly after a step, before the heat has spread as far as in the steady
  state, are resolved only as finely as that mesh allows.
- Backward Euler is first-order accurate in time; the step length is chosen by the
  schedule above, not by an error estimate.
- With radiation, the faces' slopes are those of the peak state: the linearisation is
  consistent but not exact at other temperatures, an error of the order of the step.
- A heat sink body is isothermal, and the interface is uniform; contact pressure, voids and
  spreading within the sink are not modelled.
- The sheets resolve the conduction through the thickness of the dielectric, not of the
  copper: within a copper layer the temperature is uniform through the thickness, and the
  thickness of a copper layer adds no vertical resistance (also where resin fills it).
- Components are areas of heat on the board; their bodies, and the heat they conduct or
  radiate, are not modelled. A component that dissipates mostly through its body to the
  air is better represented by the part of its power that enters the board.
- Convection is a given coefficient per face, uniform over the face.
- Radiation is to surroundings at the ambient temperature, without view factors.
