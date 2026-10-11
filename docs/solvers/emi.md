# EMI estimates

First-order estimates of the board's contribution to electromagnetic interference:

- the inductance of switching loops and their magnetic moment;
- the radiated field of a loop's current harmonics;
- traces whose return path in a reference plane is interrupted;
- the capacitance of switch nodes to a nearby metal plate.

They rank layouts and loops and show where to look; they do not predict the result of a
compliance measurement. Validation: `docs/validation/emi.md`, and against full-wave
simulations of the exported loop model (`docs/solvers/openems.md`) `docs/validation/emi.md`.

## Loop inductance

### Model

A loop is a closed chain of elements, each a part that carries the loop current from one
of its pads to another: the input capacitor, the switches of a converter, a diode. Between
consecutive elements the current flows in the copper of the net that joins them. The
inductance is computed by partial inductances (the PEEC method, Ruehli 1972) in the
quasi-static high-frequency limit, where the current distributes itself to store the least
magnetic energy: at the switching edges of power converters (tens of megahertz and above)
the skin depth is below the copper thickness and the sheet resistance is small against the
reactance.

- **Copper** is a set of thin sheets at the layers' mid-planes, discretised by triangles
  with the lowest-order divergence-conforming (Rao-Wilton-Glisson) basis: one unknown per
  interior edge, the current across it.
- **Plated holes** are tubes between the layers they join; their current enters and leaves
  the copper uniformly over the hole's area on each layer.
- **Elements** are straight round conductors between the centres of their pads, at their
  radius above (or below) the copper; their current enters and leaves the copper uniformly
  over each pad, as a lead soldered over the pad. A package inductance given per element
  (the equivalent series inductance inside the part) adds in series.
- **Other copper** in the window around the loop (planes, pours) carries the currents the
  loop induces in it: a plane under the loop holds the image current that lowers the
  loop's inductance and cancels most of its moment.

With $I$ the vector of all branch currents and $L$ the partial inductance matrix, the
current minimises $\tfrac12 I^\top L I$ subject to Kirchhoff's current law at every triangle
and a current of 1 A in the first element; the minimum is half the loop inductance. The
entries of $L$ are

$$L_{mn} = \frac{\mu_0}{4\pi} \iint \frac{\mathbf f_m(\mathbf r) \cdot \mathbf f_n(\mathbf r')}{|\mathbf r - \mathbf r'|}\, dA\, dA',$$

integrated over the supports of the basis functions $\mathbf f$: the inner integral in closed
form (Wilton et al. 1984; Graglia 1993), the outer by a 7-point rule, or 28 points for
neighbouring triangles. Straight conductors use the closed forms for parallel filaments
(Grover) and a tube for their own partial inductance, and Gauss quadrature between
conductors at an angle. The constrained minimum is solved in a basis of loops, which
satisfies the current law by construction: a spanning tree of the RWG branches in each
connected sheet, and one of the conductors between the sheets, carry the source's ampere
back to it ($\mathbf x_p$), and every other branch with the tree path that closes it is a
loop (a column of $M$). The branch currents $\mathbf x = \mathbf x_p + M \mathbf y$ minimise
$\mathbf x^T L \mathbf x$ over $\mathbf y$,

$$M^T L M\, \mathbf y = -M^T L\, \mathbf x_p,$$

a dense system of the loops alone, about a third of the branches, solved by Cholesky
factorisation; the inductance is $\mathbf x^T L \mathbf x$. Where a conductor's end spans
two sheets, the branches and the current law are solved together instead (Cholesky
factorisation of $L$ and of the Schur complement of the current law); both give the same
solution.

Every solution is an upper bound of the model's inductance: it minimises the energy over a
subset of the admissible currents.

The dense matrix $L$ of all branches grows with the square of their number, 18 GB at
32 000 triangles. Above `dense_matrix_mb` it is not kept: in the basis of loops only
$M^T L M$, $M^T L\, \mathbf x_p$ and $\mathbf x_p^T L\, \mathbf x_p$ are needed, and they
are linear in $L$, so the rows that a block of triangles contributes to are formed,
multiplied by $M$ and dropped. The inductance is then
$\mathbf x_p^T L\, \mathbf x_p + \mathbf y^T M^T L\, \mathbf x_p$. This takes about
twice the kernel evaluations, as the symmetry of $L$ is not used across blocks, and the
memory of the loops' matrix and one block; the solution is the same.

### Discretisation

- **Window.** A DC solution of 1 A through each copper connection of the loop marks the
  corridor in which the loop current flows; the copper within `window_margin_mm` of the
  corridor and of the pads is included, on the layers of the loop's pads and the next
  `depth_layers` ones. Other nets' islands smaller than `passive_min_area_mm2` and plated
  holes farther than `via_radius_mm` from the corridor are left out.
  Over a plane the moment depends on how much of the plane the window holds: the margin
  doubles (up to `max_window_doublings` times) until the moment changes by less than
  `moment_tolerance` or the window holds the board; an unsettled moment is reported.
- **Element size.** It grows from a minimum at the corridor (half the smallest pad, between
  0.1 and 0.5 mm) by 0.35 per unit distance. Away from the corridor, holes smaller than half
  the local size (antipads in planes) are closed and outlines simplified to 0.15 of it.
- **Edges.** In the high-frequency limit the sheet current grows as $d^{-1/2}$ towards a
  copper edge. Rows of vertices along the edges near the corridor, 20 µm inside and then
  three times as far each, up to half the local element size, give elements thin across the
  edge and long along it.
- **Refinement.** The first mesh uses about a quarter of the triangle budget; each further
  pass refines where the sheet current recovered at the vertices differs most from the
  computed one, until the budget (`max_triangles`) is used. The inductance of the previous
  pass is reported with the result; the remaining error is at least the difference, and at
  small budgets several times it (below). The element size of the first mesh at the loop is
  reported too; where the budget forces the coarsest elements allowed (four times the
  minimum), the result warns that the inductance is likely several per cent high.

### Accuracy and cost

Against closed forms (`docs/validation/emi.md`) the inductance converges from above: within
2 % for parallel strips, microstrip and a rectangular loop at the densities used there, and
within 0.5 % of image theory over a plane. The cost grows with the cube of the number of
loops, about a third of the branches: seconds up to 5000 triangles (measured times in the
validation note); the default budget is 4000.

The convergence is slow: the current across a trace and its return under it in the plane,
with the current singular at every copper edge, are resolved only gradually. A larger
budget always gives a lower, better inductance; with the default the result is typically 2
to 5 % high. The inductance (nH) against the budget, for loops over a plane:

| budget (`max_triangles`) | 4000 | 8000 | 16 000 | 32 000 | 64 000 | reference |
| --- | --- | --- | --- | --- | --- | --- |
| 30 mm trace over a slotted plane, 0.5 mm (`docs/validation/emi.md`) | 21.54 | 21.39 | 21.07 | 20.91 | 20.85 | about 20.82 (limit) |
| ring 20 × 10 mm, 1.5 mm over a plane | 36.87 | 36.59 | 36.20 | 35.96 | | 35.66 (openEMS) |
| the same ring, 0.2 mm over a plane | 18.00 | 17.65 | 17.37 | 16.98 | | |
| time (slotted plane) | 1.4 s | 4.4 s | 20 s | 140 s | 633 s | |

Below about 10 000 triangles such loops change little from pass to pass while they are
still several per cent high, so the difference to the previous pass understates the error
by up to tenfold there; doubling the budget and comparing is the more reliable check. Above
`dense_matrix_mb` the memory grows linearly with the triangles (1.1, 2.0 and 5.3 GB at
16 000, 32 000 and 64 000 above).

### Resistance and inductance over frequency

With `frequencies_hz` given, the loop's impedance $Z(f) = R + j\omega L$ is solved at each
frequency on the last mesh of the settled window. The loop basis gives

$$\bigl(j\omega M^T L M + M^T R(f) M\bigr)\, \mathbf y = -\bigl(j\omega M^T L + M^T R(f)\bigr)\, \mathbf x_p, \qquad Z = \mathbf x^T \bigl(j\omega L + R(f)\bigr)\, \mathbf x,$$

solved by an $LDL^T$ factorisation of the complex symmetric matrix. $R(f)$ is assembled per
layer from the overlaps $\int \mathbf f_m \cdot \mathbf f_n \, dA$ of the basis functions,
times the surface impedance of a sheet of thickness $t$ carrying current on both faces,

$$Z_s = \frac{1 + j}{2 \sigma \delta} \coth \frac{(1 + j)\, t}{2 \delta}, \qquad \delta = \frac{1}{\sqrt{\pi f \mu_0 \sigma}},$$

which is $1 / (\sigma t)$ at DC. A zero-thickness sheet's current grows without bound
towards its edges, and with it the resistance as the mesh refines; copper of thickness $t$
spreads it over about $t / 2$. A triangle nearer an edge than $t / 2$ therefore counts with
the share $d / (t / 2)$ of its overlap, blended in by
$\chi = |Z_s - 1/(\sigma t)| / |Z_s|$, from 0 at DC to 1 in the skin-effect limit. Plated
holes add the impedance of their wall, of thickness `plating_um`, on the circumference.
At low frequency the current spreads to minimise resistance and the inductance rises
above the high-frequency limit; at high frequency it tends to that limit.

### Stitches

A stitch is a part that offers the return current a path outside the loop's own chain: a
capacitor across a slot or between two planes, or a strap. It is a straight conductor
between its pad centres, like an element, with a series capacitance, resistance and
inductance (`stitches` of a loop); its nets' copper in the window is modelled whatever
its area. In the high-frequency limit a stitch is its inductance alone, and the share of
the loop current it carries is reported. With `frequencies_hz`, its series impedance
$R + j\omega L + 1/(j\omega C)$ joins the loops' impedance in the solve at each
frequency, which gives the loop's resistance and inductance with the stitch and the
stitch's current, in magnitude and phase. Below the series resonance of a capacitor with
the loop around it the stitch carries little; at its resonance it carries more than the
loop current, circulating through the detour; above it, it tends to its high-frequency
share. The return-path check of the next section finds the crossings; the stitches show
whether a part placed at one carries the return.

## Magnetic moment and radiated field

The loop's magnetic moment per ampere is

$$\mathbf m = \frac12 \int \mathbf r \times \mathbf J\, dV,$$

over all computed currents, the induced ones included. A loop over a plane therefore has a
small moment: the plane's image current cancels it. A part's conductor runs from pad centre
to pad centre one radius above the copper, while its current enters the sheet spread over
the pad; a short piece joins each end to the centroid of that spread, so that the current
is closed and the moment does not depend on the origin. With a waveform given, the loop current
is a periodic trapezoid (amplitude, frequency, duty cycle at half amplitude, rise and fall
times) whose harmonics follow exactly from its piecewise-linear shape.

Each harmonic radiates with the distribution of the solved currents: every triangle's
sheet current times its area, and every conductor in pieces of at most 1 mm, is a current
element $I\,\mathbf{dl}$, with the joining pieces. The current is closed, so the scalar
potential vanishes and

$$\mathbf E = -j\omega \mathbf A, \qquad \mathbf A(\mathbf r) = \frac{\mu_0}{4\pi} \sum_i I_i\, \mathbf{dl}_i\, \frac{e^{-jkR_i}}{R_i}, \quad R_i = |\mathbf r - \mathbf r_i|, \quad k = \frac{2\pi f}{c_0},$$

near field included and with the phase of every element, so loops comparable with the
wavelength are not reduced to their moment. A closed current has no symmetric first moment,
$\int (r_i J_j + r_j J_i)\, dV = -\int r_i r_j\, \nabla \cdot \mathbf J\, dV = 0$; pads and
triangles taken as points leave one of the order of their size squared, comparable with
the moment of a loop over a plane, and its field is subtracted. The field per ampere is computed on a
logarithmic grid of 24 frequencies per decade and interpolated to the harmonics. Three
sites are available (`field`):

- **free** (default): free space, the largest field at distance $r$ over directions 5° apart
  in elevation and 10° in azimuth around the elements' centre;
- **ground**: a measurement site. The board lies horizontal at `board_height_m` over a
  reflecting ground, which adds the image currents (horizontal components reversed); the
  antenna is at distance $r$ from the board's centre, its height scanned over
  `antenna_heights_m` in steps of 0.25 m, the board turned in steps of 10°, and the larger
  of the horizontal and vertical polarisations is taken;
- **small_loop**: the far field of the moment alone, in the loop's plane,

  $$E = \frac{\eta_0 k^2 |\mathbf m| I_n}{4 \pi r},$$

  which is $1.316 \times 10^{-14}\, f^2 A I / r$ V/m for a loop of area $A$ (Ott,
  *Electromagnetic Compatibility Engineering*).

The currents are those of the quasi-static solution: uniform around the loop. Near the
loop's self-resonance the current is no longer uniform and the field per ampere at the
source is higher (`docs/validation/emi.md`). Common-mode currents on cables, usually the
larger source, are not included.

## Return paths

A trace's return current flows in the nearest plane beneath it. Reference nets are those
given, or those whose copper fills a quarter of the board on some layer; a plane layer is
one on which they do. For every checked net (all nets but the references, unless given),
on every layer, the trace's copper without its pads and holes (and `pad_clearance_mm`
around them, where the trace changes layer or meets a part) is compared with the
reference copper of the nearest plane layer above and below. Each piece of trace over no
reference copper is a finding:

- **crossing**: the piece separates two parts of the reference under the trace: a slot, a
  row of merged antipads;
- **split**: the same, with different reference nets on either side;
- **partial**: the reference under the trace remains connected around the piece: an
  antipad under the trace's edge, the end of a plane.

The span is the piece's area divided by the trace's width, that is its extent along the
trace. Crossings and splits are listed first, by span.

## Switch-node coupling

The capacitance of a net (a switch node) to a metal plate on a face of the board, and to
the other nets around it, by the method of moments in the stratified medium of the stack:

- **Medium.** The plate is a conducting plane; between it and the board lie the insulation
  (`gap_mm`, `epsilon_r`), then the board's layers: dielectrics, mask and air beyond. A
  copper layer is a sheet at its middle, half its thickness joined to the slab on each
  side.
- **Conductors.** The node's copper over the plate and the copper of every other net within
  a window around it (`window_margin_mm`), each net one conductor at 0 V; the node at 1 V.
  Away from the node the elements grow by 0.35 per unit distance, other copper is
  simplified to the local size (antipads smaller than half of it closed), and a piece gives
  up what larger pieces on its layer already cover.
- **Green's function.** The potential of a unit charge in the medium follows from its Hankel
  transform, which a transmission-line recursion through the slabs gives in closed form:

  $$G(\rho; z, z') = \frac{1}{2\pi} \int_0^\infty \tilde\varphi(k; z, z')\, J_0(k\rho)\, k\, dk.$$

  The quasi-static direct term $C / (2\pi R)$ and, within a layer, the images of its
  nearest interface at depths $2nd$ are extracted in closed form; the remainder is
  integrated numerically and tabulated against $\rho$.
- **Solution.** Uniform charge on each triangle, Galerkin weighting: the direct and image
  terms over each triangle in closed form (Wilton et al. 1984), the remainder by quadrature,
  and Cholesky factorisation. The capacitance to each net is minus its charge, and to the
  plate the rest of the node's charge. The charges converge from below, in proportion to
  the element size at the edges; a second solution with elements twice as large
  extrapolates them, and both are reported.

With a rate of change of the node's voltage, the common-mode current $C\, dV/dt$ follows.
The parallel-plate capacitance of the node's area facing the plate, not covered on the
layers between by copper of any net,

$$C_\text{pp} = \varepsilon_0 \sum_\text{layers} \frac{A}{\sum_k t_k / \varepsilon_k},$$

and the parallel-plate overlaps with other nets on the neighbouring layers are reported
alongside. Relative permittivities are those of the stackup, or 4.3 for dielectric and
3.8 for mask where it gives none.

## Project file

```yaml
board: ../fab/
emi:
  loops:
    - name: BUCK_INPUT
      elements:                                   # in the direction of the current
        - {from: Q1.D, to: Q2.S, esl_nh: 0.5}     # the switches
        - {from: C5.2, to: C5.1, esl_nh: 0.4, radius_mm: 0.3}
      waveform: {amplitude_a: 5, frequency_hz: 500e3, duty: 0.3, rise_ns: 4, fall_ns: 6}
  loop_settings: {max_triangles: 4000, distance_m: 3}
  return_path: true                               # or {nets: [...], reference_nets: [GND]}
  plates:
    - {name: CASE, face: bottom, area: [{rect_mm: [10, 10, 60, 40]}], gap_mm: 0.5, epsilon_r: 4}
  switch_nodes:
    - {net: SW, plate: CASE, dv_dt_v_per_ns: 20}
```

| key under `emi` | default | meaning |
| --- | --- | --- |
| `loops` | | `name`, `elements` (`from`, `to` as `REF.PIN`; `radius_mm` 0.2; `esl_nh` 0), `stitches` (`from`, `to`; `c_nf` 0 for none; `esr_mohm` 0; `esl_nh` 0; `radius_mm` 0.2), `waveform` (`amplitude_a`, `frequency_hz`, `duty` 0.5, `rise_ns`, `fall_ns` = `rise_ns`) |
| `loop_settings` | | `window_margin_mm` (10), `depth_layers` (2), `via_radius_mm` (3), `passive_min_area_mm2` (10), `min_size_mm` (from the pads), `max_triangles` (4000), `estimate_error` (true), `distance_m` (3), `max_frequency_hz` ($10^9$), `moment_tolerance` (0.05; 0: one window), `max_window_doublings` (3), `dense_matrix_mb` (2048), `frequencies_hz` (none), `plating_um` (25), `field` (`free`; `ground`, `small_loop`), `board_height_m` (0.8), `antenna_heights_m` ([1, 4]) |
| `return_path` | | `true`, or `nets`, `reference_nets`, `plane_fraction` (0.25), `min_area_mm2` (0.002), `pad_clearance_mm` (0.3) |
| `plates` | | `name`, `face` (`bottom`), `area` (the whole face), `gap_mm` (0), `epsilon_r` (1) |
| `switch_nodes` | | `net`, `plate`, `dv_dt_v_per_ns` |
| `coupling_settings` | | `window_margin_mm` (three times the node's height over the plate, 2 to 30 mm), `max_panels` (4000), `min_size_mm` (a quarter of the node's distance to the plate or the nearest other copper layer) |
| `dielectric_epsilon_r`, `mask_epsilon_r` | 4.3, 3.8 | where the stackup gives none |
| `openems` | | export of a loop to openEMS (`docs/solvers/openems.md`) |

The elements of a loop must close: the second pad of each element and the first pad of the
next are on one net. Unknown keys are errors.

## Use

```
kairoq emi project.yaml [<input>...] [--view]
```

prints per loop the inductance (with that of the previous refinement pass and the part
inside the packages), with `frequencies_hz` its resistance and inductance at each frequency,
the magnetic moment, the largest radiated field overall and in the
bands 30-230 MHz and 230-1000 MHz; the return-path findings; and the switch nodes'
capacitances. `--view` shows the loop currents and the findings in the viewer.

```python
board = kairoq.open("fab/")
e = kairoq.analyse_emi(board, {"loops": [...], "return_path": True})   # or a project file
e.loops["BUCK_INPUT"]["inductance_nh"]
e.loops["BUCK_INPUT"]["harmonic_field_v_per_m"]
e.loops["BUCK_INPUT"]["impedance"]        # (frequency_hz, resistance_ohm, inductance_nh)
e.loops["BUCK_INPUT"]["field_per_ampere"] # (frequency_hz, V/m per ampere) at the site
e.loops["BUCK_INPUT"]["current_positions_m"], e.loops["BUCK_INPUT"]["current_moments_a_m"]
e.return_path            # list of findings
e.coupling               # list of switch nodes
```

## Limits

- The reported loop inductance is the high-frequency limit; at frequencies where the sheet
  resistance is comparable with the reactance (low megahertz on thin traces) the current
  spreads more and the inductance is higher. `frequencies_hz` gives it at given frequencies.
- Copper is thin: the current on the edges' side walls is not modelled. In the skin-effect
  limit the resistance is that of zero-thickness sheets; for strips 1 mm wide and 35 µm
  thick at 1 GHz it is 33 % below the cross-section solver's and 61 % below Wheeler's rule
  with Gupta's thickness correction. Elements are straight conductors, their internal
  geometry is the given package inductance, and their resistance is not included.
- The window, the layer depth and the via radius limit the model to the loop's
  surroundings; planes cut off by the window hold less image current (`docs/validation/emi.md`).
- The radiated field is that of the board's loop currents alone, without cables or an
  enclosure, and of quasi-static currents, uniform around the loop; the ground site is an
  ideal reflecting plane and a single antenna position per height and angle, not a
  calibrated site.
- The return-path check is geometric: whether a stitching capacitor or a strap carries the
  return across a gap follows from a loop with that part as a stitch. The loop's solution
  is quasi-static: the copper's capacitance, and with it the line's own resonances, is not
  part of it.
- The switch-node capacitance takes the plate as infinite, the layers as extending beyond
  the board, the mask as continuous, and every other net as held at 0 V; plated holes
  are not part of the conductors. The window cuts other nets' copper: the node's field
  then reaches the plate around the cut edges, and for a node that planes shield from
  the plate, where the plate's share is a small difference of large charges, the value is
  an upper estimate that falls as the margin grows.
