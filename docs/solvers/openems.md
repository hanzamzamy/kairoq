# openEMS export

A loop of the EMI estimates (`docs/solvers/emi.md`) can be exported to openEMS, an
open-source finite-difference time-domain (FDTD) solver (<https://openems.de>), to check
its inductance and radiated field with a full-wave simulation of the same geometry.
KairoQ writes the model and a script; openEMS, installed separately, runs it.
Validation: `docs/validation/emi.md`.

## Model

The export writes the geometry of the loop model, as the loop solver builds it:

- **Window.** The copper within the window around the loop, on the loop's layers and the
  next `depth_layers` ones, with the plated holes near the loop
  (`docs/solvers/emi.md`, Discretisation). Away from the loop the outlines are simplified
  and small holes closed, as in the estimate.
- **Copper** is perfectly conducting sheets at the layers' mid-planes. CSXCAD polygons have
  no holes; a polygon with holes is cut into strips by lines through the holes.
- **Plated holes** are perfectly conducting cylinders of the drill's radius between the
  outermost layers they join.
- **Dielectric** fills the window over the stack's dielectric layers, with their
  permittivity and loss tangent (the conductivity at half the highest frequency); copper
  layers take the permittivity of the dielectric next to them towards the inside. Solder
  mask is left out.
- **Elements** are perfectly conducting bars of square section, side twice the element's
  radius, standing on the copper from pad centre to pad centre; their centre line lies
  where the estimate's conductor does. A package inductance given for an element is not
  modelled, and the estimate written with the model leaves it out.
- **Port.** The first element has a gap of one minimum cell at its middle, bridged on
  its axis by a lumped port with a resistance (50 Ω by default): one line of mesh edges,
  so that the port's current probe encloses all of its current. The first element must
  lie along x or y; other elements that do not are thin wires on the mesh lines.
- **Stitches** of the loop are elements too; one with a capacitance has a gap like the
  port's, bridged by a lumped capacitance, and must lie along x or y. Its series
  resistance and inductance are not part of the model.

The mesh has lines at the copper edges and one minimum cell to either side, at the planes
of the sheets, at the faces of the elements and at the port. Lines closer than a minimum
cell are merged, except the sheets' planes and the port's lines, which the mesh holds
exactly (a sheet off its line by rounding is left out of the grid). Between them the cells grow
by at most 1.4 from one to the next, up to `board_cell_mm` within the window and a
twentieth of the shortest wavelength in the air around it. The air reaches a quarter of
that wavelength beyond the window; the boundary is an absorbing layer of eight cells.

The excitation is a Gaussian pulse covering 0 to `f_max_hz`. The simulation stops when the
energy in the domain has fallen by `end_criterion_db`.

## Results

`simulate.py` writes `result.json`:

- the port impedance $Z(f) = V(f) / I(f)$ at 100 frequencies from 1 % of `f_max_hz`, and the
  inductance $L(f) = \operatorname{Im} Z / (2 \pi f)$;
- at each far-field frequency, the largest field over a sphere of radius `distance_m`
  (on a 3° grid), per ampere of port current, from a near-field box ten cells inside the
  boundary;
- the equivalent moment, the moment a small loop would need to radiate that field,
  $m_\text{eq} = 4 \pi r E / (\eta_0 k^2 I)$, to compare with the estimate's moment;
- the estimate: the inductance without package inductances, the moment, and the far field
  of the solved currents (`docs/solvers/emi.md`) at the far-field frequencies, computed as
  the near-to-far-field transformation gives it (the far-field term scaled to the distance).

Below the loop's self-resonance, $L(f)$ is flat and compares with the estimate's
high-frequency limit; towards the resonance it rises. The far field includes what the
estimate leaves out: the electric dipole of the charges on the loop and the board, and
the current's variation around the loop towards its resonance.

## Project file

```yaml
emi:
  loops:
    - name: BUCK_INPUT
      elements: [...]
  openems: {f_max_hz: 1e9, min_cell_mm: 0.1}
```

| key under `emi.openems` | default | meaning |
| --- | --- | --- |
| `f_max_hz` | $10^9$ | highest frequency of the excitation |
| `min_cell_mm` | 0.1 | cell at copper edges, sheets and elements; the port gap |
| `board_cell_mm` | 1 | largest cell within the window |
| `air_margin_mm` | a quarter wavelength at `f_max_hz` | air around the window |
| `port_resistance_ohm` | 50 | resistance of the lumped port |
| `far_field_hz` | `f_max_hz` / 10, / 4 and / 2 | frequencies of the far field |
| `end_criterion_db` | −40 | energy decay at which the simulation stops |
| `max_timesteps` | $10^6$ | |
| `estimate` | true | include the loop estimate in the model |

The window, the layers and the plated holes follow `loop_settings`.

## Use

```
kairoq openems project.yaml [<input>...] -o loop_dir [--loop NAME]
python loop_dir/simulate.py [--threads N] [--setup-only]
```

The first command writes `model.json` (the geometry in metres, the settings and the
estimate) and `simulate.py`, and prints a summary; without `--loop` the first loop is
exported. The script needs openEMS with its Python modules (CSXCAD, openEMS), numpy and
h5py; `--setup-only` builds the model and the mesh and prints the cell count without
running. The simulation files go to `loop_dir/sim/`.

```python
model = kairoq.export_openems(board, "project.yaml", "loop_dir", loop="BUCK_INPUT")
model["estimate"]["inductance_nh"]
```

## Cost

The cell count and the number of time steps both grow as the minimum cell shrinks: the
time step is proportional to the smallest cell, and a loop at 1 GHz with 0.1 mm cells
needs some $10^5$ steps. Measured times are in `docs/validation/emi.md`. A coarser
`min_cell_mm` or a smaller window (`window_margin_mm`) reduces the cost.
