# PDN impedance

The impedance a load sees into a power distribution network, from DC to a few GHz: the
pairs of power and ground planes with their resonances, the decoupling capacitors and the
regulator. Validation: `docs/validation/pdn.md`.

## Model

- **Plane pairs.** A layer on which the power net covers more than `plane_fraction` of
  the board is a power plane; with the nearest ground plane above and below it forms a
  pair. The pair is the region both planes cover; holes smaller than 2.5 mm across
  (antipads, clearances) are filled, larger ones (splits, slots) are walls.
- **The cavity.** Between the planes, $h$ apart, the voltage $V(x, y)$ obeys

  $$\nabla^2 V + k^2 V = 0, \qquad k^2 = -z\,y, \quad z = j\omega\mu_0 h + 2 Z_s, \quad y = \frac{j\omega\varepsilon_0\varepsilon_r(\omega)}{h},$$

  with an open (magnetic) wall at the region's outline. $\varepsilon_r$ is the series
  permittivity of the dielectrics between, $Z_s$ the surface impedance of the plates (their
  sheet resistance at low frequency, the skin effect above). A via joining a branch to the
  power plane is a port: a circle of its radius through which its current enters.
- **Solution.** Linear finite elements on the region give the modes
  $K\varphi = \lambda M\varphi$ up to three times the highest wavenumber of the sweep, by
  shift-invert Lanczos with an exact (envelope Cholesky) factorisation. The impedance
  between ports $p$ and $q$ is

  $$Z_{pq} = z \left( \sum_m \frac{\bar\varphi_m(p)\,\bar\varphi_m(q)}{\lambda_m - k^2} + R_{pq} \right),$$

  with $\bar\varphi_m(p)$ the mode's mean over the port's circle and $R$ the static sum of
  the modes not kept. The constant mode gives the plates' capacitance
  $\varepsilon A / h$ at low frequency.
- **Branches.** Each load, capacitor and the regulator joins the power plane at the via of
  the power net nearest one of its pads (within `via_reach_mm`), or at its pin on the
  plane layer. A capacitor is $R + j\omega(L + L_m) + 1/(j\omega C)$, the regulator
  $R + j\omega(L + L_m)$, both from the port to the planes. The mounting inductance $L_m$,
  unless given, is that of the path from the pair's nearer plane up the power via, along
  the track to the power pad, across to the ground pad, along its track and down the
  ground via, over that plane: half the inductance of the path closed by its image
  (Neumann's formula over its straight pieces, vias of their radius, pads and tracks of a
  quarter of the pad's width). The part's own inductance is its ESL.
- **Pairs in parallel.** The pairs' admittance matrices add at the ports that reach them;
  the branches add on the diagonal; the loads' impedances are the inverse's diagonal.

## Project file

```yaml
pdn:
  power_net: +5V
  ground_net: GND
  f_min_hz: 1e5
  f_max_hz: 2e9
  target_ohm: 0.1
  loads: [U1]
  regulator: {component: J101, r_mohm: 20, l_nh: 2}
  capacitors:
    - {component: C12, c_uf: 22, esr_mohm: 3, esl_nh: 0.6}
  other_capacitors: {c_uf: 10, esr_mohm: 5, esl_nh: 0.5}
```

| key under `pdn` | default | meaning |
| --- | --- | --- |
| `power_net`, `ground_net` | | the two nets |
| `f_min_hz`, `f_max_hz`, `points` | $10^5$, $10^9$, 300 | logarithmic sweep |
| `target_ohm` | | target impedance; the summary lists the bands above it |
| `loads` | | components at which $Z(f)$ is reported |
| `capacitors` | | `component`, `c_uf`, `esr_mohm`, `esl_nh`, `mount_nh` (estimated when absent) |
| `other_capacitors`, `capacitor_prefix` | , `C` | values for every other two-pad part between the nets whose reference starts with the prefix |
| `regulator` | | `component`, `r_mohm` (1), `l_nh` (1) |
| `plane_fraction` | 0.1 | share of the board that makes a net's copper a plane |
| `via_reach_mm` | 2 | distance from a pad to the via that joins it to the planes |
| `dielectric` | 4.3, 0.02 | `epsilon_r`, `loss_tangent`, `frequency_hz`, `wideband` where the stackup gives none |

```
kairoq pdn project.yaml [<input>...] [-o dir]
```

prints the pairs, the branches with their mounting inductances, each load's impedance, its
resonances and the bands above the target, and with `-o` writes `<load>_z.csv`.

## Limits

- The region's outline is an open wall: radiation from the board's edges is not modelled,
  so resonances of thick pairs on small boards are sharper than in reality.
- A branch joins the planes at one via; the planes' spreading between several vias of one
  part is not resolved, and traces between a pad and its via are part of $L_m$ only through
  the via's height.
- Vias of other nets that cross the pair are filled in.
- Capacitor values are given, not read from a bill of materials; DC bias derating of
  ceramic capacitors is the user's.
