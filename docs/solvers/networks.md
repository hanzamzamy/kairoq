# Networks: S-parameters, mixed-mode parameters, TDR

A channel is a set of pads of a board taken as ports, with the copper of their nets
between them: tracks as transmission lines, vias, and the parts and loads given. Its
S-parameters over frequency, the mixed-mode parameters of its pairs, and the
time-domain reflectometry (TDR) of its ports follow, and are written as Touchstone files.
Validation: `docs/validation/networks.md`.

## Model

- **Junctions.** Track ends within 1 µm of each other on a layer are one node. A track that
  ends on the inside of another of its net and layer (a stub routed onto a line) meets it
  there: the line is split at that point.
- **Tracks.** The nets of the channel are cut as in the signal-integrity analysis
  (`docs/solvers/si.md`); each sample stands for an equal share of its track's length as a
  uniform line with the per-unit-length parameters of its cross-section at each frequency.
  A track without any usable cross-section (no reference along its whole length) is a
  series inductance of 1 nH/mm; a sample without one takes its nearest neighbour's.
- **Coupling.** Where a sample's cut holds copper of other nets of the channel (a pair's
  partner, an aggressor beside a victim), those nets' samples nearest to that copper are
  matched with it, and together they become one section of several conductors with the
  coupled parameters of the cut. Where their lengths differ, the remainder of each longer
  one is a single line. Unmatched samples are single lines.
- **Vias.** A plated hole of a channel net has a node on each copper layer it spans; track
  ends and pads on that layer at the hole join it. All reference planes are the ground of
  the matrix, and the barrel is a ladder between them (the physics-based via model):
  - nodes at each copper layer and at slices of longer gaps (a twentieth of the shortest
    wavelength);
  - each node's capacitance to the planes: its share of the barrel's charge, with its
    pads, in an axisymmetric field solution of the barrel, the pads and the planes with
    their antipads, all in their dielectrics;
  - each segment between two planes carries the pair's fundamental (TEM) mode: the voltage
    it adds is the pair's impedance between the barrels times their currents, each
    weighted by its share of the pair's height. On a pair whose modes up to three times
    the sweep's wavenumber are few (300), the impedance is that of the pair's real outline
    (`docs/solvers/pdn.md`), with the reference vias that join it short-circuited; on a
    larger pair, that of infinite planes, $\frac{\omega\mu h}{4} H_0^{(2)}(k\rho)$ between
    barrels at distance $\rho$, with the eight reference vias nearest each barrel
    short-circuited;
  - each segment's own inductance from the pair's higher modes,
    $L_{ij} = \frac{\mu}{\pi h} \sum_{n \ge 1} c_{ni} c_{nj} K_0(g_n a) I_0(g_n a)$,
    $g_n = n\pi/h$, $c_{ni}$ the integral of $\cos(n\pi z/h)$ over segment $i$;
  - outside the pairs, the wire's partial inductance with its image in the nearest plane.

  A change of return plane is the voltage the barrel's current leaves in the pairs between.
  `via_model: coaxial` takes instead, between consecutive layers, a coaxial line of
  $Z_v = \frac{60}{\sqrt{\varepsilon_r}} \ln\frac{D}{d}$ ($d$ the drill, $D$ the antipad's
  diameter of equal area, three drills where no plane has one), the barrel beyond the
  connected layers an open stub; vias without a plane around them always take it.
- **Pads** join the track ends inside them; their capacitance is not modelled.
- **Parts** join two pads: a series resistance, inductance and capacitance, or a two-port
  Touchstone file. **Loads** connect a pad to the reference: a resistance and a
  capacitance in parallel.

Each line section of $N$ conductors with per-unit-length matrices $Z = R + j\omega L$ and
$Y = G + j\omega C$ has the transfer matrix

$$\begin{pmatrix} V(0) \\ I(0) \end{pmatrix} = \exp\!\left( \ell \begin{pmatrix} 0 & Z \\ Y & 0 \end{pmatrix} \right) \begin{pmatrix} V(\ell) \\ I(\ell) \end{pmatrix},$$

in closed form for one conductor ($\cosh\gamma\ell$, $Z_c \sinh\gamma\ell$,
$\sinh\gamma\ell / Z_c$) and by scaling and squaring of a Padé approximant for two. Its
admittance blocks are stamped into the nodal admittance matrix of the channel; with unit
currents into the ports the matrix gives the port impedance matrix $\mathbf Z$, and

$$\mathbf S = (\mathbf Z_n - \mathbf I)(\mathbf Z_n + \mathbf I)^{-1}, \qquad \mathbf Z_n = \mathbf Z / Z_0,$$

with the reference impedance $Z_0$ of the single-ended ports (50 Ω by default).

**Mixed-mode parameters.** For port pairs $(p, n)$ the differential and common modes are
$(a_p - a_n)/\sqrt2$ and $(a_p + a_n)/\sqrt2$; the mixed-mode matrix lists the
differential ports first, then the common ones. Its references are $2Z_0$ and $Z_0/2$
(100 Ω and 25 Ω for 50 Ω ports).

**Step responses.** For each port in `drive`, the waveform at every port for a step into
it (the same band-limited step), as a fraction of the step: the near- and far-end
crosstalk at a victim's ends, the received edge at the far end of a line.

**TDR.** The step response of each port's reflection (each differential port's) is
computed from the sweep: a Gaussian step of the given rise time (at least $0.9/f_\text{max}$,
where its spectrum has fallen to 9 %, so that the truncated sweep does not ring), the
lowest frequency's reflection extended to DC. The impedance is
$Z(t) = Z_\text{ref}\,(1 + \rho(t))/(1 - \rho(t))$; time runs from the step at the port, a
discontinuity at delay $\tau$ appears at $2\tau$.

## Project file

```yaml
si:
  nets: [...]                               # as for the analysis; channels choose their own nets
  channels:
    - name: PCIE_TX
      ports: [U1.B5, U1.B6, J2.A12, J2.A13]  # pads, one single-ended port each
      differential: [[1, 2], [3, 4]]        # port numbers of the pairs
      impedance_ohm: 50
      parts:
        - {from: C5.1, to: C5.2, c_nf: 100, l_nh: 0.4}
        - {from: L3.1, to: L3.2, touchstone: ferrite.s2p}
      loads: [{pad: U4.7, r_ohm: 50, c_pf: 1.5}]
  network: {f_max_hz: 20e9, points: 400, rise_ps: 35}
```

| key | default | meaning |
| --- | --- | --- |
| `channels[]` | | `name`, `ports` (pads, or points `NAME@X,Y[,LAYER]` in mm: a port on the end of the nearest track on that layer, within 2 mm), `differential` (port pairs, from 1), `impedance_ohm` (50), `parts` (`from`, `to`, `r_ohm`, `l_nh`, `c_nf` in series, or `touchstone`), `loads` (`pad`, `r_ohm`, `c_pf`), `drive` (ports, from 1, to step) |
| `network` | | `f_max_hz` ($2 \times 10^{10}$), `points` (400: a linear sweep from `f_max_hz`/`points`), `vias` (true), `via_model` (`physical`, or `coaxial`), `rise_ps` (35), `tdr_ns` (twice the summed delay of the channel's nets, at least 1 ns) |

The nets of a channel are those of its ports, and the nets its parts lead to; other
pads of these nets are open unless a load is given.

## Use

```
kairoq sparams project.yaml [<input>...] [-o dir]
```

prints per channel its size, the return and insertion loss (single-ended, or differential
and the conversion to common mode) at 1, 2.5, 5, 10 and 20 GHz, and the range of its TDR;
with `-o` it writes `NAME.sNp` and, for pairs, `NAME_mixed.sMp`.

```python
r = kairoq.analyse_channels(board, {"channels": [...], "network": {...}})[0]
r["frequency_hz"], r["s"], r["mixed"], r["tdr"], r["touchstone"]
t = kairoq.read_touchstone("part.s2p")
```

## Limits

- The lines are the quasi-TEM cross-sections of the signal-integrity analysis, with
  their limits; bends, the ends of tracks and pads are not discontinuities of their own.
- The via model (`docs/validation/networks.md`): against full-wave solutions it places the
  stub notch of a via through thin pairs of planes with return vias 7 % high, and that of
  a stub in a single 2.8 mm dielectric on a small board without return vias 19 % low
  (the coaxial model: 56 % and 20 % high); radiation from the board's edges is not in it.
  The channel's export to openEMS (`docs/solvers/rf.md`) gives the full-wave answer.
- Nets are coupled only where their cuts find each other on one layer; where the routes
  separate, each is a single line, and coupling between layers is not modelled.
- The time responses come from a sampled spectrum and repeat with the period
  `points`/`f_max_hz` (20 ns by default); a longer window is reported in a warning and
  needs more points.
- Touchstone input: version 1, S-parameters, two-ports.
- Copper drawn as areas (zones, polygons) is cut only where tracks run over it as a
  reference: structures drawn as copper areas, such as some filters, are not lines of the
  model.
