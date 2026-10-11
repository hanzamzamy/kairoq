# RF structures

Passive RF analysis of a board: de-embedding of S-parameters, full-wave S-parameters of
a channel's structure with openEMS, and closed-form estimates of microstrip patch
antennas. Validation: `docs/validation/rf.md`.

## De-embedding

S-parameters measured or simulated at one reference plane are moved to another. The
operations apply in order to a Touchstone file.

**Reference impedance.** With every port's real reference changed from $Z_0$ to $Z_1$,

$$S' = (I - \rho S)^{-1} (S - \rho I), \qquad \rho = \frac{Z_1 - Z_0}{Z_1 + Z_0}.$$

**A two-port in front of a port.** A fixture $X$ whose port 2 joins port $p$ of the
network and whose port 1 becomes the new port $p$ gives, with
$\Delta = 1 - X_{22} S_{pp}$,

$$S'_{pp} = X_{11} + \frac{X_{12} X_{21} S_{pp}}{\Delta}, \quad S'_{pj} = \frac{X_{12} S_{pj}}{\Delta}, \quad S'_{ip} = \frac{S_{ip} X_{21}}{\Delta}, \quad S'_{ij} = S_{ij} + \frac{S_{ip} X_{22} S_{pj}}{\Delta}.$$

Removing a fixture of ABCD matrix $A$ connects the two-port of $A^{-1}$.

**Lines.** A uniform quasi-TEM line of impedance $Z_c$, effective permittivity
$\varepsilon_\text{eff}$, attenuation $\alpha_1$ at 1 GHz scaling with $\sqrt{f}$ (conductor
loss) and length $\ell$ has the ABCD matrix of $\gamma = \alpha + j\omega\sqrt{\varepsilon_\text{eff}}/c_0$;
a negative length removes it (the inverse of a line is the line of opposite length). The
line's parameters come from the cross-section analysis of the board
(`docs/solvers/si.md`).

**2x-thru.** Two copies of a symmetric fixture joined back to back have the ABCD matrix
$M = A^2$. The fixture is the square root

$$A = \frac{M + sI}{t}, \qquad s^2 = \det M, \quad t^2 = \operatorname{tr} M + 2s,$$

of the four roots the one nearest the previous frequency's, starting from the identity at
zero frequency; the branch is followed through electrical lengths of several
wavelengths. The thru is referred to the network's reference and interpolated to its
frequencies. The fixture must be symmetric (its two ends alike); an asymmetric fixture
needs another method (IEEE 370 describes several).

```yaml
input: measured.s2p
output: dut.s2p
operations:
  - line: {port: 1, length_mm: -5, z_ohm: 50, epsilon_eff: 3.3, loss_db_per_m_1ghz: 3}
  - thru: {port: 2, file: thru.s2p}     # port 1 of the thru at the measured port
  - renormalise_ohm: 100
```

```
kairoq deembed spec.yaml
```

```python
r = kairoq.deembed({"input": "measured.s2p", "operations": [...]})
r["frequency_hz"], r["s"], r["reference_ohm"]
```

## Full-wave S-parameters with openEMS

The channel solution (`docs/solvers/networks.md`) is quasi-TEM: lines, vias and parts
joined as a network. Where that model does not hold (discontinuities large against the
wavelength, radiation, structures that are not lines), the channel's structure is exported
to openEMS, an open-source FDTD solver (`docs/solvers/openems.md`), and solved full-wave.
The two paths answer the same question: the native one in seconds, within its model; the
full-wave one in minutes to hours, with the geometry as drawn.

**The model** holds, within a margin around the copper of the ports' nets (clipped to the
board outline):

- all copper, of every net, as sheets: outer layers on their face towards the board,
  inner layers at mid-plane; polygons without holes (cut into strips) and without
  collinear vertices; conducting sheets of the layer's conductivity and thickness
  (openEMS's model of a thin lossy conductor);
- the plated holes as cylinders between their layers;
- the dielectric layers with the SI analysis's materials; the script fits the
  Djordjevic-Sarkar model of each to six Debye poles from $f_\text{max}/200$ to
  $f_\text{max}$ (within about $10^{-3}$) and runs openEMS from a file with these Debye
  materials, as its Python interface has none; the outer copper layers as air;
- per port, where the track that ends at the pad runs straight along x or y, a microstrip
  port over its first part (up to ten widths, at least 2 mm): fed at the pad through the
  channel's reference resistance, its voltage and current measured halfway along and
  moved back to the pad with the line's measured propagation constant, attenuation
  included; elsewhere a lumped port of the reference resistance, a vertical box from the
  pad to the nearest copper of another net beneath (or above) its centre, square and as
  wide as the track (`--port-kind lumped` gives lumped ports everywhere).

**The script** builds the structure with CSXCAD, meshes it (lines on every sheet plane and
port face exactly, at copper edges with a minimum cell to either side, tracing the
circles of the pads and antipads around the channel's vias (at four angles) and halfway
between them, graded to a
twentieth of the wavelength in the densest dielectric, at most four minimum cells, inside
the model and to a twentieth of the free-space wavelength in the air margin), excites one
port at a time with a Gaussian pulse from 0 to $f_\text{max}$ and computes

$$S_{ij} = \frac{b_i}{a_j}, \qquad a = \frac{V + Z_0 I}{2}, \quad b = \frac{V - Z_0 I}{2}$$

from the ports' voltages and currents. Ports not excited (`--ports`) have their columns
from reciprocity. `--losses simple` gives perfect copper and, for each dielectric, the
conductivity of its loss tangent at the band's middle (exact there only). With far-field frequencies, port 1's radiation follows from a
near-field box: the directivity $D$, the radiation efficiency
$\eta = P_\text{rad}/P_\text{acc}$, the gain $\eta D$ and the realised gain
$D P_\text{rad}/P_\text{inc}$.

```yaml
si:
  channels: [{name: FEED, ports: [J1.1]}]
  openems: {f_max_hz: 4e9, points: 400, margin_mm: 20, min_cell_mm: 0.3, far_field_hz: [2.45e9]}
```

| key of `si.openems` | default | meaning |
| --- | --- | --- |
| `f_max_hz`, `points` | $2 \times 10^{10}$, 201 | band of the excitation; frequencies of the result |
| `margin_mm` | 3 | the model's extent around the ports' nets' copper |
| `min_cell_mm` | a quarter of the thinnest dielectric, at most 0.1 | cell at copper edges and ports |
| `board_cell_mm` | see above | largest cell inside the model |
| `air_margin_mm` | a quarter wavelength at $f_\text{max}$ | air between the model and the absorbing boundary |
| `far_field_hz` | none | frequencies of port 1's radiation |
| `end_criterion_db`, `max_timesteps` | −40, 500 000 | when the simulation stops |

```
kairoq openems project.yaml --channel FEED -o dir
python dir/simulate.py [--threads N] [--ports 1,2]
```

writes `dir/result.sNp` and `dir/result.json`. The lumped ports include the vertical path
from the reference to the pad: their inductance and the pad's capacitance are part of the
result, and can be removed by de-embedding a 2x-thru of the same ports.

## Microstrip patch antennas

A rectangular patch of width $W$ (the radiating edges) and length $L$ on a substrate of
$\varepsilon_r$ and height $h$, by the transmission-line model (Balanis, *Antenna Theory*,
chapter 14):

$$\varepsilon_\text{eff} = \frac{\varepsilon_r + 1}{2} + \frac{\varepsilon_r - 1}{2}\left(1 + 12\frac{h}{W}\right)^{-1/2}, \qquad \frac{\Delta L}{h} = 0.412\,\frac{(\varepsilon_\text{eff} + 0.3)(W/h + 0.264)}{(\varepsilon_\text{eff} - 0.258)(W/h + 0.8)},$$

$$f_r = \frac{c_0}{2 (L + 2\Delta L)\sqrt{\varepsilon_\text{eff}}}.$$

For a design frequency, $W = \frac{c_0}{2 f_r}\sqrt{2/(\varepsilon_r + 1)}$ and $L$ follows from
$f_r$. Each radiating edge is a slot of conductance $G_1$; the two couple through $G_{12}$:

$$G_1 = \frac{1}{120\pi^2}\int_0^\pi \left[\frac{\sin\left(\frac{k_0 W}{2}\cos\theta\right)}{\cos\theta}\right]^2 \sin^3\theta\, d\theta, \qquad G_{12} = \frac{1}{120\pi^2}\int_0^\pi \left[\frac{\sin\left(\frac{k_0 W}{2}\cos\theta\right)}{\cos\theta}\right]^2 J_0(k_0 L \sin\theta) \sin^3\theta\, d\theta,$$

the input resistance at an edge is $R_\text{in} = 1/(2(G_1 + G_{12}))$, and a feed inset by
$y_0$ sees $R_\text{in}\cos^2(\pi y_0/L)$. The directivity is that of the two slots, $W$ long
and $L + 2\Delta L$ apart, over an infinite ground plane, by integration over the upper
half-space; the bandwidth for a VSWR below 2 is
$3.771\,\frac{\varepsilon_r - 1}{\varepsilon_r^2}\,\frac{W}{L}\,\frac{h}{\lambda_0}$ (Jackson and
Alexopoulos, 1991), for thin substrates.

```python
kairoq.patch_antenna(epsilon_r=4.5, height_mm=1.535, frequency_hz=2.45e9)
kairoq.patch_antenna(4.5, 1.535, width_mm=37.0, length_mm=28.5, feed_ohm=50)
```

The full-wave path for an antenna on the board is the openEMS export of a channel with
one port at the feed and far-field frequencies.

## Limits

- De-embedding: real reference impedances; the 2x-thru method needs symmetric fixtures.
- openEMS models: copper without roughness; parts and loads of the channel are not in
  the model. A microstrip port needs a straight track along x or y at the pad; the
  patch's full-wave results (`docs/validation/rf.md`) are with a lumped port, and its
  directivity differs by 1.1 dB with a microstrip port on the feed (open).
- The patch model is the transmission-line model: thin substrates ($h$ below about
  $0.05\,\lambda$ in the dielectric, reported), the fundamental mode only, no feed
  reactance, an infinite ground plane.
- Passive structures only; no measured validation.
