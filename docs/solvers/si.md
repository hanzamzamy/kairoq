# Signal integrity: transmission lines

The tracks of a board as transmission lines: the per-unit-length resistance, inductance,
conductance and capacitance of each cross-section over frequency, and from them the
impedance profile, delay and loss of each net, the differential impedance of pairs, and
the places where a line leaves its target impedance or loses its reference plane.
Validation: `docs/validation/si.md`.

This is the first part of the high-speed analysis. S-parameters of nets and channels,
crosstalk, time-domain response (reflection, ringing) and eye diagrams build on these
per-unit-length parameters and follow as further parts.

## Scope and limits of the high-speed analysis

| area | covered | not covered |
| --- | --- | --- |
| transmission lines | quasi-TEM cross-sections of any number of conductors in layered dielectrics, frequency-dependent losses, microstrip dispersion | other non-TEM effects: valid while the cross-section is small against the wavelength (tens of GHz for printed-board geometry) |
| discontinuities | via barrels and their stubs (`docs/solvers/networks.md`); full-wave by the openEMS export | bends, pads and track ends as lines |
| materials | wideband dielectric model, conductor roughness | results are as good as the permittivity, loss tangent and roughness given; defaults are reported |
| validation | closed forms, published data, full-wave simulation | measured hardware |

## Cross-sections

### Model

A cross-section is a two-dimensional cut: conductors (copper of any shape) in regions of
dielectric, the rest vacuum, inside a grounded box well beyond the geometry. One
conductor, or a group, is the reference (the return); the others are signal conductors.
In the quasi-TEM approximation the per-unit-length parameters follow from two
electrostatic problems.

With each signal conductor $i$ in turn at 1 V and all others at 0 V, the potential
$\varphi_i$ solves $\nabla \cdot (\varepsilon \nabla \varphi_i) = 0$. The Maxwell
capacitance matrix is the energy

$$C_{ij} = \varepsilon_0 \int \varepsilon_r \nabla\varphi_i \cdot \nabla\varphi_j \, dA,$$

and the same with $\varepsilon_r = 1$ everywhere gives $C^{(0)}$, from which the external
inductance is

$$L = \mu_0 \varepsilon_0 \left(C^{(0)}\right)^{-1}.$$

The potentials are linear finite elements on a triangle mesh of the cut, solved by
conjugate gradients with algebraic multigrid. The mesh is graded from the signal
conductors and refined where a recovered gradient differs most from the computed one,
until the estimated relative energy error meets a target: the capacitance error is about
its square (0.08 gives about 0.2 % in impedance, 0.02 about 0.01 %). The capacitance of a
conforming solution with set potentials converges from above.

### Dielectric losses

Each dielectric has a complex permittivity
$\varepsilon_r(f) = \varepsilon'(f) - j\varepsilon''(f)$, by default the
Djordjevic-Sarkar model through the given permittivity and loss tangent at a frequency:

$$\varepsilon_r(\omega) = \varepsilon_\infty + \frac{\Delta\varepsilon}{\ln(\omega_2/\omega_1)}\,
\ln\frac{\omega_2 + j\omega}{\omega_1 + j\omega},$$

with $\omega_1 = 10^4$ and $\omega_2 = 10^{12}$ rad/s. It is causal (its real and imaginary
parts satisfy the Kramers-Kronig relations) and gives a nearly constant loss tangent
across the band, with the permittivity falling slowly with frequency. A constant model is
available. $C(f)$ is solved with $\varepsilon'(f)$, and the conductance follows to first
order in the loss:

$$G_{ij}(f) = \omega \varepsilon_0 \int \varepsilon''(f)\, \nabla\varphi_i \cdot \nabla\varphi_j \, dA.$$

Over a frequency range the field of the solution at the first frequency is kept: the
energy of each dielectric region, $\varepsilon_0 \int_\text{region} \nabla\varphi_i \cdot \nabla\varphi_j\, dA$,
is weighted with the region's $\varepsilon'(f)$ and $\varepsilon''(f)$ at each frequency. The
line calculator's option `exact_dispersion` solves instead at two frequencies per decade and
interpolates in $\log f$; the two agree within $2 \times 10^{-4}$ in capacitance on a masked
microstrip from 1 MHz to 20 GHz (`docs/validation/networks.md`). The derivation is in
`docs/methods/si.md`.

### Microstrip dispersion

A microstrip's field moves into the substrate as the frequency rises, which raises its
effective permittivity and its impedance; a quasi-TEM cross-section has no such effect.
Where a cut's track lies on an outer layer with one reference layer, the analysis applies
the formulas of Kirschning and Jansen (1982) for the effective permittivity and of Jansen
and Kirschning (1983) for the impedance of a microstrip of the track's width $w$ over the
dielectric between, of height $h$ and of the series permittivity $\varepsilon_r$ of its
layers:

$$\varepsilon_\text{eff}(f) = \varepsilon_r - \frac{\varepsilon_r - \varepsilon_\text{eff}(0)}{1 + P(w/h, \varepsilon_r, f h)},$$

with $\varepsilon_\text{eff}(0)$ the cross-section's own value. The capacitance matrix
moves accordingly from $C$ towards that of a uniform fill, $\varepsilon_r C_0$:

$$C' = \varepsilon_r C_0 - \frac{\varepsilon_r C_0 - C}{1 + P},$$

which applies to a pair and its neighbours alike and leaves a line in one dielectric
unchanged. $C'$ and the external inductance are then scaled by the same factor in opposite
directions, so that the track's impedance follows the formula of the impedance. The
formulas are fits to full-wave solutions of the bare microstrip for
$0.1 \le w/h \le 100$ and $1 \le \varepsilon_r \le 20$; above $f h = 20$ GHz mm the
correction is held at its value there, and channels whose sweep goes beyond it say so.

### Conductor losses

At high frequency the current flows on the conductor surfaces with the distribution of the
charge of the vacuum solution. With the surface impedance
$Z_s = (1+j)/(\sigma\delta)$, $\delta = 1/\sqrt{\pi f \mu_0 \sigma}$, the series impedance of
conductor $k$ is $Z_s \oint_k J_i J_j \, dl$ for unit currents $i$ and $j$. The surface
integral equals Wheeler's incremental inductance rule,

$$\oint_k J_i J_j \, dl = \frac{1}{\mu_0} \frac{\partial L_{ij}}{\partial n_k},$$

the change of inductance when the surface of conductor $k$ recedes. It is evaluated in the
domain form of the shape derivative: an area integral over the elements near the
conductor, with a velocity field that moves its surface outward by unit distance and
vanishes before the other conductors. Unlike the surface charge itself, which is singular
at corners, the area integral converges as the energy does.

Between the direct-current resistance and the skin-effect limit each conductor's surface
impedance is that of a slab of effective thickness $t_\text{eff}$,

$$Z_s(f) = K(f)\, \frac{1+j}{\sigma\delta}\, \coth\frac{(1+j)\,t_\text{eff}}{\delta},$$

with $t_\text{eff}$ chosen so that the low-frequency limit is the conductor's
direct-current resistance; this reproduces both limits exactly and the transition
approximately. The imaginary part is the internal inductance. $K(f)$ is the roughness
factor: none; Hammerstad, $1 + \tfrac{2}{\pi}\arctan\!\left(1.4\,(\Delta/\delta)^2\right)$
for an RMS height $\Delta$; or Huray, $1 + r / (1 + \delta/a + \delta^2/2a^2)$ for spheres of
radius $a$ and a ratio $r$ of sphere surface to flat area.

### Line quantities

For line $i$ (the others at rest):
$Z_0 = \sqrt{(R + j\omega L)/(G + j\omega C)}$,
$\gamma = \alpha + j\beta = \sqrt{(R + j\omega L)(G + j\omega C)}$,
$\varepsilon_\text{eff} = (\beta c_0/\omega)^2$, delay $\beta/\omega$, loss $8.686\,\alpha$ dB/m,
of which $R/2Z_0$ in the conductors and $GZ_0/2$ in the dielectric. For a pair, with the
self terms averaged and $C_m$ the (negative) mutual Maxwell capacitance:
$Z_\text{odd} = \sqrt{(L_s - L_m)/(C_s - C_m)}$,
$Z_\text{even} = \sqrt{(L_s + L_m)/(C_s + C_m)}$,
$Z_\text{diff} = 2 Z_\text{odd}$ and $Z_\text{common} = Z_\text{even}/2$.

## Board analysis

### Cuts

Each track of the analysed nets is cut across at intervals (at most `step_mm`, clear of
its ends by a track width, where bends and pads are). A cut takes the copper of every layer
along a line across the track:

- **shields**: the nearest layers above and below with reference copper across the track
  and a track width to either side; layers beyond them are left out;
- **reference**: the reference nets (given, or those whose copper fills a quarter of a
  layer) and any other copper at least 1 mm and five track widths wide;
- **partner**: for a pair, the partner's copper on the same layer;
- **neighbours**: other copper on the track's layer within three times the distance to the
  farther shield, each a quiet line; other copper in the cut is held at the reference
  potential, which leaves the track's capacitance unchanged and its inductance nearly so.

The cut reaches twenty times the distance to the farther shield beyond the track and its
partner: a microstrip's field falls as the inverse square of the lateral distance, and
there the impedance is within about 0.1 % of an unbounded cut's (`docs/validation/networks.md`).
The dielectric extends over the board only, to its outline. The dielectric layers keep
their materials, inner copper layers are filled with the material next to them, and the
solder mask covers the outer layers conformally. Etching, when given, narrows each
conductor's face away from its laminate (on inner layers, the upper face). Cuts with the same geometry (track and partner
to 1 µm, other copper to 5 µm) share one solution, and the distinct cross-sections are
solved in parallel.

### Results

Per net: its length, delay and loss at the reference frequency, the range and
length-weighted mean of its impedance and, for a pair, of its differential impedance where
the partner is coupled, and the coupled length. Distances along a net run from its first
pad through its tracks, vias and pads. Findings, each a run of consecutive samples of at
least `min_finding_mm`:

| finding | where |
| --- | --- |
| impedance, differential impedance | outside the target and its tolerance |
| no reference | no reference copper in the cut (a gap, a slot, a board without planes) |
| reference change | the shield layers change, for example through a via |
| uncoupled | a pair member without its partner in the cut |
| crosstalk | a victim's near- or far-end crosstalk from an aggressor above the threshold |

### Crosstalk

Where a cut holds a quiet line of another net beside the track (on its layer, within three
times the distance to the farther shield), the coupled matrices of the cut give the
backward and forward coupling coefficients of weak coupling,

$$K_b = \frac14\left(\frac{L_m}{L_s} + \frac{C_m}{C_s}\right), \qquad
K_f = \frac12\left(\frac{C_m}{C_s} - \frac{L_m}{L_s}\right)\sqrt{L_s C_s},$$

with $L_s$, $C_s$ the track's own and $L_m$, $C_m$ the mutual inductance and (the magnitude
of the) mutual capacitance per length. Summed over the length the two nets run together,
for a step of rise time $t_r$ into the aggressor and matched ends:

- near-end crosstalk $K_b \min(1, 2\tau/t_r)$, with $\tau$ the delay of the coupled length:
  it saturates once the round trip outlasts the rise;
- far-end crosstalk $\lvert K_f\rvert\, \ell / t_r$; it vanishes in a homogeneous dielectric
  (stripline), where $L_m/L_s = C_m/C_s$.

Each victim lists its aggressors by their crosstalk; above `crosstalk.threshold` a finding
is made. A pair's partner is not an aggressor. The waveforms with the actual ends follow
from a channel of the nets with steps (`docs/solvers/networks.md`).

Pairs are given, or found by name: two nets of the analysed set that differ only in a last
`P`/`N`, `p`/`n` or `+`/`-` (before trailing digits).

## Project file

```yaml
board: ../fab/
si:
  nets: ["*_P", "*_N", CLK*]       # patterns; default every net with tracks but the references
  reference_nets: [GND]            # default: nets filling a quarter of a layer
  pairs: [{p: USB_DP, n: USB_DM}]  # besides those found by name
  targets:
    - {nets: ["*USB*"], ohm: 90, differential: true}
    - {ohm: 100, differential: true, tolerance: 0.1}
    - {ohm: 50}
  materials:
    default: {epsilon_r: 4.3, loss_tangent: 0.02, frequency_hz: 1e9}
    DIELECTRIC_2: {epsilon_r: 4.6, loss_tangent: 0.018}
    mask: {epsilon_r: 3.8, loss_tangent: 0.025}
  copper: {roughness: {model: huray, radius_um: 0.5, ratio: 2.0}}
```

| key under `si` | default | meaning |
| --- | --- | --- |
| `nets` | all with tracks | net name patterns (`*`, `?`) |
| `reference_nets`, `plane_fraction` | auto, 0.25 | reference nets |
| `pairs`, `pair_by_name` | , true | differential pairs |
| `targets` | | `nets` (patterns), `ohm`, `differential` (false), `tolerance` (0.1) |
| `reference_hz` | $10^9$ | frequency of the reported quantities |
| `dispersion` | true | microstrip dispersion of tracks on outer layers |
| `materials` | stackup | per stackup layer name, `default` and `mask`: `epsilon_r`, `loss_tangent`, `frequency_hz` ($10^9$), `wideband` (true) |
| `copper` | | `conductivity_s_per_m` (stackup or $5.8 \times 10^7$), `roughness`: `model` (`none`, `hammerstad` with `rms_um`, `huray` with `radius_um` and `ratio`) |
| `etch_um`, `etch_face` | 0, | narrowing of each side of the far face: outer layers' outer face; inner layers' top face, or the face named per layer (`{In1.Cu: bottom}`, copper or stackup layer names) |
| `cut` | | `step_mm` (1), `lateral_mm` (25, the largest half width; a cut reaches twenty times the distance to its shields), `neighbour_mm` (1), `neighbour_heights` (3: quiet lines within this many shield distances), `large_copper_mm` (1), `target_error` (0.08), `reference_merge` (0.1: reference pieces of a layer join across gaps below this fraction of their distance to the track; 0: never) |
| `crosstalk` | | `rise_ps` (100), `threshold` (0.05) |
| `min_finding_mm` | 0.5 | shorter runs are not reported |
| `findings` | | `rise_ps` (100), `fraction` (0.2): runs whose delay is below this fraction of the rise time are electrically short, counted per net in the summary |

## Use

```
kairoq si project.yaml [<input>...] [--view]
```

prints the nets with layers, length, delay, impedance and differential impedance ranges
and loss, then the findings; `--view` shows the impedance on the tracks
(`docs/viewer.md`).

```python
r = kairoq.analyse_si(board, {"nets": ["*_P", "*_N"]})   # or a project file
r.nets, r.samples, r.findings
kairoq.serve(board, si=r)

# One cross-section: shapes in millimetres (rect_mm, polygon_mm, circle_mm, ring_mm).
x = kairoq.cross_section({
    "conductors": [{"name": "s", "rect_mm": [-0.1, 0.2, 0.1, 0.235]},
                   {"name": "g", "rect_mm": [-3, -0.035, 3, 0], "reference": True}],
    "dielectrics": [{"rect_mm": [-3, 0, 3, 0.2], "epsilon_r": 4.3, "loss_tangent": 0.02}],
    "frequencies_hz": [1e8, 1e9, 1e10]})
x["frequencies"][1]["lines"][0]["z0"]   # [real, imaginary] ohm
```

The cross-section returns, per frequency, the matrices `r`, `l`, `g`, `c` (row major, SI
units), the quantities of each line and, for two signal conductors, of the pair.

## Limits

- Quasi-TEM, with the dispersion of a bare microstrip on outer layers: the effect of
  solder mask on the dispersion, the different dispersion of a pair's even and odd modes,
  and dispersion of lines between dielectrics of different permittivity are not modelled.
- The surface impedance model is exact in its low- and high-frequency limits and
  approximate where the skin depth is comparable with the copper thickness.
- A cut is two-dimensional: bends, pads, vias and the ends of tracks are not lines and
  are excluded from the cuts; their effect is a discontinuity, treated by the network
  analysis and the openEMS export.
- Quiet neighbours beyond three times the distance to the shield are held at the reference
  potential, and are not aggressors; crosstalk between layers (broadside) is not scanned.
- The crosstalk estimates are those of weak coupling and matched ends; a channel gives the
  waveforms with the actual ends.
- The reference is the copper given or found as such; a return through a capacitor or a
  far via is not seen.
