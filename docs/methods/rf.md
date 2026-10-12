# RF structures: method

This document derives the native RF analyses: the de-embedding operations on
S-parameters (change of reference, connection and removal of a two-port, lines, the
2x-thru) and the transmission-line model of a rectangular microstrip patch. The full-wave
export to openEMS is in `docs/methods/external.md`. The user guide is `docs/solvers/rf.md`;
the validation, `docs/validation/rf.md`.

## 1. Symbols

| symbol | meaning | unit |
| --- | --- | --- |
| $\mathbf{S}$ | scattering matrix of an $n$-port with the real reference $Z_0$ at every port | |
| $X_{ij}$ | scattering parameters of a two-port fixture at the same reference | |
| $\mathbf{T} = \begin{pmatrix} A & B \\ C & D \end{pmatrix}$ | transfer (ABCD) matrix of a two-port | |
| $\gamma = \alpha + j\beta$ | propagation constant | 1/m |
| $W, L, h$ | patch width (radiating edges), length, substrate height | m |
| $k_0 = 2\pi/\lambda_0$ | free-space wavenumber | 1/m |

## 2. De-embedding

### 2.1 Change of reference impedance

With the power waves of Kurokawa (1965) and a real reference $Z_0$ at every port,
$\mathbf{Z} = Z_0(\mathbf{I} + \mathbf{S})(\mathbf{I} - \mathbf{S})^{-1}$. For the new reference $Z_1$,
$\mathbf{S}' = (\mathbf{Z} - Z_1\mathbf{I})(\mathbf{Z} + Z_1\mathbf{I})^{-1}$; substituting and
multiplying through by $(\mathbf{I} - \mathbf{S})$,

$$\mathbf{S}' = \big[(Z_0 - Z_1)\mathbf{I} + (Z_0 + Z_1)\mathbf{S}\big]\big[(Z_0 + Z_1)\mathbf{I} + (Z_0 - Z_1)\mathbf{S}\big]^{-1} = (\mathbf{I} - \rho\mathbf{S})^{-1}(\mathbf{S} - \rho\mathbf{I}), \qquad \rho = \frac{Z_1 - Z_0}{Z_1 + Z_0}, \tag{1}$$

the two factors commuting as functions of $\mathbf{S}$.

### 2.2 A two-port in front of a port

Let a fixture $X$ have its port 2 joined to port $p$ of the network and its port 1 become
the new port $p$, with incident wave $a'$ and reflected $b'$. At the junction the wave
leaving the network enters the fixture and the reverse: $a_{X2} = b_p$, $a_p = b_{X2}$. From
$b_{X2} = X_{21}a' + X_{22}b_p$ and $b_p = S_{pp}a_p + \sum_{j \neq p} S_{pj}a_j$,

$$a_p = \frac{X_{21}a' + X_{22}\sum_{j \neq p} S_{pj}a_j}{\Delta}, \qquad \Delta = 1 - X_{22}S_{pp}.$$

Substituting into $b' = X_{11}a' + X_{12}b_p$ and $b_i = S_{ip}a_p + \sum_{j \neq p} S_{ij}a_j$,

$$S'_{pp} = X_{11} + \frac{X_{12}X_{21}S_{pp}}{\Delta}, \quad S'_{pj} = \frac{X_{12}S_{pj}}{\Delta}, \quad S'_{ip} = \frac{S_{ip}X_{21}}{\Delta}, \quad S'_{ij} = S_{ij} + \frac{S_{ip}X_{22}S_{pj}}{\Delta}, \tag{2}$$

for $i, j \neq p$, using $X_{12}(S_{pp}X_{22} + \Delta) = X_{12}$ for $S'_{pj}$. Fixture and network
must share the reference; a thru is renormalised to the network's by (1) first.

### 2.3 Removing a fixture

In transfer matrices, a cascade multiplies: a fixture of ABCD matrix $\mathbf{T}_f$ (its port 1
outward) in front of the device $\mathbf{T}_d$ is measured as $\mathbf{T}_f\mathbf{T}_d$. Placing the
two-port of $\mathbf{T}_f^{-1}$ in front of the measured port by (2) gives
$\mathbf{T}_f^{-1}\mathbf{T}_f\mathbf{T}_d = \mathbf{T}_d$: the fixture is removed. The conversions
between S and ABCD at the reference $Z_0$ are the standard ones (Pozar 2012, §4.4):

$$A = \frac{(1 + S_{11})(1 - S_{22}) + S_{12}S_{21}}{2S_{21}}, \quad B = Z_0\frac{(1 + S_{11})(1 + S_{22}) - S_{12}S_{21}}{2S_{21}},$$
$$C = \frac{(1 - S_{11})(1 - S_{22}) - S_{12}S_{21}}{2S_{21}Z_0}, \quad D = \frac{(1 - S_{11})(1 + S_{22}) + S_{12}S_{21}}{2S_{21}},$$

and their inverse.

### 2.4 Lines

A uniform line of impedance $Z_c$, effective permittivity $\varepsilon_\text{eff}$ and length
$\ell$ has $\mathbf{T} = \exp(\ell\mathbf{M})$ in the closed form of `docs/methods/networks.md`,
section 3.2, with

$$\gamma = \alpha_1\sqrt{f / 1\ \text{GHz}} + j\,\frac{2\pi f \sqrt{\varepsilon_\text{eff}}}{c_0},$$

$\alpha_1$ the attenuation at 1 GHz in Np/m (the given dB/m divided by $20/\ln 10$), scaled
as the skin effect. Since $\exp(\ell\mathbf{M})^{-1} = \exp(-\ell\mathbf{M})$, a negative length
removes a line exactly, loss included.

### 2.5 The 2x-thru

Two copies of a fixture joined back to back measure $\mathbf{M} = \mathbf{T}_f\mathbf{T}_f^\text{rev}$,
with $\mathbf{T}^\text{rev} = \begin{pmatrix} D & B \\ C & A \end{pmatrix}$ the fixture seen from
its other end (reciprocal, $AD - BC = 1$). For a symmetric fixture, $A = D$ and
$\mathbf{T}_f^\text{rev} = \mathbf{T}_f$, so $\mathbf{M} = \mathbf{T}_f^2$ and the fixture is a square root
of $\mathbf{M}$.

**Square roots of a 2×2 matrix.** By the Cayley–Hamilton theorem
$\mathbf{M}^2 = (\operatorname{tr}\mathbf{M})\,\mathbf{M} - (\det\mathbf{M})\,\mathbf{I}$. For
$\mathbf{R} = (\mathbf{M} + s\mathbf{I})/t$ with $s^2 = \det\mathbf{M}$,

$$\mathbf{R}^2 = \frac{\mathbf{M}^2 + 2s\mathbf{M} + s^2\mathbf{I}}{t^2} = \frac{(\operatorname{tr}\mathbf{M} + 2s)\,\mathbf{M}}{t^2} = \mathbf{M} \quad \text{if } t^2 = \operatorname{tr}\mathbf{M} + 2s. \tag{3}$$

The signs of $s$ and $t$ give four roots (where $t \neq 0$). The physical one is continuous in
frequency and is the identity at zero frequency: at each frequency the root nearest the
previous one is taken, with $B/Z_0$ and $CZ_0$ making the entries dimensionless. For a
lossless line of electrical length $\theta$ per half, $\det\mathbf{M} = 1$,
$\operatorname{tr}\mathbf{M} = 2\cos 2\theta$, and (3) gives $\pm\mathbf{T}(\theta)$ for $s = 1$ and
$\pm j\,\mathbf{T}(\theta - \pi/2)$ for $s = -1$, $\mathbf{T}(\theta)$ the line's matrix. In the
distance used (the sum of the moduli of the dimensionless entries) each spurious root lies
at least 2 from the physical one at every $\theta$, so following the branch from zero frequency picks $\mathbf{T}(\theta)$ through any
number of wavelengths, provided a step of the sweep moves $\theta$ by well under a radian. The thru's
data are interpolated to the network's frequencies linearly in magnitude and in unwrapped
phase. The method needs a symmetric fixture; IEEE 370 (2020) describes methods for
others.

## 3. Rectangular microstrip patch

### 3.1 Transmission-line model

The patch is a section of wide microstrip, open at both ends, resonating where its
effective length is half a wavelength. Its effective permittivity and the extension of
each open end by the fringing field are the closed forms (Balanis 2016, ch. 14, after
Hammerstad 1975 and Hammerstad and Jensen 1980):

$$\varepsilon_\text{eff} = \frac{\varepsilon_r + 1}{2} + \frac{\varepsilon_r - 1}{2}\Big(1 + 12\frac{h}{W}\Big)^{-1/2}, \qquad \frac{\Delta L}{h} = 0.412\,\frac{(\varepsilon_\text{eff} + 0.3)(W/h + 0.264)}{(\varepsilon_\text{eff} - 0.258)(W/h + 0.8)},$$

$$f_r = \frac{c_0}{2(L + 2\Delta L)\sqrt{\varepsilon_\text{eff}}}.$$

For a design frequency the width $W = \frac{c_0}{2f_r}\sqrt{2/(\varepsilon_r + 1)}$ is the usual
choice for efficient radiation, and $L$ follows from $f_r$ with $\Delta L$.

### 3.2 Radiating edges

Each radiating edge is a slot of length $W$ whose aperture field is that of the open end.
Its radiated power gives the conductance

$$G_1 = \frac{1}{120\pi^2}\int_0^\pi \Big[\frac{\sin\big(\tfrac{k_0W}{2}\cos\theta\big)}{\cos\theta}\Big]^2\sin^3\theta\, d\theta,$$

and the two edges, $L$ apart and excited in antiphase by the fundamental mode, couple
through

$$G_{12} = \frac{1}{120\pi^2}\int_0^\pi \Big[\frac{\sin\big(\tfrac{k_0W}{2}\cos\theta\big)}{\cos\theta}\Big]^2 J_0(k_0L\sin\theta)\sin^3\theta\, d\theta.$$

At resonance the susceptances cancel and the input resistance at an edge is
$R_\text{in} = 1/\big(2(G_1 + G_{12})\big)$ (Balanis 2016, ch. 14). The voltage of the
fundamental mode varies as $\cos(\pi y/L)$ along the patch, and the current as its
derivative, so a feed inset by $y_0$ sees $R_\text{in}\cos^2(\pi y_0/L)$, and the inset
matching a feed resistance $R_f \leq R_\text{in}$ is $y_0 = \frac{L}{\pi}\arccos\sqrt{R_f/R_\text{in}}$.
The integrals are evaluated by the composite Simpson rule with 4000 intervals, with the
limit $\tfrac{k_0W}{2}$ of the integrand's bracket at $\cos\theta = 0$.

### 3.3 Directivity

The two slots radiate over an infinite ground plane in $z = 0$ as magnetic currents along
$y$ at $x = \pm s/2$, $s = L + 2\Delta L$. Their radiation intensity, normalised to 1 at
broadside, is

$$U(\theta, \phi) = \operatorname{sinc}^2\Big(\frac{k_0W}{2}\sin\theta\sin\phi\Big)\big(1 - \sin^2\theta\sin^2\phi\big)\cos^2\Big(\frac{k_0s}{2}\sin\theta\cos\phi\Big):$$

the element factor of a uniform line source of length $W$, the pattern of a magnetic
dipole along $y$ ($|\hat{\mathbf{y}} \times \hat{\mathbf{r}}|^2$), and the array factor of the
two slots in phase. The directivity is $D = 4\pi / \int_\text{upper half-space} U\, d\Omega$,
by Simpson's rule with 180 intervals in $\theta$ over $[0, \pi/2]$ and 360 in $\phi$.

### 3.4 Bandwidth

The bandwidth for a VSWR below 2 is the closed form of Jackson and Alexopoulos (1991) for
thin substrates,

$$\text{BW} = 3.771\,\frac{\varepsilon_r - 1}{\varepsilon_r^2}\,\frac{W}{L}\,\frac{h}{\lambda_0}.$$

## 4. Error

| part | error |
| --- | --- |
| de-embedding (1)–(3) | exact algebra; real, equal references at all ports |
| 2x-thru | exact for a symmetric fixture; an asymmetric one is split wrongly; branch tracking needs sweep steps of well under a radian of the half's electrical length |
| thru interpolation | linear in magnitude and phase between the thru's samples |
| patch | the transmission-line model: thin substrates ($h$ below about $0.05$ wavelengths in the dielectric, reported), fundamental mode, no feed reactance, infinite ground plane; against openEMS, `docs/validation/rf.md` |

## 5. Implementation

| step | where |
| --- | --- |
| reference change (1) | `src/core/src/deembed.cpp`, `renormalise` |
| connection (2) | `deembed.cpp`, `connect_fixture` |
| S and ABCD | `deembed.cpp`, `abcd_from_s`; `src/core/src/network.cpp`, `s_from_abcd` |
| lines 2.4 | `deembed.cpp`, `add_line`; `network.cpp`, `line_abcd` |
| 2x-thru (3) | `deembed.cpp`, `half_thru_abcd`, `deembed_thru` |
| interpolation | `network.cpp`, `Touchstone::at` |
| patch 3.1–3.4 | `src/core/src/rf.cpp`, `patch_antenna`, `two_slot_directivity` |

## References

- Balanis, C. A. (2016). *Antenna Theory: Analysis and Design*, 4th ed. Wiley.
- Hammerstad, E. O. (1975). Equations for microstrip circuit design. *Proceedings of the
  5th European Microwave Conference*, 268–272. doi:10.1109/EUMA.1975.332206.
- Hammerstad, E., and Jensen, Ø. (1980). Accurate models for microstrip computer-aided
  design. *IEEE MTT-S International Microwave Symposium Digest*, 407–409.
  doi:10.1109/MWSYM.1980.1124303.
- IEEE (2020). *IEEE Standard for Electrical Characterization of Printed Circuit Board and
  Related Interconnects at Frequencies up to 50 GHz*, IEEE Std 370-2020.
  doi:10.1109/IEEESTD.2021.9316329.
- Jackson, D. R., and Alexopoulos, N. G. (1991). Simple approximate formulas for input
  resistance, bandwidth, and efficiency of a resonant rectangular patch. *IEEE
  Transactions on Antennas and Propagation* 39(3), 407–410. doi:10.1109/8.76341.
- Kurokawa, K. (1965). Power waves and the scattering matrix. *IEEE Transactions on
  Microwave Theory and Techniques* 13(2), 194–202. doi:10.1109/TMTT.1965.1125964.
- Pozar, D. M. (2012). *Microwave Engineering*, 4th ed. Wiley.
