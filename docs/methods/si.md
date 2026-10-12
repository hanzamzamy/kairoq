# Transmission lines: method

This document derives the per-unit-length parameters of a transmission-line cross-section:
capacitance and external inductance from two electrostatic problems, the dielectric model,
microstrip dispersion, the conductor losses by Wheeler's rule in the form of a shape
derivative, and the line quantities. The user guide is `docs/solvers/si.md`; the
validation, `docs/validation/si.md` and `docs/validation/networks.md`.

## 1. Symbols

| symbol | meaning | unit |
| --- | --- | --- |
| $\varphi_i$ | potential with conductor $i$ at 1 V, the others at 0 V | V |
| $C$, $C^{(0)}$ | Maxwell capacitance matrix with the dielectrics, in vacuum | F/m |
| $L$ | inductance matrix | H/m |
| $R$, $G$ | resistance and conductance matrices | Ω/m, S/m |
| $\varepsilon_r = \varepsilon' - j\varepsilon''$ | complex relative permittivity | – |
| $\delta$ | skin depth $1/\sqrt{\pi f \mu_0 \sigma}$ | m |
| $\mathbf{V}$ | velocity field of a shape perturbation | – |

## 2. Quasi-TEM model

A line uniform along $z$ whose cross-section is small against the wavelength carries,
to first order, transverse fields that satisfy the two-dimensional Laplace equation (Paul
2008, ch. 3). The electric field is that of the electrostatic potentials $\varphi_i$:

$$\nabla \cdot (\varepsilon_r \nabla \varphi_i) = 0 \text{ outside the conductors}, \qquad \varphi_i = \delta_{ik} \text{ on conductor } k, \qquad \varphi_i = 0 \text{ on the reference and the box}. \tag{1}$$

The capacitance matrix is the energy of these potentials,

$$C_{ij} = \varepsilon_0 \int \varepsilon_r\, \nabla \varphi_i \cdot \nabla \varphi_j\, dA, \tag{2}$$

which equals the charge form $C_{ij} = q_j$ on conductor $i$ by Green's identity. The magnetic
field of TEM currents is that of the vacuum electrostatic problem rotated by 90°, so in
vacuum $L C^{(0)} = \mu_0 \varepsilon_0 I$; the magnetic field does not see the dielectrics, so
the external inductance is

$$L = \mu_0 \varepsilon_0 \big(C^{(0)}\big)^{-1}. \tag{3}$$

The domain is bounded by a grounded box well beyond the geometry. For a via in the
cross-section solver's axisymmetric mode the cut's horizontal axis is the radius of a body
of revolution: each element is weighted by $2\pi r$ at its centroid and the axis is left
free (natural condition).

## 3. Finite elements and error

(1) is discretised with linear triangles as in the DC solver (`docs/methods/dc.md`, §4): the
element matrix $\varepsilon_r (b_i b_j + c_i c_j)/(4A)$, the conductor nodes fixed, one
right-hand side per signal conductor, solved by conjugate gradients with algebraic
multigrid to $10^{-11}$. (2) is evaluated exactly per triangle.

**The computed capacitance is an upper bound.** By the Dirichlet principle the exact
$\varphi_i$ minimises the energy (2) among potentials with the conductors' values; the
discrete one minimises it over a subspace. For any vector of conductor voltages
$\mathbf{v}$, $\mathbf{v}^\mathsf{T} C_h \mathbf{v} \geq \mathbf{v}^\mathsf{T} C \mathbf{v}$, and the
difference is the squared energy error (as in `docs/methods/dc.md`, §6.1). In particular the
self capacitances converge from above, the external inductances from below.

The error is estimated by gradient recovery per material (`docs/methods/thermal.md`, §7),
summed over the solutions, each relative to its own energy; refinement continues until the
estimate meets the target. The relative error of $C$ is about the square of the estimate,
and that of an impedance $\sqrt{L/C}$ about half of it.

## 4. Dielectric model

By default each dielectric follows the Djordjevic–Sarkar model (Djordjevic et al. 2001), a
continuous distribution of relaxation frequencies between $\omega_1 = 10^4$ and
$\omega_2 = 10^{12}$ rad/s:

$$\varepsilon_r(\omega) = \varepsilon_\infty + \frac{\Delta}{\ln(\omega_2/\omega_1)} \ln \frac{\omega_2 + j\omega}{\omega_1 + j\omega}. \tag{4}$$

It is analytic in the lower half of the complex $\omega$ plane, so its real and imaginary
parts satisfy the Kramers–Kronig relations: it is causal. Its two parameters follow from the
permittivity $\varepsilon_r$ and loss tangent $\tan\delta$ given at a frequency $\omega_0$.
With $S = \ln(\omega_2/\omega_1)$, the imaginary part of (4) is
$-\frac{\Delta}{S}\big(\arctan\frac{\omega}{\omega_1} - \arctan\frac{\omega}{\omega_2}\big)$, so
$\varepsilon''(\omega_0) = \varepsilon_r \tan\delta$ gives

$$\Delta = \frac{\varepsilon_r \tan\delta\, S}{\arctan(\omega_0/\omega_1) - \arctan(\omega_0/\omega_2)},$$

and the real part, $\varepsilon_\infty + \frac{\Delta}{2S} \ln \frac{\omega_2^2 + \omega^2}{\omega_1^2 + \omega^2}$,
equal to $\varepsilon_r$ at $\omega_0$, gives
$\varepsilon_\infty = \varepsilon_r - \frac{\Delta}{2S} \ln \frac{\omega_2^2 + \omega_0^2}{\omega_1^2 + \omega_0^2}$.
Between $\omega_1$ and $\omega_2$ the loss tangent is nearly constant and the permittivity
falls slowly, about linearly in $\log\omega$. A constant model is available.

## 5. Frequency dependence

The potentials are solved once, at the first frequency. The energy of each dielectric region
$r$ is kept, $E^{(r)}_{ij} = \varepsilon_0 \int_r \nabla\varphi_i \cdot \nabla\varphi_j\, dA$, and at
any frequency

$$C_{ij}(f) = \sum_r \varepsilon'_r(f)\, E^{(r)}_{ij}, \qquad G_{ij}(f) = \omega \sum_r \varepsilon''_r(f)\, E^{(r)}_{ij}. \tag{5}$$

For $C$ this is the energy of the old potentials in the new medium; by the Dirichlet
principle it bounds the new capacitance from above, and the difference is of second order in
the change of the field, hence in the change of the permittivities' ratios. For $G$, (5) is
the first-order perturbation of the capacitance by the imaginary part of the permittivity.
The option that solves at two frequencies per decade and interpolates in $\log f$ agrees with
(5) within $2 \times 10^{-4}$ in capacitance on a masked microstrip from 1 MHz to 20 GHz
(`docs/validation/networks.md`).

## 6. Microstrip dispersion

A microstrip's field moves into the substrate as the frequency rises; the quasi-TEM model
has no such effect. For a track on an outer layer with one reference layer, the effective
permittivity follows Kirschning and Jansen (1982),

$$\varepsilon_\text{eff}(f) = \varepsilon_r - \frac{\varepsilon_r - \varepsilon_\text{eff}(0)}{1 + P(f)}, \qquad P = P_1 P_2 \big((0.1844 + P_3 P_4) f_n\big)^{1.5763},$$

with $f_n = f h$ in GHz mm, $u = w/h$ and

$$P_1 = 0.27488 + \Big(0.6315 + \frac{0.525}{(1 + 0.0157 f_n)^{20}}\Big) u - 0.065683\, e^{-8.7513 u}, \qquad P_2 = 0.33622 \big(1 - e^{-0.03442 \varepsilon_r}\big),$$

$$P_3 = 0.0363\, e^{-4.6 u} \big(1 - e^{-(f_n/38.7)^{4.97}}\big), \qquad P_4 = 1 + 2.751 \big(1 - e^{-(\varepsilon_r/15.916)^8}\big),$$

$\varepsilon_\text{eff}(0)$ being the cross-section's own value and $\varepsilon_r$ the series
permittivity of the layers under the track. The impedance ratio $Z(f)/Z(0)$ follows Jansen
and Kirschning (1983) (coefficients $R_1 \dots R_{17}$, implemented as published).

The capacitance matrix moves accordingly from $C$ towards that of a uniform fill,
$\varepsilon_r C^{(0)}$:

$$C' = \varepsilon_r C^{(0)} - \frac{\varepsilon_r C^{(0)} - C}{1 + P}.$$

For the track itself, $C'_{ss}/C^{(0)}_{ss} = \varepsilon_r - (\varepsilon_r - \varepsilon_\text{eff}(0))/(1+P) = \varepsilon_\text{eff}(f)$,
the formula's value; for a line in one dielectric, $C = \varepsilon_r C^{(0)}$ and $C' = C$.
Then $C'$ is divided and $L$ multiplied by $s = (Z(f)/Z(0)) \sqrt{\varepsilon_\text{eff}(f)/\varepsilon_\text{eff}(0)}$:
the product $LC'$, hence the velocity, is unchanged, and the impedance
$\sqrt{L s / (C'/s)} = s\, Z(0) \sqrt{\varepsilon_\text{eff}(0)/\varepsilon_\text{eff}(f)} = Z(f)$
is the formula's. The formulas fit full-wave solutions for $0.1 \leq w/h \leq 100$ and
$1 \leq \varepsilon_r \leq 20$; above $f_n = 20$ the correction is held.

## 7. Conductor losses

### 7.1 Surface currents and Wheeler's rule

In the skin-effect limit the current flows on the conductor surfaces, and its distribution
is that of the vacuum solution's surface charge (the tangential magnetic field of the TEM
mode equals the normal electric field divided by the impedance of free space). With the
surface impedance $Z_s$ the dissipated power per length gives the series impedance

$$Z_{ij} = Z_s \sum_k \oint_k J_i J_j\, dl, \tag{6}$$

with $J_i$ the surface current of unit current on line $i$. Wheeler (1942) observed that the
surface integral is the change of inductance when the conductor's surface recedes:

$$\oint_k J_i J_j\, dl = \frac{1}{\mu_0} \frac{\partial L_{ij}}{\partial n_k}.$$

### 7.2 The shape derivative in domain form

The surface charge, and with it the integrand of (6), is singular at conductor corners, so
the surface integral converges poorly. The derivative is evaluated instead as the shape
derivative of the energy (2) in vacuum, in its domain form (Delfour and Zolésio 2011): when
the surface of conductor $k$ moves with a velocity field $\mathbf{V}$ (outward unit normal
on the surface, extended smoothly into the field and vanishing before any other conductor),

$$\frac{\partial C^{(0)}_{ab}}{\partial n_k} = \varepsilon_0 \int \Big( \nabla\varphi_a \cdot \nabla\varphi_b \, \nabla \cdot \mathbf{V} - \nabla\varphi_a \cdot (D\mathbf{V} + D\mathbf{V}^\mathsf{T}) \nabla\varphi_b \Big)\, dA, \tag{7}$$

an area integral near the conductor that converges as the energy does. $\mathbf{V}$ is the
direction away from the nearest surface point, of length 1 within $a = 0.45\, d_\text{min}$
and falling linearly to 0 at $2a$, with $d_\text{min}$ the distance to the nearest other
conductor or the box. Differentiating (3),
$\partial L / \partial n = -\mu_0\varepsilon_0\, (C^{(0)})^{-1} \frac{\partial C^{(0)}}{\partial n} (C^{(0)})^{-1}$,
which gives the matrices $M_k$ of (6) per conductor.

### 7.3 Between DC and the skin effect

Each conductor's surface impedance is taken as that of a slab of effective thickness
$t_\text{eff}$ carrying current on one face (`docs/methods/emi.md`, §5):

$$Z_s(f) = K(f)\, \frac{1 + j}{\sigma\delta} \coth \frac{(1 + j)\, t_\text{eff}}{\delta} \ \xrightarrow{\ f \to 0\ }\ \frac{1}{\sigma t_\text{eff}}.$$

$t_\text{eff}$ is chosen so that the DC limit is the conductor's resistance: with $A_k$ its
cross-section, $\frac{1}{\sigma t_\text{eff}} (M_k)_{kk} = \frac{1}{\sigma A_k}$, so
$t_\text{eff} = A_k (M_k)_{kk}$ (for the reference, the mean over the lines). Both limits are
exact; the transition is approximate. The imaginary part of (6), divided by $\omega$, is the
internal inductance, added to (3).

The roughness factor $K(f)$ is 1, or Hammerstad's
$1 + \frac{2}{\pi} \arctan\big(1.4 (\Delta/\delta)^2\big)$ for an RMS height $\Delta$
(Hammerstad and Jensen 1980), or Huray's $1 + r / (1 + \delta/a + \delta^2/(2a^2))$ for spheres
of radius $a$ and surface ratio $r$ (Huray 2009).

## 8. Line quantities

For line $i$ with the others at rest,
$Z_0 = \sqrt{(R_{ii} + j\omega L_{ii})/(G_{ii} + j\omega C_{ii})}$ and
$\gamma = \alpha + j\beta = \sqrt{(R_{ii} + j\omega L_{ii})(G_{ii} + j\omega C_{ii})}$; the delay
is $\beta/\omega$, the effective permittivity $(\beta c_0/\omega)^2$, the loss $8.686\,\alpha$
dB/m, of which $R/(2Z_0)$ in the conductors and $G Z_0/2$ in the dielectric to first order
(Pozar 2012, §2.7). For a symmetric pair, with $C_m < 0$ the mutual Maxwell capacitance,
$Z_\text{odd} = \sqrt{(L_s - L_m)/(C_s - C_m)}$, $Z_\text{even} = \sqrt{(L_s + L_m)/(C_s + C_m)}$,
$Z_\text{diff} = 2 Z_\text{odd}$ and $Z_\text{common} = Z_\text{even}/2$.

## 9. Board analysis

Tracks are cut across at intervals; a cut takes the copper of every layer along a line
across the track, out to twenty times the distance to the nearest plane (at most 25 mm)
beyond the track and its partner, within the board's outline. Cuts with the same geometry
are solved once. Reference pieces of a layer are joined across gaps smaller than a tenth of
their distance to the track, which changes the fields by about the square of that ratio. Each net's
impedance profile, delay and loss follow from its cuts weighted by their length.

## 10. Implementation

| step | where |
| --- | --- |
| domain, box, mesh, axisymmetric weight | `src/core/src/cross_section.cpp`, `Solver::build_domain`, `Solver::mesh` |
| potentials (1), capacitance (2), estimate | `cross_section.cpp`, `Solver::solve`, `Solver::capacitance`, `Solver::estimate` |
| shape derivative (7), effective thickness | `cross_section.cpp`, `Solver::surface_products`, `Solver::run` |
| frequency model (5), losses, dispersion | `cross_section.cpp`, `CsModel::at`, `kirschning_jansen` |
| dielectric (4), roughness | `cross_section.cpp`, `DielectricModel::at`, `ConductorModel::roughness_factor` |
| cuts, nets, pairs | `src/core/src/si.cpp` |

## References

- Delfour, M. C., and Zolésio, J.-P. (2011). *Shapes and Geometries: Metrics, Analysis,
  Differential Calculus, and Optimization*, 2nd ed. SIAM. doi:10.1137/1.9780898719826.
- Djordjevic, A. R., Biljić, R. M., Likar-Smiljanić, V. D., and Sarkar, T. K. (2001).
  Wideband frequency-domain characterization of FR-4 and time-domain causality. *IEEE
  Transactions on Electromagnetic Compatibility* 43(4), 662–667. doi:10.1109/15.974647.
- Hammerstad, E., and Jensen, Ø. (1980). Accurate models for microstrip computer-aided
  design. *IEEE MTT-S International Microwave Symposium Digest*, 407–409.
  doi:10.1109/MWSYM.1980.1124303.
- Huray, P. G. (2009). *The Foundations of Signal Integrity*. Wiley–IEEE Press.
  doi:10.1002/9780470543481.
- Jansen, R. H., and Kirschning, M. (1983). Arguments and an accurate model for the
  power-current formulation of microstrip characteristic impedance. *Archiv für Elektronik
  und Übertragungstechnik (AEÜ)* 37(3/4), 108–112.
- Kirschning, M., and Jansen, R. H. (1982). Accurate model for effective dielectric constant
  of microstrip with validity up to millimetre-wave frequencies. *Electronics Letters* 18(6),
  272–273. doi:10.1049/el:19820186.
- Paul, C. R. (2008). *Analysis of Multiconductor Transmission Lines*, 2nd ed. Wiley.
  doi:10.1002/9780470131541.
- Pozar, D. M. (2012). *Microwave Engineering*, 4th ed. Wiley.
- Wheeler, H. A. (1942). Formulas for the skin effect. *Proceedings of the IRE* 30(9),
  412–424. doi:10.1109/JRPROC.1942.232015.
