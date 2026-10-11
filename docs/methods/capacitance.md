# Switch-node capacitance: method

This document derives the capacitance of a net (a switch node) to a metal plate on a face of
the board and to the other nets around it: the electrostatic Green's function of the
stratified stack, the method of moments that uses it, and the error. The user guide is
`docs/solvers/emi.md`, "Switch-node coupling"; the validation, `docs/validation/emi.md`,
"Switch-node capacitance with fringing".

## 1. Symbols

| symbol | meaning | unit |
| --- | --- | --- |
| $\varepsilon(z)$ | permittivity of the stack at height $z$ | F/m |
| $G(\rho; z, z')$ | potential at $(\rho, z)$ of a unit charge at $(0, z')$ | V/C |
| $\tilde\varphi(k; z, z')$ | its Hankel transform | V m/C |
| $\sigma_s$ | surface charge density on a conductor | C/m² |
| $P$ | potential coefficient matrix | V/C |
| $q_i$ | charge of panel $i$ | C |

## 2. Physical model

The plate is a perfectly conducting plane. Between it and the board lie the insulation
(gap and permittivity given) and then the board's layers, with air beyond the board. Each
copper layer is a conducting sheet at its mid-plane; half its thickness is joined to the slab
on each side, so that the conductors lie on interfaces of the stratified medium. The node is
at 1 V; every other net within a window around it, and the plate, at 0 V. The potential
satisfies

$$\nabla \cdot (\varepsilon(z) \nabla \phi) = 0$$

outside the conductors, $\phi$ is constant on each conductor, and the charge on a conductor
is the jump of $\varepsilon\, \partial \phi / \partial n$ across it. The capacitance of the node
to a net is minus the charge on that net; to the plate, the rest of the node's charge (the
plate's charge, since the total is zero).

## 3. The Green's function of the stack

### 3.1 Spectral form

The potential of a unit point charge at height $z'$ satisfies
$\nabla \cdot (\varepsilon \nabla G) = -\delta(\mathbf{r} - \mathbf{r}')$. It is radially
symmetric; with the Hankel transform

$$G(\rho; z, z') = \frac{1}{2\pi} \int_0^\infty \tilde\varphi(k; z, z')\, J_0(k\rho)\, k\, dk, \tag{1}$$

the transform satisfies in each slab $\tilde\varphi'' = k^2 \tilde\varphi$, with
$\tilde\varphi$ and $D = \varepsilon \tilde\varphi'$ continuous across interfaces, $D$ jumping by
$-1$ at the source, $\tilde\varphi = 0$ on the plate and decay in the half-space above.

Define at any height the ratio $Y_\downarrow = D / \tilde\varphi$ seen looking down and
$Y_\uparrow = -D / \tilde\varphi$ looking up. Across a slab of permittivity $\varepsilon$ and
thickness $t$, with the normalised ratio $y = Y / (\varepsilon k)$ at its far side, the
solutions $e^{\pm kz}$ transform it as a transmission line transforms an admittance:

$$y_\text{near} = \frac{\tanh kt + y_\text{far}}{1 + y_\text{far} \tanh kt}, \tag{2}$$

starting from $y = \infty$ on the plate and $y = 1$ in a half-space. At the source the jump
condition gives

$$\tilde\varphi(k; z', z') = \frac{1}{Y_\downarrow + Y_\uparrow},$$

and between levels the transform decays across each slab segment of thickness $t$ by
$1 / (\cosh kt + y \sinh kt)$, with $y$ the normalised ratio at the segment's far end looking
away from the source. In a homogeneous medium $Y_\downarrow = Y_\uparrow = \varepsilon k$, so
$\tilde\varphi = e^{-k|z - z'|} / (2\varepsilon k)$ and (1) gives $1/(4\pi\varepsilon R)$, as it
must.

### 3.2 Extraction of the singular terms

For large $k$ the transform behaves as $C\, e^{-k|z - z'|} / k$ with
$C = \lim_{k \to \infty} k\, \tilde\varphi\, e^{k|z - z'|}$ ($1/(\varepsilon_1 + \varepsilon_2)$ for a
source on an interface). Since $\int_0^\infty e^{-k|\Delta z|} J_0(k\rho)\, dk = 1/R$, this
part of (1) is the quasi-static direct term $C / (2\pi R)$, extracted in closed form.

When source and field lie in one plane, an interface at a small distance $d$ makes the
transform vary on the scale $1/d$, which numerical integration resolves only with effort.
Its effect is extracted too. Let the source plane separate a medium $\varepsilon_o$ from a slab
$\varepsilon_m$ of thickness $d$, beyond which lies $\varepsilon_f$ (or the plate). With
$K = (\varepsilon_m - \varepsilon_o)/(\varepsilon_m + \varepsilon_o)$,
$\Gamma = (\varepsilon_m - \varepsilon_f)/(\varepsilon_m + \varepsilon_f)$ ($-1$ for the plate) and
$u = e^{-2kd}$, the transform of this three-medium problem is

$$k\tilde\varphi = C\, \frac{1 + \Gamma u}{1 - \Gamma K u} = C \Big(1 + \sum_{n \geq 1} (1 + K)\, \Gamma^n K^{n-1} u^n\Big),$$

by expanding the geometric series. Each power $u^n$ transforms back to a point at the extra
depth $2nd$, so the nearest interface contributes the images

$$G_\text{images}(\rho) = \sum_{n \geq 1} \frac{C (1 + K) \Gamma^n K^{n-1}}{2\pi \sqrt{\rho^2 + (2nd)^2}}, \tag{3}$$

kept while their coefficients exceed $10^{-7}\,C$. The other interfaces are left to the
remainder.

### 3.3 The remainder

The remainder is (1) with the direct term and the images (3) subtracted from the integrand.
It varies on scales larger than the slabs and is integrated numerically: 8-point
Gauss–Legendre panels in $k$, no wider than $\pi/\rho$ (the oscillation of $J_0$) and growing
where the integrand no longer varies on their scale, until the integrand has stayed below
$10^{-12}\,C$ over four panels in a row. It is tabulated at 32 points per decade of $\rho$, as
$\text{rem}(\rho) \sqrt{\rho^2 + a^2}$ ($a$ the stack's height), which tends to constants at
both ends, and interpolated by cubic (Catmull–Rom) splines in $\log\rho$. Below the first
point it is interpolated linearly to its value at $\rho = 0$.

The split is exact: the direct term and images are subtracted from the transform before the
remainder is integrated, so they only move the singular parts to closed forms. The spectral
approach and the extraction of quasi-static terms are those of Michalski and Mosig (1997).

## 4. Method of moments

### 4.1 Discretisation

The copper of the node and of the nets within the window is meshed into triangles
(panels), finer at the node and growing away from it; other copper is simplified to the
local size, and a piece of copper gives up what larger pieces on its layer already cover.
Each panel carries a uniform charge density, $\sigma_s = q_i / A_i$ on panel $i$.

### 4.2 Galerkin system

Requiring the mean potential over each panel to equal its conductor's potential (Galerkin
testing with the same pulse functions; Harrington 1968) gives

$$\sum_j P_{ij} q_j = V_i, \qquad P_{ij} = \frac{1}{A_i A_j} \int_{T_i} \int_{T_j} G\, dA'\, dA. \tag{4}$$

The inner integral of the direct term and of each image over panel $j$ is the potential
integral $\int_{T_j} dA'/R$ in closed form (`docs/methods/emi.md`, equation 4); the remainder
by a 7-point rule over $j$. The outer mean over $i$ uses the 7-point rule, or the 28-point
rule where the panels are large against their vertical distance. Panels farther apart than
2.5 times the sum of their sizes use $G$ between centroids.

$P$ is symmetric (reciprocity) and positive definite: $\tfrac{1}{2} \mathbf{q}^\mathsf{T} P \mathbf{q}$
is the electrostatic energy of the charges $\mathbf{q}$, positive for any non-zero
distribution. It is factorised by Cholesky, and $\mathbf{q} = P^{-1} \mathbf{V}$ with
$V = 1$ on the node's panels and 0 elsewhere.

### 4.3 The computed capacitance is a lower bound

For one conductor at potential $V$ and charge $Q$, the true surface charge minimises the
energy $W = \tfrac{1}{2} \iint \sigma_s G \sigma_s$ among all distributions of total charge $Q$
(Thomson's theorem), and $W = Q^2 / (2C)$. The Galerkin solution minimises the same energy
over the piecewise-constant distributions, a subset, so its energy is larger and its
capacitance smaller: $C_h \leq C$. The same holds for the coefficients of a system of
conductors in the sense of the energy. Refining increases the computed charges towards the
true ones.

### 4.4 Extrapolation

Near the edges of thin conductors the charge density grows as $d^{-1/2}$, which uniform
panels resolve to first order in the element size $h$ there: $Q_h \approx Q - c\, h$. A second
solution on panels twice as large gives $Q_{2h} \approx Q - 2ch$, and the extrapolated charge
is

$$Q \approx 2 Q_h - Q_{2h} \tag{5}$$

(Richardson 1911). Both the extrapolated and the computed values are reported.

## 5. Parallel-plate estimate

Alongside, the area $A$ of the node facing the plate, not covered by copper of any net on
the layers between, gives the parallel-plate capacitance through the dielectrics in series:

$$C_\text{pp} = \varepsilon_0 \sum \frac{A}{\sum_k t_k / \varepsilon_{r,k}} .$$

## 6. Implementation

| step | where |
| --- | --- |
| stack, panels, window | `src/core/src/emi_capacitance.cpp`, `node_capacitance` |
| spectral transform (2), direct term, images (3), remainder | `src/core/src/layered_green.cpp`, `spectral_scaled`, `LayeredGreen` |
| Galerkin system (4), Cholesky | `emi_capacitance.cpp`, `charges` |
| potential integrals, quadrature | `src/core/src/peec.cpp` |
| extrapolation (5) | `emi_capacitance.cpp`, after `charges` |
| parallel-plate estimate | `src/core/src/emi.cpp` |

## References

- Harrington, R. F. (1968). *Field Computation by Moment Methods*. Macmillan; reprinted by
  IEEE Press (1993). doi:10.1109/9780470544631.
- Michalski, K. A., and Mosig, J. R. (1997). Multilayered media Green's functions in
  integral equation formulations. *IEEE Transactions on Antennas and Propagation* 45(3),
  508–519. doi:10.1109/8.558666.
- Richardson, L. F. (1911). The approximate arithmetical solution by finite differences of
  physical problems involving differential equations, with an application to the stresses
  in a masonry dam. *Philosophical Transactions of the Royal Society A* 210, 307–357.
  doi:10.1098/rsta.1911.0009.
- The potential integrals and quadrature rules: `docs/methods/emi.md`.
