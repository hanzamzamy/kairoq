# EMI loops: method

This document derives the loop solver of the EMI estimates: the inductance and the current
distribution of a switching loop by partial inductances, its resistance and inductance over
frequency, its magnetic moment and the field it radiates. The switch-node capacitance is
derived in `docs/methods/capacitance.md`. The user guide is `docs/solvers/emi.md`; the
validation, `docs/validation/emi.md`.

## 1. Symbols

| symbol | meaning | unit |
| --- | --- | --- |
| $\mathbf{J}$, $\mathbf{K}$ | volume and sheet current density | A/m², A/m |
| $\mathbf{A}$ | magnetic vector potential | V s/m |
| $\mu_0$ | permeability of free space, $4\pi \times 10^{-7}$ | H/m |
| $\mathbf{f}_n$ | basis function of branch $n$ | 1/m |
| $x_n$ | current of branch $n$ | A |
| $L$ | partial inductance matrix of the branches | H |
| $D$ | incidence matrix of the current law | – |
| $\mathbf{m}$ | magnetic moment | A m² |

## 2. Physical model

### 2.1 The inductive limit

At the frequencies of switching edges, tens of megahertz and more, the skin depth in copper
is a fraction of the layer thickness and the sheet resistance is small against $\omega L$.
In the limit $R / \omega L \to 0$ the current on the conductors takes the distribution that
stores the least magnetic energy among all distributions that satisfy charge conservation
and carry the loop current:

$$W = \frac{1}{2} \int \mathbf{A} \cdot \mathbf{J}\, dV = \frac{\mu_0}{8 \pi} \iint \frac{\mathbf{J}(\mathbf{r}) \cdot \mathbf{J}(\mathbf{r}')}{|\mathbf{r} - \mathbf{r}'|}\, dV\, dV', \qquad \nabla \cdot \mathbf{J} = 0, \tag{1}$$

with the loop's inductance $L_\text{loop} = 2W / I^2$ (Ruehli et al. 2017). The
retardation $e^{-jkR}$ is neglected inside the loop: the loop is small against the
wavelength (quasi-static).

### 2.2 Conductors

- **Copper** is a set of current sheets at the layers' mid-planes; their sheet current
  replaces $\mathbf{J}\, dV$ by $\mathbf{K}\, dA$ in (1). Towards an edge of a zero-thickness
  sheet the current grows as $d^{-1/2}$; real copper spreads it over about $t/2$.
- **Elements** (the parts of the loop) are straight round conductors between the centres of
  their pads, at their radius above the copper; their current enters and leaves the copper
  spread over each pad. An equivalent series inductance inside the part adds in series.
- **Plated holes** are vertical tubes between the layers they join.
- **Other copper** in a window around the loop carries the currents the loop induces in it.

## 3. Partial inductances

### 3.1 Branches

The sheet current on each mesh is expanded in Rao–Wilton–Glisson (RWG) functions (Rao,
Wilton and Glisson 1982). Interior edge $n$ is shared by triangles $T^+$ and $T^-$ with
opposite vertices $\mathbf{p}^\pm$:

$$\mathbf{f}_n(\mathbf{r}) = \pm \frac{\mathbf{r} - \mathbf{p}^\pm}{2 A^\pm}, \quad \mathbf{r} \in T^\pm. \tag{2}$$

On $T^\pm$ the normal component of $\mathbf{f}_n$ across the edge is $1/\ell_n$ ($\ell_n$ its
length) and zero across the other edges, so $x_n$ is the total current across edge $n$ and
the normal current is continuous: no charge accumulates on the edges. The divergence is
constant on each triangle, $\nabla \cdot \mathbf{f}_n = \pm 1/A^\pm$, so the current law is
one equation per triangle: the currents of its edges and of any conductor ending on it sum
to zero. Each element and hole is a branch of its own, $\mathbf{f} = \hat{\mathbf{l}}\, \delta$
along its axis, whose ends take their current from sets of triangles with weights summing
to one (the pad's or the hole's area).

With $\mathbf{K} = \sum_n x_n \mathbf{f}_n$, (1) becomes $W = \tfrac{1}{2} \mathbf{x}^\mathsf{T} L\, \mathbf{x}$ with

$$L_{mn} = \frac{\mu_0}{4\pi} \iint \frac{\mathbf{f}_m(\mathbf{r}) \cdot \mathbf{f}_n(\mathbf{r}')}{|\mathbf{r} - \mathbf{r}'|}\, dA\, dA', \tag{3}$$

the partial inductances of Ruehli (1972). $L$ is symmetric and positive definite: (3) is the
energy of a current distribution, positive for any non-zero one.

### 3.2 Sheets with sheets

For a pair of triangles $T$ (vertices $\mathbf{p}_i$) and $U$ (vertices $\mathbf{p}'_j$), (2) and
(3) give contributions $\frac{\mu_0}{4\pi} \frac{s_i s_j}{4 A_T A_U} K_{ij}$ with signs $s = \pm 1$
and

$$K_{ij} = \int_T \int_U \frac{(\mathbf{r} - \mathbf{p}_i) \cdot (\mathbf{r}' - \mathbf{p}'_j)}{R}\, dA'\, dA, \qquad R = |\mathbf{r} - \mathbf{r}'|.$$

The inner integral is singular when the triangles touch and is evaluated in closed form.
Let $\boldsymbol{\rho}_0$ be the projection of $\mathbf{r}$ onto the plane of $U$ and $d$ its
height above it. With $I_0 = \int_U dA'/R$ and $\mathbf{I}_1 = \int_U (\boldsymbol{\rho}' - \boldsymbol{\rho}_0)/R\, dA'$,

$$\int_U \frac{\mathbf{r}' - \mathbf{p}'_j}{R}\, dA' = \mathbf{I}_1 + (\boldsymbol{\rho}_0 - \mathbf{p}'_j)\, I_0 .$$

Both reduce to sums over the edges of $U$ (Wilton et al. 1984; Graglia 1993). For an edge
with unit tangent $\hat{\mathbf{l}}$, outward normal $\hat{\mathbf{u}}$, end coordinates
$s^\mp$ along $\hat{\mathbf{l}}$ measured from the foot of $\boldsymbol{\rho}_0$, signed
distance $t_0$ of $\boldsymbol{\rho}_0$ from its line, $R_0^2 = t_0^2 + d^2$ and
$R^\pm = \sqrt{R_0^2 + (s^\pm)^2}$:

$$f_2 = \ln \frac{R^+ + s^+}{R^- + s^-}, \qquad \beta = \arctan \frac{t_0 s^+}{R_0^2 + |d| R^+} - \arctan \frac{t_0 s^-}{R_0^2 + |d| R^-},$$

$$I_0 = \sum_\text{edges} \big(t_0 f_2 - |d|\, \beta\big), \qquad \mathbf{I}_1 = \sum_\text{edges} \frac{\hat{\mathbf{u}}}{2} \big(R_0^2 f_2 + s^+ R^+ - s^- R^-\big). \tag{4}$$

(The logarithm is evaluated in the equivalent forms that avoid cancellation when the edge's
line passes near $\boldsymbol{\rho}_0$, and $f_2 = 0$ on the line.) The outer integral over
$T$ is smooth enough for a symmetric 7-point rule of degree 5 (Dunavant 1985); for triangles
closer than 2.5 element sizes, where $I_0$ and $\mathbf{I}_1$ vary sharply, the same rule on
the four midpoint sub-triangles (28 points). For $T = U$ the $3 \times 3$ kernel is
symmetrised, $K \leftarrow (K + K^\mathsf{T})/2$, as (3) is symmetric.

### 3.3 Straight conductors

The partial mutual inductance of two straight filaments is Neumann's integral
$\frac{\mu_0}{4\pi} \iint d\mathbf{l} \cdot d\mathbf{l}' / R$:

- **Parallel**, at distance $\rho$, spanning $[a_1, b_1]$ and $[a_2, b_2]$ along a common
  axis, it integrates in closed form (Grover 1946):
  $$M = \frac{\mu_0}{4\pi} \big[G(b_1 - a_2) - G(a_1 - a_2) - G(b_1 - b_2) + G(a_1 - b_2)\big], \qquad G(u) = u \operatorname{arsinh}\frac{u}{\rho} - \sqrt{u^2 + \rho^2}. \tag{5}$$
  With $\rho$ the radius and equal spans, (5) is the partial self-inductance of a tube
  carrying its current on its surface.
- **At an angle**, by 8-point Gauss–Legendre quadrature on both, each subdivided so that
  the pieces are shorter than half their distance.
- **Perpendicular** conductors, and vertical holes against horizontal sheets, do not couple:
  $d\mathbf{l} \cdot d\mathbf{l}' = 0$.

A horizontal conductor and a triangle couple through the inner integral (4) at Gauss points
along the conductor, subdivided by their distance.

## 4. The minimum-energy problem

With the first element's current fixed at 1 A and $D$ the current law (one row per
triangle, one equation dropped in each connected group, where it is redundant), the branch
currents solve

$$\min_\mathbf{x}\ \tfrac{1}{2} \mathbf{x}^\mathsf{T} L\, \mathbf{x} \quad \text{subject to } D \mathbf{x} = \mathbf{b}, \tag{6}$$

and the inductance is $\mathbf{x}^\mathsf{T} L \mathbf{x}$, twice the minimum (1 A).

**The computed inductance is an upper bound.** Every discrete current that satisfies the
current law is an admissible current of the continuous model, so the discrete minimum of
(6) is the minimum of (1) over a subset, never below the minimum over all; refining the
mesh enlarges the subset, and the computed inductance decreases towards the model's from
above (`docs/validation/emi.md`, "Convergence of the loop inductance").

### 4.1 Solution in a basis of loops

A spanning forest of the graph whose nodes are the triangles and whose edges are the RWG
branches is grown breadth-first in each connected sheet, and one conductor between the
sheets closes the forests. The source's ampere returns to it through the trees: a
particular solution $\mathbf{x}_p$ with $D \mathbf{x}_p = \mathbf{b}$. Every branch outside
the trees closes a cycle with the tree path between its ends; these cycles are the columns
of $M$, and $D M = 0$. Every admissible current is $\mathbf{x} = \mathbf{x}_p + M \mathbf{y}$,
and (6) becomes unconstrained:

$$M^\mathsf{T} L M\, \mathbf{y} = -M^\mathsf{T} L\, \mathbf{x}_p, \qquad L_\text{loop} = \mathbf{x}_p^\mathsf{T} L\, \mathbf{x}_p + \mathbf{y}^\mathsf{T} M^\mathsf{T} L\, \mathbf{x}_p, \tag{7}$$

the second form following from the first by substitution. $M^\mathsf{T} L M$ is symmetric
positive definite (the columns of $M$ are independent) and is factorised by Cholesky
(Golub and Van Loan 2013). Its size, the number of cycles, is about a third of the branches.

When a conductor's end spans several sheets the forests cannot be closed this way, and (6)
is solved by the range-space method (Nocedal and Wright 2006, §16.2): with $L_F$ the free
branches' block, $\mathbf{l}$ their coupling to the source and $D_F$ their part of the current
law, the optimality conditions $L_F \mathbf{x} + \mathbf{l} + D_F^\mathsf{T} \boldsymbol{\lambda} = 0$,
$D_F \mathbf{x} = -\mathbf{d}_s$ give, with $L_F = C C^\mathsf{T}$,

$$(D_F L_F^{-1} D_F^\mathsf{T})\, \boldsymbol{\lambda} = \mathbf{d}_s - D_F L_F^{-1} \mathbf{l}, \qquad \mathbf{x} = -L_F^{-1}(\mathbf{l} + D_F^\mathsf{T} \boldsymbol{\lambda}),$$

and $L_\text{loop} = l_{ss} + 2 \mathbf{l}^\mathsf{T}\mathbf{x} + \mathbf{x}^\mathsf{T} L_F \mathbf{x} = l_{ss} + \mathbf{l}^\mathsf{T}\mathbf{x} + \mathbf{d}_s^\mathsf{T} \boldsymbol{\lambda}$
by the first condition. Both methods give the same solution.

### 4.2 Without the full matrix

$M^\mathsf{T} L M$, $M^\mathsf{T} L \mathbf{x}_p$ and $\mathbf{x}_p^\mathsf{T} L \mathbf{x}_p$ are
linear in $L$. Above a memory budget, the rows of $L$ that a block of triangles contributes
are formed, multiplied by $M$ and $\mathbf{x}_p$ and dropped; the result equals that of the
full matrix to rounding.

## 5. Resistance and inductance over frequency

At a frequency $f$ the branch currents minimise the complex power, and with the copper's
loss the system of (7) becomes

$$M^\mathsf{T} (j\omega L + R(f)) M\, \mathbf{y} = -M^\mathsf{T} (j\omega L + R(f))\, \mathbf{x}_p, \qquad Z = \mathbf{x}^\mathsf{T} (j\omega L + R(f))\, \mathbf{x},$$

complex symmetric, solved by an $LDL^\mathsf{T}$ factorisation. $R$ is assembled per layer
from the overlaps of the basis functions,
$\int \mathbf{f}_i \cdot \mathbf{f}_j\, dA = \frac{s_i s_j}{4A}\big(\overline{|\mathbf{r}|^2} - (\mathbf{p}_i + \mathbf{p}_j) \cdot \mathbf{c} + \mathbf{p}_i \cdot \mathbf{p}_j\big)$
on a triangle with centroid $\mathbf{c}$ and $\overline{|\mathbf{r}|^2} = \tfrac{1}{6}\big(\sum_a |\mathbf{v}_a|^2 + \sum_{a<b} \mathbf{v}_a \cdot \mathbf{v}_b\big)$,
times the surface impedance of the sheet.

**Surface impedance of a sheet.** A copper layer of thickness $t$ carrying current on both
faces (symmetric field) satisfies $d^2 E / dz^2 = j\omega\mu_0\sigma E$, so
$E(z) = E_s \cosh(\gamma z) / \cosh(\gamma t/2)$ with $\gamma = (1+j)/\delta$ and
$\delta = 1/\sqrt{\pi f \mu_0 \sigma}$. The sheet current is
$K = \int \sigma E\, dz = (2\sigma E_s / \gamma) \tanh(\gamma t / 2)$, and

$$Z_s = \frac{E_s}{K} = \frac{1 + j}{2\sigma\delta} \coth\frac{(1+j)\, t}{2\delta} \ \xrightarrow{\ f \to 0\ }\ \frac{1}{\sigma t}. \tag{8}$$

A plated hole's wall carries current on its inner face only; the same argument for one
face of a wall of thickness $t_p$ gives $\frac{1+j}{\sigma\delta} \coth \frac{(1+j) t_p}{\delta}$,
divided by the circumference $2\pi r$ for the impedance per length.

**Edges.** The sheet model concentrates current at the edges as $d^{-1/2}$, whose
resistance grows without bound as the mesh refines; copper of thickness $t$ spreads it over
about $t/2$. In the skin-effect limit a triangle nearer an edge than $t/2$ counts with the
share $d / (t/2)$ of its overlap, and the two overlaps are blended by
$\chi = |Z_s - 1/(\sigma t)| / |Z_s|$, from 0 at DC to 1 in the skin-effect limit. This is
an approximation, checked against references in `docs/validation/emi.md`; at high
frequency the resistance is a lower bound.

## 6. Moment and radiated field

**Current elements.** The solved currents are represented by elements $I\, d\mathbf{l}$ at
points: each triangle's current $\int_T \mathbf{K}\, dA$ at its centroid, each conductor in
pieces of at most 1 mm, and short joining pieces from each conductor end to the centroid of
the triangles that feed it, so that the represented current is closed exactly.

**Moment.** For a closed current the magnetic moment
$\mathbf{m} = \tfrac{1}{2} \sum \mathbf{r} \times I\, d\mathbf{l}$ does not depend on the
origin. A triangle's linear current contributes nothing beyond its centroid's term.

**Harmonics.** A trapezoidal current (amplitude $a$, period $T$, rise $t_r$, fall $t_f$,
width $\tau$ at half amplitude) is continuous and piecewise linear. Integrating its Fourier
coefficient twice by parts leaves only the jumps $\Delta_k$ of its slope at the corners
$t_k$:

$$c_n = -\frac{1}{T (n\omega_0)^2} \sum_k \Delta_k\, e^{-j n \omega_0 t_k}, \qquad \omega_0 = \frac{2\pi}{T},$$

and the amplitude of harmonic $n$ is $2|c_n|$ (Paul 2006, ch. 3).

**Field.** Outside the sources the closed current has no charge, so the scalar potential
vanishes and $\mathbf{E} = -j\omega \mathbf{A}$ with the retarded potential

$$\mathbf{A}(\mathbf{r}) = \frac{\mu_0}{4\pi} \sum I\, d\mathbf{l}\, \frac{e^{-jkR}}{R}, \qquad k = \omega / c_0,$$

near field included. Point elements of finite size leave a symmetric first moment
$Q = \tfrac{1}{2} \sum \big((\mathbf{x} - \mathbf{c}) I\, d\mathbf{l}^\mathsf{T} + I\, d\mathbf{l}\,(\mathbf{x} - \mathbf{c})^\mathsf{T}\big)$,
which a closed current does not have ($\int (r_i J_j + r_j J_i)\, dV = -\int r_i r_j \nabla \cdot \mathbf{J}\, dV = 0$);
its field, $\frac{\mu_0}{4\pi} (jk + 1/R) \frac{e^{-jkR}}{R} Q \hat{\mathbf{R}}$, is
subtracted. Over a conducting ground the image of each element (mirrored in the ground,
horizontal components reversed) is added. The reported field is the largest over directions
at the given distance in free space (polar angle in 5° and azimuth in 10° steps), or, over a
ground, the largest of the
horizontal and vertical components over azimuths in 10° steps and antenna heights in
0.25 m steps, as a radiated-emission site scans them. Far from the loop the result tends to
the small-loop field $E = \eta_0 k^2 m / (4\pi r)$ (Balanis 2016, ch. 5), the third option.

## 7. Mesh and refinement

A DC solution of 1 A through each copper connection of the loop marks the corridor of the
current; the copper within a margin of it is meshed (`docs/methods/dc.md`, §4.1), with an
element size growing from the corridor and rows of vertices along copper edges near it. Each
refinement pass recovers the sheet current at the vertices (area-weighted averages) and
refines where the triangle's current differs most from the recovered one, until the triangle
budget is used. Over a plane, the window grows until the moment settles.

The convergence in the triangle budget is slow and from above (§4), and the difference to
the previous pass understates the error at small budgets (`docs/solvers/emi.md`,
"Accuracy and cost").

## 8. Implementation

| step | where |
| --- | --- |
| window, corridor, mesh, conductors | `src/core/src/emi_loop.cpp`, `LoopSolver::setup`, `initial_seeds`, `discretise` |
| closed forms (4), quadrature rules | `src/core/src/peec.cpp`, `triangle_potential`, `triangle_rule_7`, `triangle_rule_28` |
| triangle pairs (3) | `emi_loop.cpp`, `pair_kernel` |
| conductors (5) | `peec.cpp`, `parallel_filaments`, `segment_mutual`; `emi_loop.cpp`, `filament_blocks` |
| current law, range-space solution | `emi_loop.cpp`, `LoopSolver::solve` |
| basis of loops (7), streamed form | `emi_loop.cpp`, `solve_by_loops`, `loop_system_streamed` |
| impedance (8) | `emi_loop.cpp`, `sheet_impedance`, `wall_impedance`, `frequency_sweep` |
| current elements, moment | `emi_loop.cpp`, end of `LoopSolver::solve` |
| harmonics, field | `emi_loop.cpp`, `trapezoid_harmonics`, `loop_field_per_ampere`, `radiate` |
| refinement passes, window growth | `emi_loop.cpp`, `LoopSolver::run`, `solve_loop` |

## References

- Balanis, C. A. (2016). *Antenna Theory: Analysis and Design*, 4th ed. Wiley.
- Dunavant, D. A. (1985). High degree efficient symmetrical Gaussian quadrature rules for
  the triangle. *International Journal for Numerical Methods in Engineering* 21(6),
  1129–1148. doi:10.1002/nme.1620210612.
- Golub, G. H., and Van Loan, C. F. (2013). *Matrix Computations*, 4th ed. Johns Hopkins
  University Press.
- Graglia, R. D. (1993). On the numerical integration of the linear shape functions times
  the 3-D Green's function or its gradient on a plane triangle. *IEEE Transactions on
  Antennas and Propagation* 41(10), 1448–1455. doi:10.1109/8.247786.
- Grover, F. W. (1946). *Inductance Calculations: Working Formulas and Tables*. Van
  Nostrand; reprinted by Dover (2004).
- Nocedal, J., and Wright, S. J. (2006). *Numerical Optimization*, 2nd ed. Springer.
  doi:10.1007/978-0-387-40065-5.
- Paul, C. R. (2006). *Introduction to Electromagnetic Compatibility*, 2nd ed. Wiley.
  doi:10.1002/0471758159.
- Rao, S. M., Wilton, D. R., and Glisson, A. W. (1982). Electromagnetic scattering by
  surfaces of arbitrary shape. *IEEE Transactions on Antennas and Propagation* 30(3),
  409–418. doi:10.1109/TAP.1982.1142818.
- Ruehli, A. E. (1972). Inductance calculations in a complex integrated circuit
  environment. *IBM Journal of Research and Development* 16(5), 470–481.
  doi:10.1147/rd.165.0470.
- Ruehli, A. E., Antonini, G., and Jiang, L. (2017). *Circuit Oriented Electromagnetic
  Modeling Using the PEEC Techniques*. Wiley–IEEE Press. doi:10.1002/9781119078388.
- Wilton, D. R., Rao, S. M., Glisson, A. W., Schaubert, D. H., Al-Bundak, O. M., and Butler,
  C. M. (1984). Potential integrals for uniform and linear source distributions on
  polygonal and polyhedral domains. *IEEE Transactions on Antennas and Propagation* 32(3),
  276–281. doi:10.1109/TAP.1984.1143304.
