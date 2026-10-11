# Thermal: method

This document derives the thermal solver: the reduction of three-dimensional conduction in
the board to a stack of sheets, the boundary conditions, the discretisation on non-matching
meshes, the nonlinear and coupled iterations, the transient scheme and the error estimate.
The user guide is `docs/solvers/thermal.md`; the validation, `docs/validation/thermal.md` and
`docs/validation/thermal-transients.md`. The DC solver it couples with is derived in
`docs/methods/dc.md`, whose notation it follows.

## 1. Symbols

| symbol | meaning | unit |
| --- | --- | --- |
| $T$ | temperature | °C (K for radiation) |
| $k_{xy}, k_z$ | in-plane and through-thickness conductivity | W/(m K) |
| $z_i$ | height of sheet $i$; $T_i(x, y)$ its temperature | m |
| $\tau_i = z_i - z_{i+1}$ | thickness of the slice between sheets $i$ and $i+1$ | m |
| $G_i$ | in-plane conductance of sheet $i$ | W/K |
| $c_i$ | conductance per area between sheets $i$ and $i+1$ | W/(m² K) |
| $h, \varepsilon$ | convection coefficient, emissivity of a face | W/(m² K), – |
| $r$ | series resistance per area (mask, interface) | m² K/W |
| $C$ | heat capacity | J/K |

## 2. Physical model

### 2.1 Conduction in the board

The board is a stack of copper and dielectric layers over its outline. Steady conduction
with heat input $Q$ per volume satisfies

$$-\nabla \cdot (\mathbf{k} \nabla T) = Q, \qquad \mathbf{k} = \mathrm{diag}(k_{xy}, k_{xy}, k_z), \tag{1}$$

with the conductivities of each layer: copper is isotropic; a glass-reinforced dielectric
conducts better along its plies than across them. The sides of the board are adiabatic; the
faces exchange heat with the surroundings (§2.4) or with heat sinks (§2.6).

### 2.2 Sheets: layerwise discretisation through the thickness

(1) is the stationarity condition of the energy

$$\mathcal{E}(T) = \frac{1}{2} \int_\text{board} \Big( k_{xy} |\nabla_{xy} T|^2 + k_z \Big(\frac{\partial T}{\partial z}\Big)^2 \Big) dV - \int_\text{board} Q\, T\, dV + (\text{face terms}). \tag{2}$$

Through the thickness the temperature is approximated by functions that are linear between
chosen heights $z_1 > z_2 > \dots > z_n$, the sheets:

$$T(x, y, z) = T_i(x, y)\, \frac{z - z_{i+1}}{\tau_i} + T_{i+1}(x, y)\, \frac{z_i - z}{\tau_i}, \qquad z_{i+1} \leq z \leq z_i. \tag{3}$$

This is a Galerkin method in $z$ with piecewise-linear functions, the layerwise approach of
laminate theory (Reddy 2004). Inserting (3) into (2), slice by slice:

- **Through-thickness conduction.** $\partial T / \partial z = (T_i - T_{i+1}) / \tau_i$ is
  constant in the slice, so its energy is exactly
  $\tfrac{1}{2} c_i \int (T_i - T_{i+1})^2 dA$ with $c_i = k_z / \tau_i$.
- **In-plane conduction.** The exact integral over the slice is
  $\tfrac{1}{2} k_{xy} \tau_i \int \tfrac{1}{3}\big(|\nabla T_i|^2 + \nabla T_i \cdot \nabla T_{i+1} + |\nabla T_{i+1}|^2\big) dA$.
  It is evaluated with the trapezoidal rule in $z$, which gives
  $\tfrac{1}{2} k_{xy} \tfrac{\tau_i}{2} \int \big(|\nabla T_i|^2 + |\nabla T_{i+1}|^2\big) dA$:
  each slice gives half its in-plane conductance $k_{xy} \tau_i$ to each of its two sheets.
  The rule is exact for temperatures that vary linearly through the slice, consistent to
  second order in $\tau_i$, and removes every in-plane coupling between different sheets, so
  that each sheet can be meshed on its own.

A copper layer is one sheet at its mid-plane, in which the copper conducts $k_\text{Cu} t$
where there is copper. Between the copper of an inner layer is resin, taken as the
dielectric below the layer, conducting $k_{xy} t$; between the copper of an outer layer is
air (nothing) or, optionally, solder mask filling the gaps, $k_\text{mask} t$. The copper is
taken as isothermal through its thickness.

Collecting the terms, each sheet satisfies

$$-\nabla \cdot (G_i \nabla T_i) + c_{i-1} (T_i - T_{i-1}) + c_i (T_i - T_{i+1}) = q_i, \tag{4}$$

with $G_i$ the sum of the in-plane shares and the copper of sheet $i$ and $q_i$ the heat
input per area. With only one sheet per copper layer, (4) is the classical layered (2.5D)
model; intermediate sheets refine (3) towards the three-dimensional solution.

**Placing the sheets.** The heat from a small source spreads through the dielectric within
a distance of the order of its size, so the slices are thin next to the copper. In a
dielectric of thickness $t_d$, slices of thickness $\tau_1, 2\tau_1, 4\tau_1, \dots$ (capped at
$\tau_\text{max}$) are laid from both faces inwards while more than twice the next slice
remains; the rest is divided into the fewest equal slices no thicker than $\tau_\text{max}$.
The defaults, $\tau_1 = 0.025$ mm and $\tau_\text{max} = 0.2$ mm, were chosen by the
convergence study of `docs/validation/thermal.md`.

### 2.3 Plated holes

A plated hole conducts along its wall between the layers whose copper it touches, as in the
DC model (`docs/methods/dc.md`, equation 2), with the thermal conductivity of copper:
$G = k_\text{Cu}\, \pi t_p (d + t_p) / L$, the wall isothermal around its circumference on
each layer.

### 2.4 Faces: convection and radiation through the mask

The outermost sheet of each face loses heat through the solder mask, a series resistance
per area $r = t_\text{mask} / k_\text{mask}$, to a surface at temperature $T_s$ cooled by
convection and radiation:

$$\frac{T - T_s}{r} = h (T_s - T_a) + \varepsilon \sigma_\text{SB} \big((T_s + 273.15)^4 - (T_a + 273.15)^4\big) \equiv q_\text{out}(T_s). \tag{5}$$

The right-hand side increases monotonically with $T_s$, so (5) has one root between $T_a$
and $T$; it is found by bisection (60 halvings), and the flux leaving the sheet is
$q(T) = (T - T_s)/r$. Radiation is to surroundings at the ambient temperature, without view
factors. Where a face is dielectric under air gaps of the outer copper, its sheet carries
the face (§2.2).

### 2.5 Sources

A source is an area of copper. A set power $P$ is distributed over the nodes of the area in
proportion to their lumped areas (§4.1). With a junction-to-case resistance $R_{jc}$, the
power enters a junction node instead, joined to each node $v$ of the area by the conductance
$(A_v / A) / R_{jc}$, so that the parallel sum is $1/R_{jc}$. A set temperature ties the
area's nodes to one unknown at that temperature; the heat it takes out is the residual of
its row.

### 2.6 Heat sinks

Under a heat sink, the face exchanges heat with the sink body through the mask and the
interface material: per area, $(T - T_\text{sink}) / (r_\text{mask} + r_\text{TIM})$ with
$r_\text{TIM} = t_\text{TIM} / k_\text{TIM}$. The body has one temperature: either set, or
joined to ambient by its resistance $R_{sa}$ with $T_\text{sink}$ an unknown.

## 3. Variational form

Collecting (2) with (3), the faces, sinks, junctions and barrels, the discrete-in-$z$
problem minimises

$$J(T) = \sum_i \frac{1}{2}\int G_i |\nabla T_i|^2 + \sum_i \frac{1}{2}\int c_i (T_i - T_{i+1})^2 + \sum_\text{barrels} \frac{G}{2} (T_a - T_b)^2 + (\text{junctions, sinks}) + \int_\text{faces} \Phi(T) - \sum_i \int q_i T_i, \tag{6}$$

with $\Phi' = q$ the face flux of §2.4 (convex, because $q$ increases with $T$). Every term is
convex, and strictly so once heat can leave (a face with $h > 0$ or $\varepsilon > 0$, a set
temperature, or a sink): the solution is unique.

## 4. Discretisation

### 4.1 Sheets

Each sheet is triangulated over the board outline (without non-plated holes), with the
copper outlines, hole walls and source areas of its layer as constraints (`docs/methods/dc.md`,
§4.1), and $T_i$ is piecewise linear. The in-plane term gives the element matrix of the DC
solver, $K^e_{ij} = G\,(b_i b_j + c_i c_j)/(4A)$, with $G$ constant per triangle (copper
or fill). The lumped area of a node is a third of the area of each triangle around it,
$A_v = \sum_{T \ni v} A_T / 3$; it carries the sources, faces and capacities.

### 4.2 Coupling between sheets on non-matching meshes

The sheets have different meshes, so the term $\tfrac{1}{2} \int c (T_i - T_{i+1})^2 dA$ is
evaluated by nodal quadrature from both sides: half of it on the nodes of sheet $i$, half
on those of sheet $i+1$,

$$\frac{1}{2} \int c (T_i - T_{i+1})^2 \approx \frac{1}{2} \sum_{v \in i} \frac{c A_v}{2} \big(T_v - T_{i+1}(x_v)\big)^2 + \frac{1}{2} \sum_{w \in i+1} \frac{c A_w}{2} \big(T_w - T_i(x_w)\big)^2, \tag{7}$$

where $T_{i+1}(x_v)$ interpolates the neighbouring sheet at the point below node $v$:
$\sum_k \lambda_k T_k$ over the triangle containing it, with barycentric weights $\lambda_k$.
Each term is a square of a linear form in four unknowns, so it contributes the symmetric
positive semi-definite matrix $\tfrac{1}{2} c A_v\, \mathbf{u}\mathbf{u}^\mathsf{T}$ with
$\mathbf{u} = (1, -\lambda_1, -\lambda_2, -\lambda_3)$. The scheme is exact for equal constant
temperatures and consistent for smooth fields; symmetric splitting makes the global matrix
symmetric.

### 4.3 Faces, sinks, junctions

Faces and sink contacts are lumped at the nodes of the outer sheets with the same areas
$A_v$. A sink contact adds the conductance $A_v / (r_\text{mask} + r_\text{TIM})$ between the
node and the sink unknown; a sink with $R_{sa}$ adds $1/R_{sa}$ to its diagonal and
$T_a / R_{sa}$ to its load. The barrel and junction conductances are two-node stamps, as in
the DC solver.

### 4.4 Linearisation of the faces

The face flux is nonlinear (radiation, and the mask in series with it). At a state $T^0$
it is replaced by its tangent, $q(T) \approx q(T^0) + g\,(T - T^0)$, with the slope
$g = dq/dT$ taken by a central difference of ±0.01 K. Per node this adds $A_v g$ to the
diagonal and $A_v (g T^0_v - q(T^0_v))$ to the load: one Newton step. The linear system,
symmetric positive definite, is solved as in the DC solver (`docs/methods/dc.md`, §4.4), and
the steps repeat until the temperature settles (§5). With convection alone and no mask the
face term is linear and one step is exact.

## 5. Electrothermal coupling

With DC ports, the copper's conductivity depends on the temperature
(`docs/methods/dc.md`, §2.2) and the Joule heat is a source. The coupled state is a fixed
point of $T \mapsto \mathcal{G}(T)$: solve DC with $\sigma(T)$, deposit its Joule heat, solve
the thermal problem. The Joule heat of each DC triangle is placed at its centroid on the
thermal sheet of its layer by barycentric weights, and that of each barrel half at each
of its two ends.

The plain iteration $T^{k+1} = \mathcal{G}(T^k)$ contracts by about the feedback gain
$a \alpha$, the fractional rise of the losses per kelvin times the temperature rise per watt;
near thermal runaway ($a \alpha \to 1$) it slows without bound, and beyond it there is no
steady state. Anderson acceleration (Anderson 1965; Walker and Ni 2011) combines the last
$m \leq 3$ iterates. With $f_k = \mathcal{G}(T^k) - T^k$, it solves the least-squares
problem

$$\gamma = \arg\min_\gamma \Big\| f_k - \sum_{j} \gamma_j (f_{j+1} - f_j) \Big\|_2$$

over the differences of consecutive residuals (normal equations, regularised by $10^{-10}$
of their trace), and continues from

$$T^{k+1} = \mathcal{G}(T^k) - \sum_j \gamma_j \big(\mathcal{G}(T^{j+1}) - \mathcal{G}(T^j)\big).$$

The iteration stops when the largest change of a nodal temperature is below the tolerance
(0.02 K). It reports thermal runaway when the change has grown three times in a row, or
when a temperature falls below $20 - 1/\alpha$ °C, where the linear law would make the
resistivity vanish.

## 6. Transients

With $C$ the heat capacity of each unknown (lumped: the copper, the resin between the copper
of inner layers, half of each adjacent slice and, on the outer sheets, the mask, each by its
volumetric capacity times the node's share of the volume; plus junction and sink bodies),
the semi-discrete system $C \dot{T} + K T = f(t)$ is integrated by the backward Euler
method,

$$\Big(\frac{C}{\Delta t} + K\Big) T^{n+1} = \frac{C}{\Delta t} T^n + f(t^{n+1}), \tag{8}$$

which is unconditionally stable and first-order accurate in time (Thomée 2006, ch. 1). The
faces keep the slopes $g$ of the steady solution at the peak power, so that the matrix of
(8) depends on $\Delta t$ only and one multigrid hierarchy serves every step of a given
length; with convection alone this is exact. Steps start at `step_s`, double every
`steps_per_doubling` steps up to `max_step_s`, and land on every corner of a power schedule,
after which they start again small. Each step solves for the change $T^{n+1} - T^n$, to a
tolerance relative to the residual of the previous state.

## 7. Error estimate

The estimate in the energy norm has two parts per triangle:

- **In-plane conduction**, by gradient recovery as in the DC solver
  (`docs/methods/dc.md`, §6.2), with the gradient recovered within copper and within the
  dielectric separately: across a copper edge $G$ jumps and so does the exact gradient, and
  averaging over the jump would report an error that refinement cannot remove.
- **Vertical coupling.** The coupling (7) lumps $\int c\, w^2 dA$, with $w$ the difference
  of the temperatures of two sheets, at the nodes. For a linear $w$ on a triangle with
  vertex values $w_1, w_2, w_3$,
  $$\frac{A}{3} \sum_i w_i^2 - \int_T w^2 dA = \frac{A}{3} \sum_i w_i^2 - \frac{A}{6} \Big(\sum_i w_i^2 + \sum_{i<j} w_i w_j\Big) = \frac{A}{12} \sum_{i<j} (w_i - w_j)^2,$$
  so the lumping error of the term grows with the variation of $w$ across the triangle. Per
  triangle and per neighbouring sheet, with the half coupling assembled from each side, the
  estimate adds
  $$\eta_{V,T}^2 = \frac{c}{2} \cdot \frac{A}{12} \sum_{i<j} (w_i - w_j)^2 .$$
  It dominates near small sources, where the thin slices next to the copper, stiffly
  coupled and nearly without in-plane conduction, change steeply under coarse elements
  (`docs/validation/thermal.md`, "Real board against Elmer").

The relative error is $e = \sqrt{\sum \eta^2 / E}$ with $E = T^\mathsf{T} (K + \text{faces}) T$,
temperatures taken above ambient; refinement marks as in the DC solver with the margin
$m = 0.8$. The error of the lumped face terms, of order $(h_e / \ell)^2$ for elements of size
$h_e$ and a thermal length $\ell$, is not estimated.

## 8. Reported quantities

- **Source temperatures.** Mean over the area, weighted by the lumped areas; maximum over
  its nodes; the junction unknown when there is one.
- **Heat balance.** Heat in (sources, Joule heat) against heat out through the faces, the
  sinks and the fixed temperatures (each the residual of its row), and, in a transient, the
  rate of change of the stored heat $C (T^{n+1} - T^n)/\Delta t$.

## 9. Implementation

| step | where |
| --- | --- |
| sheets, slices, fills, barrels (§2.2, §2.3) | `src/core/src/thermal.cpp`, `build_model`, `slice_thicknesses` |
| meshes of the sheets (§4.1) | `thermal.cpp`, `Run::mesh_pieces`; `mesh2d.cpp` |
| unknowns, lumped areas, faces, sinks, sources (§2.5, §2.6, §4.3) | `thermal.cpp`, `Run::make_layout` |
| matrix: in-plane, vertical coupling (7), barrels, junctions, sinks | `thermal.cpp`, `Run::make_system` |
| face flux (5) and its slope | `thermal.cpp`, `surface_flux`, `surface_slope`, `Run::face_terms` |
| steady passes, Newton, coupling, Anderson (§4.4, §5) | `thermal.cpp`, `Run::steady` |
| Joule heat transfer | `thermal.cpp`, `Run::joule_load` |
| transient (8) | `thermal.cpp`, `Run::transient` |
| error estimate (§7) | `fem_util.cpp`, `zz_estimate`; `thermal.cpp`, `Run::add_vertical_estimate` |

## References

- Anderson, D. G. (1965). Iterative procedures for nonlinear integral equations. *Journal
  of the ACM* 12(4), 547–560. doi:10.1145/321296.321305.
- Reddy, J. N. (2004). *Mechanics of Laminated Composite Plates and Shells: Theory and
  Analysis*, 2nd ed. CRC Press.
- Thomée, V. (2006). *Galerkin Finite Element Methods for Parabolic Problems*, 2nd ed.
  Springer. doi:10.1007/3-540-33122-0.
- Walker, H. F., and Ni, P. (2011). Anderson acceleration for fixed-point iterations. *SIAM
  Journal on Numerical Analysis* 49(4), 1715–1735. doi:10.1137/10078356X.
- Further references of the shared steps (meshing, multigrid, gradient recovery):
  `docs/methods/dc.md`.
