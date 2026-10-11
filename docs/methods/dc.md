# DC current: method

This document derives the DC solver from the conduction equation to the numbers it reports,
and states where each step is implemented. The user-facing description, the project-file
keys and the limits are in `docs/solvers/dc.md`; the measured accuracy is in
`docs/validation/dc.md`. References are listed at the end and in `docs/references.bib`.

## 1. Symbols

| symbol | meaning | unit |
| --- | --- | --- |
| $V$ | electric potential | V |
| $\sigma$ | conductivity of copper | S/m |
| $t$ | thickness of a copper layer | m |
| $g = \sigma t$ | sheet conductance | S |
| $\mathbf{K} = -g \nabla V$ | sheet current density | A/m |
| $\Omega_\ell$ | copper of layer $\ell$, a planar region | |
| $\Gamma_p$ | area of port $p$ (a union of triangles) | |
| $a(\cdot,\cdot)$ | energy (bilinear) form | W/V² |
| $P$ | dissipated power | W |

## 2. Physical model

### 2.1 From three-dimensional conduction to sheets

In a conductor at DC, charge conservation $\nabla \cdot \mathbf{J} = 0$ and Ohm's law
$\mathbf{J} = -\sigma \nabla \varphi$ give $\nabla \cdot (\sigma \nabla \varphi) = 0$. A copper
layer occupies $\Omega_\ell \times [z_0, z_0 + t]$, with no current across its top and bottom
faces except where it meets a hole wall. If the potential does not vary through the
thickness, $\varphi(x, y, z) = V(x, y)$, integrating the conservation law over $z$ gives

$$\int_{z_0}^{z_0 + t} \nabla \cdot (\sigma \nabla \varphi)\, dz = \nabla_{xy} \cdot (\sigma t\, \nabla_{xy} V) + \Big[\sigma \frac{\partial \varphi}{\partial z}\Big]_{z_0}^{z_0+t} = 0,$$

and the bracket vanishes on insulated faces. So in each layer

$$\nabla \cdot (g \nabla V) = 0 \quad \text{in } \Omega_\ell, \qquad g \nabla V \cdot \mathbf{n} = 0 \quad \text{on the copper edges}. \tag{1}$$

The assumption $\partial \varphi / \partial z = 0$ holds when $t$ is small against the distance
over which the current turns. It fails where current enters the copper from a hole wall or a
pad: there the current crowds within about $t$ of the contact. This constriction is not
modelled; against a three-dimensional solution it amounts to about 0.1 mΩ per via
(`docs/validation/dc.md`, section "3D reference (Elmer)").

### 2.2 Temperature

The resistivity of copper rises linearly with temperature near room temperature:

$$\sigma(T) = \frac{\sigma_{20}}{1 + \alpha (T - 20\ \mathrm{^\circ C})}, \qquad \alpha = 0.00393\ \mathrm{K^{-1}},$$

with $\sigma_{20} = 5.8 \times 10^7$ S/m for annealed copper (IEC 60028) unless the stackup
states another value. With a temperature field (coupled thermal solution), $\sigma$ is
evaluated at the centroid of each triangle, so $g$ is constant per triangle.

### 2.3 Plated holes

A plated hole of finished diameter $d$ has a copper wall of thickness $t_p$ outside that
diameter. Its cross-section is the annulus between the radii $d/2$ and $d/2 + t_p$:

$$A = \pi \Big(\frac{d}{2} + t_p\Big)^2 - \pi \Big(\frac{d}{2}\Big)^2 = \pi t_p (d + t_p).$$

A slot whose centres are $\ell$ apart adds two straight walls, $2 t_p \ell$. Between two
consecutive layers $a$ and $b$ that the hole connects, the current flows along the wall over
the distance $L$ between the layers' mid-planes, where the sheets of (1) lie, so the wall is
a conductance

$$G = \frac{\sigma_p A}{L}, \tag{2}$$

with $\sigma_p$ the plating conductivity at the mean temperature of the two layers. The wall
is taken as equipotential around its circumference on each layer: the current enters it
from all sides of the hole.

### 2.4 Ports

A port is an area $\Gamma_p$ of copper held at one potential, as a soldered pin holds it.
Either its potential is set, $V = U_p$ on $\Gamma_p$, or the current it injects is set:
$V = V_p$ on $\Gamma_p$ with $V_p$ unknown and

$$\oint_{\partial \Gamma_p} \mathbf{K} \cdot \mathbf{n}\, ds = I_p, \tag{3}$$

the outward sheet current around the port equal to the set current. Copper connected only
to set-current ports has no potential reference; one of its ports is taken as 0 V, and the
set currents must sum to zero there, since no current can leave.

## 3. Variational form

Let $\mathcal{V}$ be the potentials that are piecewise smooth on the copper, constant on each
port area and on each hole wall of each layer. For $u, w \in \mathcal{V}$ define

$$a(u, w) = \sum_\ell \int_{\Omega_\ell} g\, \nabla u \cdot \nabla w \, dA + \sum_{\text{barrels}} G\, (u_a - u_b)(w_a - w_b), \tag{4}$$

where $u_a$ is the value of $u$ on the hole wall at layer $a$. $a(u, u)$ is the power
dissipated by the potential $u$: the Joule heat $\int g |\nabla u|^2$ in the sheets and
$G (u_a - u_b)^2$ in each barrel.

Multiplying (1) by a test function $w \in \mathcal{V}$, integrating by parts, and using the
insulated edges, the equipotential walls (where the barrel currents enter) and (3):

$$a(V, w) = \sum_{p \in \mathcal{I}} I_p\, w_p \quad \text{for all } w \in \mathcal{V}_0, \tag{5}$$

where $\mathcal{I}$ are the set-current ports, $w_p$ is the value of $w$ on $\Gamma_p$, and
$\mathcal{V}_0 \subset \mathcal{V}$ are the potentials that vanish on the set-voltage ports
and on the reference ports. $V$ itself takes the set voltages on their ports. Equation (5)
is the stationarity condition of the energy

$$J(v) = \tfrac{1}{2} a(v, v) - \sum_{p \in \mathcal{I}} I_p\, v_p \tag{6}$$

over the potentials with the set voltages. $a$ is symmetric, and positive on $\mathcal{V}_0$
for every connected group of copper that holds at least one fixed potential, so $J$ is
strictly convex there and (5) has exactly one solution (Brenner and Scott 2008, ch. 2).

## 4. Discretisation

### 4.1 Mesh

Each copper island connected to a port is triangulated by a constrained Delaunay
triangulation (CDT library) whose constraints are the copper outline, the hole walls and the
port outlines, so that each equipotential area is an exact union of triangles. Triangles
are refined by Delaunay refinement (insertion of circumcentres, with encroached boundary
edges split first) until their smallest angle is about 20°, the method of Ruppert (1995);
see Shewchuk (2002) for its analysis. Outlines are simplified within 5 µm and hole walls are
polygons of at least 16 sides of equal area: these are the only geometric approximations.

### 4.2 Linear elements

$V$ is approximated by $V_h = \sum_i V_i N_i$, with $N_i$ the piecewise linear function
that is 1 at node $i$ and 0 at the others. On a triangle with vertices
$(x_1, y_1), (x_2, y_2), (x_3, y_3)$ in counter-clockwise order and area $A$, the gradients
of the three shape functions are constant:

$$\nabla N_i = \frac{1}{2A} \begin{pmatrix} b_i \\ c_i \end{pmatrix}, \qquad b_i = y_j - y_k, \quad c_i = x_k - x_j,$$

for $(i, j, k)$ cyclic. With $g$ constant on the triangle, the element matrix is

$$K^e_{ij} = \int_T g\, \nabla N_i \cdot \nabla N_j \, dA = g\, \frac{b_i b_j + c_i c_j}{4A}. \tag{7}$$

The code computes $2A$ as `area2` and adds `gs * (bb[i]*bb[j] + cc[i]*cc[j]) / (2*area2)`,
which is (7). The quadrature is exact, since the integrand is constant.

Each barrel adds the matrix $G \begin{pmatrix} 1 & -1 \\ -1 & 1 \end{pmatrix}$ between the
unknowns of its two hole walls, which is the barrel term of (4).

### 4.3 Equipotential areas

The nodes inside a port area or on a hole wall of one layer are tied to one unknown: the
discrete space keeps $V_h$ constant there, as $\mathcal{V}$ requires. In the matrix this
sums the rows and columns of the tied nodes. A set-voltage port's unknown is fixed; a
set-current port's unknown stays free, and its right-hand side is $I_p$, which is (5) with
$w = N_p$, the function that is 1 on the whole port area. One solution therefore gives all
ports, whatever their number.

Writing the unknowns as free (f) and fixed (c), the system is

$$K_{ff} V_f = F_f - K_{fc} V_c, \tag{8}$$

with $F_f$ the set currents. $K_{ff}$ is symmetric; it is positive definite because every
connected group of copper holds a fixed unknown (a set voltage, the reference of a group of
set-current ports, or 0 V for copper without ports, which carries no current).

### 4.4 Solution

(8) is solved by the conjugate gradient method (Hestenes and Stiefel 1952) preconditioned
with algebraic multigrid: classical Ruge–Stüben coarsening with direct interpolation
(Stüben 1999) and a sparse approximate inverse smoother, SPAI-0 (Bröker and Grote 2002), as
implemented in AMGCL (Demidov 2019). Iteration stops at a relative residual of
$10^{-10}$. On the two benchmark boards, classical coarsening needs 24 and 30 iterations
against 31 and 30 for smoothed aggregation, with the same solution
(`docs/validation/dc.md`, section "Performance").

## 5. Reported quantities

- **Sheet current and dissipation.** On each triangle, $\nabla V_h$ is constant:
  $\nabla V_h = \frac{1}{2A} \sum_i V_i (b_i, c_i)$. The sheet current is
  $\mathbf{K} = -g \nabla V_h$, the volume current density $\mathbf{K}/t$, and the power
  $g |\nabla V_h|^2 A$.
- **Barrel currents.** $I_b = G (V_a - V_b)$, power $I_b^2 / G$.
- **Port currents.** The current into the copper at an unknown is the row of the full
  system, $(K V)_d$. A set-current port reports its set current; the set-voltage ports on
  one unknown share what the set-current ports there do not supply.
- **Power balance.** With $V$ the discrete solution, $V^\mathsf{T} K V = \sum_d V_d (K V)_d$
  is both the sum of the element and barrel dissipations and the power the ports deliver,
  $\sum_p U_p I_p$. The two are reported separately; their agreement checks the assembly.

## 6. Error

### 6.1 The energy error and the resistance

Let $V$ be the exact solution of (5) and $V_h$ the discrete one. Because the discrete space
is a subspace of $\mathcal{V}$ (linear elements on triangles that fit the equipotential areas
exactly) and the set voltages are represented exactly, the error $e = V - V_h$ satisfies
Galerkin orthogonality, $a(e, w_h) = 0$ for all discrete $w_h$ that vanish on the fixed
unknowns. Two cases follow (compare Brenner and Scott 2008, §2.5).

**Set voltages only.** $V - V_h \in \mathcal{V}_0$, so with $w = V_h - V$ in (5) and the
symmetry of $a$,

$$a(V_h, V_h) = a(V, V) + a(e, e) \quad \Rightarrow \quad P_h = P + \|e\|_a^2 \geq P.$$

For two ports at a voltage difference $U$, $R = U^2 / P$, so $R_h \leq R$.

**Set currents with one reference.** $P = a(V, V) = \sum_p I_p V_p$ by (5) with $w = V$,
and the minimum of (6) is $J(V) = -\tfrac{1}{2} P$. The discrete solution minimises $J$ over
a subspace, so $J(V_h) \geq J(V)$; with Galerkin orthogonality,

$$P_h = P - \|e\|_a^2 \leq P.$$

For two ports carrying $I$, $R = P / I^2$, so again $R_h \leq R$.

In both cases the discrete copper conducts too well: **the computed resistance between two
ports converges from below**, and its relative error is

$$\frac{R - R_h}{R} \approx \frac{\|e\|_a^2}{a(V, V)}. \tag{9}$$

The validation cases converge from below, as (9) requires (`docs/validation/dc.md`). The
statement covers the discretisation error only: the simplified outlines and polygonal hole
walls change the geometry itself, by up to 5 µm.

### 6.2 Gradient recovery

The exact gradient is continuous where $g$ is; the discrete one jumps between triangles.
Zienkiewicz and Zhu (1987) recover a continuous gradient and take its difference from the
discrete one as the error. At each node $i$, the recovered gradient is the area-weighted
mean of the gradients of the triangles around it,

$$\mathbf{G}^*_i = \frac{\sum_{T \ni i} A_T\, \nabla V_h|_T}{\sum_{T \ni i} A_T},$$

leaving out triangles inside equipotential areas (they are boundary conditions, not field).
The estimate on a triangle evaluates $\int_T g |\mathbf{G}^* - \nabla V_h|^2 dA$ with the
vertex quadrature rule:

$$\eta_T^2 = g\, A_T\, \frac{1}{3} \sum_{i \in T} \big| \mathbf{G}^*_i - \nabla V_h|_T \big|^2, \tag{10}$$

and the relative error in the energy norm is estimated as

$$e = \sqrt{\frac{\sum_T \eta_T^2}{\sum_T g\, A_T\, |\nabla V_h|_T|^2}}.$$

By (9), $e^2$ estimates the relative error of a two-port resistance and of the total
dissipation. Recovery estimators are asymptotically exact for smooth solutions on regular
meshes and underestimate near singularities (Ainsworth and Oden 2000, ch. 4): at re-entrant
corners of copper, where the gradient is unbounded, the measured error is up to twice $e^2$
(`docs/validation/dc.md`).

### 6.3 Refinement

While $e$ exceeds the target $\tau$, triangles are refined. For linear elements and a smooth
solution, $\eta_T^2$ scales with $h_T^2$ times the integral of the second derivatives over
$T$; splitting a triangle into $n^2$ children of size $h_T / n$ therefore divides its
contribution by about $n^2$. The goal is

$$\sum_T \eta_T^2 \leq (m \tau)^2 \sum_T g A_T |\nabla V_h|^2, \qquad m = 0.9,$$

the margin $m$ aiming a little below the target. The split factor follows from the overall
reduction needed, $n = \min(3, \max(2, \lceil \sqrt{\eta^2_\text{total} / \text{goal}} \rceil))$;
triangles are marked in descending order of $\eta_T^2$, each mark lowering the predicted total
by $\eta_T^2 (1 - 1/n^2)$, until the prediction meets the goal or the triangle budget is
reached. Each marked triangle receives the $n(n+1)/2$ points of a barycentric lattice
shifted off its edges, $\lambda = ((i + \tfrac{1}{3})/n, (j + \tfrac{1}{3})/n)$, and the
island is triangulated again from all previous vertices and the new ones (§4.1).

## 7. Implementation

| step | where |
| --- | --- |
| islands, ports, hole links, the passes | `src/core/src/dc.cpp`, `solve_dc` |
| mesh (§4.1) | `src/core/src/mesh2d.cpp`, `mesh_region`; CDT |
| conductivity at temperature (§2.2) | `dc.cpp`, `sigma_at` |
| element matrix (7), barrels (2) | `dc.cpp`, "Element and barrel conductances" |
| tied nodes, reference potentials, system (8) | `dc.cpp`, "Unknowns" and "Boundary conditions" |
| CG with AMG (§4.4) | `src/core/src/fem_util.cpp`, `SpdSolver`, `solve_spd`; AMGCL |
| recovery estimate (10) | `fem_util.cpp`, `zz_estimate` |
| marking and lattice points (§6.3) | `fem_util.cpp`, `mark_and_seed` |
| port currents, fields, barrel currents (§5) | `dc.cpp`, after the solution |

## References

- Ainsworth, M., and Oden, J. T. (2000). *A Posteriori Error Estimation in Finite Element
  Analysis*. Wiley. doi:10.1002/9781118032824.
- Brenner, S. C., and Scott, L. R. (2008). *The Mathematical Theory of Finite Element
  Methods*, 3rd ed. Springer. doi:10.1007/978-0-387-75934-0.
- Bröker, O., and Grote, M. J. (2002). Sparse approximate inverse smoothers for geometric
  and algebraic multigrid. *Applied Numerical Mathematics* 41(1), 61–80.
  doi:10.1016/S0168-9274(01)00110-6.
- Demidov, D. (2019). AMGCL: an efficient, flexible, and extensible algebraic multigrid
  implementation. *Lobachevskii Journal of Mathematics* 40(5), 535–546.
  doi:10.1134/S1995080219050056.
- Hestenes, M. R., and Stiefel, E. (1952). Methods of conjugate gradients for solving linear
  systems. *Journal of Research of the National Bureau of Standards* 49(6), 409–436.
  doi:10.6028/jres.049.044.
- IEC 60028 (1925). International standard of resistance for copper.
- Ruppert, J. (1995). A Delaunay refinement algorithm for quality 2-dimensional mesh
  generation. *Journal of Algorithms* 18(3), 548–585. doi:10.1006/jagm.1995.1021.
- Shewchuk, J. R. (2002). Delaunay refinement algorithms for triangular mesh generation.
  *Computational Geometry* 22(1–3), 21–74. doi:10.1016/S0925-7721(01)00047-5.
- Stüben, K. (1999). *Algebraic Multigrid (AMG): An Introduction with Applications*. GMD
  Report 70, GMD-Forschungszentrum Informationstechnik; also in Trottenberg, U., Oosterlee,
  C., and Schüller, A. (2001), *Multigrid*, Academic Press, appendix A.
- Zienkiewicz, O. C., and Zhu, J. Z. (1987). A simple error estimator and adaptive procedure
  for practical engineering analysis. *International Journal for Numerical Methods in
  Engineering* 24(2), 337–357. doi:10.1002/nme.1620240206.
