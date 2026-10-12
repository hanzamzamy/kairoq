# PDN impedance and plane-pair cavities: method

This document derives the impedance of a power distribution network: the plane pair as a
two-dimensional cavity, its finite-element modes and their computation, the impedance
between ports from the modes, the branches of capacitors and regulators with their
mounting inductance, and the error. The same cavity gives the fundamental mode of the via
model of a channel (`docs/methods/networks.md`, section 4.3). The user guide is
`docs/solvers/pdn.md`; the validation, `docs/validation/pdn.md`.

## 1. Symbols

| symbol | meaning | unit |
| --- | --- | --- |
| $\Omega$ | the region both planes of a pair cover, without the port circles | |
| $h$ | separation of the planes | m |
| $V(x, y)$ | voltage between the planes | V |
| $\mathbf{K}_s$ | surface current density of the pair (current per width) | A/m |
| $z, y$ | series impedance and shunt admittance of the pair per square | Ω, S |
| $Z_s$ | surface impedance of one plate | Ω |
| $k^2 = -zy$ | squared wavenumber | 1/m² |
| $\mathbf{K}, \mathbf{M}$ | stiffness and mass matrices of linear elements on $\Omega$ | , m² |
| $\lambda_m, \boldsymbol\varphi_m$ | eigenvalues and $\mathbf{M}$-orthonormal eigenvectors | 1/m², |
| $\mathbf{f}_p$ | port $p$'s averaging vector over its circle | |

## 2. Model

### 2.1 The pair as a two-dimensional line

Between two parallel plates $h$ apart, much less than a wavelength, the field is the
pair's fundamental mode: $E_z$ uniform across the height and the currents on the plates
equal and opposite. With the voltage $V = -\int E_z\, dz$ and the surface current
$\mathbf{K}_s$ on the upper plate, the telegrapher's equations in two dimensions are

$$\nabla V = -z\, \mathbf{K}_s, \qquad \nabla \cdot \mathbf{K}_s = -y\, V, \tag{1}$$

$$z = j\omega\mu_0 h + 2 Z_s, \qquad y = \frac{j\omega\varepsilon_0\varepsilon_r(\omega)}{h},$$

so that $\nabla^2 V + k^2 V = 0$ with $k^2 = -zy$. The series impedance holds the magnetic
field between the plates and the two plates' surface impedance; the shunt admittance the
dielectric, complex with its loss. Several dielectrics between the plates are in series
for the uniform $E_z$: $\varepsilon_r = t / \sum_i (t_i/\varepsilon_{r,i}(f))$, each from the
wideband model (`docs/methods/si.md`, section 4).

**Surface impedance.** A plate of conductivity $\sigma$ and thickness $t$ carrying current on
the face towards the pair, with the other face free, has

$$Z_s = \frac{1 + j}{\sigma\delta} \coth\frac{(1 + j)\, t}{\delta}, \qquad \delta = \frac{1}{\sqrt{\pi f \mu_0 \sigma}},$$

the ratio $E/H$ at the face for the field $E \propto \cosh(\gamma(t - z))$ in the plate,
$\gamma = (1 + j)/\delta$, whose magnetic field vanishes at the free face. It tends to the
sheet resistance $1/(\sigma t)$ for $t \ll \delta$ and to $(1 + j)/(\sigma\delta)$ for
$t \gg \delta$.

### 2.2 Boundary conditions

- **Outline.** At the edge of the region the current normal to it vanishes: an open
  (magnetic) wall, $\partial V/\partial n = 0$. Radiation from the edge is neglected.
- **Holes.** Holes in the region smaller than 2.5 mm across (antipads, clearances) are
  filled; larger ones (splits, slots) are open walls of their own.
- **Ports.** A barrel of radius $a$ that carries the current $I$ into the pair is a circle
  removed from the region. The current spreads radially, $\mathbf{K}_s \cdot \hat\rho = I/(2\pi a)$
  on the circle; with the normal $\mathbf{n}$ pointing out of $\Omega$ (into the circle,
  $\mathbf{n} = -\hat\rho$), (1) gives

  $$\frac{\partial V}{\partial n} = -\nabla V \cdot \hat\rho = z\, \frac{I}{2\pi a}. \tag{2}$$

  The port's voltage is the mean of $V$ over its circle.
- **Shorted ports.** A barrel of the reference net joining the two plates holds $V = 0$ on
  its circle.

## 3. Finite elements

### 3.1 Weak form

Multiplying $\nabla^2 V + k^2 V = 0$ by a test function $v$ vanishing on the shorted
circles and integrating by parts, with (2) on the open ports,

$$\int_\Omega \nabla V \cdot \nabla v\, dA - k^2 \int_\Omega V v\, dA = z \sum_p I_p\, \frac{1}{2\pi a_p}\oint_{C_p} v\, ds. \tag{3}$$

With linear elements on triangles, $V = \sum_i V_i N_i$, this is

$$(\mathbf{K} - k^2 \mathbf{M})\, \mathbf{V} = z \sum_p I_p\, \mathbf{f}_p, \qquad K_{ij} = \int \nabla N_i \cdot \nabla N_j, \quad M_{ij} = \int N_i N_j, \tag{4}$$

with $(\mathbf{f}_p)_i = \oint_{C_p} N_i\, ds / |C_p|$ by the trapezoidal rule on the circle's
edges, divided by the discrete circumference $|C_p|$. The same vector gives the port's
mean voltage, $\mathbf{f}_p^\mathsf{T}\mathbf{V}$, so the impedance matrix of the ports is

$$\mathbf{Z} = z\, \mathbf{F}^\mathsf{T} (\mathbf{K} - k^2 \mathbf{M})^{-1} \mathbf{F}, \qquad \mathbf{F} = [\mathbf{f}_1 \cdots \mathbf{f}_P], \tag{5}$$

symmetric (reciprocal). On a triangle with area $A$ and the coefficients
$b_i = y_{i+1} - y_{i+2}$, $c_i = x_{i+2} - x_{i+1}$, $K_{ij} = (b_i b_j + c_i c_j)/(4A)$ and
$M_{ij} = A(1 + \delta_{ij})/12$. The nodes on shorted circles are removed.

### 3.2 Modal expansion

$\mathbf{K}$ is symmetric positive semi-definite (singular with the constant vector when no
port is shorted), $\mathbf{M}$ symmetric positive definite. The generalised eigenproblem
$\mathbf{K}\boldsymbol\varphi = \lambda\mathbf{M}\boldsymbol\varphi$ has real eigenvalues
$0 \leq \lambda_1 \leq \lambda_2 \leq \ldots$ and eigenvectors with
$\boldsymbol\varphi_m^\mathsf{T}\mathbf{M}\boldsymbol\varphi_n = \delta_{mn}$, so
$(\mathbf{K} - k^2\mathbf{M})^{-1} = \sum_m \boldsymbol\varphi_m\boldsymbol\varphi_m^\mathsf{T}/(\lambda_m - k^2)$ and

$$Z_{pq} = z \sum_m \frac{\bar\varphi_m(p)\, \bar\varphi_m(q)}{\lambda_m - k^2}, \qquad \bar\varphi_m(p) = \mathbf{f}_p^\mathsf{T}\boldsymbol\varphi_m. \tag{6}$$

The frequency enters only through $z$ and $k^2$: the modes are computed once per pair. The
constant mode $\lambda_1 = 0$, $\bar\varphi_1 = 1/\sqrt{A}$, gives at low frequency
$Z \approx z/(-k^2 A) = 1/(yA) = h/(j\omega\varepsilon A)$, the capacitance of the plates.
Loss makes $k^2$ complex and bounds each resonance at $\lambda_m = \operatorname{Re} k^2$.

### 3.3 Modes kept and the static remainder

Modes with $\lambda_m \leq (3 k_\text{max})^2$ are computed, $k_\text{max}$ the wavenumber at the
highest frequency in the most permittive dielectric. The others enter by their static sum:
with a shift $s > 0$,

$$\mathbf{R} = \mathbf{F}^\mathsf{T}(\mathbf{K} + s\mathbf{M})^{-1}\mathbf{F} - \sum_{m\ \text{kept}} \frac{\bar\varphi_m \bar\varphi_m^\mathsf{T}}{\lambda_m + s} = \sum_{m\ \text{not kept}} \frac{\bar\varphi_m \bar\varphi_m^\mathsf{T}}{\lambda_m + s}, \tag{7}$$

computed exactly from one solve per port, and

$$\mathbf{Z} = z\Big(\sum_{m\ \text{kept}} \frac{\bar\varphi_m \bar\varphi_m^\mathsf{T}}{\lambda_m - k^2} + \mathbf{R}\Big).$$

Each omitted term is used with $\lambda_m + s$ in place of $\lambda_m - k^2$, a relative error of
$|k^2 + s|/|\lambda_m - k^2|$: at most about $1/8$ for the first omitted mode at the top of
the sweep, and falling as $1/\lambda_m$. The remainder holds above all the high modes that
make up the logarithmic field near each port, for which the error is far smaller.

### 3.4 Mesh

The mesh (`docs/methods/dc.md`, section 4.1, for the generator) has edges no longer than a
twentieth of the wavelength at $2k_\text{max}$, a quarter of that along the outline, and is
graded to half the port's radius at each circle. Conforming elements give eigenvalues that
are upper bounds of the true ones (the min-max principle restricted to a subspace); for
linear elements on a uniform one-dimensional mesh,

$$\frac{\lambda_h}{\lambda} = 1 + \frac{(\kappa h)^2}{12} + O\big((\kappa h)^4\big)$$

for a mode of wavenumber $\kappa$ (Strang and Fix 1973), and triangles behave alike to a
factor of order one. At $\kappa = k_\text{max}$, $\kappa h = \pi/20$ and the eigenvalue is about
0.2 % high, the resonance frequency 0.1 %.

## 4. Eigenvalues

### 4.1 Shift-invert Lanczos

With the shift $s = 0.1\,(\pi/\sqrt{A})^2$, below the first non-zero mode of a region of area
$A$, $\mathbf{K} + s\mathbf{M}$ is positive definite. The operator
$\mathbf{T} = (\mathbf{K} + s\mathbf{M})^{-1}\mathbf{M}$ is self-adjoint in the $\mathbf{M}$ inner
product with the eigenvalues $\theta_m = 1/(\lambda_m + s)$: the lowest modes become the
largest and best separated, which Lanczos finds first (Ericsson and Ruhe 1980). From a
random start, the Lanczos recursion in the $\mathbf{M}$ inner product,

$$\mathbf{w} = \mathbf{T}\mathbf{q}_j, \quad \alpha_j = \mathbf{q}_j^\mathsf{T}\mathbf{M}\mathbf{w}, \quad \mathbf{w} \leftarrow \mathbf{w} - \sum_{i \leq j} (\mathbf{q}_i^\mathsf{T}\mathbf{M}\mathbf{w})\, \mathbf{q}_i, \quad \beta_j = \lVert\mathbf{w}\rVert_\mathbf{M}, \quad \mathbf{q}_{j+1} = \mathbf{w}/\beta_j,$$

builds the tridiagonal matrix of the $\alpha_j$ and $\beta_j$; the orthogonalisation against
all previous vectors is done twice, which keeps them orthogonal to rounding. Every ten
steps the tridiagonal matrix is diagonalised by the implicit QL method (Bowdler et al.
1968). A Ritz pair $(\theta, \mathbf{y})$ of the $j$ steps has the residual
$\lVert\mathbf{T}\mathbf{x} - \theta\mathbf{x}\rVert_\mathbf{M} = |\beta_j\, y_j|$, and its eigenvalue is in
error by about the residual squared over the gap to the next (Parlett 1998). The
iteration stops when every Ritz value with $\lambda \leq (3k_\text{max})^2$ has a residual below
$10^{-7}$ of the largest Ritz value and three more lie beyond with residuals below
$10^{-4}$ of it, so that no mode inside the target is missing between converged ones.

### 4.2 Envelope Cholesky

The solves with $\mathbf{K} + s\mathbf{M}$ use its exact Cholesky factor in envelope (profile)
storage after reverse Cuthill–McKee ordering (Cuthill and McKee 1969; George and Liu 1981):
on two-dimensional meshes the envelope is narrow, about the square root of the number of
nodes, and the factor is computed once and reused for every Lanczos step and for (7).

## 5. Branches

### 5.1 Capacitors and the regulator

Each load, capacitor and the regulator joins the power plane at the via of the power net
nearest one of its pads within reach, or at its pin on the plane layer; that via is a port
of every pair whose power plane it reaches. A capacitor is the series impedance
$R + j\omega(L + L_m) + 1/(j\omega C)$ from the port to the planes, the regulator
$R + j\omega(L + L_m)$, with the part's ESR and ESL and the mounting inductance $L_m$.

### 5.2 Mounting inductance

The mounting path leaves the face of the nearer plane, runs up the power via to the pad
layer, along to the power pad, across to the ground pad, to the ground via and down it to
the same face. The plane is taken as perfectly conducting, so the field of the path above
it is that of the path closed by its image below the face, and the inductance of the path
is half that of the closed loop of path and image. The loop's inductance is Neumann's
formula over its straight pieces,

$$L = \sum_i L_{p,i} + \frac{\mu_0}{4\pi} \sum_{i \neq j} \hat{\mathbf{u}}_i \cdot \hat{\mathbf{u}}_j \int_{\ell_i}\int_{\ell_j} \frac{ds\, ds'}{R},$$

the partial self-inductance $L_{p,i}$ of each piece as a tube of its radius
(`docs/methods/networks.md`, section 4.5) and the mutual terms by 8-point Gauss–Legendre
on each half of each piece, with $R$ held at least the pieces' mean radius (Grover 1946).
Vias take their radius, pads and tracks a quarter of the pad's width. Collinear
consecutive pieces of one radius are merged first, so that their mutual term, singular
where they meet, is part of the closed-form self term.

## 6. Assembly

At each frequency of the logarithmic sweep, each pair's impedance matrix (5) over its
ports is inverted, and the admittance matrices of all pairs add at the ports they share
(a via reaching several pairs). The branches' admittances $1/Z_\text{branch}$ add on the
diagonal. The impedance a load sees is the diagonal entry of the inverse at its port.
Resonances are the local maxima of $|Z|$ on the sweep, and the bands above the target are
those of the sampled $|Z|$.

## 7. Error

| source | effect |
| --- | --- |
| model | the fundamental mode only ($h \ll$ wavelength); no radiation at the edges, so resonances of thick pairs on small boards are sharper than measured |
| mesh | eigenvalues upper bounds, about 0.2 % high at $k_\text{max}$; against the closed form of a rectangle, the first ten within $2 \times 10^{-3}$ (`docs/validation/pdn.md`) |
| static remainder | relative error of each omitted mode $\lvert k^2 + s\rvert/\lvert\lambda_m - k^2\rvert$, at most about 1/8 at the top of the sweep |
| Lanczos | residuals below $10^{-7}$ of the scale; eigenvalue error about their square over the gap |
| near resonances | a 0.1 % error in a mode's frequency moves $Z$ there by tens of times as much: 4 % against the double modal sum of a rectangle |
| mounting inductance | ideal plane and thin wires: +1.5 % and +0.8 % against Grover's rectangle with its image |
| branches | one via per branch; spreading between several vias of one part is not resolved |

## 8. Implementation

| step | where |
| --- | --- |
| region, filled holes, pairs | `src/core/src/pdn.cpp`, `analyse_pdn` |
| mesh, $\mathbf{K}$, $\mathbf{M}$, $\mathbf{F}$ (4) | `src/core/src/cavity.cpp`, `Cavity::Cavity` |
| shift-invert Lanczos 4.1, QL | `cavity.cpp`, `Cavity::Cavity`, `tridiagonal_eigen` |
| envelope Cholesky 4.2 | `cavity.cpp`, `EnvelopeCholesky` |
| static remainder (7) | `cavity.cpp`, end of `Cavity::Cavity` |
| impedance (5), (6) at a frequency, $z$, $k^2$ | `cavity.cpp`, `Cavity::impedance` |
| surface impedance | `pdn.cpp`, `surface_impedance` |
| mounting inductance 5.2 | `pdn.cpp`, `image_loop_inductance` |
| assembly 6 | `pdn.cpp`, `analyse_pdn`, the frequency loop |

## References

- Bowdler, H., Martin, R. S., Reinsch, C., and Wilkinson, J. H. (1968). The QR and QL
  algorithms for symmetric matrices. *Numerische Mathematik* 11, 293–306.
  doi:10.1007/BF02166681.
- Cuthill, E., and McKee, J. (1969). Reducing the bandwidth of sparse symmetric matrices.
  *Proceedings of the 24th National Conference of the ACM*, 157–172.
  doi:10.1145/800195.805928.
- Ericsson, T., and Ruhe, A. (1980). The spectral transformation Lanczos method for the
  numerical solution of large sparse generalized symmetric eigenvalue problems.
  *Mathematics of Computation* 35(152), 1251–1268. doi:10.1090/S0025-5718-1980-0583502-2.
- George, A., and Liu, J. W. H. (1981). *Computer Solution of Large Sparse Positive
  Definite Systems*. Prentice-Hall.
- Grover, F. W. (1946). *Inductance Calculations: Working Formulas and Tables*. Van Nostrand.
- Parlett, B. N. (1998). *The Symmetric Eigenvalue Problem*. SIAM Classics in Applied
  Mathematics 20. doi:10.1137/1.9781611971163.
- Strang, G., and Fix, G. J. (1973). *An Analysis of the Finite Element Method*.
  Prentice-Hall.
