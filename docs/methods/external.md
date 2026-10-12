# External reference solvers: method

KairoQ exports four problems to open-source solvers that model more of the physics than
the native analyses, so that the native results can be checked against a solution of
the full problem on the same geometry: EMI loops and channels to openEMS (full-wave,
finite-difference time-domain), DC current and heat to Elmer (three-dimensional finite
elements). This document states the problem each export poses, the quantities KairoQ
derives from the solvers' output, and where the reference itself is approximate. The user
guides are `docs/solvers/openems.md`, `docs/solvers/rf.md` and `docs/solvers/elmer.md`.

## 1. openEMS

### 1.1 The solver

openEMS (Liebig et al. 2013) integrates Maxwell's curl equations on a rectilinear Yee grid
(Yee 1966): electric and magnetic field components staggered by half a cell in space and
half a step in time, second-order accurate in both, with the time step at the Courant
limit of the smallest cell. Conductors are perfect conductors on the grid's edges or, for
the channel export, conducting sheets of a layer's conductivity and thickness (openEMS's
model of a thin lossy conductor). The domain ends in eight cells of perfectly matched
layer. A Gaussian pulse covers the band from 0 to $f_\text{max}$; the run stops when the
energy in the domain has fallen by the end criterion. The frequency responses are the
discrete Fourier transforms of the recorded port voltages and currents.

**Errors of the reference.** The grid's dispersion error is second order in the cell
against the wavelength (cells of at most a twentieth of the wavelength in the densest
dielectric: a few tenths of a percent in phase velocity; Taflove and Hagness 2005). Circles (barrels,
pads, antipads) are staircased: grid lines trace them at four angles, which holds the
staircase within a few percent of the radius. Copper edges carry lines at the edge and one
minimum cell to either side, which places the edge singularity of the current on the
grid. The truncation by the end criterion leaves a residual of about the criterion's
amplitude in the transforms.

### 1.2 Geometry

Copper is a set of sheets: in the loop export perfect conductors at the layers'
mid-planes; in the channel export outer layers on their face towards the board, inner
layers at mid-plane. CSXCAD polygons have no holes, so a polygon with holes is cut into
strips by lines through its holes. Plated holes are cylinders between their outermost
layers. The model's extent is a window or margin around the copper of interest, clipped to
the board, with air of a quarter wavelength around it.

### 1.3 Dielectrics

The cross-sections use the Djordjevic–Sarkar model (`docs/methods/si.md`, section 4).
openEMS's dispersive materials are Lorentz and Debye media; the channel export fits

$$\varepsilon(f) \approx \varepsilon_\infty + \sum_{k=1}^{6} \frac{\Delta\varepsilon_k}{1 + j 2\pi f \tau_k} \tag{1}$$

to it, with the relaxation times $\tau_k = 1/(2\pi f_k)$ at six frequencies spaced
logarithmically from $f_\text{max}/600$ to $3 f_\text{max}$, and $\varepsilon_\infty$ and the
$\Delta\varepsilon_k$ from a least-squares fit of the real and imaginary parts at 80
logarithmically spaced frequencies from $f_\text{max}/200$ to $f_\text{max}$, under the constraint
that all coefficients are non-negative (Lawson and Hanson 1995). Non-negative
$\Delta\varepsilon_k$ make every term passive ($\varepsilon'' \geq 0$ at every frequency) and causal,
so the fitted medium is too. The fit reproduces the model within $2 \times 10^{-4}$
(`docs/validation/rf.md`). As openEMS's Python interface has no Debye material, the script
writes the structure to its XML form, replaces the dielectrics there, and runs openEMS on
the file.

With `--losses simple` (and always in the loop export) a dielectric is instead
non-dispersive with the conductivity $\kappa = 2\pi f_c\, \varepsilon_0 \varepsilon_r \tan\delta$
of its loss at the band's centre $f_c$, exact there only.

### 1.4 Ports and S-parameters

A **lumped port** is a resistance $Z_0$ in series with a source across a gap or a vertical
box, with its voltage along the box and its current from a loop around it. A **microstrip
port** feeds the first part of a straight track at the pad and measures voltage and
current halfway along; openEMS gives the line's measured propagation constant
$\gamma = j\beta$ ($\beta$ complex, its imaginary part the loss) and impedance $Z_L$. The
voltage and current are moved back to the pad over the distance $d$ by the transfer matrix
of the line (`docs/methods/networks.md`, equation 2):

$$V_\text{pad} = V\cosh\gamma d + Z_L I \sinh\gamma d, \qquad I_\text{pad} = I\cosh\gamma d + \frac{V}{Z_L}\sinh\gamma d,$$

with the attenuation, which openEMS's own shift of the reference plane leaves out.

At each port's reference plane, $a = \tfrac12(V + Z_0 I)$ and $b = \tfrac12(V - Z_0 I)$; with
port $j$ excited, $S_{ij} = b_i/a_j$, the ratio of power waves (the common factor
$1/\sqrt{Z_0}$ cancels). Columns of ports not excited are taken from reciprocity,
$S_{ij} = S_{ji}$, where the rows are known.

### 1.5 Loop impedance and inductance

The loop export excites one lumped port in a gap of the loop's first element; its impedance
is $Z(f) = V/I$ and its inductance $L(f) = \operatorname{Im} Z / (2\pi f)$. Below the loop's
self-resonance $L(f)$ is flat and compares with the native estimate's high-frequency limit
(`docs/methods/emi.md`, section 5).

### 1.6 Far field

The near-to-far-field transformation integrates the equivalent surface currents on a box
inside the boundary. For the loop, the largest field magnitude $E$ over a sphere of radius
$r$ per ampere of port current gives the equivalent moment: a small loop of moment $m$ has
the far field $E = \eta_0 k^2 m \sin\theta / (4\pi r)$ (Balanis 2016, ch. 5), so

$$m_\text{eq} = \frac{4\pi r E}{\eta_0 k^2 I},$$

the moment a small loop would need for the same largest field, to compare with the
estimate's moment (`docs/methods/emi.md`, section 6). For a channel's port 1 the figures of
merit are, with the accepted power $P_\text{acc} = \tfrac12\operatorname{Re}(V I^*)$, the incident
power $P_\text{inc} = |a|^2/(2Z_0)$ and the radiated power $P_\text{rad}$ from the box:

$$\eta = \frac{P_\text{rad}}{P_\text{acc}}, \qquad G = \eta D, \qquad G_\text{realised} = D\frac{P_\text{rad}}{P_\text{inc}} = G\,(1 - |\Gamma|^2),$$

$D$ the directivity from the transformation.

## 2. Elmer

### 2.1 The solver

Elmer (CSC – IT Center for Science) solves the exported problems by Galerkin finite
elements on the exported mesh: linear six-node wedges, all layers sharing one
triangulation, so the volumes join conformally. The meshes and the reasons for their
shape are in `docs/solvers/elmer.md`.

### 2.2 DC current

The potential satisfies $\nabla \cdot (\boldsymbol\sigma \nabla\phi) = 0$ in the copper, with
the copper's conductivity at the DC temperature, the outer faces of the copper inside each
port's area as Dirichlet boundaries, and no current through every other surface. A plated
hole is a column of wedges over its span with an anisotropic conductivity: along the hole,
the plating's conductance spread over the hole's section,

$$\sigma_z = \sigma_p \frac{\pi t_p (d + t_p)}{\pi d^2 / 4},$$

the annulus of thickness $t_p$ around the drill $d$ over the drill's area, so that the
column has the plating's resistance per length; across it, copper, as the wall conducts
around the hole.

For $n$ ports Elmer solves $n$ cases, port $k$ at 1 V and the others at 0 V, by conjugate
gradients (Hestenes and Stiefel 1952) with an incomplete LU preconditioner to a relative
residual of $10^{-12}$. The **nodal loads**, the residual $\mathbf{r} = \mathbf{K}\boldsymbol\phi - \mathbf{f}$
of the discrete equations at the Dirichlet nodes, are the discrete currents leaving
through them: summed over a port's nodes, they are the current into the port, exactly
consistent with the discrete solution (it is the discrete form of
$\oint \boldsymbol\sigma\nabla\phi \cdot \mathbf{n}\, dA$, without differentiating the
field). Case $k$ gives column $k$ of the conductance matrix $\mathbf{G}$ of the
equipotential ports. The project's conditions then give the port voltages from
$\mathbf{I} = \mathbf{G}\mathbf{V}$: set voltages are moved to the right side, set currents
remain unknowns, and a group of ports connected through $\mathbf{G}$ without a set voltage
takes its first port as 0 V, as in the native solver (`docs/methods/dc.md`, section 2.4);
the system is solved by Gaussian elimination with partial pivoting.

### 2.3 Heat

The temperature satisfies $\nabla \cdot (\mathbf{k}\nabla T) + q = 0$ through the whole
board: copper, the dielectrics with their in-plane and through-thickness conductivities,
the gaps of the copper layers filled with the dielectric below, mask or nothing (air), and
plated holes as columns built as in 2.2 with copper's thermal conductivity in place of
the electrical one. Sources heat their copper uniformly by volume at
the set power. The top and bottom faces lose heat by convection through the solder mask
as a series resistance,

$$-\mathbf{k}\nabla T \cdot \mathbf{n} = h_\text{eff}(T - T_\text{amb}), \qquad h_\text{eff} = \frac{1}{1/h + t_\text{mask}/k_\text{mask}},$$

as in the layered model (`docs/methods/thermal.md`); the sides, including the cut edges of
a crop, are adiabatic. Radiation is left out (set `emissivity: 0` in the native model to
compare). The linear system is solved by IDR(s) (Sonneveld and van Gijzen 2008) with an
ILU(1) preconditioner to a relative residual of $10^{-12}$: with conductivities a thousand
times apart and thin copper the system is badly conditioned, and on these meshes IDR(s)
needs about a fifteenth of the time of conjugate gradients.

### 2.4 Errors of the reference

For smooth solutions, linear elements give the error in the energy norm first order in
the element size and the conductance, an energy, second order; the conductance computed from the nodal loads
is that of the discrete solution, an upper bound of the exported model's conductance (the
discrete potential minimises the dissipation over a subset of the admissible fields, so
its dissipation at 1 V is larger; the same argument as the lower bound on resistance in
`docs/methods/dc.md`, section 6.1). The element sizes are set by hand; there is no adaptive refinement, so the
reference's own error is judged by refining it (`docs/validation/dc.md`,
`docs/validation/thermal.md`).

## 3. Implementation

| step | where |
| --- | --- |
| loop model, script | `src/core/src/openems.cpp`; `src/core/resources/openems/simulate.py` |
| channel model, script | `openems.cpp`; `src/core/resources/openems/simulate_ports.py` |
| Debye fit (1) | `simulate_ports.py`, `djordjevic_sarkar`, `nnls`, `debye_poles`, `run_debye` |
| microstrip port moved to the pad | `simulate_ports.py`, `move_to_pad` |
| S-parameters, reciprocity, far field | `simulate_ports.py`, `main` |
| loop impedance, moment | `simulate.py`, `main` |
| DC mesh, conductivities, cases | `src/core/src/elmer.cpp`, DC export |
| conductance matrix, port conditions | `elmer.cpp`, reading of the results |
| thermal mesh, materials, boundaries | `elmer.cpp`, thermal export |

## References

- Balanis, C. A. (2016). *Antenna Theory: Analysis and Design*, 4th ed. Wiley.
- CSC – IT Center for Science. *ElmerSolver Manual* and *Elmer Models Manual*.
  https://www.elmerfem.org.
- Hestenes, M. R., and Stiefel, E. (1952). Methods of conjugate gradients for solving
  linear systems. *Journal of Research of the National Bureau of Standards* 49(6),
  409–436. doi:10.6028/jres.049.044.
- Lawson, C. L., and Hanson, R. J. (1995). *Solving Least Squares Problems*. SIAM Classics
  in Applied Mathematics 15. doi:10.1137/1.9781611971217.
- Liebig, T., Rennings, A., Held, S., and Erni, D. (2013). openEMS – a free and open
  source equivalent-circuit (EC) FDTD simulation platform supporting cylindrical
  coordinates suitable for the analysis of traveling wave MRI applications.
  *International Journal of Numerical Modelling* 26(6), 680–696. doi:10.1002/jnm.1875.
- Sonneveld, P., and van Gijzen, M. B. (2008). IDR(s): a family of simple and fast
  algorithms for solving large nonsymmetric systems of linear equations. *SIAM Journal on
  Scientific Computing* 31(2), 1035–1062. doi:10.1137/070685804.
- Taflove, A., and Hagness, S. C. (2005). *Computational Electrodynamics: The
  Finite-Difference Time-Domain Method*, 3rd ed. Artech House.
- Yee, K. S. (1966). Numerical solution of initial boundary value problems involving
  Maxwell's equations in isotropic media. *IEEE Transactions on Antennas and Propagation*
  14(3), 302–307. doi:10.1109/TAP.1966.1138693.
