# Networks: method

This document derives the network analysis of a channel: the circuit that stands for its
copper, the transmission-line sections and the via model it is built from, the nodal
solution, the S-parameters and mixed-mode parameters, and the time-domain responses. The
user guide is `docs/solvers/networks.md`; the validation, `docs/validation/networks.md`. The
per-unit-length parameters of the lines are those of `docs/methods/si.md`.

## 1. Symbols

| symbol | meaning | unit |
| --- | --- | --- |
| $\omega = 2\pi f$ | angular frequency | rad/s |
| $\mathbf{Z}', \mathbf{Y}'$ | per-unit-length series impedance and shunt admittance matrices | Ω/m, S/m |
| $\ell$ | length of a line section | m |
| $\mathbf{T}$ | transfer (ABCD) matrix of a section | |
| $\mathbf{Y}$ | nodal admittance matrix of the channel | S |
| $\mathbf{Z}$ | impedance matrix of the ports | Ω |
| $Z_0$ | reference impedance of the single-ended ports | Ω |
| $h$ | separation of a pair of planes | m |
| $a$ | radius of a via barrel | m |
| $k$ | wavenumber in the dielectric of a pair, $\omega \sqrt{\varepsilon_r} / c_0$ | 1/m |
| $H_n^{(2)}, K_0, I_0$ | Hankel function of the second kind, modified Bessel functions | |

The time dependence is $e^{j\omega t}$.

## 2. The channel as a circuit

All copper of the reference nets is one node, the reference of the circuit (ground). The
copper of the channel's nets is reduced to nodes: track ends (ends within 1 µm on a layer
are one node), pads, the ends of the samples of each track, and the nodes of each via. The
elements between them are line sections (section 3), vias (section 4), and parts and loads
(section 5). Every element is linear, so at each frequency it is described by the currents
it draws from its nodes as a linear function of their voltages, and its contributions add
into the nodal admittance matrix $\mathbf{Y}$ (section 6).

The model is quasi-TEM throughout: the lines carry their transverse modes only, and the
fields of different elements interact only through shared nodes, through the coupled
sections of section 3.4, and through the plane pairs of section 4.

## 3. Line sections

### 3.1 The telegrapher's equations

On a uniform section of $N$ conductors over the reference, the vectors of voltage
$\mathbf{V}(z)$ and current $\mathbf{I}(z)$ (in the direction of $z$) satisfy

$$\frac{d\mathbf{V}}{dz} = -\mathbf{Z}' \mathbf{I}, \qquad \frac{d\mathbf{I}}{dz} = -\mathbf{Y}' \mathbf{V}, \tag{1}$$

with $\mathbf{Z}' = \mathbf{R} + j\omega\mathbf{L}$ and $\mathbf{Y}' = \mathbf{G} + j\omega\mathbf{C}$
from the cross-section at the frequency (Paul 2008, ch. 3). Writing
$\mathbf{M} = \begin{pmatrix} 0 & \mathbf{Z}' \\ \mathbf{Y}' & 0 \end{pmatrix}$, (1) is
$\frac{d}{dz}(\mathbf{V}, \mathbf{I}) = -\mathbf{M} (\mathbf{V}, \mathbf{I})$, whose solution
from the far end $z = \ell$ back to the near end is

$$\begin{pmatrix} \mathbf{V}(0) \\ \mathbf{I}(0) \end{pmatrix} = \mathbf{T} \begin{pmatrix} \mathbf{V}(\ell) \\ \mathbf{I}(\ell) \end{pmatrix}, \qquad \mathbf{T} = \exp(\ell \mathbf{M}) = \begin{pmatrix} \mathbf{A} & \mathbf{B} \\ \mathbf{C} & \mathbf{D} \end{pmatrix}. \tag{2}$$

### 3.2 One conductor

For $N = 1$, $\mathbf{M}^2 = Z'Y'\,\mathbf{I}_2 = \gamma^2 \mathbf{I}_2$, so the even and odd powers of
the exponential series sum separately:

$$\mathbf{T} = \cosh(\gamma\ell)\, \mathbf{I}_2 + \frac{\sinh(\gamma\ell)}{\gamma}\, \mathbf{M} = \begin{pmatrix} \cosh\gamma\ell & Z_c \sinh\gamma\ell \\ \sinh(\gamma\ell)/Z_c & \cosh\gamma\ell \end{pmatrix},$$

with $\gamma = \sqrt{Z'Y'}$ and $Z_c = \sqrt{Z'/Y'}$, the principal roots ($\operatorname{Re}\gamma \geq 0$
for a passive line). This form is evaluated directly.

### 3.3 Several conductors: scaling and squaring

For $N \geq 2$, the exponential (2) is computed by scaling and squaring with the diagonal
Padé approximant of degree $(6, 6)$ (Moler and Van Loan 2003; Golub and Van Loan 2013,
§9.3). With $s = \max(0, \lceil \log_2 (2 \lVert \ell\mathbf{M} \rVert_\infty) \rceil)$ and
$\mathbf{X} = \ell\mathbf{M} / 2^s$, so that $\lVert\mathbf{X}\rVert_\infty \leq \tfrac12$,

$$R_{66}(\mathbf{X}) = Q(-\mathbf{X})^{-1} Q(\mathbf{X}), \qquad Q(\mathbf{X}) = \sum_{j=0}^{6} c_j \mathbf{X}^j, \qquad c_j = \frac{(12 - j)!\, 6!}{12!\, j!\, (6 - j)!},$$

that is $c = (1, \tfrac12, \tfrac{5}{44}, \tfrac{1}{66}, \tfrac{1}{792}, \tfrac{1}{15840}, \tfrac{1}{665280})$,
and $\mathbf{T} = R_{66}(\mathbf{X})^{2^s}$ by $s$ squarings. For $\lVert\mathbf{X}\rVert_\infty \leq \tfrac12$,
$R_{pq}(\mathbf{X}) = e^{\mathbf{X} + \mathbf{E}}$ with $\mathbf{E}$ commuting with $\mathbf{X}$ and

$$\frac{\lVert \mathbf{E} \rVert_\infty}{\lVert \mathbf{X} \rVert_\infty} \leq \varepsilon_{pq} = 2^{3 - (p+q)} \frac{p!\, q!}{(p+q)!\, (p+q+1)!}, \qquad \varepsilon_{66} \approx 3.4 \times 10^{-16}$$

(Moler and Van Loan 2003, §3; Golub and Van Loan 2013, §9.3). Since $\mathbf{E}$ commutes with $\mathbf{X}$, the
squarings give exactly $e^{\ell\mathbf{M} + 2^s\mathbf{E}}$: the computed $\mathbf{T}$ is the
exact exponential of a matrix within a relative distance $\varepsilon_{66}$ of $\ell\mathbf{M}$, at
the level of rounding. Rounding in the squarings themselves is not covered by the bound;
as the sections are short (a sample of a track, a fraction of a wavelength at the highest
frequency), $s$ is small.

### 3.4 Samples and coupled sections

Each track is cut at intervals as in the signal-integrity analysis
(`docs/methods/si.md`, section 9). Sample $i$ of a track of length $L$ with $n$ samples
stands for a uniform section of length $L/n$ with the parameters of its cut. A non-uniform
line is thereby replaced by a cascade of uniform sections; the difference is first order
in the change of the parameters over a section.

Where the cut of a sample holds copper of other nets of the channel, the samples of those
nets nearest to that copper (on the same layer, within 0.6 times the longer sample's
length plus 0.1 mm) form with it one section of several conductors, of the common length
$\ell_c$, the shortest of their lengths (their mean where they differ by less than
10 µm). The rest of each longer sample is a single line in series. A member running against
the sample's direction is connected with its ends exchanged.

A track without any cut with a reference is a series inductance of 1 nH/mm; a sample
without a usable cut takes the cross-section of its nearest neighbour on the track.

### 3.5 Admittance of a section

With $\mathbf{I}_1 = \mathbf{I}(0)$ into the section at its near end and
$\mathbf{I}_2 = -\mathbf{I}(\ell)$ into it at its far end, (2) gives
$\mathbf{V}_1 = \mathbf{A}\mathbf{V}_2 - \mathbf{B}\mathbf{I}_2$, so
$\mathbf{I}_2 = \mathbf{B}^{-1}(\mathbf{A}\mathbf{V}_2 - \mathbf{V}_1)$, and
$\mathbf{I}_1 = \mathbf{C}\mathbf{V}_2 - \mathbf{D}\mathbf{I}_2$:

$$\begin{pmatrix} \mathbf{I}_1 \\ \mathbf{I}_2 \end{pmatrix} = \begin{pmatrix} \mathbf{D}\mathbf{B}^{-1} & \mathbf{C} - \mathbf{D}\mathbf{B}^{-1}\mathbf{A} \\ -\mathbf{B}^{-1} & \mathbf{B}^{-1}\mathbf{A} \end{pmatrix} \begin{pmatrix} \mathbf{V}_1 \\ \mathbf{V}_2 \end{pmatrix}. \tag{3}$$

The four blocks are added into $\mathbf{Y}$ at the rows and columns of the section's near
and far nodes. $\mathbf{B} = \ell\mathbf{Z}' + O(\ell^3)$ is invertible for any section of
non-zero length with $\mathbf{Z}'$ invertible.

## 4. Vias

### 4.1 Structure

A plated hole of a channel net has a node on each copper layer of its barrel. Where
reference planes surround it, the barrel is described by the physics-based via model
(Rimolo-Donadio et al. 2009): the region between two consecutive reference planes is a
parallel-plate waveguide, the barrel a current source in it, and the barrel's charge a
capacitance to the planes. The barrel's segments, between nodes at each copper layer and
at slices that keep each segment below a twentieth of the shortest wavelength in the stack,
carry the currents of the ladder.

### 4.2 Fields of a barrel between two planes

Between perfectly conducting planes at $z = 0$ and $z = h$, filled with a dielectric of
wavenumber $k$, the vector potential of a vertical current on the barrel's surface
satisfies the Neumann condition $\partial A_z / \partial z = 0$ at the planes (the image of a
vertical current in a conducting plane has the same direction). With
$\delta(z - z') = \tfrac1h + \tfrac2h \sum_{n \geq 1} \cos\frac{n\pi z}{h} \cos\frac{n\pi z'}{h}$,
the field separates into the modes $\cos(n\pi z/h)$:

- $n = 0$, the fundamental (TEM) mode: uniform in $z$, a radial wave with the
  two-dimensional Green's function $-\tfrac{j}{4} H_0^{(2)}(k\rho)$;
- $n \geq 1$: radial dependence $K_0(g_n\rho)$ with $g_n^2 = (n\pi/h)^2 - k^2$; for the
  sweeps considered $n\pi/h \gg k$, and $g_n = n\pi/h$ is used, which is the quasi-static
  limit. These modes decay within a few $h$ of the barrel.

The fundamental mode carries the coupling between barrels and to the board; the higher
modes store magnetic energy near each barrel and give its own inductance.

### 4.3 The fundamental mode

A barrel of radius $a$ carrying the total current $I$, uniform in $z$, excites
$E_z = A\, H_0^{(2)}(k\rho)$ for $\rho \geq a$. From
$H_\phi = \frac{1}{j\omega\mu}\frac{\partial E_z}{\partial\rho} = -\frac{A k}{j\omega\mu} H_1^{(2)}(k\rho)$
and $2\pi a H_\phi(a) = I$,

$$A = -\frac{j\omega\mu I}{2\pi k a\, H_1^{(2)}(ka)},$$

and the voltage between the planes at the barrel and at a distance $\rho$ are

$$Z_{ii} = \frac{j\omega\mu h}{2\pi k a} \frac{H_0^{(2)}(ka)}{H_1^{(2)}(ka)}, \qquad Z_{ij} = \frac{\omega\mu h}{4} H_0^{(2)}(k\rho_{ij}) \tag{4}$$

(the second using $H_1^{(2)}(ka) \approx \tfrac{2j}{\pi ka}$ for $ka \ll 1$, so that the source is
a line current; Rimolo-Donadio et al. 2009). For $k\rho \to 0$ both tend to
$\frac{j\omega\mu h}{2\pi}\ln\frac{1}{k\rho}$ plus constants, an inductance, as they must; at
larger $k\rho$ the outward wave carries power, $\operatorname{Re} Z > 0$. The distance is
taken no smaller than $a_i + a_j$.

The dielectrics between the planes are in series for the uniform $E_z$ of this mode
($D_z$ is continuous), so the pair takes the permittivity
$\varepsilon = t / \sum_k (t_k / \varepsilon_k(f))$, with each layer's complex wideband
permittivity (`docs/methods/si.md`, section 4).

**Infinite planes.** On a pair whose outline holds more than 300 modes up to three times the
sweep's wavenumber, the planes are taken as infinite. The barrels of the reference nets
that join the two planes short the pair where they stand: the eight nearest each channel
barrel within $\max(5\ \text{mm}, 10h)$ are included as ports of zero voltage. With the
impedance matrix (4) over all ports, ordered open then shorted, and its inverse
$\mathbf{Y}_\text{all} = \mathbf{Z}_\text{all}^{-1}$, setting the shorted voltages to zero leaves

$$\mathbf{Z}_v = \big( (\mathbf{Y}_\text{all})_{oo} \big)^{-1}, \tag{5}$$

the impedance between the open barrels with the shorts in place.

**Finite pairs.** Otherwise the impedance between the barrels is that of the pair's real
outline, a cavity with open edges, from its modes with the joining reference barrels as
shorted ports (`docs/methods/pdn.md`).

**Height weights.** A segment of the barrel covering the height $\ell_i$ of a pair carries
its share of the pair's mode: the mode's voltage is driven by
$\sum_j w_j I_j$ and appears across segment $i$ as the fraction $w_i$, with $w_i = \ell_i / h$.
This is the $n = 0$ term of the expansion of section 4.2 tested on the segments: the
coupling is $w_i (\mathbf{Z}_v)_{v(i) v(j)} w_j$, symmetric.

### 4.4 The higher modes: inductance of a barrel's segments

For $n \geq 1$, the potential of the surface current of segment $j$ (a ring current density
on the cylinder $\rho = a$), averaged over the circumference of segment $i$ on the same
barrel, is given by the ring-to-ring kernel of the modified Helmholtz operator: by the
addition theorem, the mean of $K_0(g|\mathbf{r} - \mathbf{r}'|)$ over two circles of radius $a$
about one axis is $K_0(ga) I_0(ga)$. With the expansion of $\delta(z - z')$ and the
two-dimensional Green's function $K_0(g\rho)/(2\pi)$ of $-\nabla^2 + g^2$,

$$L_{ij} = \frac{\mu}{\pi h} \sum_{n \geq 1} c_{ni}\, c_{nj}\, K_0(g_n a)\, I_0(g_n a), \qquad c_{ni} = \int_{\text{segment } i} \cos\frac{n\pi z}{h}\, dz, \tag{6}$$

taken between segments of one barrel in one pair. The series is summed to $n = 4000$.
With $|c_{ni}| \leq 2/g_n$ and $K_0(x)I_0(x) \leq 1.07/(2x)$ (the maximum of
$2x K_0(x) I_0(x)$ is 1.067, near $x = 1.08$), the terms beyond $N$ are bounded by

$$\frac{\mu}{\pi h} \sum_{n > N} \frac{4}{g_n^2} \frac{1.07}{2 g_n a} \leq \frac{1.07\, \mu h^2}{\pi^4 a N^2},$$

which for $N = 4000$, $h \leq 2$ mm and $h/a \leq 100$ is below $10^{-15}$ H, negligible
against the segments' inductance. $K_0 I_0$ is taken from the standard library below $x = 20$ and from
its asymptotic series
$K_0(x)I_0(x) \approx \frac{1}{2x}\big(1 + \frac{1}{8x^2} + \frac{27}{128x^4}\big)$ above.

### 4.5 Segments outside the pairs

A segment above the top plane or below the bottom one (a stub towards an outer layer)
takes the partial self-inductance of a straight tube of length $\ell$ and radius $a$,

$$L = \frac{\mu_0}{2\pi}\Big(\ell \operatorname{arsinh}\frac{\ell}{a} - \sqrt{\ell^2 + a^2} + a\Big),$$

plus its mutual inductance with its image in the nearest plane face, a collinear segment
of equal length at the gap $s = 2d$ ($d$ the distance of its nearer end to the face; the
image of a vertical current has the same direction):

$$M = \frac{\mu_0}{4\pi}\big(F(2\ell + s) - 2F(\ell + s) + F(s)\big), \qquad F(x) = x \ln x,$$

both from Neumann's formula (Grover 1946).

### 4.6 Capacitance of the nodes

The barrel is cut midway between its nodes; each piece, with the pad of its layer, is a
conductor of an axisymmetric electrostatic problem with the planes and their antipads,
all grounded, in the stack's dielectrics. It is solved by the cross-section solver in its
axisymmetric mode (`docs/methods/si.md`, sections 2 and 3), which gives the Maxwell matrix
$C_{ij}$ of the pieces.

The pieces touch, so the field between neighbours is singular where they meet and their
mutual capacitances grow without bound as the mesh refines. The row sum
$C_i = \sum_j C_{ij}$, the charge on piece $i$ when the whole barrel is at 1 V, is the
capacitance of piece $i$ to the planes, and converges. Each node takes $C_i$ to ground, and
the mutual terms are dropped. This assumes neighbouring nodes at nearly equal potential,
which holds while the segments are short against the wavelength. The capacitive part of
the via is the role of the higher modes in Zhang et al. (2008), here computed numerically
for the actual pads and antipads.

### 4.7 Assembly

Per pair of planes, the branches are the segments of the channel's barrels in it; their
impedance matrix is

$$(\mathbf{Z}_b)_{ij} = w_i\, (\mathbf{Z}_v)_{v(i)\, v(j)}\, w_j + j\omega L_{ij},$$

with $v(i)$ the barrel of branch $i$ and $L_{ij}$ from (6) for branches of one barrel, zero
otherwise. $\mathbf{Y}_b = \mathbf{Z}_b^{-1}$ gives the current of each branch from the voltages
across all of them; branch $i$ draws $(\mathbf{Y}_b\mathbf{v})_i$ from its upper node and
returns it to its lower one, with $v_j$ the voltage from upper to lower node of branch $j$.
The node capacitances (4.6) and the outside segments (4.5) are added as shunt and series
elements. A change of return plane is not a separate element: it is the voltage that the
barrel's current leaves in the pairs between the two layers.

### 4.8 Coaxial model

Without planes around the barrel, or on request, each stretch of barrel between
consecutive copper layers is a coaxial line of impedance and propagation constant

$$Z_v = \frac{\eta_0}{2\pi\sqrt{\varepsilon_r}}\ln\frac{D}{d} = \frac{60\ \Omega}{\sqrt{\varepsilon_r}}\ln\frac{D}{d}, \qquad \gamma = \frac{j\omega\sqrt{\varepsilon_r}}{c_0},$$

with $d$ the drill, $D$ the mean equal-area diameter of the antipads in the reference copper
it passes (three drills where there is none), and $\varepsilon_r$ the thickness-weighted mean
of the dielectrics between the two layers (Pozar 2012, §2.2, for the coaxial line). The
barrel's end beyond its last connected layer is an open stub. The sections are stamped by
(3) with the closed form of 3.2.

### 4.9 Special functions

$H_0^{(2)}$ and $H_1^{(2)}$ are $J_n - jY_n$ from their power series for $|z| < 12$
(Abramowitz and Stegun 1964, 9.1.10, 9.1.11, 9.1.13), summed until the terms fall below
$10^{-17}$ of the sums, and from the Hankel asymptotic expansion (9.2.7–9.2.8) beyond, summed
to its smallest term. At $|z| = 12$ the largest term of the series is about $4 \times 10^3$
and the result about $5 \times 10^{-2}$, so about $10^{-11}$ relative accuracy remains.

## 5. Parts and loads

A part given by values is the series impedance $R + j\omega L + 1/(j\omega C)$ between its
two pads; a load is the admittance $1/R + j\omega C$ from a pad to ground. A two-port
Touchstone part with reference resistance $R_0$ enters through
$\mathbf{Y}_p = \frac{1}{R_0}(\mathbf{I} - \mathbf{S}_p)(\mathbf{I} + \mathbf{S}_p)^{-1}$, which
follows from the definition of the scattering matrix (section 6.2) with $\mathbf{V} = \mathbf{Z}\mathbf{I}$
solved for $\mathbf{Y} = \mathbf{Z}^{-1}$; its ports are referred to ground.

## 6. Solution

### 6.1 Nodal analysis

With a branch admittance $y$ between nodes $a$ and $b$ added as $+y$ at $(a, a)$ and
$(b, b)$ and $-y$ at $(a, b)$ and $(b, a)$ (omitting ground), and the multi-port elements by
their blocks, Kirchhoff's current law at every node reads $\mathbf{Y}\mathbf{V} = \mathbf{J}$
with $\mathbf{J}$ the injected currents. Injecting a unit current at port $q$ and none
elsewhere gives the voltages at the port nodes as column $q$ of the port impedance matrix:

$$\mathbf{Z} = \mathbf{P}\, \mathbf{Y}^{-1} \mathbf{P}^\mathsf{T},$$

with $\mathbf{P}$ selecting the port nodes. $\mathbf{Y}$ is dense-factorised by Gaussian
elimination with partial pivoting, for all ports at once; the frequencies are independent
and run in parallel.

### 6.2 Scattering parameters

With the power waves $\mathbf{a} = (\mathbf{V} + Z_0\mathbf{I})/(2\sqrt{Z_0})$ and
$\mathbf{b} = (\mathbf{V} - Z_0\mathbf{I})/(2\sqrt{Z_0})$ of ports with the real reference
$Z_0$ (Kurokawa 1965), and $\mathbf{V} = \mathbf{Z}\mathbf{I}$,

$$\mathbf{S} = (\mathbf{Z} - Z_0\mathbf{I})(\mathbf{Z} + Z_0\mathbf{I})^{-1} = (\mathbf{Z}_n - \mathbf{I})(\mathbf{Z}_n + \mathbf{I})^{-1}, \qquad \mathbf{Z}_n = \mathbf{Z}/Z_0. \tag{7}$$

**Reciprocity.** Every element above except a non-reciprocal Touchstone part has a
symmetric admittance, so $\mathbf{Y}$, $\mathbf{Z}$ and $\mathbf{S}$ are symmetric.

**Passivity.** If every element is passive, $\operatorname{Re}(\mathbf{x}^\mathsf{H}\mathbf{Y}\mathbf{x}) \geq 0$
for all $\mathbf{x}$ (the power it absorbs). For port currents $\mathbf{i}$ and
$\mathbf{y} = \mathbf{Y}^{-1}\mathbf{P}^\mathsf{T}\mathbf{i}$,
$\mathbf{i}^\mathsf{H}\mathbf{Z}\mathbf{i} = (\mathbf{Y}\mathbf{y})^\mathsf{H}\mathbf{y}$, whose real
part is $\operatorname{Re}(\mathbf{y}^\mathsf{H}\mathbf{Y}\mathbf{y}) \geq 0$. Then from (7),

$$\mathbf{I} - \mathbf{S}^\mathsf{H}\mathbf{S} = 2\,(\mathbf{Z}_n + \mathbf{I})^{-\mathsf{H}}\, (\mathbf{Z}_n + \mathbf{Z}_n^\mathsf{H})\, (\mathbf{Z}_n + \mathbf{I})^{-1} \succeq 0,$$

so no singular value of $\mathbf{S}$ exceeds one. A computed $|S_{ij}| > 1$ therefore points
to a non-passive element: a Touchstone part that is not passive, or a line whose fitted
$\mathbf{R}$ or $\mathbf{G}$ is not positive semi-definite.

### 6.3 Mixed-mode parameters

For a pair of ports $(p, n)$, the differential and common voltages and currents
$V_d = V_p - V_n$, $I_d = (I_p - I_n)/2$, $V_c = (V_p + V_n)/2$, $I_c = I_p + I_n$ with the
references $2Z_0$ and $Z_0/2$ have the waves

$$a_d = \frac{V_d + 2Z_0 I_d}{2\sqrt{2Z_0}} = \frac{a_p - a_n}{\sqrt2}, \qquad a_c = \frac{V_c + \tfrac{Z_0}{2} I_c}{2\sqrt{Z_0/2}} = \frac{a_p + a_n}{\sqrt2},$$

and likewise for $b$ (Bockelman and Eisenstadt 1995). With $\mathbf{M}$ the orthogonal matrix
whose rows are these combinations (differential ports first, then common, then the ports
in no pair), $\mathbf{a}_m = \mathbf{M}\mathbf{a}$ and

$$\mathbf{S}_m = \mathbf{M}\mathbf{S}\mathbf{M}^\mathsf{T}, \tag{8}$$

since $\mathbf{M}^{-1} = \mathbf{M}^\mathsf{T}$. The block of the paired ports is reported. As
$\mathbf{M}$ is orthogonal, $\mathbf{S}_m$ is reciprocal and passive when $\mathbf{S}$ is.

## 7. Time-domain responses

### 7.1 Step response from the sweep

The sweep has the frequencies $f_k = k\,\Delta f$, $k = 1 \ldots K$, $\Delta f = f_\text{max}/K$.
A parameter $\Gamma(f)$ (a reflection or a transmission) is extended to DC by the real part
of its lowest sample, $\Gamma_0 = \operatorname{Re}\Gamma(f_1)$, and to negative frequencies by
$\Gamma(-f) = \overline{\Gamma(f)}$. The excitation is a Gaussian step: its derivative, the
Gaussian pulse of standard deviation $\sigma$, has the spectrum
$H(f) = \exp\big(-\tfrac12(2\pi f\sigma)^2\big)$, and its 10–90 % rise time is
$t_r = 2 \cdot 1.2816\,\sigma = 2.563\,\sigma$. The pulse response is the inverse Fourier
integral, approximated by the rectangle rule over the sweep,

$$h(t) = \Delta f \Big( \Gamma_0 + 2 \sum_{k=1}^{K} \operatorname{Re}\big[ \Gamma(f_k)\, H(f_k)\, e^{j 2\pi f_k (t - t_0)} \big] \Big),$$

delayed by $t_0 = 4\sigma$ so that it starts from rest, and the step response is its running
integral $s(t) = \int_0^t h\,dt'$ by the rectangle rule with step
$\Delta t = \min(t_r/10, 1/(4 f_\text{max}))$. Times are reported from $t_0$, the middle of the
incident step.

### 7.2 Properties and limits

- **Band limit.** The sweep truncates the spectrum at $f_\text{max}$, which produces ringing
  of the order of $H(f_\text{max})$. The rise time is therefore held at no less than
  $0.9/f_\text{max}$, where $2\pi f_\text{max}\sigma = 2.21$ and $H(f_\text{max}) = 0.088$.
- **Periodicity.** A spectrum sampled at the spacing $\Delta f$ is the transform of a function
  of period $1/\Delta f = K/f_\text{max}$: the computed $h$ is the true response folded with
  that period. The window $t_0$ plus the duration must stay below it, and the true
  response must have died out by then; a longer window is reported.
- **Final value.** Over a period only the DC term integrates to a non-zero value, so the
  step response settles to $\Gamma_0$: the reflection at the lowest frequency stands for
  the one at DC.

### 7.3 TDR

For each port (each differential port, with $S_{dd,kk}$ and the reference $2Z_0$), the step
response of the reflection is $\rho(t)$, and the impedance is

$$Z(t) = Z_\text{ref}\,\frac{1 + \rho(t)}{1 - \rho(t)}, \tag{9}$$

the impedance of a lumped load with the reflection $\rho$. On a distributed line it reads
the impedance at the distance $v t/2$ to first order: reflections of reflections and the
loss of the line in between are not removed. The step responses at every port for a step
into a driven port are formed the same way from the column $S_{jd}$.

## 8. Error

| part | error |
| --- | --- |
| line parameters | as in `docs/methods/si.md` |
| samples as uniform sections | first order in the change of the cross-section along a sample |
| transfer matrix | rounding level (section 3.3) |
| via model | against full-wave solutions: the stub notch of a via through thin pairs of planes with return vias 7 % high, that of a stub in a single 2.8 mm dielectric on a small board without return vias 19 % low (`docs/validation/networks.md`); no radiation from the board's edges |
| coaxial via model | the same cases 56 % and 20 % high |
| higher-mode series | below $10^{-15}$ H (section 4.4) |
| TDR | band limit and periodicity (section 7.2); first-order reading of (9) |

## 9. Implementation

| step | where |
| --- | --- |
| nodes, junctions, samples, coupled sections | `src/core/src/network.cpp`, `Channel::run` |
| transfer matrix (2), closed form, Padé (3.3) | `network.cpp`, `section_abcd`, `line_abcd`, `expm` |
| admittance blocks (3) | `network.cpp`, `stamp_abcd` |
| via ladder: nodes, weights, (6), 4.5, 4.6 | `src/core/src/via_model.cpp`, `via_ladder` |
| radial impedance (4), Hankel functions | `via_model.cpp`, `radial_self`, `radial_transfer`, `hankel2` |
| shorted barrels (5), assembly 4.7 | `network.cpp`, `Channel::run`, case `Element::Cavity` |
| finite pairs | `src/core/src/cavity.cpp`, `Cavity::impedance` |
| coaxial model 4.8 | `network.cpp`, `Channel::run`, `Element::Barrel` |
| parts and loads | `network.cpp`, cases `Series`, `Shunt`, `Part` |
| nodal solution, (7) | `network.cpp`, `solve`, after the element loop |
| mixed mode (8) | `network.cpp`, after the frequency loop |
| step responses, TDR (9) | `network.cpp`, `step_response`, `tdr` |

## References

- Abramowitz, M., and Stegun, I. A. (1964). *Handbook of Mathematical Functions*. National
  Bureau of Standards, Applied Mathematics Series 55.
- Bockelman, D. E., and Eisenstadt, W. R. (1995). Combined differential and common-mode
  scattering parameters: theory and simulation. *IEEE Transactions on Microwave Theory and
  Techniques* 43(7), 1530–1539. doi:10.1109/22.392911.
- Golub, G. H., and Van Loan, C. F. (2013). *Matrix Computations*, 4th ed. Johns Hopkins
  University Press.
- Grover, F. W. (1946). *Inductance Calculations: Working Formulas and Tables*. Van Nostrand.
- Kurokawa, K. (1965). Power waves and the scattering matrix. *IEEE Transactions on
  Microwave Theory and Techniques* 13(2), 194–202. doi:10.1109/TMTT.1965.1125964.
- Moler, C., and Van Loan, C. (2003). Nineteen dubious ways to compute the exponential of a
  matrix, twenty-five years later. *SIAM Review* 45(1), 3–49. doi:10.1137/S00361445024180.
- Paul, C. R. (2008). *Analysis of Multiconductor Transmission Lines*, 2nd ed. Wiley.
  doi:10.1002/9780470131541.
- Pozar, D. M. (2012). *Microwave Engineering*, 4th ed. Wiley.
- Rimolo-Donadio, R., Gu, X., Kwark, Y. H., Ritter, M. B., Archambeault, B., de Paulis, F.,
  Zhang, Y., Fan, J., Brüns, H.-D., and Schuster, C. (2009). Physics-based via and trace
  models for efficient link simulation on multilayer structures up to 40 GHz. *IEEE
  Transactions on Microwave Theory and Techniques* 57(8), 2072–2083.
  doi:10.1109/TMTT.2009.2025470.
- Zhang, Y., Fan, J., Selli, G., Cocchini, M., and de Paulis, F. (2008). Analytical
  evaluation of via-plate capacitance for multilayer printed circuit boards and packages.
  *IEEE Transactions on Microwave Theory and Techniques* 56(9), 2118–2128.
  doi:10.1109/TMTT.2008.2002237.
