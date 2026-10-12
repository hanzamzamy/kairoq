# Time-domain response: method

This document derives the time-domain simulation of a channel with its drivers and
receivers: the impulse responses from the swept scattering matrix, the convolution in
wave variables, the terminations and their integration, the switching coefficients of IBIS
buffers, and the error of each step. The user guide is `docs/solvers/timedomain.md`; the
validation, `docs/validation/timedomain.md`. The channel's sweep is that of
`docs/methods/networks.md`.

## 1. Symbols

| symbol | meaning | unit |
| --- | --- | --- |
| $S_{ij}(f_k)$ | the channel's scattering matrix at $f_k = k\,\Delta f$, $k = 1 \ldots K$ | |
| $f_\text{max} = K \Delta f$ | highest frequency of the sweep | Hz |
| $Z_0$ | the channel's reference impedance | Ω |
| $a_i, b_i$ | waves entering and leaving port $i$ | $\sqrt{\text{W}}$ |
| $s_{ij}[m]$ | sampled impulse response times the step | |
| $\Delta t = 1/(2 n f_\text{max})$ | time step, $n$ the oversampling | s |
| $k_u, k_d$ | switching coefficients of an IBIS output's pullup and pulldown | |

## 2. The channel in the time domain

### 2.1 Wave variables

With $a_i = (V_i + Z_0 I_i)/(2\sqrt{Z_0})$ and $b_i = (V_i - Z_0 I_i)/(2\sqrt{Z_0})$, $I_i$ into
the channel (Kurokawa 1965), the port's voltage and current are

$$V_i = \sqrt{Z_0}\,(a_i + b_i), \qquad I_i = \frac{a_i - b_i}{\sqrt{Z_0}}, \tag{1}$$

and in the frequency domain $\mathbf{b} = \mathbf{S}\mathbf{a}$. The channel is linear and time
invariant, so in the time domain $\mathbf{b} = \mathbf{s} * \mathbf{a}$ with the impulse
responses $\mathbf{s}(t)$, the inverse Fourier transforms of $\mathbf{S}(f)$. Wave variables are
used because $\mathbf{S}$ is bounded for every channel, including open and short ends, where
$\mathbf{Z}$ or $\mathbf{Y}$ would be singular.

### 2.2 Impulse responses from the sweep

$\mathbf{S}$ is known on the grid $f_k$; it is extended to DC by the real part of its lowest
sample and to negative frequencies by $S(-f) = \overline{S(f)}$ (the response is real). The
inverse transform is approximated by the rectangle rule over the grid, tapered by a raised
cosine over the last quarter of the band,

$$W(f) = \begin{cases} 1, & f \leq \tfrac34 f_\text{max} \\ \tfrac12\Big(1 + \cos\dfrac{\pi (f - \tfrac34 f_\text{max})}{\tfrac14 f_\text{max}}\Big), & f > \tfrac34 f_\text{max} \end{cases}$$

and multiplied by $\Delta t$ so that the convolution of section 2.3 is a sum:

$$s_{ij}[m] = \Delta f\, \Delta t \Big( S_{ij,0} + 2 \sum_{k=1}^{K} W(f_k) \operatorname{Re}\big[ S_{ij}(f_k)\, e^{j 2\pi f_k m \Delta t} \big] \Big). \tag{2}$$

The taper removes the jump of the truncated spectrum at $f_\text{max}$, which would make the
responses ring at $f_\text{max}$ with an amplitude that decays only as $1/t$ (the Gibbs
phenomenon); in exchange the band is reduced, to about $0.875 f_\text{max}$ by its mean.

**Periodicity.** A spectrum sampled at the spacing $\Delta f$ is the transform of a function
of period $T = 1/\Delta f$; (2) is periodic in $m$ with $T/\Delta t$ samples. The simulation
lasts at most $T/2$, and the true responses must have decayed within $T$.

**Causality.** A band-limited spectrum cannot hold a response that starts abruptly at
$t = 0$: an instantaneous part of $\mathbf{s}$, such as the reflection at a port's own
junction, spreads to both sides of $t = 0$, and the part before it appears, wrapped, at the
end of the period. The samples $m = -w \ldots -1$ with $w = \min(T/(4\Delta t), 3/(f_\text{max}\Delta t))$
are added to $m = 0$. Over a period the sum of (2) is exactly the DC value $S_{ij,0}$
(only the $k = 0$ term survives the sum over $m$), so moving the precursor to $t = 0$ keeps
the response causal and its DC value exact; without it a step into a mismatched port
settles at the wrong level (`docs/validation/timedomain.md`, "A defect found").

### 2.3 Convolution

The waves leaving the channel at step $n$ are

$$b_i[n] = \sum_j s_{ij}[0]\, a_j[n] + h_i[n], \qquad h_i[n] = \sum_j \sum_{m=1}^{n} s_{ij}[m]\, a_j[n - m], \tag{3}$$

the rectangle rule for $\int_0^t s(\tau) a(t - \tau)\, d\tau$. The history $h_i[n]$ depends only
on earlier steps; the instantaneous part $s_{ij}[0]$ couples the unknowns $a_j[n]$ to the
terminations. The cost is $O(P^2 N^2)$ for $P$ ports and $N$ steps.

## 3. Terminations

Each port is closed by a termination that gives the current $I_i^\text{t}(V_i, \mathbf{x}_i)$
it drives into the channel at the port voltage $V_i$, with internal unknowns $\mathbf{x}_i$ and
equations $\mathbf{r}_i(V_i, \mathbf{x}_i) = 0$. Capacitors and inductors are integrated by
backward Euler: $C\,dV/dt \to C (V^n - V^{n-1}) / \Delta t$ and likewise for $L\,dI/dt$.

- **Open:** $I^\text{t} = 0$.
- **Linear driver:** $I^\text{t} = (e(t) - V)/R$, with $e(t)$ moving linearly between its levels
  at each transition over $1.25\,t_{10\text{–}90}$ (a linear ramp spends 80 % of its duration
  between 10 % and 90 %).
- **Linear receiver:** $I^\text{t} = -(V - V_\text{term})/R - C (V - V^{n-1})/\Delta t$, with the
  resistance absent when it is not given.
- **IBIS buffer** (section 4): the pin at $V$, the package ($R_\text{pkg}$, $L_\text{pkg}$ in
  series, $C_\text{pkg}$ at the pin), the die at $V_d$ with $C_\text{comp}$ and the buffer's
  current. With the package current $i_L$ from pin to die as the second unknown,

$$\begin{aligned}
r_1 &= i_L - I_\text{buf}(V_d, t) - C_\text{comp}\frac{V_d - V_d^{n-1}}{\Delta t} = 0, \\
r_2 &= V - V_d - R_\text{pkg}\, i_L - L_\text{pkg}\frac{i_L - i_L^{n-1}}{\Delta t} = 0, \\
I^\text{t} &= -i_L - C_\text{pkg}\frac{V - V^{n-1}}{\Delta t}.
\end{aligned}$$

### 3.1 The system at each step

The unknowns are the entering waves $a_i[n]$ and the terminations' own unknowns. With
$b_i$ from (3) and $V_i$, $I_i$ from (1), the equations are

$$I_i - I_i^\text{t}(V_i, \mathbf{x}_i) = 0, \qquad \mathbf{r}_i(V_i, \mathbf{x}_i) = 0 \quad \text{for every port } i,$$

solved by Newton's method (Nocedal and Wright 2006, ch. 11) from the previous step's
values, with the Jacobian by forward differences (step $10^{-7}\max(1, |x|)$) and Gaussian
elimination with partial pivoting, until the largest residual is below $10^{-12}$ or the
largest update below $10^{-12}$, at most 30 iterations. Steps that do not converge are
counted and reported. For linear terminations the system is linear and one iteration
solves it.

### 3.2 Backward Euler

Backward Euler is first-order accurate and L-stable: it does not oscillate at
discontinuities of the sources and damps the stiff modes of the package. It also damps
oscillations: applied to $\dot y = j\omega y$ it multiplies $y$ by
$(1 - j\omega\Delta t)^{-1}$ per step, of modulus $(1 + \omega^2\Delta t^2)^{-1/2}$, which over
time is the decay $e^{-\omega^2\Delta t\, t / 2}$, an added damping ratio of

$$\zeta_\text{num} \approx \frac{\omega\Delta t}{2} \tag{4}$$

(Hairer and Wanner 1996). For a resonance of the package or the loads at 2 GHz
with $\Delta t = 6.25$ ps (sweep to 20 GHz, $n = 4$), $\zeta_\text{num} \approx 0.04$: ringing at
the die is damped more than a low-loss package would damp it. Raising the oversampling $n$
reduces $\zeta_\text{num}$ in proportion. The channel itself is not integrated, so its
reflections and ringing are not affected.

## 4. IBIS buffers

### 4.1 Current of the buffer

An IBIS model (IBIS Open Forum 2023) tabulates the currents into the pin of its stages
against voltage, per corner: the pullup against $V_\text{pu} - V$, the pulldown against
$V - V_\text{pd}$, the ground and power clamps against $V - V_\text{gc}$ and $V_\text{pc} - V$.
They are interpolated linearly, and extrapolated linearly from the last two points of a
table. The current at the die is

$$I_\text{buf}(V, t) = k_u(t)\, I_\text{pu}(V_\text{pu} - V) + k_d(t)\, I_\text{pd}(V - V_\text{pd}) + I_\text{gc}(V - V_\text{gc}) + I_\text{pc}(V_\text{pc} - V), \tag{5}$$

with $k_u = 0$, $k_d = 1$ for an output at rest low; an input has $k_u = k_d = 0$.

### 4.2 Switching coefficients

A V-t waveform of an edge gives the die voltage $V(t)$ of the buffer driving a fixture
$R_f$ to $V_f$ with $C_f$. At the pin the buffer's current balances the fixture's and the
capacitances':

$$I_\text{buf}(V, t) = -\frac{V - V_f}{R_f} - (C_f + C_\text{comp}) \frac{dV}{dt}, \tag{6}$$

with $dV/dt$ by a central difference of the interpolated waveform. With (5) this is one
linear equation in $k_u(t)$ and $k_d(t)$ at each time.

- **Two waveforms** (fixtures to different voltages): at each time of either waveform,
  the two equations
  $$\begin{pmatrix} p_a & q_a \\ p_b & q_b \end{pmatrix} \begin{pmatrix} k_u \\ k_d \end{pmatrix} = \begin{pmatrix} r_a \\ r_b \end{pmatrix},$$
  $p$ and $q$ the pullup and pulldown currents at the waveform's voltage and $r$ the right
  side of (6) less the clamps, are solved by Cramer's rule; times where the determinant
  vanishes are skipped.
- **One waveform:** the stages are taken to switch complementarily, $k_d = 1 - k_u$, and
  $k_u = (r - q)/(p - q)$.
- **No waveform:** a linear transition of $k_u$ from 0 to 1 (and $k_d$ from 1 to 0) over the
  time the `[Ramp]` slope of 20–80 % takes for the full swing $V_\text{pu} - V_\text{pd}$.

During a pattern, each transition applies the rising or falling coefficients from its
time, interpolated linearly; the last transition that has begun decides.

## 5. Patterns

A PRBS of order $n$ is the output of the Fibonacci shift register of the polynomial
$x^n + x^m + 1$ from the all-ones state, with $(n, m)$ = (7, 6), (9, 5), (15, 14), (23, 18),
(31, 28), each a primitive trinomial, so the sequence has the maximal period $2^n - 1$
with $2^{n-1}$ ones per period (Golomb 1967). The transitions of a pattern are its bit
boundaries where the level changes, at the bit rate after the delay.

## 6. Reported quantities

Per port: the voltage over time, its extremes and final value; for single edges, the
overshoot $\max(0, V_\text{max} - \max(V_\text{final}, V_\text{initial}))$ and the undershoot
$\max(0, \min(V_\text{final}, V_\text{initial}) - V_\text{min})$.

## 7. Error

| source | effect |
| --- | --- |
| band limit | details finer than about $1/f_\text{max}$ are lost; edges faster than $0.9/f_\text{max}$ are reported |
| taper | bandwidth reduced to about $0.875 f_\text{max}$ |
| DC extension | the lowest sample stands for DC: on the validation channel the final value is within 0.1 % |
| periodicity | at most half the period is simulated |
| rectangle rule (2), (3) | with the sampling at $2n f_\text{max}$ the spectrum is resolved; the differences against an independent frequency-domain solution are below 0.6 % of the swing (`docs/validation/timedomain.md`) |
| backward Euler | first order; damping ratio (4) added to resonances of the terminations |
| switching coefficients | from interpolated, differentiated waveforms: 0.011 against a known ramp |

## 8. Implementation

| step | where |
| --- | --- |
| impulse responses (2), taper, precursor | `src/core/src/timedomain.cpp`, `simulate` |
| convolution (3), Newton step of 3.1 | `timedomain.cpp`, `simulate`, the time loop |
| terminations of section 3 | `timedomain.cpp`, `Open`, `Linear`, `Ibis` |
| buffer current (5), tables | `src/core/src/ibis.cpp`, `IbisModel::current`, `IbisTable::at` |
| switching coefficients 4.2 | `ibis.cpp`, `switching_coefficients` |
| patterns | `timedomain.cpp`, `prbs`, `Stimulus::sequence` |

## References

- Golomb, S. W. (1967). *Shift Register Sequences*. Holden-Day.
- Hairer, E., and Wanner, G. (1996). *Solving Ordinary Differential Equations II: Stiff and
  Differential-Algebraic Problems*, 2nd ed. Springer. doi:10.1007/978-3-642-05221-7.
- IBIS Open Forum (2023). *I/O Buffer Information Specification*, version 7.2.
  https://ibis.org/ver7.2/
- Kurokawa, K. (1965). Power waves and the scattering matrix. *IEEE Transactions on
  Microwave Theory and Techniques* 13(2), 194–202. doi:10.1109/TMTT.1965.1125964.
- Nocedal, J., and Wright, S. J. (2006). *Numerical Optimization*, 2nd ed. Springer.
  doi:10.1007/978-0-387-40065-5.
