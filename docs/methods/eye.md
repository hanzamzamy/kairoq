# Eye diagrams: method

This document derives the eye analysis of a channel between a linear driver and a linear
receiver: the transfer with the terminations, the single-bit response and its
equalisation, the eye of a bit pattern by superposition, and the statistical eye at a bit
error ratio with noise and jitter. The user guide is `docs/solvers/eye.md`; the validation,
`docs/validation/eye.md`. The channel's sweep is that of `docs/methods/networks.md`.

## 1. Symbols

| symbol | meaning | unit |
| --- | --- | --- |
| $\mathbf{S}, \mathbf{Z}$ | the channel's scattering and impedance matrices at a frequency | , Ω |
| $\mathbf{Y}_T$ | diagonal matrix of the terminations' admittances to the reference | S |
| $H(f)$ | transfer from a unit source to the receiver's voltage | |
| $T$ | unit interval (UI), $1/\text{bit rate}$ | s |
| $p(t)$ | single-bit response for a unit source swing | V/V |
| $A$ | half the source's swing | V |
| $c_k$ | cursors: $p$ at the sampling phase $\varphi$ plus $k$ UI | V/V |
| $\sigma$ | rms of the Gaussian noise at the sampler | V |
| $\Phi$ | standard normal distribution function | |

## 2. Transfer

### 2.1 Impedance matrix

From $\mathbf{b} = \mathbf{S}\mathbf{a}$ and the wave definitions (`docs/methods/timedomain.md`,
equation 1), $\mathbf{V} = \sqrt{Z_0}(\mathbf{I} + \mathbf{S})\mathbf{a}$ and
$\mathbf{I}_p = (\mathbf{I} - \mathbf{S})\mathbf{a}/\sqrt{Z_0}$, so

$$\mathbf{Z} = Z_0 (\mathbf{I} + \mathbf{S})(\mathbf{I} - \mathbf{S})^{-1},$$

where $\mathbf{I} - \mathbf{S}$ is invertible for a channel that is not a short circuit at a port.
As both factors are functions of $\mathbf{S}$ they commute; $\mathbf{Z}$ is computed row by row
from $\mathbf{z}_i^\mathsf{T} (\mathbf{I} - \mathbf{S}) = Z_0 (\mathbf{e}_i + \mathbf{S}_i)^\mathsf{T}$.

### 2.2 Terminations

Each port $i$ is closed by a source $V_{s,i}$ behind the admittance $y_i$ (a load has
$V_{s,i} = 0$; an open port $y_i = 0$). The current into the channel is
$I_i = y_i (V_{s,i} - V_i)$, and with $\mathbf{V} = \mathbf{Z}\mathbf{I}$:

$$(\mathbf{I} + \mathbf{Z}\mathbf{Y}_T)\, \mathbf{V} = \mathbf{Z}\mathbf{Y}_T \mathbf{V}_s. \tag{1}$$

The driver's legs have $y = 1/R_s$ and $V_s = 1$ (single-ended) or $\pm\tfrac12$
(differential); the receiver's $y = 1/R_t + j\omega C_t$. $H(f)$ is the receiver's voltage, or
the difference of its legs, solved from (1) by Gaussian elimination with partial pivoting.

### 2.3 CTLE

A continuous-time linear equaliser multiplies $H$ by

$$H_\text{CTLE}(f) = A_\text{DC}\,\frac{1 + jf/f_z}{(1 + jf/f_{p1})(1 + jf/f_{p2})},$$

a pole absent when its frequency is not given. At DC it is exactly $A_\text{DC}$.

## 3. Single-bit response

### 3.1 The bit

An edge from 0 to 1 starting at $t_s$ and rising linearly over $\tau$ is the integral of a
rectangular pulse of height $1/\tau$; its transform is

$$E_{t_s}(f) = \frac{e^{-s t_s} - e^{-s (t_s + \tau)}}{s^2 \tau}, \qquad s = j 2\pi f,$$

with the limit $e^{-s t_s}/s$ for $\tau \to 0$. One bit is a rising edge at 0 and a falling
one at $T$:

$$B(f) = E_0(f) - E_T(f) = \frac{(1 - e^{-sT})(1 - e^{-s\tau})}{s^2 \tau}, \qquad B(0) = T. \tag{2}$$

The ramp lasts $\tau = 1.25\, t_r$ for the 10–90 % time $t_r$.

### 3.2 Inverse transform

As in the time-domain response, the product is tapered by the raised cosine $W(f)$ over the
last quarter of the sweep (`docs/methods/timedomain.md`, section 2.2), and

$$p(t) = \Delta f \Big( H_0 T + 2 \sum_{k=1}^{K} \operatorname{Re}\big[ W(f_k)\, H(f_k)\, B(f_k)\, e^{j 2\pi f_k t} \big] \Big), \tag{3}$$

with $H_0$ the real part of the lowest frequency's transfer times the CTLE's $A_\text{DC}$. It
is sampled at $\Delta t = T/n$ ($n$ samples per UI) over one period $1/\Delta f$, starting a
tenth of the period before the bit, since the band-limited edges spread to both sides of
it. The period must hold at least 8 UI, and the response must decay within it: its span is
from the first to the last sample above $10^{-4}$ of its peak, and a span that reaches the
ends of the period is reported.

### 3.3 Transmitter FFE

A feed-forward equaliser with taps $w_i$ and main tap $m$ sends
$\sum_i w_i\, \text{bit}(t - (i - m)T)$ for each bit; by linearity the response becomes

$$p_\text{FFE}(t) = \sum_i w_i\, p\big(t - (i - m) T\big).$$

## 4. Levels

The symbols are $s_k = \pm 1$, so a waveform is

$$v(t) = v_\text{mid} + A \sum_k s_k\, p(t - kT). \tag{4}$$

For a differential driver $v_\text{mid} = 0$. For a single-ended driver with levels 0 and
$2A$, a bit $b_k \in \{0, 1\}$ is $A(s_k + 1)$; the sum of $A\, p(t - kT)$ over all $k$ is the
response to a constant source $A$, which is $A H_0 \sum_i w_i$, so
$v_\text{mid} = A H_0 \sum_i w_i$. Since the distribution of the $s_k$ is symmetric, the
levels of a zero are those of a one reflected about $v_\text{mid}$: only the ones are
computed.

## 5. Pattern eye

The waveform (4) of the pattern (a PRBS of `docs/methods/timedomain.md`, section 5, or a
given string, repeated to the length asked for) is formed by superposition of the
response's span and folded over two UI centred on the best phase of the statistical eye
(section 6), leaving out the bits within one span of either end. Its density counts, per
time column and voltage bin (128 between the extremes, with 5 % margins), the fraction of
the traces that pass: each trace marks the bins from its value to the midpoints towards
its neighbouring samples.

## 6. Statistical eye

### 6.1 Cursors

At the sampling sample $\varphi$ the cursors are $c_k = p(\varphi + kT)$ for the bits $k$ UI
earlier ($k > 0$, post-cursors) and later ($k < 0$, pre-cursors) within the response's span,
and the level of a one is

$$v_1 = A c_0 + Y, \qquad Y = A \sum_{k \neq 0} s_k c_k \tag{5}$$

(Casper et al. 2002). A decision-feedback equaliser of $N$ taps subtracts $d_k$ from
$c_1 \ldots c_N$, with $d_k$ the post-cursors at the response's peak, each held within
$\pm\lambda |c_0|$ ($\lambda$ the limit); the decisions are taken as correct.

### 6.2 Distribution of the interference

With independent, equiprobable symbols, $Y$ is the sum of independent two-point variables
$\pm A c_k$, so its distribution is the successive convolution of
$\tfrac12(\delta_{A c_k} + \delta_{-A c_k})$. It is computed on a histogram of $M = 2048$ bins
of width $h = 2\sum|A c_k| / (M - 1)$ over $\pm\sum |A c_k|$: starting from a unit mass at 0,
each cursor moves every bin's mass half to $+A c_k$ and half to $-A c_k$. A mass at a value
between two bins is split between them in proportion to its distance from each, which
keeps the mass and the mean exactly and adds a variance of at most $h^2/4$ per cursor: after
$N_c$ cursors the histogram is smeared by at most $\sqrt{N_c}\, h/2$ in standard deviation.

### 6.3 Jitter

Jitter moves the sampling instant. Random jitter of rms $\sigma_j$ is the set of offsets
$o\,\Delta t$, $|o| \leq \lceil 9\sigma_j/\Delta t \rceil$, each with the probability of its
interval, $\Phi\big((o + \tfrac12)\Delta t/\sigma_j\big) - \Phi\big((o - \tfrac12)\Delta t/\sigma_j\big)$;
the mass beyond nine deviations, $2 \times 10^{-19}$, is below any bit error ratio of
interest. Deterministic jitter of peak-to-peak $J$ is the dual-Dirac model, offsets of
$\pm J/2$ with probability $\tfrac12$ each (INCITS 2004); the two combine by convolution.
Offsets of probability below $10^{-3}$ times the BER are dropped. The level of a one is then
a mixture over the offsets $j$ with probabilities $w_j$:

$$P(v_1 < x) = \sum_j w_j \sum_i h_{ij}\, \Phi\Big(\frac{x - A c_{0j} - y_i}{\sigma}\Big), \tag{6}$$

with $h_{ij}$ the histogram of $Y$ at the offset sample.

### 6.4 Inner edge at a BER

The mixture (6) is gathered on a grid of 8192 points, split linearly as in 6.2, spanning
the extremes of $v_1$ widened by $10\sigma$. With noise, the inner edge $x(\varphi)$ of the ones
is the solution of $P(v_1 < x) = \text{BER}$ by bisection (60 halvings, to below
$10^{-18}$ of the span), with $\Phi$ evaluated through the complementary error function so
that tail probabilities keep their relative accuracy. Without noise it is read from the
cumulative distribution, each grid point's mass spread uniformly over its interval. The
eye's inner edges at $\varphi$ are $v_\text{mid} \pm x(\varphi)$.

### 6.5 Reported quantities

The phases cover two UI centred on the response's peak, at $T/n$. The eye height is the
largest opening $2x(\varphi)$, at the best phase; the eye width the run of phases around it
where the eye is open, at most one UI. The worst case is the peak-distortion opening

$$2\Big(A c_0 - \sum_{k \neq 0} A |c_k|\Big),$$

the largest over the phases: the opening for the worst sign combination, without noise or
jitter. That combination has the probability $2^{-N_c}$, so without noise and jitter the
statistical eye at a BER below $2^{-N_c}$ is exactly the worst case.

## 7. Error

| source | effect |
| --- | --- |
| channel sweep | as for the time-domain response: band limit, taper, DC from the lowest frequency, periodicity |
| sampling at $T/n$ | phases and jitter offsets resolved to $1/n$ UI: the eye width to about one sample |
| histogram 6.2 | spread at most $\sqrt{N_c}\, h/2$; against an enumeration of all sign combinations, $5 \times 10^{-7}$ V on a 0.37 V eye (`docs/validation/eye.md`) |
| span | cursors below $10^{-4}$ of the peak beyond the span are left out |
| jitter tails | below $10^{-3}$ BER per dropped offset; beyond $9\sigma_j$ negligible |
| DFE | fixed at the peak; error propagation not modelled |
| symbols | independent and equiprobable; codings with run-length limits are not modelled |

## 8. Implementation

| step | where |
| --- | --- |
| impedance matrix, terminations (1), CTLE | `src/core/src/eye.cpp`, `eye`, the transfer loop |
| bit spectrum (2), response (3), span | `eye.cpp`, `eye`, the single-bit response |
| FFE, levels (4) | `eye.cpp`, `eye` |
| cursors, DFE (5) | `eye.cpp`, `cursors_at`, `equalised` |
| interference histogram 6.2 | `eye.cpp`, `isi_at`, `deposit` |
| jitter 6.3, mixture (6), bisection 6.4 | `eye.cpp`, the phase loop; `normal_mass`, `tail` |
| height, width, worst case | `eye.cpp`, after the phase loop |
| pattern eye | `eye.cpp`, the pattern block |

## References

- Casper, B. K., Haycock, M., and Mooney, R. (2002). An accurate and efficient analysis
  method for multi-Gb/s chip-to-chip signaling schemes. *2002 Symposium on VLSI Circuits,
  Digest of Technical Papers*, 54–57. doi:10.1109/VLSIC.2002.1015043.
- INCITS (2004). *Information Technology – Fibre Channel – Methodologies for Jitter and
  Signal Quality Specification (FC-MJSQ)*. Technical Report INCITS TR-35-2004.
