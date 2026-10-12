# Crosstalk: method

This document derives the crosstalk scan of the signal-integrity analysis: the backward and
forward coupling coefficients of a cut, their sum along a coupled run, and the near- and
far-end crosstalk estimated from them. The waveforms with the actual ends come from a
channel (`docs/methods/networks.md`). The user guide is `docs/solvers/si.md`, "Crosstalk";
the validation, `docs/validation/crosstalk.md`.

## 1. Symbols

| symbol | meaning | unit |
| --- | --- | --- |
| $L_s, C_s$ | the victim's own inductance and Maxwell capacitance per length | H/m, F/m |
| $L_m$ | mutual inductance per length | H/m |
| $C_m$ | magnitude of the mutual Maxwell capacitance per length ($C_m = -C_{12} > 0$) | F/m |
| $Z_c = \sqrt{L_s/C_s}$, $v = 1/\sqrt{L_s C_s}$ | impedance and velocity of the victim | Ω, m/s |
| $f(t)$ | the aggressor's incident wave at its near end | V |
| $\ell$, $\tau = \ell/v$ | coupled length and its delay | m, s |
| $t_r$ | 10–90 % rise time of the step | s |
| $K_b$, $K_f$ | backward (dimensionless) and forward (s/m) coupling coefficients | |

## 2. Model and assumptions

Two lines run together over a length $\ell$: the aggressor (1), driven at its near end
($z = 0$) by a wave $f(t)$, and the victim (2), quiet. Both are lossless and quasi-TEM with
per-unit-length matrices from the cut (`docs/methods/si.md`) at the reference frequency of
the analysis. The assumptions are:

1. **Weak coupling.** The victim's current and voltage do not act back on the aggressor:
   terms of second order in $L_m/L_s$ and $C_m/C_s$ are dropped. The aggressor then
   carries the undisturbed wave $V_1 = f(t - z/v)$, $I_1 = V_1/Z_c$.
2. **Equal lines.** Aggressor and victim have one impedance $Z_c$ and velocity $v$; those of
   the victim are used.
3. **Matched ends.** No wave enters the victim at either end, and none is reflected.

## 3. Derivation

The victim obeys the coupled telegrapher's equations with the aggressor's wave as source:

$$\frac{\partial V_2}{\partial z} + L_s \frac{\partial I_2}{\partial t} = -L_m \frac{\partial I_1}{\partial t}, \qquad \frac{\partial I_2}{\partial z} + C_s \frac{\partial V_2}{\partial t} = C_m \frac{\partial V_1}{\partial t}. \tag{1}$$

Define its forward and backward waves $F = \tfrac12(V_2 + Z_c I_2)$ and
$B = \tfrac12(V_2 - Z_c I_2)$. Using $L_s = Z_c/v$ and $Z_c C_s = 1/v$, the first equation of
(1) plus $Z_c$ times the second, and the first minus $Z_c$ times the second, give

$$\Big(\frac{\partial}{\partial z} + \frac1v \frac{\partial}{\partial t}\Big) F = \frac12 \Big(Z_c C_m - \frac{L_m}{Z_c}\Big) f'\Big(t - \frac zv\Big), \qquad \Big(\frac{\partial}{\partial z} - \frac1v \frac{\partial}{\partial t}\Big) B = -\frac12 \Big(Z_c C_m + \frac{L_m}{Z_c}\Big) f'\Big(t - \frac zv\Big). \tag{2}$$

**Far end.** Along the forward characteristic $t - z/v = \text{const}$ the source of the first
equation is constant. With $F(0, t) = 0$ (matched near end),

$$V_2(\ell, t) = F(\ell, t) = K_f\, \ell\, f'(t - \tau), \qquad K_f = \frac12\Big(Z_c C_m - \frac{L_m}{Z_c}\Big) = \frac12 \Big(\frac{C_m}{C_s} - \frac{L_m}{L_s}\Big)\sqrt{L_s C_s}. \tag{3}$$

**Near end.** Along the backward characteristic $t + z/v = \text{const}$, from $B(\ell, t) = 0$
(matched far end) to $z = 0$, the source is met at the retarded time $t - 2z/v$:

$$V_2(0, t) = B(0, t) = \frac12 \Big(Z_c C_m + \frac{L_m}{Z_c}\Big) \int_0^\ell f'\Big(t - \frac{2z}{v}\Big) dz = K_b \big[ f(t) - f(t - 2\tau) \big], \qquad K_b = \frac14 \Big(\frac{L_m}{L_s} + \frac{C_m}{C_s}\Big), \tag{4}$$

using $\tfrac v2 \cdot \tfrac{L_m}{Z_c} = \tfrac12 \tfrac{L_m}{L_s}$ and $\tfrac v2 \cdot Z_c C_m = \tfrac12 \tfrac{C_m}{C_s}$.
These are the classical weak-coupling results (Paul 2008).

**Homogeneous dielectric.** In one dielectric $\mathbf{L}\mathbf{C} = \mu\varepsilon\, \mathbf{I}$ for
the TEM lines, so $L_m/L_s = C_m/C_s$ to first order and $K_f = 0$: no far-end crosstalk.
On microstrip the field partly in air makes $L_m/L_s > C_m/C_s$ and $K_f < 0$.

## 4. Estimates for a step

For a step of amplitude 1 and rise time $t_r$:

- **Near end.** $f(t) - f(t - 2\tau)$ peaks at $\min(1, 2\tau/t_r)$ for a linear ramp: the
  near-end crosstalk is $K_b \min(1, 2\tau/t_r)$. It saturates at $K_b$ once the round trip
  of the coupled run outlasts the rise.
- **Far end.** The largest slope of the ramp is $1/t_r$, so the far-end crosstalk is
  $\lvert K_f \rvert \ell / t_r$, held at most 1 (beyond which the weak-coupling assumption
  has failed). For the Gaussian step of the channel's waveforms the largest slope is
  $2.563/(\sqrt{2\pi}\, t_r) = 1.02/t_r$.

## 5. Along a board

### 5.1 Coefficients of a cut

A cut of the victim's track holds, besides it, the copper of other nets on its layer within
three times the distance to the farther shield, at most 1 mm, as quiet lines
(`docs/methods/si.md`, section 9). The
coupled matrices of the cut give $L_s = L_{11}$, $C_s = C_{11}$ and, for each other net $o$,
$L_m = L_{1o}$ and $C_m = -C_{1o}$, and with them $K_b$ and $K_f$ by (3) and (4). A pair's
partner is not an aggressor. Pairs are those given, or two nets that differ only in a last
`P`/`N`, `p`/`n` or `+`/`-` before trailing digits.

### 5.2 Sum over the run

With $v$ uniform along the run, the sources of (2) add linearly, and coefficients that
vary along $z$ give

$$V_2(\ell, t) = \Big(\int_0^\ell K_f(z)\, dz\Big) f'(t - \tau), \qquad V_2(0, t) = \frac2v \int_0^\ell K_b(z)\, f'\Big(t - \frac{2z}{v}\Big) dz.$$

The far-end coefficient is therefore the signed sum $\sum_k K_{f,k} w_k$ over the victim's
samples where the aggressor lies in the cut, $w_k$ their lengths: exact under the
assumptions of section 2. The near end is a weighted mean of $K_b$ over a window of length
$v t_r / 2$ that slides along the run. When the run is short against that window
($2\tau < t_r$) it equals $\big(\sum_k K_{b,k} w_k\big) \cdot \tfrac{2}{v t_r}$, which is the
estimate with the length-weighted mean $\bar K_b$; on a longer run with varying coupling
the peak lies between $\bar K_b$ and the largest local $K_b$, both of which are reported,
together with the place of the largest. The delay $\tau$ is the sum of the samples'
delays.

A victim and an aggressor whose larger estimate exceeds the threshold make a finding.

## 6. Error

| source | effect |
| --- | --- |
| weak coupling | second order in the coupling; the directional coupler, exact for TEM lines, agrees with the channel to 0.5 % at a coupling of 0.1 (`docs/validation/crosstalk.md`) |
| matched ends | ends that are not matched change the waveforms: at 50 Ω on 63 Ω lines the far end is −6.3 % instead of −5.5 % |
| equal lines | unequal aggressor and victim change both coefficients to first order in their difference |
| near end on long, varying runs | between $\bar K_b$ and the largest $K_b$ (section 5.2) |
| separate runs | runs between the same nets at different places are summed as one |
| coverage | aggressors on the track's own layer only; broadside coupling is not scanned |

## 7. Implementation

| step | where |
| --- | --- |
| quiet lines in the cut | `src/core/src/si.cpp`, `Analysis::make_cut` |
| coefficients (3), (4) per sample, sums of 5.2 | `si.cpp`, `Analysis::run`, the crosstalk block |
| estimates of section 4, findings | `si.cpp`, same block |
| pairs by name | `si.cpp`, `name_pair` |
| waveforms with actual ends | `src/core/src/network.cpp` (`docs/methods/networks.md`) |

## References

- Oliver, B. M. (1954). Directional electromagnetic couplers. *Proceedings of the IRE*
  42(11), 1686–1692. doi:10.1109/JRPROC.1954.274664.
- Paul, C. R. (2008). *Analysis of Multiconductor Transmission Lines*, 2nd ed. Wiley.
  doi:10.1002/9780470131541.
