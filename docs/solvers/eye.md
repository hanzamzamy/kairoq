# Eye diagrams

The eye of a channel (`docs/solvers/networks.md`) between a linear driver and a linear
receiver: the response to a single bit, the eye of a bit pattern, and the statistical eye
at a bit error ratio (BER) with noise, jitter and equalisation. Validation:
`docs/validation/eye.md`.

## Method

**Transfer.** At each frequency of the channel's sweep the impedance matrix
$Z = Z_0 (I + S)(I - S)^{-1}$ and the terminations give the port voltages

$$(I + Z Y_T)\, V = Z Y_T V_s,$$

with $Y_T$ the terminations' admittances to the reference (the source resistance at the
driver's ports, the receiver's resistance and capacitance at its ports, none elsewhere)
and $V_s$ the source: 1 at a single-ended driver, $\pm\tfrac12$ on the legs of a
differential one. The transfer $H(f)$ is the receiver's voltage (differential: the legs'
difference) for this unit source, multiplied by the CTLE

$$H_\text{CTLE}(f) = A_\text{DC}\,\frac{1 + jf/f_z}{(1 + jf/f_{p1})(1 + jf/f_{p2})}.$$

**Single-bit response.** One bit of duration $T$ with linear edges over
$\tau = 1.25\, t_r$ (the 10–90 % time $t_r$) has the spectrum

$$B(f) = \frac{\left(1 - e^{-j 2\pi f T}\right)\left(1 - e^{-j 2\pi f \tau}\right)}{(j 2\pi f)^2\, \tau}.$$

The response $p(t)$ is the inverse Fourier series of $H B$ over the sweep's period
$1/\Delta f$, tapered over the sweep's last quarter as for the time-domain response, with
the DC term $H(0)\, T$ from the channel's lowest frequency and the CTLE's DC gain. It is
sampled at $T/n$ ($n$ samples per UI) from a tenth of the period before the bit. A
transmitter FFE with taps $w_i$ and main tap $m$ makes it
$\sum_i w_i\, p(t - (i - m) T)$.

**Pattern eye.** The symbols are $s_k = \pm 1$ times half the swing $A$ about the middle
level, which is 0 for a differential driver and $A\, H(0) \sum_i w_i$ for a single-ended
one (levels 0 and the swing). The waveform of a pattern is the superposition

$$v(t) = v_\text{mid} + A \sum_k s_k\, p(t - kT),$$

folded over two UI centred on the eye's best phase, after the response's span has filled.
Its density counts, per time and voltage bin, the fraction of the bits whose trace passes.

**Statistical eye.** At a sampling phase $\varphi$ the cursors are
$c_k = A\, p(\varphi + kT)$, for the bits $k$ UI earlier ($k > 0$) and later ($k < 0$)
within the response's span (the samples above $10^{-4}$ of its peak). A DFE of $N$ taps
subtracts $d_k$ from the post-cursors $c_1 \ldots c_N$: the post-cursors at the
response's peak, each limited to a fraction of the main cursor. With equiprobable symbols
the intersymbol interference $Y = \sum_{k \ne 0} s_k c_k$ has the distribution of
successive two-point convolutions, built on a histogram of 2048 bins between
$\pm \sum |c_k|$, each value split linearly between its neighbouring bins. The level of a
one is $c_0 + Y$ plus Gaussian noise of deviation $\sigma$:

$$P(\text{one below } x) = \sum_j w_j \sum_i h_{ij}\, \Phi\!\left(\frac{x - c_{0j} - y_i}{\sigma}\right),$$

a mixture over the sampling offsets $j$ of the jitter with probabilities $w_j$: Gaussian
random jitter, each offset of one sample taking the probability of its interval to
$\pm 9$ deviations, and deterministic jitter as a dual Dirac of $\pm\tfrac12$ its
peak-to-peak value. The inner edge of the ones is the $x$ where this probability is the
BER (by bisection; without noise, from the cumulative histogram). Zeros mirror ones about
the middle level, so the eye's inner edges are $v_\text{mid} \pm x(\varphi)$.

The phases cover two UI centred on the response's peak. The **eye height** is the largest
opening $2x$, at the **best phase**; the **eye width** is the run of open phases around it,
at most one UI. The **worst case** is the peak-distortion opening
$2\left(c_0 - \sum_{k \ne 0} |c_k|\right)$ without noise or jitter: at a BER below
$2^{-N}$ for $N$ cursors the statistical eye reaches it.

## Project file

```yaml
si:
  channels: [{name: PCIE_TX, ports: [U1.A1, U1.A2, J1.3, J1.4]}]
  network: {f_max_hz: 20e9, points: 400}
  eyes:
    - name: TX_8G
      channel: PCIE_TX
      tx: [1, 2]            # differential: ports of the positive and negative legs
      rx: [3, 4]
      bit_rate: 8e9
      swing_v: 0.8
      noise_mv: 2
      rj_ps: 1
      dj_ps: 10
      ffe: [0.8, -0.2]
      ctle: {dc_gain_db: -6, f_zero_hz: 2e9, f_pole1_hz: 8e9, f_pole2_hz: 20e9}
      dfe: {taps: 2}
```

| key | default | meaning |
| --- | --- | --- |
| `name`, `channel` | | the eye, and a channel of `si.channels` |
| `tx`, `rx` | | `[p]` single-ended or `[p, n]` differential, ports from 1; other ports open |
| `bit_rate` | $5 \times 10^9$ | bit/s |
| `swing_v` | 1 | the source's open-circuit peak-to-peak swing (differential: of the legs' difference) |
| `r_source_ohm`, `rise_ps` | 50, 30 | per leg; the 10–90 % time of the source's edges |
| `r_term_ohm`, `c_term_pf` | 50, 0 | the receiver per leg to the reference |
| `pattern`, `bits` | `prbs7`, 4096 | the pattern's eye: a PRBS or a string of 0 and 1, repeated to `bits` |
| `samples_per_ui` | 32 | time resolution of the response and of the phases |
| `noise_mv`, `rj_ps`, `dj_ps` | 0 | Gaussian noise (rms), random jitter (rms), deterministic jitter (peak-to-peak) |
| `ber` | $10^{-12}$ | the bit error ratio of the statistical eye |
| `ffe`, `ffe_main` | none, 1 | transmitter taps in order; the main tap, from 1 |
| `ctle` | none | `dc_gain_db`, `f_zero_hz`, `f_pole1_hz`, `f_pole2_hz` (0: none) |
| `dfe` | none | `taps`, `limit` (1: the largest tap, as a fraction of the main cursor) |

## Use

```
kairoq eye project.yaml [<input>...] [-o dir]
```

prints per eye the height and width at the BER, the worst case, the main cursor, the
number of pre- and post-cursors and the DFE's taps. With `-o`, `NAME.svg` is the pattern's
eye with the statistical eye's inner edges, `NAME_pulse.csv` the single-bit response and
`NAME_contour.csv` the inner edges per phase.

```python
r = kairoq.eye(board, {"channels": [...], "network": {...}, "eyes": [...]})[0]
r["eye_height_v"], r["eye_width_ui"], r["cursors"], r["pulse_t"], r["pulse_v"], r["svg"]
```

## Limits

- Linear drivers and receivers only; the buffers' nonlinearity (IBIS) is in the
  time-domain response, not in the eye. No IBIS-AMI models: the equalisers are the
  parametric FFE, CTLE and DFE above, not a device's adaptation.
- No compliance metrics or masks of any standard.
- The DFE's taps are fixed at the response's peak; decision errors are not propagated.
- Symbols are independent and equiprobable: the statistical eye does not model a coding's
  run-length limits.
- Jitter is applied at the sampler (the sampling phase), not to the transmitter's edges;
  crosstalk from other aggressors is not added to the statistical eye.
- The channel is linear and band-limited: edges faster than $0.9/f_\text{max}$ and bit rates
  above $0.75\, f_\text{max}$ are reported. The response must decay within the sweep's
  period, which is reported otherwise.
