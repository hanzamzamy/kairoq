# Validation: eye diagrams

## Scope

The eye diagrams (`docs/solvers/eye.md`) are checked against:

1. the time-domain response (`docs/solvers/timedomain.md`) of the same channel and
   terminations, single-ended and differential;
2. an enumeration of all signs of the cursors with Gaussian noise, and the closed form of
   a noise-limited eye;
3. peak-distortion analysis, and the eye of a PRBS7 pattern;
4. the equalisers' definitions (DFE, CTLE, FFE) and the jitter's effect on the width.

The tests are in `tests/python/test_eye.py`. The channels are the synthetic microstrip
(20 mm, 81.8 Ω, 118 ps) and the synthetic coupled microstrip pair (20 mm), swept to
40 GHz in 800 points (period 20 ns), between 50 Ω terminations.

## Results

| case | quantity | reference | result | deviation |
| --- | --- | --- | --- | --- |
| microstrip, 2.5 Gbit/s, 16 bits, 0 to 1 V behind 50 Ω, 30 ps edges | far-end waveform | time-domain response, both at 3.125 ps | | 0.63 mV at most |
| coupled pair, the same, complementary legs of 0.5 V | differential far-end waveform | time-domain response | | 0.72 mV at most |
| microstrip, 10 Gbit/s, noise 5 mV | eye height at $10^{-12}$ | enumeration of the cursors' signs | 0.365649 V | $5 \times 10^{-7}$ V |
| the same, noise 20 mV | eye height at $10^{-12}$ | enumeration | 0.161416 V | $-7 \times 10^{-8}$ V |
| microstrip, 1 Gbit/s, noise 10 mV | eye height at $10^{-12}$ | $2(A - Q\sigma)$, $Q = 7.034$ | 0.358739 V | $-1.6 \times 10^{-5}$ V |
| microstrip, 10 Gbit/s, no noise | eye height at $10^{-12}$ | peak distortion $c_0 - \sum \lvert c_k \rvert$ | 0.430719 V | $-5 \times 10^{-6}$ V |
| the same | worst case | peak distortion | 0.430723 V | 0 |
| the same, PRBS7 of 2048 bits | inner opening of the pattern's eye at its centre | peak distortion | 0.4345 V | +3.8 mV (within 1.5 voltage bins) |
| microstrip, 10 Gbit/s, DFE of 3 taps | taps | the single-bit response 1, 2, 3 UI after its peak | | 0 |
| CTLE, zero and pole at 3 GHz, −6 dB | single-bit response | the response without CTLE times $10^{-6/20}$ | | $7 \times 10^{-16}$ V |
| CTLE, −6 dB, zero 1 GHz, poles 5 and 20 GHz | area of the single-bit response, relative | the DC gain, 0.501187 | 0.501187 | $< 10^{-6}$ |
| FFE $[-0.1, 0.8, -0.1]$ | single-bit response | the shifted, weighted responses | | $< 10^{-12}$ V |
| microstrip, 10 Gbit/s, DJ 12.5 ps (4 samples) | eye width | without jitter, less the DJ: 0.84375 UI | 0.84375 UI | 0 |
| microstrip, 10 Gbit/s, RJ 1.5 ps | eye width | without jitter, less $2 Q \sigma$: 0.7577 UI | 0.7813 UI | +0.75 samples |

- **The single-bit response** and the convolution engine of the time-domain response
  share only the channel's sweep and its taper; the 0.6 mV difference is the time sampling
  of the edges (sampling the response at 32 points per UI and interpolating it gives
  14 mV).
- **The enumeration** takes the cursors at the best phase and the noise, and finds the
  level where the probability of a one below it is the BER by bisection over all
  $2^{N}$ sign combinations, independently of the histograms and the grid of the solver.
- **The closed form** holds where the intersymbol interference is negligible: at
  1 Gbit/s the line's reflections settle within the bit (residual ISI below 0.2 % of the
  main level, used as the tolerance).
- **Peak distortion.** Without noise or jitter, the worst sign combination of $N$ cursors
  has probability $2^{-N}$, above $10^{-12}$ for fewer than 40 cursors, so the
  statistical eye's opening is the worst case; the 5 µV difference is the histogram's
  linear splitting. PRBS7 holds every pattern of 7 bits, and the response spans fewer
  than 7 UI, so the pattern's eye reaches the worst case to within its voltage bins.
- **Jitter.** A dual Dirac of $\pm 2$ samples closes every phase within 2 samples of a
  closed one: the width falls by exactly 4 samples. Gaussian jitter closes a phase whose
  distance to the nearest closed one is below $Q\sigma$, to within the sampling of the
  offsets.

## Defects found

- The CTLE's DC gain was first taken from the channel's lowest frequency (50 MHz), where
  a CTLE with its zero at 1 GHz already adds 0.2 %; the DC term now uses the CTLE's DC
  gain exactly.
- The statistical eye first scanned one UI of phases centred on the response's peak; the
  eye's opening lies elsewhere when the response has a flat top, and its width was cut
  at the window's edge. The phases now cover two UI.

## Real board

`examples/cm5_minima_si.yaml` computes the eyes of the PCIe transmit pair (CM5 connector
to the M.2 connector) at the Gen2 and Gen3 rates: 0.8 V of differential swing behind
50 Ω per leg, 50 ps edges, 50 Ω per leg at the far end, noise 2 mV, RJ 1 ps, DJ 10 ps, at
$10^{-12}$:

| eye | height | width | worst case | main cursor |
| --- | --- | --- | --- | --- |
| 5 Gbit/s | 0.350 V | 0.84 UI | 0.378 V | 0.389 V |
| 8 Gbit/s, de-emphasis $[0.8, -0.2]$ | 0.208 V | 0.72 UI | 0.234 V | 0.306 V |

The channel's sweep takes most of the 14 s; each eye takes 0.1 s.

## Not validated

No eye is compared with a commercial statistical simulator or a measured eye; the
receiver models are linear (no IBIS input or AMI model in the eye).
