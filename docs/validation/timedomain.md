# Validation: time-domain response

## Scope

The time-domain response (`docs/solvers/timedomain.md`) is checked against:

1. the lattice diagram of a step into a line with an open end;
2. an independent frequency-domain solution of the same linear network;
3. IBIS buffers against equivalent linear sources and loads;
4. switching coefficients recovered from V-t waveforms computed for a known buffer;
5. the IBIS reader and the bit patterns.

The tests are in `tests/python/test_timedomain.py`. The channel throughout is the
synthetic microstrip (20 mm, 81.8 Ω, 118 ps) swept to 40 GHz in 800 points; time step
3.1 ps.

## Results

| case | quantity | reference | result | deviation |
| --- | --- | --- | --- | --- |
| 1 V behind 50 Ω, 50 ps edge, far end open | near-end step | lattice: $Z/(Z + 50)$ = 0.6207 V | 0.6214 V | +0.11 % |
| the same | far-end step | lattice, lossless: twice that, 1.2414 V | 1.2281 V | −1.1 % |
| the same | final value | 1 V | 1.0009 V | +0.09 % |
| 50 Ω source, far end 100 Ω and 2 pF | far-end waveform over 4 ns | periodic solution in the frequency domain (numpy) | | 6 mV at most (0.6 % of the swing) |
| IBIS buffer: 50 Ω stages, 250 ps ramp, no package or $C_\text{comp}$ | waveforms at both ends | linear source of 0 to 1.8 V behind 50 Ω, the same edge | | $6 \times 10^{-12}$ V |
| IBIS input: 2 pF, clamps off | far-end waveform | a 2 pF load | | $9 \times 10^{-12}$ V |
| IBIS buffer: 40 Ω stages, 1 pF, 200 ps ramp; its rising waveforms into 50 Ω to ground and to $V_{CC}$ (integrated at 0.01 ps) | recovered $k_u$, $k_d$ | the ramp | | 0.011 at most |

- **The lattice diagram** is for a lossless line; the line's conductor and dielectric loss
  lower the step arriving at the far end by about a percent.
- **The frequency-domain solution** is independent of the time stepping: from the
  channel's impedance matrix $Z = Z_0 (I + S)(I - S)^{-1}$, the terminations and the Fourier
  coefficients of a square wave of the sweep's period, written in the test. It shares
  only the channel's sweep and its taper. The 6 mV difference combines the extension of
  the sweep to DC, the folding of the impulse responses' precursor and the time step.
- **Equivalent sources and loads.** A buffer whose stages are linear resistances is a
  Thevenin source; without its clamps conducting, an input is its die capacitance. The
  agreement to $10^{-11}$ V checks the IBIS path (tables, switching from the ramp, the
  Newton step with internal unknowns) against the linear one.
- **Switching coefficients.** The two waveforms are computed in the test by integrating
  the buffer (stages, $C_\text{comp}$, fixture) with the known ramp; the reader recovers the
  ramp within 0.011, an error consistent with interpolating the tabulated waveforms
  (10 ps apart) and differentiating them for the capacitance's current.
- **Reader and patterns.** The comment character, scaled units (`1.5pF`, `250.0p`),
  `NA` values, `[Ramp]` fractions and unread keywords are read as specified; PRBS7 and
  PRBS9 have their maximal periods (127, 511) and their $2^{n-1}$ ones.

## A defect found

The first impulse responses lost the negative-time half of a port's instantaneous
reflection (wrapped to the end of the sweep's period): a step into the 81.8 Ω line behind
50 Ω settled at 0.904 V instead of 1 V. The part before $t = 0$ is now added at $t = 0$.

## Real board

`examples/cm5_minima_si.yaml` sends 20 bits of PRBS7 at 2.5 Gbit/s, 0.8 V behind 50 Ω,
into the positive member of the PCIe transmit pair, all other ends of the pair
terminated with 50 Ω: the far end swings 0 to 0.397 V; the near and far ends of the
negative member see ±0.026 V and ±0.055 V of crosstalk. The four-port of 238 nodes and
the 1601 steps take well under a second after the channel's sweep.

## Not validated

No IBIS model of a real device is in the tests (vendors' models carry their own terms of
use); the package's R, L and C and the supply-referenced tables are implemented as
specified but not compared with a reference simulator.
