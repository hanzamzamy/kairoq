# Validation: crosstalk

## Scope

The crosstalk analysis (`docs/solvers/si.md`, Crosstalk; `docs/solvers/networks.md`,
Coupling and Step responses) is checked against:

1. the exact theory of the TEM directional coupler: coupled stripline in one dielectric;
2. the cancellation of far-end crosstalk in a homogeneous dielectric;
3. the weak-coupling estimates of near- and far-end crosstalk, against the waveforms of
   a channel with matched ends;
4. a real board's pairs.

The tests are in `tests/python/test_crosstalk.py`.

## Coupled lines

Synthetic coupled lines, tracks of 0.2 mm 0.2 mm apart, lossless dielectric of
$\varepsilon_r = 4.5$ (the theory is for TEM lines without loss).

| case | quantity | reference | result | deviation |
| --- | --- | --- | --- | --- |
| stripline between planes 0.2 mm above and below, 7 mm, ports of $\sqrt{Z_e Z_o}$ | largest coupling $\lvert S_{21} \rvert$ | $(Z_e - Z_o)/(Z_e + Z_o)$ of the cross-section: 0.09949 | 0.09999 | +0.5 % |
| the same | its frequency | quarter wave: 5.047 GHz | 5.05 GHz | +0.05 % |
| the same | isolated port, reflection | 0 | below 0.01 | |
| stripline, 20 mm | far-end crosstalk, 100 ps | 0 | $2.6 \times 10^{-5}$ | |
| microstrip on 0.2 mm, 20 mm, ports matched to the lines, 100 ps | near-end waveform, peak | $K_b$ = 0.0884 | 0.0890 | +0.7 % |
| the same | far-end waveform, peak | $1.02\, K_f \ell / t_r$ = −0.0553 | −0.0553 | +0.1 % |

- **The directional coupler** (Oliver 1954) is exact for coupled TEM lines in one
  dielectric with ports of $\sqrt{Z_e Z_o}$: the coupled port receives
  $(Z_e - Z_o)/(Z_e + Z_o)$ where the coupled length is a quarter wave, the isolated port
  nothing, and the ports are matched. It checks the coupled sections of a channel, their
  matrix exponential and the assembly of two nets into one section. $Z_e$ and $Z_o$ come
  from a cross-section built separately; its difference from the board's cut is within
  the deviation.
- **Far-end crosstalk** in a homogeneous dielectric cancels exactly
  ($L_m/L_s = C_m/C_s$); the scan finds $2.6 \times 10^{-5}$ of the step.
- **Weak coupling.** For a step into one microstrip with all ends matched, the near end of
  the other sees $K_b$ while the step travels and its reflection returns (the coupled
  run's round trip, 243 ps, outlasts the rise), the far end the derivative of the step
  scaled by $K_f \ell$; for the Gaussian step, whose largest slope is
  $2.563/(\sqrt{2\pi}\, t_r)$, that is $1.02\, K_f \ell / t_r$. The scan's coefficients and
  the channel's waveforms agree within 0.7 %. With the ends at 50 Ω instead (the lines are
  63 Ω) the far end is −6.3 % instead of −5.5 %: the estimate assumes matched ends.

## A real board

The HDMI data pairs of `cm5_minima`: the scan finds coupling between neighbouring pairs
and none attributed to a pair's own partner. The pairs' differential impedance, with
quiet lines now included out to three times the distance to the shield (instead of two),
is unchanged (`docs/validation/networks.md`, within 0.02 Ω).

## Limits

- The scan finds aggressors on the track's own layer only; coupling between layers
  (broadside) is not scanned, and in channels not modelled.
- The scan's crosstalk is that of weak coupling with matched ends, summed over coupled
  runs as if they were one; the channel gives the waveforms of the actual geometry and
  ends.
