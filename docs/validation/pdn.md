# Validation: power distribution

## Scope

The plane-pair cavity and the impedance of a power distribution network (`docs/solvers/pdn.md`)
against the closed-form modal sum, the static limit and closed forms of mounting inductance.
The tests are in `tests/python/test_pdn.py` and `tests/cpp/test_cavity.cpp`.

## Plane-pair cavity and PDN impedance

| case | quantity | reference | deviation |
| --- | --- | --- | --- |
| rectangle 40 × 10 mm, two ports of 0.15 mm | first 10 modes | $(m\pi/a)^2 + (n\pi/b)^2$ | $2 \times 10^{-3}$ at most |
| the same, 1 MHz | capacitance between the planes | $\varepsilon A / h$ | below $5 \times 10^{-3}$ |
| the same, loss tangent 0.02, 0.5 to 9.5 GHz | transfer impedance | double sum of the modes | 4 % at most, near resonances |
| plane pair 40 × 30 mm, 0.2 mm, a load and a 0.1 µF capacitor (ESR, ESL, mounting), lossy plates and dielectric, 1 MHz to 4 GHz | load impedance | double sum of the modes with the same ports and branch | 2.5 % of $\lvert Z \rvert$ + 0.3 Ω |
| the same, 1 MHz | $\lvert Z \rvert$ | capacitor and plates in parallel: 1.588 Ω | −0.55 % |
| capacitor with via-in-pad, pitch 2 and 3 mm, 0.52 and 1.02 mm over the plane, radius 0.05 mm | mounting inductance | half the rectangle with its image (Grover): 1.542, 3.133 nH | +1.5 %, +0.8 % |

- The largest differences lie just below resonances, where an error of 0.1 % in a mode's
  frequency moves the impedance by tens of times as much, and at the top of the band, where
  the static sum stands in for the modes not kept.
- **Real board.** `examples/cm5_minima_pdn.yaml`: the +5 V plane of `cm5_minima` with the
  ground planes on either side, two pairs of 0.13 and 0.55 mm (20 and 24 cm², 97 and 30
  modes to 2 GHz), the CM5 module as the load, the USB-C input as the source and seven
  capacitors. Below 23 MHz the impedance stays under 0.1 Ω; above, the pairs' resonances
  (largest at 522 MHz, 24 Ω) with the capacitors' mounting inductances of 1.2 to 1.6 nH
  (pads, tracks to their vias, vias) and 6.6 nH at the USB-C connector. 12 s.

## Not validated

The PDN impedance is not compared with a full-wave solution or a measured board.
