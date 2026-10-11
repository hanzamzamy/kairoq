# Validation: transmission lines

## Scope

The signal-integrity analysis (`docs/solvers/si.md`) is checked for:

1. the cross-section solver against exact and published closed forms: coaxial line,
   stripline, microstrip, edge-coupled stripline;
2. the loss models: dielectric loss and effective permittivity of a homogeneous line,
   conductor loss against an exact case and against an independent finite difference, the
   dielectric and roughness models at their defining points;
3. convergence and cost;
4. the board analysis: tracks from the three importers, synthetic boards, and the
   differential pairs of a real design against its net classes;
5. microstrip dispersion against closed forms and full-wave solutions.

The tests are in `tests/python/test_si.py` and, for the tracks, `tests/python/test_kicad.py`.

## Cross-sections against closed forms

Strips of a thousandth of their width in thickness against the zero-thickness forms;
target error 0.02.

| case | quantity | reference | result | deviation |
| --- | --- | --- | --- | --- |
| coaxial line, $a = 0.5$, $b = 1.75$ mm, $\varepsilon_r = 2.25$ | capacitance | $2\pi\varepsilon/\ln(b/a)$: 99.918 pF/m | 99.926 pF/m | $+7.9 \times 10^{-5}$ |
| the same | external inductance | $\mu_0 \ln(b/a)/2\pi$: 250.553 nH/m | 250.533 nH/m | $-7.9 \times 10^{-5}$ |
| stripline $w = 0.2$, $b = 0.5$ mm, $\varepsilon_r = 4$ | impedance | conformal mapping: 56.463 Ω | 56.321 Ω | $-0.25\ \%$ |
| microstrip $w/h = 0.5$, $\varepsilon_r = 4.3$ | impedance, effective permittivity | Hammerstad and Jensen: 96.404 Ω, 2.9901 | 96.198 Ω, 2.9890 | $-0.21\ \%$, $-0.04\ \%$ |
| microstrip $w/h = 1.5$ | the same | 58.154 Ω, 3.1985 | 58.099 Ω, 3.1965 | $-0.10\ \%$, $-0.06\ \%$ |
| microstrip $w/h = 5$ | the same | 26.167 Ω, 3.5595 | 26.126 Ω, 3.5642 | $-0.16\ \%$, $+0.13\ \%$ |
| edge-coupled stripline $w = s = 0.15$, $b = 0.5$ mm | odd, even mode impedance | Cohn, exact: 53.930, 74.634 Ω | 53.831, 74.602 Ω | $-0.18\ \%$, $-0.04\ \%$ |

- **Thickness.** The closed forms are for strips of zero thickness. A finite thickness $t$
  widens a strip electrically by $\Delta w = (t/\pi)(1 + \ln(4\pi w/t))$ (Wheeler), which
  lowers its impedance. On the stripline, with $t$ from 1.6 to 0.1 µm, the deviation
  follows this correction (a slope of $-0.51$ against $\Delta w/w$); extrapolated to zero
  thickness it is $-0.08\ \%$. The tolerances of the tests (0.4 %) allow for the
  thickness of a thousandth of the width.
- **Reference accuracy.** Hammerstad and Jensen state 0.01 % for the impedance in air and
  0.2 % for the effective permittivity; the conformal-mapping and Cohn forms are exact.

## Losses

| case | quantity | reference | result | deviation |
| --- | --- | --- | --- | --- |
| coaxial line as above, copper | resistance at 1 GHz (skin depth 2.1 µm) | $R_s (1/a + 1/b) / 2\pi$: 3.3765 Ω/m | 3.3760 Ω/m | $-1.4 \times 10^{-4}$ |
| the same | resistance at DC | $1/\sigma A$ of both conductors | | within $2 \times 10^{-3}$ |
| stripline as above, loss tangent 0.02 | effective permittivity | $\varepsilon_r$ | 4.0000 | $2 \times 10^{-15}$ |
| the same | dielectric loss at 1 GHz | $\pi f \sqrt{\varepsilon_r} \tan\delta / c_0$: 3.6409 dB/m | 3.6409 dB/m | $10^{-14}$ |
| microstrip $w = h = 0.2$, $t = 0.035$ mm, air | resistance at 100 GHz | Wheeler's rule by finite difference (faces moved by ±2 µm): 312.63 Ω/m | 312.79 Ω/m | $+5.2 \times 10^{-4}$ |

- The first form of the conductor loss integrated the surface charge at the nodes. At the
  corners of rectangular conductors the charge density is singular, and the resistance of
  the 35 µm microstrip at 10 GHz still changed by 1.8 % between the finest meshes. The
  domain form of the shape derivative replaced it: the same resistance changes by 0.1 %
  between the same meshes, and it agrees with the independent finite difference to
  $5 \times 10^{-4}$.
- The dielectric model reproduces the given permittivity and loss tangent at their
  frequency (to $10^{-3}$ with the mesh), its permittivity falls with frequency, and the
  loss tangent stays within 15 % from 1 MHz to 10 GHz. The roughness factors equal their
  definitions at the defining points (Hammerstad where the skin depth equals the RMS
  height, Huray at a skin depth of 0.02 µm) to $2 \times 10^{-3}$.

## Microstrip dispersion

| case | quantity | reference | deviation |
| --- | --- | --- | --- |
| microstrip 0.35 mm on 0.2 mm, $\varepsilon_r$ 4.5, 0.1 to 40 GHz | effective permittivity | scikit-rf's Kirschning-Jansen formulas with the cross-section's static value | $2.2 \times 10^{-4}$ at most |
| the same | impedance | scikit-rf's Jansen-Kirschning formula | $5 \times 10^{-5}$ at most |
| the microstrip of `docs/validation/rf.md`, 20 mm, ports at its ends | group delay of $S_{21}$, 1–10 GHz | openEMS, cells of 0.025 mm: 122.43 ps | 121.06 ps, −1.1 % (without dispersion 120.67 ps, −1.4 %) |
| the same | phase of $S_{21}$ at 20 GHz | openEMS: −884.9° | −873.4° (without dispersion −866.9°) |

- The formulas are implemented as published and agree with an independent
  implementation; the remaining difference is the cross-section's internal inductance,
  which its effective permittivity includes and the formulas' static input does not.
- Against the full-wave solution the dispersion closes about a third of the phase
  difference at 20 GHz ($f h = 4$ GHz mm). The remaining 11.5° is not attributed; the
  strip's thickness, which the formulas do not have, and the full-wave model's lumped
  ports are candidates.
- On `cm5_minima`'s PCIe pair (dielectric about 0.1 mm) the group delay changes by 0.2 ps.

## Convergence and cost

A microstrip of 0.2 mm on 0.2 mm, 35 µm thick, with a plane 10 mm wide, on a 20-thread x86-64
computer:

| target error | triangles | time | impedance | resistance at 10 GHz |
| --- | --- | --- | --- | --- |
| 0.05 | 40 000 | 1.6 s | 66.840 Ω | 99.28 Ω/m |
| 0.02 | 242 000 | 8.5 s | 66.886 Ω | 98.77 Ω/m |
| 0.01 | 300 000 (budget) | 30 s | 66.892 Ω | 98.86 Ω/m |

The impedance converges from below (the capacitance from above), 0.08 % low at 0.05. On
the board the default target is 0.08: on the HDMI pairs of `cm5_minima` it gives
impedances 0.17 % below those at 0.05, eight times faster.

## Board analysis

- **Tracks.** All three importers give the same tracks on the three KiCad boards (59,
  1580 and 2084, of them 200 arcs on `cm5_minima`), with ends within 2 µm (the Gerber
  coordinates' resolution).
- **Synthetic microstrip.** Every cut along a straight track is one cross-section; the
  delay is the length times the delay per length; distances run from the first pad.
- **Synthetic pair** ($w = s = h = 0.2$ mm, with mask): its pair is found by name, and its
  differential impedance on the board, 104.73 Ω, agrees with a cross-section of the same
  geometry built by hand, 104.92 Ω, to 0.19 %.
- **Pair over a slot** (1 mm across a plane): a "no reference" finding at the slot,
  1.0 ± 0.3 mm long, on each member; no change of reference.

### A real design

`cm5_minima` (six layers, the stackup and permittivities of its IPC-2581 export), all 36
nets of its pairs: 2026 cuts, 1471 distinct cross-sections. The design assigns its
pairs to net classes of 100 Ω (Ethernet, HDMI, MIPI D-PHY) and 90 Ω (USB, PCIe):

| pairs | nets | coupled length | differential impedance, mean | against the class | samples within 10 % |
| --- | --- | --- | --- | --- | --- |
| Ethernet | 8 | 332 mm | 97.8 Ω | $-2.2\ \%$ | 84 % |
| HDMI | 8 | 125 mm | 93.2 Ω | $-6.8\ \%$ | 50 % |
| MIPI D-PHY | 10 | 333 mm | 99.6 Ω | $-0.4\ \%$ | 89 % |
| USB | 4 | 65 mm | 87.4 Ω | $-2.9\ \%$ | 59 % |
| PCIe | 6 | 319 mm | 91.3 Ω | $+1.4\ \%$ | 94 % |

The means agree with the classes within 7 %; the classes' geometry itself agrees with the
stackup within 2.7 % (`docs/validation/networks.md`, "Cross-sections of channels"), and the
HDMI pairs, the furthest, run with a ground strip between their members. The samples
outside 10 % lie on short
sections at connectors, fan-outs and pads, where the pair spreads or meets wider copper;
on HDMI, whose routes are short and end at a connector, they are half of the samples.
The findings also show where the Ethernet pairs change layer through vias between the
connector J1 and the protection device U402, without a plane under them on either side
for about a millimetre. The classes are the designer's targets for the fabricator's stackup, not
measurements; the agreement shows the analysis consistent with the design's own intent.

## Defects found

- **Inside test of the mesher.** In the segment grid that decides whether a point lies
  inside a region, a crossing on the edge of a grid cell could be counted in no cell: the
  cell bounds were computed otherwise than the cells in which the segments are registered.
  A region was then taken as inside a hole, and a conductor of a cross-section was meshed
  over. The crossing is now counted in the cell found as the segments were registered.
  The same test serves the DC, thermal and EMI meshes.
- **Tracks of IPC-2581.** KiCad writes tracks as `Line` features inside `UserSpecial`;
  they are read there as well.

## Performance

Cross-sections of a track over a plane: 0.5 to 1.5 s at the default target; the board
analysis solves the distinct cross-sections in parallel (`cm5_minima`, 1471 of them: 58 s
on 20 threads). Cuts through fan-outs with dozens of copper pieces take up
to 7 s.
In the Python tests the analysis of `cm5_minima`'s pairs takes 73 s.

Reference pieces of a layer join across gaps below a tenth of their distance to the
track (`reference_merge`; a gap $g$ at distance $D$ changes the fields by about $(g/D)^2$). On
the USB and PCIe pairs of `cm5_minima` this joins 749 cross-sections into 729 distinct ones,
halves the conductors of the slowest (from 66 to 146, to 24 to 62) and reduces the
cross-sections' time from 621 to 563 s of processor time (33 to 30 s); the mean
impedances change by 0.006 % at most, the differential ones by 0.0007 %. The slowest
cross-sections are still those through fan-outs, whose copper lies near the track.
