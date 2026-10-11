# Validation: EMI estimates

## Scope

The EMI estimates (`docs/solvers/emi.md`) are checked for:

1. loop inductance against closed forms: parallel strips, microstrip, a rectangular loop,
   and a loop over a plane by image theory;
2. the magnetic moment: its value, its independence of the origin, and its cancellation
   over a plane;
3. the harmonics of the loop current and the radiated field;
4. the return-path check on boards with a slot, a split plane, merged antipads and a
   single antipad;
5. the switch-node capacitance with and without shielding, and with fringing by the method
   of moments;
6. the openEMS export of a loop, and the inductance and far field against full-wave
   simulations of it;
7. the resistance and inductance of loops over frequency, and the radiated field of the
   loop currents;
8. the return current through stitches across a slot, and the convergence of the loop
   inductance with the triangle budget;
9. a real board, and time and memory.

The tests are in `tests/python/test_emi.py`; they write their compared quantities to
`build/validation/report.md`. Kernels (the triangle integrals against brute-force
quadrature, the filament formulas, the dense factorisation) are checked in
`tests/cpp/test_peec.cpp`, the layered Green's function of the method of moments in
`tests/cpp/test_layered_green.cpp`. The full-wave comparisons are in
`tests/python/test_openems.py`; they run when `KAIROQ_OPENEMS_PYTHON` names a Python
interpreter with openEMS (`CONTRIBUTING.md`).

## Loop inductance

| case | quantity | reference | result | deviation |
| --- | --- | --- | --- | --- |
| two parallel strips, $w = s = 1$ mm, 8 and 16 mm long | inductance per length | $\mu_0 K(k)/K(k')$, $k = s/(s+2w)$: 0.8038 nH/mm | 0.8157 nH/mm | +1.5 % |
| microstrip, $w = 1$ mm, $h = 1.535$ mm, 8 and 16 mm long | inductance per length | $Z_0/c_0$, Hammerstad and Jensen in air: 0.5042 nH/mm | 0.5140 nH/mm | +1.9 % |
| rectangular loop 20 × 10 mm of a 0.5 mm strip | inductance | partial inductances of four sides, radius $w/4$: 46.96 nH | 47.77 nH | +1.7 % |
| the same loop 1.535 mm over a plane | inductance | $L_\text{free} - M(\text{loop}, \text{image at } 2h)$: 36.77 nH | 36.93 nH | +0.4 % |

- **References.** The per-length inductances of parallel strips (coplanar strips, by
  conformal mapping; Simons, *Coplanar Waveguide Circuits*) and of microstrip (Hammerstad and
  Jensen, accurate to 0.01 % in air) are exact for thin perfect conductors. They are compared
  with the difference of two lengths, which cancels the ends, each solved at the same density
  of about 200 triangles per millimetre of length. The rectangular loop's reference sums the
  partial self and mutual inductances of its four sides (Grover), each a round conductor of
  the equivalent radius of a thin strip in the high-frequency limit, $w/4$; it neglects the
  corners and is accurate to about a percent itself. The image reference subtracts the
  mutual inductance between the loop and its mirror image, for an infinite plane.
- **Convergence.** Every solution is an upper bound of the model's value. On the parallel
  strips the per-length inductance fell from +4.8 % to +2.0 % as the density doubled from
  100 to 200 triangles per millimetre, on the microstrip from +3.4 % to +2.4 %. The rate is
  slow because the high-frequency current grows as $d^{-1/2}$ towards a copper edge; rows of
  thin elements along the edges and adaptive refinement are what make it usable.
- **The adaptive passes** agree within 3 % of each other in the tests; the inductance of
  the previous pass is reported with each result as an indication of the remaining error.
  Successive meshes are not nested, so the values need not fall strictly from pass to pass.
- **Refinement of the loop.** The rectangular loop's inductance falls by about 0.16 nH
  (0.33 %) at each doubling of the triangles from 834 to 4504 (47.94, 47.78, 47.61,
  47.46 nH; the last in 1.8 s). No power law that would justify extrapolation has set in,
  and the solution has not converged at the finest mesh; part of the 3 to 4 % between this
  estimate and openEMS (section "Full-wave comparison of the loop estimates") is this
  discretisation.

## Magnetic moment

| case | quantity | reference | result |
| --- | --- | --- | --- |
| rectangular loop 20 × 10 mm | moment per ampere | area of the centre line, 200 mm² | 199.45 mm² (−0.3 %) |
| the same loop 30 mm farther from the origin | moment per ampere | the same | within $2 \times 10^{-3}$ |
| the same loop over a plane reaching 20 mm beyond it | moment left | 0 for an infinite plane | 6.9 % of the free loop's |

The moment lies slightly under the centre line's area because the high-frequency current
shifts to the inner edge of the strip. It does not depend on the origin because the computed
current satisfies the current law everywhere.

This check found a defect in the first implementation. The DC solver ties the triangles of
a port into one node, which is right for an equipotential contact. In the inductance model
the same ties let the current vanish at one point of a pad and reappear at another without
storing magnetic energy. The moment then depended on the loop's position (212 against
250 mm² for the same loop moved by 30 mm). Elements and plated holes now inject their
current over the pad or hole by area, and every triangle keeps its own current law.

Over a plane, the residual moment depends on how far the plane reaches beyond the loop:
37, 20 and 14 mm² for 5, 12 and 20 mm, while the inductance changed by less than 0.3 %. The
moment, and with it the radiated estimate, is the more sensitive quantity. The window now
grows until the moment settles or the window holds the board: from margins of 5 and 20 mm
alike it reaches the whole board here (11.8 mm², 5.9 % of the free loop's), in up to
four solutions instead of one.

## Harmonics and radiated field

| case | quantity | reference | result |
| --- | --- | --- | --- |
| trapezoid 2 A, 1 MHz, $d = 0.3$, 20 ns edges | harmonic amplitudes to 50 MHz | $2 A d\, \lvert\mathrm{sinc}(n \pi d)\rvert\, \lvert\mathrm{sinc}(n \pi t_r/T)\rvert$ (Paul) | within $3 \times 10^{-16}$ A |
| the same with 50 ns fall time | harmonic amplitudes to 20 MHz | discrete Fourier transform of the sampled waveform | within $10^{-6}$ A |
| loop 20 × 10 mm, 3 A, 2 MHz, 5 ns edges, at 3 m | field of each harmonic to 300 MHz | $1.316 \times 10^{-14} f^2 A I / r$ (Ott) with the computed moment | within $6.5 \times 10^{-4}$ |

The last deviation is the rounding of Ott's constant; the field is the exact small-loop far
field $\eta_0 k^2 m I / (4 \pi r)$.

## Return paths

On a two-layer board with a 0.3 mm trace over a plane 0.2 mm below:

| plane | expected | found |
| --- | --- | --- |
| solid | nothing | nothing |
| a slot 1 mm wide across the trace | one crossing | one crossing at the slot, span 1.011 mm |
| split into GND and PWR by a 0.5 mm gap | one split between GND and PWR | one split, span 0.505 mm, nets GND and PWR |
| a row of overlapping antipads across the trace | one crossing | one crossing |
| one antipad under the trace's edge | one partial gap | one partial gap |

The spans are 1 % above the drawn gaps: the trace's width is estimated from its area and
perimeter, which includes the round ends.

## Switch-node coupling

| case | quantity | reference | deviation |
| --- | --- | --- | --- |
| patch 10 × 8 mm on the top, plate under the board, 0.5 mm gap of permittivity 4 | capacitance | layered parallel plates: 1.5065 pF | $10^{-16}$ |
| the same with a plane under half of the patch | capacitance | half of it | $10^{-16}$ |
| the same | overlap with the plane | $\varepsilon_0 \varepsilon_r A / t$ | $10^{-16}$ |

These check the bookkeeping (stack, shielding, overlap). The parallel-plate model itself
leaves out fringing, which adds to the capacitance of small areas.

## Real board

`examples/cm5_minima_emi.yaml` on the `cm5_minima` carrier (six layers):

- **Loop.** The input loop of the converter U702 (input capacitor C706, the converter between
  its input and ground pins), with 0.2 and 0.4 nH inside the parts: 1.58 nH, of which
  0.98 nH in the copper; moment 0.32 mm² per ampere, small because the inner ground plane
  carries the image current; with 2 A at 1.5 MHz, duty 0.66 and 2 ns edges, the largest field
  at 3 m is −7.3 dBµV/m (262 MHz). The window grows to 20 mm around the corridor, where the
  moment changes by 3 % from the last doubling; the copper of the top three layers in it,
  nine vias: 4635 triangles, above the budget of 4000, which the outlines alone need; 6.8 s.
- **Return paths.** 20 crossings or splits and 11 partial gaps over the board, mostly
  traces on the third layer over gaps in the +5 V plane of the fourth.
- **Switch node.** The converter's switch node, on the top layer, faces a plate under the
  board through four inner planes: its capacitance is negligible ($10^{-7}$ pF, slivers).

Developing this example found three defects, corrected:

- edge rows placed per boundary segment produced so many vertices on curved outlines that
  the matrix could not be allocated; they are placed by arc length;
- the planes' antipads and the via rings, all modelled, needed 22 000 triangles before any
  refinement; holes far from the loop are closed, outlines simplified to the local size,
  vias limited to those near the loop and layers to those near the loop's own;
- the window and the element size followed straight lines between the pads; they follow the
  DC corridor of the loop current, which also sizes the window around the actual copper
  path.

## openEMS export: the model

On the rectangular ring loop above (20 × 10 mm centre line, 0.5 mm strip, closed
by a part across a 1 mm gap):

- the exported sheet is the ring's copper with the outline simplified as in the loop
  model: within 0.045 mm (0.15 of the element size at the loop) and 0.3 % in area;
- the element runs from pad centre to pad centre, its axis one radius above the copper;
- the dielectric slabs add up to the board's copper and dielectric thickness;
- the estimate written with the model equals the loop solution;
- a plane with three square openings is written as non-overlapping hole-free pieces whose
  union is the window less the openings, to $10^{-6}$ in area.

## Full-wave comparison of the loop estimates

The ring loop free and 1.535 mm over a plane reaching 25 mm beyond it (the cases above,
element radius 0.25 mm; the exported model held the plane within the window's 10 mm
margin, and section "Radiated field of the loop currents" repeats the case over the whole
plane), exported
with the default settings (excitation to
1 GHz, cells of 0.1 mm at the copper, FR-4 of permittivity 4.5 under the window) and
simulated with openEMS 0.0.36. The estimate uses 3000 triangles.

| case | quantity | estimate | openEMS | estimate / openEMS − 1 |
| --- | --- | --- | --- | --- |
| free | inductance (openEMS at 100 MHz) | 47.66 nH | 45.65 nH | +4.4 % |
| free | far field per ampere at 3 m, 100 MHz | 8.757 mV/m | 8.674 mV/m | +1.0 % |
| free | the same, 250 MHz | 54.73 mV/m | 56.27 mV/m | −2.7 % |
| free | the same, 500 MHz | 218.9 mV/m | 257.4 mV/m | −15 % |
| over a plane | inductance (openEMS at 100 MHz) | 36.89 nH | 35.66 nH | +3.5 % |
| over a plane | far field per ampere at 3 m, 100 MHz | 1.026 mV/m | 1.014 mV/m | +1.2 % |
| over a plane | the same, 250 MHz | 6.414 mV/m | 6.523 mV/m | −1.7 % |
| over a plane | the same, 500 MHz | 25.65 mV/m | 28.82 mV/m | −11 % |

- **Radiated field.** Where the loop is small against the wavelength, the small-loop far
  field of the computed moment agrees with the full-wave field within 3 %, also over the
  plane, where the plane cancels 88 % of the moment and only the remainder radiates: the
  estimate's moment includes the image current in the finite plane correctly. At 500 MHz
  the loop's perimeter (61 mm) is a tenth of the free-space wavelength and about a sixth of
  the wavelength along a strip on the board; the current
  is no longer uniform around the loop and the field per ampere at the port exceeds the
  small-loop value by 11 to 15 %. The estimate is a small-loop estimate and should be read
  as such above about a tenth of the wavelength in perimeter.
- **Inductance.** openEMS gives a flat inductance below 120 MHz (within a ripple of about
  1 %), rising towards the loop's self-resonance above. The ripple is not the truncation
  of the time signal: ending the free loop's run at −50 dB of energy instead of −40 dB took
  3.5 times the time steps (60.8 instead of 17.5 ns) and left it at 0.90 % (0.86 %), the
  inductance at 100 MHz changing
  by 0.09 %. The estimate lies 3.5 to 4.4 % above it. The estimate converges from above
  (section "Loop inductance": 1.7 % above the partial-inductance reference for the free loop, which is
  itself accurate to about 1 %); FDTD with metal edges on mesh lines widens strips by a
  fraction of a cell and gives a lower inductance as the cell grows. With cells of
  0.05 mm, openEMS gives 45.93 nH for the free loop, 0.6 % more than with 0.1 mm (26
  minutes, 239 000 time steps); extrapolated linearly in the cell size, about 46.2 nH,
  3 % below the estimate and 1.6 % below the partial-inductance reference. The far field
  changes by less than 0.6 % at 100 and 250 MHz.
- **Time.** 4.7 and 5.7 minutes on 16 threads (0.35 and 0.66 million cells, 83 000 and
  29 000 time steps), against 3.2 and 4.5 s for the estimates.

Developing the comparison found two defects of the export, corrected:

- the sheets' planes were merged with nearby dielectric boundaries in the mesh; a sheet
  off the mesh lines is ignored by openEMS. The sheets' planes and the port's lines are
  now kept exactly;
- the port was a box of the element's section. The impedance it gave behaved as the loop
  in parallel with about 290 Ω, and the inductance fell with frequency; a lossless small
  loop has neither. With the port as one line of edges on the element's axis, the real
  part of the impedance is the small radiation resistance and the inductance is flat at
  low frequency.

The energy in the domain falls quickly to −35 dB and then slowly (a run to −50 dB was
stopped at −40 dB after 90 000 steps); the end criterion is −40 dB by default.

## Loop resistance and inductance over frequency

Two coplanar strips 1 mm wide, 1 mm apart, 35 µm thick, closed at one end by a resistor;
the values per length are the difference of loops 8 and 16 mm long.

| frequency | quantity | loop solver | reference | source |
| --- | --- | --- | --- | --- |
| 1 kHz | $R'$ | 0.9865 Ω/m | 0.9852 Ω/m | $2 / (\sigma t w)$ |
| 1 GHz | $L'$ | 0.815 nH/mm | 0.804 nH/mm | conformal mapping, zero thickness |
| 1 GHz | $R'$ | 12.02 Ω/m | 12.16 Ω/m | Wheeler's rule, zero thickness, both faces |

Without the edge weighting the resistance at 1 GHz grows with refinement (12.86 Ω/m with
the weighting over $t/4$, 13.07 Ω/m over $t/10$). The model is a sheet: copper of finite
thickness carries current on the edges' side walls, and the references that include it
are higher, 15.95 Ω/m by the cross-section solver and 19.3 Ω/m by Wheeler's rule with
Gupta's thickness correction. They differ from each other by a fifth; the loop solver's
high-frequency resistance is therefore a lower bound.

On the CM5 example's input loop the resistance rises from 4.2 mΩ at 1 kHz to 15.1 mΩ at
100 MHz and the inductance falls from 2.43 to 1.60 nH, against 1.58 nH in the
high-frequency limit.

## Radiated field of the loop currents

| case | site | reference | largest deviation |
| --- | --- | --- | --- |
| loop 20 × 10 mm, 3 m, 1 to 300 MHz | free space | small loop with its near-field term | $< 10^{-3}$ |
| loop 100 × 100 mm, 300 MHz to 1 GHz | free space, far field | radiation integral of a uniform current on the centre line, 2° grid | $< 2 \times 10^{-3}$ |
| loop 10 × 10 mm 0.8 m over ground, 3 m, 1 to 300 MHz | ground, antenna 1 to 4 m | loop and image, near and far field | $< 5 \times 10^{-4}$ |

At 1 MHz the field at 3 m ($kr = 0.06$) is 16 times the far-field term alone. The large
loop is one wavelength around at 750 MHz; there the small-loop formula is 11 % high, at
1 GHz 20 %. The 5° by 10° grid of directions loses less than $2 \times 10^{-4}$ against a 1°
grid.

Developing the field found a defect in the moment. A part's conductor ends at its pad
centres, one radius above the copper, and its current enters the sheet spread over the
pad: the path was not closed. On the free loops the omission is the conductor's raised
area, about $10^{-3}$ of the moment; on the `cm5_minima` input loop, whose moment the inner plane
cancels, it was half: 0.17 against 0.32 mm² per ampere, and the radiated estimate rises by
5.5 dB to −7.3 dBµV/m. The same open ends, offset from the centroids of their pads by
micrometres, gave the current elements a net $\sum I\,\mathbf{dl}$ of
$2.9 \times 10^{-6}$ A m per ampere, which radiated as an electric dipole eleven times the
loop's field at 37.5 MHz. Joining each end to the centroid of its spread closes the current; the
moment is then independent of the origin to rounding, and the free-space field of the
`cm5_minima` loop equals the small-loop field of its moment within $10^{-3}$ up to 250 MHz.

The full-wave comparison above holds the plane only within the window's 10 mm margin
(40 × 30 mm). Here it is repeated with the window grown until the moment settles, which
reaches the whole plane (70 × 60 mm), whose larger image current halves the moment. openEMS 0.0.36, default export
settings, 1.0 million cells:

| case | quantity | estimate | openEMS | estimate / openEMS − 1 |
| --- | --- | --- | --- | --- |
| free | far field per ampere at 3 m, 100 MHz | 8.757 mV/m | 8.674 mV/m | +1.0 % |
| free | the same, 500 MHz | 218.8 mV/m | 257.4 mV/m | −15 % |
| over the whole plane | inductance (openEMS at 100 MHz) | 36.86 nH | 35.79 nH | +3.0 % |
| over the whole plane | far field per ampere at 3 m, 100 MHz | 0.5197 mV/m | 0.5047 mV/m | +3.0 % |
| over the whole plane | the same, 250 MHz | 3.259 mV/m | 3.244 mV/m | +0.5 % |
| over the whole plane | the same, 500 MHz | 13.19 mV/m | 14.26 mV/m | −7.5 % |

The field of the currents and the small-loop field of the moment differ by less than 2 %
on these loops, which are small against the wavelength: the deficit at 500 MHz is not the
field's computation but the currents': the model's are quasi-static and uniform around the
loop, while towards the self-resonance the current along the loop exceeds the current at
the port.

## Switch-node capacitance with fringing

The potential of a unit charge in a stratified medium against closed forms: a uniform
medium, $1/(4\pi\varepsilon R)$, within $10^{-6}$; a charge on a slab of FR-4 (4.4) 1.5 mm
thick over a ground plane against its image series, from 20 µm to 30 mm, within
$2 \times 10^{-4}$, also at 30 mm, where the direct term and the remainder cancel to
$10^{-3}$; reciprocity between two heights of a five-slab stack within $10^{-6}$.

| case | method of moments | reference | difference |
| --- | --- | --- | --- |
| disc 5 mm radius, 0.503 mm over a plate, vacuum | 1.8245 pF | 1.8257 pF | −0.06 % |
| the same disc on FR-4 (4.5) of 0.5 mm, masks of 10 µm (3.8), 0.5 mm air gap | 1.8463 pF | 1.8486 pF | −0.13 % |
| strip 1 mm wide on FR-4 (4.5) 0.5 mm over a plate, per length | 0.12682 pF/mm | 0.12772 pF/mm | −0.7 % |

The discs' reference is the axisymmetric field solution of the same stack with a disc
0.5 µm thick; the strip's is the cross-section solver, its conductors of 5 and 2 µm
extrapolated to zero thickness. Refining the elements at the vacuum disc's edge from 0.16
to 0.02 mm moves its capacitance from 1.8168 to 1.8244 pF, halving the difference each
time, towards 1.8255 pF; the extrapolation from two meshes gives that value from the
coarsest ones. On a 10 × 8 mm patch over a plate, through 1.5 mm of FR-4 and 0.5 mm of
insulation, fringing adds 59 % to the parallel-plate capacitance.

Developing the solver found four defects, corrected:

- a conductor in the middle of its copper layer lies microns from two interfaces, and a
  mask puts a third 10 µm away; their images varied on a scale far below the elements,
  and the capacitance came out 18 % low. Conductors now lie on the interface between the
  slabs beside them, and the images of the nearest interface are integrated in closed form;
- a height that differed from an interface by rounding was taken as a distinct one, and
  the direct term was taken in the wrong medium: heights within $10^{-12}$ m of an
  interface are now on it;
- the size field acted on outlines only, so a plane whose outline lay at the window's edge
  was meshed coarsely under the node, and its image charge was lost (the capacitance to it
  fell by a third as the window grew); vertices at the node's scale now cover the node's
  footprint on other copper;
- large elements of planes stacked 0.1 mm apart made the system lose positive
  definiteness; the largest element is now eight times the smallest layer spacing, with a
  finer rule over such pairs.

On the `cm5_minima` switch node, which four inner planes shield from the enclosure, the
solution has 2766 panels and takes 2 s; its capacitance to the inner ground plane is
3.70 pF (3.54 pF as a parallel-plate overlap). The plate's share is 0.014 pF, between
0.002 and 0.008 pF on the two meshes and falling as the window grows: there it is the
small difference of large charges, and the window's cut edges let the field reach the
plate. It is an upper estimate. On a node over a plane wider than the window, it falls from
0.040 to 0.006 pF as the margin grows from 1 to 16 mm while the capacitance to the plane
changes by 1.5 %.

The cross-section solver, used as the strip's reference, gave 0.184 instead of
0.128 pF/mm for a conductor 1 µm thick with a margin of three times the extent: its
default outline tolerance (2 µm) exceeds such a conductor. Copper is not that thin, but
the tolerance is not scaled to the geometry.

## Return current through stitches

A trace of 30 mm on 0.5 mm of FR-4 over a plane with a 1 mm slot across its path, open to
the board's edge, so that the return detours 12 mm around the slot's end; a part across
the slot beside the trace is the stitch.

| case | loop solver | openEMS, 100 MHz |
| --- | --- | --- |
| no slot | 15.06 nH | |
| slot | 21.54 nH | 19.95 nH |
| slot with a strap across it | 16.59 nH | 14.95 nH |
| inductance the strap removes | 4.95 nH | 5.00 nH (−1 %) |

The strap carries 0.86 of the loop current in the high-frequency limit. With a 100 pF
capacitor as the stitch, the solution at each frequency agrees within 1 % in inductance
and in the stitch's current with the circuit of two meshes, loop and stitch, whose
inductances the solutions with the stitch open and strapped give; the copper's resistance,
which the circuit leaves out, makes the difference at the antiresonance (2 % at 200 MHz).
Below 30 MHz the capacitor carries less than 3 % of the current, at 300 MHz 1.4 times it,
circulating through the detour, and at 1 GHz 0.89, near the strap's share. The full-wave
model of the capacitor rings: the capacitor and the detour form a tank that the port
hardly damps, and the run had not decayed to the end criterion after 175 000 time steps;
the capacitor case is therefore checked against the circuit only.

The absolute inductances lie 8 and 11 % above openEMS, more than on the ring loops above. The loop solver has not converged on this board: doubling its triangle budget
from 4000 lowers them by 0.2 nH and again by 0.3 nH at 16 000, while the reported value with
half the triangles differs by only 0.03 nH; the return current concentrates in a band of
the order of the dielectric's thickness under the 30 mm trace. At 32 000 triangles the dense
matrix of the branches no longer fits in 10 GB.

## Convergence of the loop inductance

On the slotted board above, the loop inductance falls by 0.16, then 0.31 nH as the triangle
budget doubles from 4000 to 16 000 (21.54, 21.39, 21.07 nH; with the strap 16.59, 16.41,
16.13 nH), while the inductance with half the triangles, which the solution reports as its
indication of the error, differs by only 0.03 to 0.08 nH at 4000 and 8000. The decrements
grow, so the solution is not yet in the range where an extrapolation from the passes would
hold. Two causes were examined:

- the return current in the plane flows in a band of the order of the dielectric's
  thickness under the trace, where the plane's elements stay at 0.5 to 0.9 mm: the
  gradient recovery does not mark a region that the mesh resolves too coarsely to show a
  gradient. Keeping the a-priori size there at half the layer spacing (0.29 to 0.44 mm)
  changed the inductance by less than 0.6 %;
- the decrease follows the refinement at the trace's edges, where the current of a sheet
  without thickness grows as $d^{-1/2}$ and the energy converges slowly from above.

The dense matrix of the partial inductances of all branches limited the solution to about
16 000 triangles in 10 GB. Above `dense_matrix_mb` its rows are now formed by blocks and
multiplied out (`docs/solvers/emi.md`); the solution equals the dense one to $10^{-14}$, and
the budgets beyond that limit converge:

| triangle budget | slot | with the strap | unknowns (slot) | memory (slot) | the dense matrix | time (slot) |
| --- | --- | --- | --- | --- | --- | --- |
| 16 000 | 21.072 nH | 16.129 nH | 16 748 | 1.1 GB | 2.2 GB | 20 s |
| 32 000 | 20.905 nH | 16.012 nH | 37 858 | 2.0 GB | 11.5 GB | 140 s |
| 64 000 | 20.853 nH | 15.992 nH | 70 547 | 5.3 GB | 39.8 GB | 633 s |

Towards about 20.82 and 15.98 nH, the solutions lie 4.4 and 6.9 % above openEMS (19.95 and
14.95 nH), of the order of the ring loops above, and the inductance the strap removes,
4.86 nH, 2.8 % below openEMS's 5.00 nH.

The slow convergence is not particular to the slot. A rectangular ring over a plane loses
0.7 to 2 % per doubling of the budget up to 32 000 triangles (`docs/solvers/emi.md`,
"Accuracy and cost"). On the slotted board the refinement puts 22 times more elements on
the trace at 32 000 than at 4000 triangles, along all its length and within 0.1 mm of its
edges, and refining the trace alone or the plane alone recovers only part of the change
(21.18 and 21.39 nH at 16 000, against 21.07 nH with both): what converges slowly is the
current across the trace and its return under it. At the default budget the first mesh has
the coarsest elements allowed, 1.2 mm at the loop for a 0.5 mm trace over 0.5 mm of
dielectric. Three changes did not help: weighting the refinement indicator by the element
size, rows of vertices along the edges at a fixed spacing, and estimates of the error from
the rate of change over the passes, which erred by up to 60 times. The solution therefore
reports the element size at the loop, warns when the budget sets it, and the documentation
gives the convergence against the budget, from above in every case.

## Performance

A 20-thread x86-64 computer, Release build; the loop 20 × 10 mm
over a plane, default settings (window growth and the error estimate included):

| budget | triangles | unknowns | time | peak memory |
| --- | --- | --- | --- | --- |
| 1000 | 1196 | 1560 | 0.4 s | 125 MB |
| 2000 | 1508 | 2000 | 0.9 s | 337 MB |
| 4000 | 3229 | 4528 | 2.9 s | 312 MB |
| 6000 | 4755 | 6673 | 5.7 s | 831 MB |

The time is that of all adaptive passes and window sizes. The loop basis
(`docs/solvers/emi.md`) solves the loops alone, about a third of the branches, without the
current law; against the constrained solution it replaces, one solution of the loop over
the plane takes 0.75 s at 2473 triangles and 2.4 s at 4884, against 4.2 and 57.6 s, and the
free loop 1.8 s at 4504 triangles against 45.1 s, with the same inductance and moment to
$10^{-7}$. The constrained solution grows with the cube of the unknowns (dense
factorisation); the loop basis with fewer unknowns and a smaller dense block. The return-path check of the whole `cm5_minima` board
takes 0.6 s, the switch-node coupling 0.1 s. The EMI tests run in 12 s.
