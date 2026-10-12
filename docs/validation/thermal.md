# Validation: steady thermal and electrothermal coupling

## Scope

The thermal solver (`docs/solvers/thermal.md`) is checked for:

1. exactness and convergence against closed-form solutions;
2. the conservation of energy, alone and coupled with the DC solver;
3. the layered approximation against an independent three-dimensional solution, and where
   it breaks down, with air or mask in the gaps of thick outer copper, and on a region of a
   real board against Elmer;
4. the error estimate where the heat leaves through the thickness;
5. trace temperature rise against the IPC-2152 data;
6. time and memory.

The tests are in `tests/python/test_thermal.py`; they write their compared quantities to
`build/validation/report.md`. The three-dimensional comparison is `tools/ref3d.py`, the
IPC-2152 comparison `tools/ipc2152.py` (logs in `logs/`, not tracked).

Material values unless stated otherwise: copper $k = 385$ W/(m K), dielectric
$k_{xy} = 0.8$, $k_z = 0.3$ W/(m K), solder mask $k = 0.25$ W/(m K) and 10 µm thick,
copper 35 µm.

## Closed-form cases

| case | quantity | reference | relative error |
| --- | --- | --- | --- |
| slab 20 × 20 mm, uniform 2 W on the top copper, $h = 10$ on both faces | top and bottom temperature | series resistances | $4 \times 10^{-11}$ |
| same | heat out − heat in | 0 | $4 \times 10^{-11}$ of the input |
| same, radiation only, $\varepsilon = 0.9$ | top temperature | two-node Stefan-Boltzmann balance | $1 \times 10^{-13}$ |
| fin 5 × 60 mm held at 75 °C at one end, $h = 10$ | heat from the fixed end | $G w \theta m \tanh(m \ell)$ + pad faces | $3.8 \times 10^{-5}$ |
| same | tip temperature | $\theta / \cosh(m \ell)$ | $8.7 \times 10^{-6}$ |
| copper discs joined by one via, rims at 125 and 25 °C | thermal resistance | $\ln(r_2/r_1)/(\pi k t)$ + barrel | $-2.5 \times 10^{-4}$ |
| strip 5 × 400 mm carrying 10 A, coupled | temperature in the middle | local balance with $\rho(T)$ | $-5.8 \times 10^{-4}$ |
| `cm5_minima`, 3 A and 2 W, coupled | (heat out − heat in) / heat in | 0 | $2.5 \times 10^{-9}$ |

The slab is one-dimensional, which the sheets represent exactly. In the fin, the closed
form lumps the cross-section at one temperature; the test gives the dielectric negligible
in-plane conductivity so that the closed form and the model describe the same problem
(with the dielectric's in-plane conduction, the sheet model lies slightly below the lumped
closed form, because that conduction takes place at a lower temperature than the
copper's). Elements of at most 0.5 mm keep the error of the lumped surface term, of order
$(h_e m)^2/12$, below $10^{-4}$. The error estimate does not include that term: with the
default 2 mm elements the fin's heat is 0.05 % high, and a smaller target error does not
change it. In the electrothermal strip the temperature in the middle satisfies

$$\theta = \frac{a\,(1 + \alpha (T_a - 20))}{1 - a \alpha}, \qquad a = \frac{I^2 \rho_{20}}{w^2 t H},$$

with $H$ the conductance per area of both faces.

## Order of convergence

As for the DC solver (`docs/validation/dc.md`, "Order of convergence"): the thermal
resistance between the rims of an annulus of copper (radii 1 and 4.5 mm, copper only) on
meshes of 0.2, 0.1 and 0.05 mm converges at the observed order 1.87, against 1.945 expected
for the inner rim's polygon of 70 sides; on the finest mesh it is $1.2 \times 10^{-4}$ below
$\ln(r_2/r_1)/(2\pi k t)$. On the fin of the closed-form cases, a one-dimensional field on
unstructured triangles, the error falls faster than $h^2$ but irregularly (successive
ratios between 4 and 13), so no order is stated for it.

## Layered model against a three-dimensional solution

`tools/ref3d.py` solves steady conduction on a graded tensor-product grid in which every
layer is resolved through its thickness, with anisotropic materials, convection through the
mask on both faces and adiabatic sides. It is written independently of the library. The
board is 30 × 30 mm, $h = 10$ W/(m² K) on both faces, 35 µm copper unless stated, and
0.5 W enters the copper of a square pad or a trace on the top layer (on the middle layer for
the inner cases).

The reference itself was checked in three ways:

- on the one-dimensional slab, against the series resistances ($-2 \times 10^{-7}$);
- per solution, by the residual (below $10^{-6}$) and by the heat balance of the faces
  (within $10^{-6}$ of the input);
- for grid convergence. The edges of an isothermal pad make the heat flux singular, and the
  solution converges at first order in the cell size near the source: for the 0.5 mm pad
  the mean rise is 951.9, 908.2, 886.8 and 876.0 K at 50, 25, 12.5 and 6.25 µm. Each case is
  therefore solved at 25 and 12.5 µm next to the source (in the plane and through the 0.1 mm
  of dielectric beside it) and extrapolated as $2T_{h/2} - T_h$; the difference between the
  two resolutions is given as its uncertainty.

What fills the gaps of an outer copper layer changes the temperature of a small pad: the
heat leaves an isothermal pad mostly at its edges, where the gap material touches it. A
real outer layer has a thin solder mask on the laminate and air above it. The reference is
run at the two ends: air in the gaps (not conducting; the face is cooled on the copper and
on the dielectric below the gap) and mask to the full height of the copper. The layered
model does not distinguish the two. Gaps of an inner layer are filled with dielectric in
both.

Mean temperature rise of the source copper (K), and the error of the layered model at the
default slicing (sheets) and with one sheet per copper layer (single), against air / mask:

| case | 3D, air | 3D, mask | sheets | error, sheets | error, single |
| --- | --- | --- | --- | --- | --- |
| 0.5 mm square over a plane, 1.5 mm | 865.4 ± 2.5 % | 781.8 ± 1.5 % | 809.5 | −6.5 % / +3.5 % | −67.2 % / −63.7 % |
| 1 mm square over a plane, 1.5 mm | 423.0 ± 1.3 % | 398.6 ± 0.8 % | 407.2 | −3.7 % / +2.1 % | −52.5 % / −49.6 % |
| 2 mm square over a plane, 1.5 mm | 204.3 ± 0.6 % | 197.6 ± 0.4 % | 200.4 | −1.9 % / +1.4 % | −36.5 % / −34.3 % |
| 5 mm square over a plane, 1.5 mm | 79.5 ± 0.2 % | 78.4 ± 0.1 % | 79.0 | −0.6 % / +0.7 % | −17.0 % / −15.9 % |
| 0.5 mm square over a plane, 0.2 mm | 628.1 ± 3.1 % | 559.0 ± 1.9 % | 590.2 | −6.0 % / +5.6 % | −30.1 % / −21.4 % |
| 1 mm square over a plane, 0.2 mm | 248.2 ± 1.6 % | 232.8 ± 1.0 % | 240.1 | −3.3 % / +3.1 % | −18.1 % / −12.7 % |
| 2 mm square over a plane, 0.2 mm | 103.0 ± 0.7 % | 100.1 ± 0.4 % | 101.8 | −1.1 % / +1.7 % | −8.5 % / −5.8 % |
| trace 0.25 × 20 mm, 1.5 mm | 106.0 ± 0.8 % | 102.0 ± 0.4 % | 103.5 | −2.3 % / +1.5 % | −28.4 % / −25.6 % |
| trace 1 × 20 mm, 1.5 mm | 83.1 ± 0.3 % | 81.6 ± 0.2 % | 82.4 | −0.9 % / +0.9 % | −16.8 % / −15.3 % |
| trace 5 × 20 mm, 1.5 mm | 56.4 ± 0.1 % | 55.9 ± 0.1 % | 56.2 | −0.4 % / +0.5 % | −7.8 % / −7.0 % |
| 2 mm square, 10 oz copper, over a plane, 1.5 mm | 200.3 ± 0.7 % | 165.8 ± 0.2 % | 196.5 | −1.9 % / +18.5 % | −37.7 % / −24.8 % |
| trace 1 × 20 mm, 3 oz, inner layer | 65.8 ± 0.1 % | 65.1 ± 0.1 % | 65.6 | −0.4 % / +0.7 % | −5.2 % / −4.1 % |
| trace 1 × 20 mm, 10 oz, inner layer | 59.9 ± 0.0 % | 58.4 ± 0.0 % | 59.7 | −0.4 % / +2.2 % | −3.6 % / −1.1 % |

- The sheets lie between the two references in every case. With 35 µm copper they are
  within 2.3 % of both for traces and for sources of 2 mm and more (5 mm: within 0.7 %); for
  a 0.5 mm pad they are 6 to 6.5 % below the air reference and 3.5 to 5.6 % above the mask
  reference, which differ from each other by 11 to 12 %.
- One sheet per copper layer underestimates the temperature by two thirds for a 0.5 mm pad
  and by 16 % for a 5 mm one over 1.5 mm of dielectric; this is where the plain layered
  approximation breaks down, and why the dielectrics are sliced.
- Thick copper on an outer layer is where the gap material matters most: for a 2 mm pad of
  10 oz copper the two references differ by 21 %, and the sheets are 2 % below the air
  reference. Thick inner copper, embedded in resin, is within 2.2 % of both.

`test_against_3d_reference` keeps three of these cases (the 1 mm trace and the 1 and 2 mm
pads) in the test suite with the stored reference values, at the default target error: the
layered result must lie between the two references, widened by their uncertainties.

The slicing was chosen on the 0.5 and 1 mm pads, against the mask reference at 50 µm cells,
where the trend rather than the level is meaningful: a first
slice of 0.05, 0.025 and 0.0125 mm raised the computed rise by 4 and 2 % in turn, and a
largest slice of 0.1 instead of 0.2 mm by 1 %. The default first slice of 0.025 mm costs
about 12 % more triangles than 0.05 mm.

## Thick outer copper

Between the copper of an outer layer the model has air by default; the three-dimensional
references above differ by what fills those gaps, most for thick copper. With the
gaps filled with mask to the copper's thickness (`outer_gaps: mask`), the mean rise of the
source against the reference with the same fill, 0.5 W, 1.5 mm over a plane:

| source | air | mask |
| --- | --- | --- |
| 0.5 mm pad | 805.5 K (−6.9 %) | 745.7 K (−4.6 %) |
| 1 mm pad | 404.1 K (−4.5 %) | 386.3 K (−3.1 %) |
| 2 mm pad | 200.2 K (−2.0 %) | 195.1 K (−1.3 %) |
| 2 mm pad, 10 oz | 198.4 K (−0.9 %) | 164.3 K (−0.9 %) |

The 10 oz pad was 18.5 % above the mask reference with air in the gaps; with the matching
fill it is within 1 % of both. The difference lies in the conduction of the gap material
next to the copper, not in the copper's own vertical resistance, which at 400 W/(m K) is
negligible: sheets through the thickness of the copper were not needed.

## Real board against Elmer

Elmer, an open-source finite element package, solves the same problem in three dimensions
on a region of `cm5_minima` (`docs/solvers/elmer.md`, "Thermal reference"). The region is
22 × 22 mm (x from 114 to 136 mm, y from −83 to −61 mm; `--crop`), cut out in both models
with adiabatic edges. It holds the 2 A buck converter U702 (a 2 × 2 mm DFN with an exposed
pad, 0.5 W) and its inductor L701 (5 × 5 mm, 0.2 W), over six copper layers, mostly planes,
and their vias. Both faces have $h = 10$ W/(m² K), radiation is off, and the gaps of the
outer layers are air. One hole cut by the region's edge is left out of both.

Temperature rise above 25 °C (K):

| model | size | U702, mean | U702, largest | L701, mean | bottom layer, largest |
| --- | --- | --- | --- | --- | --- |
| Elmer, 0.2 mm | 1.49 M nodes, 2.73 M elements | 118.30 | 127.45 | 97.92 | 94.89 |
| Elmer, 0.1 mm | 1.95 M nodes, 3.61 M elements | 118.49 | 127.51 | 97.95 | 94.97 |
| Elmer, 0.05 mm | 4.12 M nodes, 7.80 M elements | 118.77 | 127.62 | 97.99 | 95.17 |
| layered, target 0.1 (default) | 0.20 M triangles | 117.13 | 128.89 | 97.63 | 93.91 |
| layered, target 0.05 | 2.60 M triangles | 119.02 | 130.00 | 98.11 | 95.98 |
| layered, target 0.03 | 4.96 M triangles (near the limit) | 119.22 | 129.95 | 98.13 | 95.93 |
| layered, target 0.1, in-plane estimate only | 0.17 M triangles | 111.12 | 121.61 | 96.36 | 93.41 |

- **The reference.** Elmer's size is the element size at the copper's edges, holes and the
  sources, and half of it the largest element layer through a dielectric. Halving it changes
  the results by 0.3 % at most; the copper's outlines, not the size, set most of the
  triangulation (0.22 to 0.36 M triangles). The same export on the 2 mm source of the
  section above gives a mean rise of 189.3, 197.3 and 201.8 K at 0.2, 0.1 and 0.05 mm; the
  first-order extrapolation, 206.3 K, is within 1 % of the finite-volume reference
  (204.3 ± 0.6 %).
- **Converged, the layered model agrees with the 3D solution** within 0.4 % at the sources'
  mean temperatures and 0.8 % on the bottom layer. The largest temperature of U702, a
  value at the hottest edge of a pad, is 1.8 % higher.
- **At the default target error** the layered model is 1.4 % low at U702 and 0.4 % at
  L701, in 6 s; at 0.05 it is within 0.2 % of its converged value, in 43 s.
- **The estimate's vertical term** (`docs/solvers/thermal.md`, "Adaptive refinement")
  decides this. Refined by the in-plane estimate alone, the mesh stops at 0.17 M triangles
  with U702 6.4 % low while the estimate reads 0.086: the error sits in the thin
  dielectric slices under the package, which conduct little in the plane and so hardly
  enter that estimate. On the whole board (`examples/cm5_minima_thermal.yaml`, 2 W in U1
  and 3 A coupled), U1's mean rise is 0.3 % below a solution at 5 M triangles with the
  vertical term, and 6.3 % below without it, at 2.1 M triangles and 88 s instead of
  1.3 M and 31 s. On the 3D comparison boards above, at a target of 0.03, the term changes
  the results by 0.25 % at most.
- Elmer took 14, 10 and 28 minutes (one thread), the layered model 6 to 107 s.

## Error estimate where the heat leaves through the thickness

The estimate covers the in-plane conduction and the lumped vertical coupling between the
sheets (`docs/solvers/thermal.md`, "Adaptive refinement"), not the exchange with the faces
and the heat sinks. Two cases where most of the heat leaves through the thickness, 0.5 W in
a 2 mm pad over a plane: on 0.2 mm of dielectric with the bottom face on a sink through
0.5 mm of a 3 W/(m K) interface, the mean rise is 67.73 K at the default target and
67.10 K at 0.03 (4.9 million triangles), 0.9 % apart; on 1.5 mm without a sink, 198.87 and
200.07 K, 0.6 % apart. A gradient-recovery estimate of the exchange with the faces and the
sinks refined both cases to about 4 million triangles at the default target, without a
gain that the references could show; it is not used.

## Trace temperature rise against IPC-2152

`tools/ipc2152.py` reproduces the IPC-2152 test fixture (IPC-TM-650 2.5.4.1A): a 12 in
(304.8 mm) trace in the middle of a 365.8 × 142.2 mm board of 0.070 in (1.78 mm) polyimide
without planes, cooled on both faces; the temperature rise is the mean over the copper of
the middle 6 in, as the fixture measures it by resistance. External traces are on the top
of a two-layer board, internal ones on the middle layer of a three-layer board. The trace
is heated by its own current with $\rho_{20} = 1.7\ \mu\Omega\,\mathrm{cm}$ and
$\alpha = 0.00393\ \mathrm{K^{-1}}$ (electrothermal coupling), $h = 10$ W/(m² K) on both
faces including radiation, ambient 25 °C.

The IPC-2152 charts enter through the curve fits of Brooks and Adam (*Trace Currents and
Temperatures Revisited*, 2015, Table 3-1 and Appendix 6), with $W$ and $T_h$ in mil and $I$
in A:

$$\Delta T_\mathrm{ext} = 215.3\, I^2 W^{-1.15} T_h^{-1}, \qquad
\Delta T_\mathrm{int} = C\, I^{a} W^{-b} T_h^{-1.52}$$

with $(C, a, b) = (200, 1.9, 1.10)$ for 1 oz, $(300.3, 2, 1.15)$ for 2 oz and
$(300, 1.9, 1.15)$ for 3 oz (225 at 10 mil). For each case the current is the one at which
the fit gives a rise of 20 or 60 °C.

The conductivity of the fixture's polyimide is stated in IPC-2152 only; the model is run at
0.2 and 0.6 W/(m K), isotropic, as the ends of the range for polyimide-glass laminates.

| layer | oz | width (mil) | current (A) | IPC rise (°C) | model, $k = 0.2$ | model, $k = 0.6$ |
| --- | --- | --- | --- | --- | --- | --- |
| external | 1 | 10 | 1.34 | 20 | 36.6 (+83 %) | 17.0 (−15 %) |
| external | 1 | 10 | 2.33 | 60 | 153.2 (+155 %) | 58.9 (−2 %) |
| external | 1 | 50 | 3.39 | 20 | 33.9 (+69 %) | 17.6 (−12 %) |
| external | 1 | 50 | 5.88 | 60 | 137.6 (+129 %) | 61.1 (+2 %) |
| external | 1 | 200 | 7.53 | 20 | 26.1 (+31 %) | 16.4 (−18 %) |
| external | 1 | 200 | 13.04 | 60 | 98.1 (+63 %) | 56.5 (−6 %) |
| external | 2 | 10 | 1.90 | 20 | 36.9 (+85 %) | 17.3 (−14 %) |
| external | 2 | 10 | 3.29 | 60 | 155.0 (+158 %) | 59.7 (0 %) |
| external | 2 | 50 | 4.80 | 20 | 34.2 (+71 %) | 17.9 (−11 %) |
| external | 2 | 50 | 8.31 | 60 | 139.4 (+132 %) | 62.1 (+4 %) |
| external | 2 | 200 | 10.65 | 20 | 26.2 (+31 %) | 16.6 (−17 %) |
| external | 2 | 200 | 18.44 | 60 | 98.1 (+63 %) | 57.1 (−5 %) |
| external | 3 | 10 | 2.33 | 20 | 37.3 (+86 %) | 17.5 (−13 %) |
| external | 3 | 10 | 4.03 | 60 | 157.0 (+162 %) | 60.6 (+1 %) |
| external | 3 | 50 | 5.88 | 20 | 34.6 (+73 %) | 18.1 (−9 %) |
| external | 3 | 50 | 10.18 | 60 | 141.3 (+135 %) | 63.2 (+5 %) |
| external | 3 | 200 | 13.04 | 20 | 26.2 (+31 %) | 16.8 (−16 %) |
| external | 3 | 200 | 22.59 | 60 | 98.1 (+64 %) | 57.6 (−4 %) |
| internal | 1 | 10 | 1.46 | 20 | 33.2 (+66 %) | 16.8 (−16 %) |
| internal | 1 | 10 | 2.60 | 60 | 146.2 (+144 %) | 62.4 (+4 %) |
| internal | 1 | 50 | 3.70 | 20 | 35.3 (+76 %) | 19.4 (−3 %) |
| internal | 1 | 50 | 6.60 | 60 | 159.3 (+166 %) | 73.7 (+23 %) |
| internal | 1 | 200 | 8.27 | 20 | 30.1 (+50 %) | 19.2 (−4 %) |
| internal | 1 | 200 | 14.74 | 60 | 127.7 (+113 %) | 72.7 (+21 %) |
| internal | 2 | 10 | 2.10 | 20 | 34.3 (+71 %) | 17.5 (−13 %) |
| internal | 2 | 10 | 3.63 | 60 | 139.6 (+133 %) | 60.7 (+1 %) |
| internal | 2 | 50 | 5.29 | 20 | 36.4 (+82 %) | 20.1 (0 %) |
| internal | 2 | 50 | 9.16 | 60 | 151.5 (+152 %) | 71.3 (+19 %) |
| internal | 2 | 200 | 11.73 | 20 | 30.3 (+52 %) | 19.5 (−2 %) |
| internal | 2 | 200 | 20.32 | 60 | 118.3 (+97 %) | 68.8 (+15 %) |
| internal | 3 | 10 | 3.51 | 20 | 72.0 (+260 %) | 35.1 (+75 %) |
| internal | 3 | 10 | 6.25 | 60 | 578.9 (+865 %) | 158.1 (+163 %) |
| internal | 3 | 50 | 7.99 | 20 | 60.3 (+201 %) | 32.4 (+62 %) |
| internal | 3 | 50 | 14.24 | 60 | 386.0 (+543 %) | 141.2 (+135 %) |
| internal | 3 | 200 | 18.48 | 20 | 54.2 (+171 %) | 34.2 (+71 %) |
| internal | 3 | 200 | 32.96 | 60 | 310.2 (+417 %) | 151.8 (+153 %) |

- The IPC data lie inside the band of the board conductivity for all external traces, and
  for internal traces of 1 and 2 oz at a 20 °C rise. At $k = 0.6$ W/(m K) the model is 9
  to 18 % below the fit at 20 °C and within 6 % at 60 °C (external); the internal traces
  of 1 and 2 oz are within −16 to +23 %, and above the fit at 60 °C.
- The dependence on the copper weight follows the fit for external traces within 3.5 %
  (nearly the same rise at the same width for 1, 2 and 3 oz), and the dependence on the
  width within about ±6 % across 10 to 200 mil: the scaling is reproduced, the level
  depends on the board's conductivity.
- The model's rise grows faster with current than the fit's. This is the constant heat
  transfer coefficient: natural convection and radiation increase with the surface
  temperature, and Brooks and Adam needed 11 to 14 W/(m² K) instead of 10 to match the
  larger rises.
- The 3 oz internal fit is an outlier: for the same rise it lets an internal trace carry
  36 % more current than an external one, against 9 to 10 % for 1 and 2 oz, which the
  physics of the fixture does not explain; Brooks and Adam attribute the spread of the
  internal constants to the control of the tests. The model is 60 to 165 % above it.
- The 10 mil traces at $k = 0.2$ W/(m K) reach the limit of 5 M triangles with an estimated
  error of 0.12 instead of 0.1.

**Lateral profile.** Jouppi (*The Value of IPC-2152*, Figure 10) gives the temperature
next to a 10 mil, 1 oz external trace at a 10 °C rise. Rise in °C:

| distance from the trace | Jouppi | model, $k = 0.2$ | $k = 0.4$ | $k = 0.6$ |
| --- | --- | --- | --- | --- |
| 0.00 in | 10.0 | 15.35 | 9.87 | 7.78 |
| 0.33 in | 4.0 | 1.49 | 1.86 | 1.97 |
| 0.60 in | 1.4 | 0.29 | 0.58 | 0.75 |
| 0.90 in | 0.5 | 0.05 | 0.16 | 0.26 |

The decay length beyond 0.33 in, 7 mm, agrees at $k = 0.6$ W/(m K); nearer the trace the
model falls off faster (at 0.33 in, 25 % of the rise at the trace, against 40 %).

**Planes.** Brooks and Adam simulated (TRM) a 6 in, 1 oz, 200 mil external trace at 15 A,
$h = 11$, with large pads at the ends: rise at the middle 94.7 °C without a plane, 47.9 °C
with a 1 oz plane on the bottom, 33.7 °C with one 10 mil below the trace. The pad size is
not given; the model uses 1 in squares.

| configuration | TRM | model, $k = 0.2$ | model, $k = 0.6$ |
| --- | --- | --- | --- |
| no plane | 94.7 | 126.0 | 72.9 |
| plane on the bottom | 47.9 (0.51) | 48.6 (0.39) | 28.0 (0.38) |
| plane 10 mil below | 33.7 (0.36) | 25.4 (0.20) | 20.2 (0.28) |

(ratio to the case without a plane in parentheses). The model predicts a larger benefit of
the planes than TRM. Both are models, with board size, pads and material not fully stated
in the source; this comparison is qualitative.

## Performance

A 20-thread x86-64 computer, Release build, default settings.

| case | size | triangles | unknowns | time | peak memory |
| --- | --- | --- | --- | --- | --- |
| `cm5_minima`, 2 W and 3 A coupled (`examples/cm5_minima_thermal.yaml`) | 61 × 64 mm, 6 layers, 24 sheets | 2.11 M | 1.01 M | 80 to 89 s | 4.9 GB |
| 3D comparison boards, target error 0.03 | 30 × 30 mm, 2 or 3 layers | 2.2 to 4.9 M | | 37 to 116 s | |
| IPC-2152 fixture, coupled | 366 × 142 mm, 0.070 in | 0.85 to 5.0 M | | 13 to 219 s | |

On `cm5_minima` the time divides into meshing 4.7 s, assembly 7.1 s, the linear solutions
59.2 s and the DC solutions 3.8 s, over three coupling iterations; the result converged to
an estimated error of 0.083 in four passes, with the temperature changing by 0.0002 K in
the last iteration. The times of the comparison runs were measured with other runs on the
same machine and are upper bounds.

The electrothermal coupling is a fixed-point iteration whose plain form slows as the
feedback of the resistivity on the heat grows: 3 to 6 iterations for rises up to 75 K in
the IPC-2152 cases, 6 to 9 up to 200 K, and 21 for a rise of 566 K, close to thermal
runaway. Anderson acceleration (depth 3) now combines the last iterations; with it the
IPC-2152 cases take 3 to 4, 4 to 5 and 5 to 8 iterations in the same ranges (8 at 579 K). On the strip of
the electrothermal test, at a feedback gain $a\alpha$ (the iteration's contraction) of
0.5, 0.8, 0.9 and 0.95, the plain iteration takes 14, 42, more than 60 and more than 60
iterations, the accelerated one 5, 8, 8 and 9, on the same temperature. Above a gain of 1
there is no steady state: the plain iteration grows, and the accelerated one can settle
below $20 - 1/\alpha$ °C, where the resistivity of the linear law vanishes; both are
reported as thermal runaway.

In the test suite (`tests/python`) the coupled `cm5_minima` case takes 76 s.

## Defects found and corrected during this validation

- The estimator stalled at the copper edges, where the in-plane conductance jumps by two
  orders of magnitude; the gradient is now recovered per material.
- Radiation was iterated as a fixed point and needed up to 11 coupled iterations; it is now
  linearised by Newton's method, and the mesh of the last adaptive pass is kept for the
  coupling iterations (3 iterations on `cm5_minima`).
- One sheet per copper layer underestimated the temperature of small sources by up to
  two thirds; dielectrics are now resolved by graded sheets.
- The three-dimensional reference first filled the gaps of the outer copper with
  dielectric, used 50 µm cells at the source and did not check its own convergence; the
  first made thick outer copper appear 39 % too hot in the layered model, the second put
  the reference 10 % high for a 0.5 mm pad. It now brackets the gap material, extrapolates
  from two resolutions and checks the residual and the heat balance of every solution.
- The resin between the copper of inner layers was not represented; it now adds its
  in-plane conductance.
