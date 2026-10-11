# Validation: import, canonical model and viewer

## Scope

The importers read fabrication data into the canonical model, and the viewer displays it.
Their results are geometry and connectivity; validation therefore asks two questions:

1. Does an import reproduce known geometry exactly, up to the stated discretisation?
2. Do independent descriptions of the same board, exported in different formats, give the
   same model?

The checks are tests in `tests/`. Running the Python tests writes the report
`build/validation/report.md` (and `report.json`) with every compared quantity, its
reference and its tolerance.

## Analytical cases

`tests/python/synthetic.py` writes Gerber X2, Excellon, job file and IPC-D-356 files
directly, without the library, for boards whose geometry is known in closed form.

| case | quantity | reference | result |
| --- | --- | --- | --- |
| straight track between two pads, widths 0.1 to 3 mm | copper area | $2A_\mathrm{pad} + (L - 2)\,w$, or the stadium $Lw + \pi w^2/4$ when the track covers the pads | relative error ≤ $4 \times 10^{-9}$ |
| 4-layer plane with a cleared square (clear polarity) | copper area | $(s - 2m)^2 - h^2$ | exact |
| 6 mm disc | copper area | $\pi r^2$ | relative error $3 \times 10^{-8}$ |
| stackup from the job file | thickness | sum of layer thicknesses | exact |
| board profile | outline area | rectangle | exact |
| via chain without X2 attributes | nets from connectivity and IPC-D-356 | one named net | as expected |

The C++ unit tests (`tests/cpp/`) cover the geometric primitives (circle and capsule areas,
arc chord deviation, the area-preserving arc construction), every aperture template and
macro primitive including rotation and exposure, clear polarity, aperture blocks, step and
repeat, inch units with trailing-zero omission, single-quadrant arcs, the Excellon number
formats and slot forms, IPC-D-356 records, job files, placement files, stackup files, the
WKB encoding and the round trip of the model through a directory and a zip archive.

Curves are replaced by polygons whose vertices are rounded to $1\ \mathrm{nm}$. Circles and
the arcs bounding filled areas keep their area by construction; for a circle of radius
$50\ \mu\mathrm{m}$ the vertex rounding alone leaves a relative area error of up to
$2 \times 10^{-5}$.

## Cross-check of the importers

Each KiCad test board (`testdata/README.md`) is exported once as a fabrication set and once
as IPC-2581. The two importers share only the geometric primitives and the final merging
step; parsing, shape construction, transforms, net assignment and drill handling are
independent. Measured agreement (kicad-cli 10.0.6 exports):

| quantity | ecc83 (2 layers) | pic_programmer (2 layers) | cm5_minima (6 layers) |
| --- | --- | --- | --- |
| stackup thickness | equal | equal | equal |
| outline area, relative difference | 0 | 0 | $3.6 \times 10^{-6}$ |
| copper differing between the imports, largest layer (XOR area / copper area) | $1.5 \times 10^{-4}$ | $2.0 \times 10^{-4}$ | $8.1 \times 10^{-6}$ (inner layers: 0) |
| nets | identical | identical | identical |
| drills (position, diameter, plating, slot length) | 33 of 33 | 251 of 251 | 482 of 482 |
| plated hole nets: derived from copper (Gerber) against stated (IPC-2581) | 33 of 33 | 244 of 244 | 468 of 468 |
| pads (component, pin, layer, net) | 66 of 66 | 474 of 474 named | 737 of 737 |
| largest pad shape difference (XOR area / pad area) | $1.2 \times 10^{-3}$ | $1.5 \times 10^{-3}$ | $1.5 \times 10^{-3}$ |
| component positions | within 1 µm | within 1 µm | within 1 µm |
| Gerber drill files against Excellon | identical | identical | identical |
| IPC-D-356 netlist against X2 pad nets | no disagreement | no disagreement | no disagreement |

The remaining copper and pad differences are discretisation: the two exports describe
circles differently (flashed apertures against IPC-2581 primitives), so the polygons place
their vertices at different angles. The areas agree.

Two differences between KiCad's exports are expected and not counted as errors: the two
mounting pads of `pic_programmer` P3 have an empty pin name in the Gerber files and the
names `PAD0`, `PAD1` in IPC-2581; and the IPC-2581 file lists components that the
placement file omits.

During validation, a first export refilled the zones for the Gerber files only
(`--check-zones`). The cross-check showed the resulting copper differences (up to 1.4 % of
the copper of an inner layer, located at the zone clearances). The export now uses the
stored fill for both formats.

## Performance

Import times on a 20-thread x86-64 computer, Release build:

| board | fabrication set | IPC-2581 |
| --- | --- | --- |
| ecc83 | 0.01 s | 0.02 s |
| pic_programmer | 0.13 s | 0.15 s (zip) |
| cm5_minima | 0.34 s | 0.64 s (zip) |

## Viewer

The viewer was checked in a headless Chromium (SwiftShader WebGL) driven by Playwright on
`cm5_minima`: all 14 meshes load, layer visibility, colouring by net, net selection by
name and by pointer, the section plane and the vertical scale work, and the page logs no
errors. This check is manual; it is not part of the automated tests.

## ODB++ importer

The three KiCad demonstration boards (`testdata/README.md`) exported as ODB++ (KiCad 10.0.6,
ODB++ 8.1, millimetres with six decimals), compared with their IPC-2581 exports:

| board | copper layers | nets | pads | holes | components (ODB++ / IPC-2581 / placement file) | copper differing (XOR / area, worst layer) | pad shape differing (worst pad) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `ecc83` | 2 | 14 | 66 | 33 | 15 / 15 / 11 | 0 | 0 |
| `pic_programmer` | 2 | 112 | 478 | 251 | 63 / 63 / 63 | $1.1 \times 10^{-4}$ | $8.8 \times 10^{-4}$ |
| `cm5_minima` | 6 | 221 | 737 | 482 | 84 / 112 / 71 | $2.9 \times 10^{-5}$ | $2.0 \times 10^{-3}$ |

- Nets, pads with component, pin and net, holes with position, diameter, plating and net,
  the outlines (equal areas), the stackup thicknesses and the dielectric permittivities
  agree exactly. The remaining copper and pad differences are the vertices of arcs, placed
  at different angles by the two exports; the tolerances in the tests are $5 \times 10^{-4}$
  for copper and $5 \times 10^{-3}$ for pads, as for the comparison of IPC-2581 with Gerber
  above.
- Component positions, sides and rotations agree for every part. The sets differ by
  design: ODB++ lists the parts with pads, IPC-2581 also the footprints without pads (logos),
  the placement file only the parts marked for assembly.
- Import times: 0.02, 0.14 and 0.42 s, as for IPC-2581.
- The same job as a directory, as tgz archives in the POSIX and the GNU tar format with a
  root name longer than a tar header's name field, and with its matrix, netlist, top
  copper and top components compressed with Unix compress gives the same board. The
  compressor of the test is checked against `gzip -d`; the compressed copper features
  (469 kB) use every code width up to 16 bits.

Developing the comparison found three defects, corrected:

- the feature count line (`F n`) was read as a feature record, which shifted every net
  assignment of `eda/data` by one feature;
- component rotations were taken as counter-clockwise; ODB++ gives them clockwise;
- the IPC-2581 importer took the rotation of mirrored (bottom) placements as it stands; it
  is given in the mirrored frame and is now converted, so that all importers give the
  rotation counter-clockwise as seen from the top (`docs/canonical-model.md`).

## Limits of this validation

- All real test boards come from one design program (KiCad). The importers implement the
  published formats, but exports of other programs may use features the KiCad exports do
  not (negative copper layers, Gerber step and repeat in copper, IPC-2581 shapes such as
  `Moire` or `Thermal`). These are implemented from the specifications and covered by unit
  tests, not by real boards.
- Nets derived from connectivity alone are tested on a synthetic board only.
- The geometry is planar per layer. Plating thickness, etching profile and hole wall
  geometry are not represented in the model; the solvers take them from their own settings.
