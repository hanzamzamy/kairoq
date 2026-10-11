# Importers

Three importers produce the canonical model (`docs/canonical-model.md`): one for a set of
fabrication files, one for IPC-2581 and one for ODB++. They read only what printed circuit
board design software exports in standard formats; nothing tool-specific is required.

## Fabrication file set

Input: files, directories (their files, not recursively) or zip archives, in any
combination.

| data | formats read |
| --- | --- |
| copper, solder mask, board profile | Gerber RS-274X with X2 attributes; RS-274X without attributes |
| drilled holes and slots | Excellon (IPC-NC-349 dialects), Gerber X2 drill files |
| stackup, file roles | Gerber job file (`.gbrjob`), or a stackup YAML file |
| nets | X2 object attributes; IPC-D-356(A) netlist; copper connectivity |
| component placement | pick-and-place CSV |

### File roles

The role of each Gerber file comes from, in order: its `.FileFunction` attribute, the
job file entry for its name, and the Protel file extension (`.gtl` top copper, `.gbl`
bottom copper, `.g<n>` inner copper layer $n+1$, `.gts`/`.gbs` solder mask, `.gko`,
`.gm1`, `.gml` profile). Other files are recognised by content: Excellon by its `M48`
header, IPC-D-356 by its parameter and test records, placement CSV by its header row.
A file named `*stackup*.yaml` is read as a stackup. Unrecognised files are listed in the
warnings.

### Gerber

The reader implements the Gerber Layer Format Specification (Ucamco, revision 2023.08):
standard and macro apertures with all macro primitives (including the deprecated 2, 6 and
22), arithmetic expressions and variables, linear and circular interpolation in multi- and
single-quadrant mode, regions, dark and clear polarity, load mirroring, rotation and
scaling, step and repeat, aperture blocks, and the attributes `TF`, `TA`, `TO` and `TD`,
also in the comment form `G04 #@! ...` that KiCad writes for compatibility. Deprecated
statements that change the image (`IPNEG`, `OF`, `SF`, `MI`, `IR`) are reported and
ignored. Incremental coordinates and trailing-zero omission are supported.

Strokes with a circular aperture become capsules (arcs: offset polylines) of the exact
width. Strokes with other apertures become the swept convex hull of the aperture; this is
exact for rectangles along straight lines and is reported otherwise. Objects are merged per
net in file order, so that clear polarity removes only what was drawn before it.

### Nets

- When every non-empty copper file carries `.N` attributes, the nets are those attributes.
  The net name `N/C` (not connected) is net 0. Holes take the net of the copper at their
  centre on the layers they span; holes touching copper of two nets are reported.
- Otherwise the nets are derived from the copper: each connected piece of copper on a layer
  is a node, plated holes join the nodes at their centres, and each connected group is a
  net. An IPC-D-356 netlist names the groups through its test points; disagreements
  (shorts, missing copper) are reported. Groups without a netlist point are named `N$<k>`.
  Generators write the netlist in an origin of their own (KiCad: another origin than its
  Gerber files): the netlist is placed by the translation, with or without a mirror of $y$,
  under which most of its holes land on plated drills, and the placement is reported. The
  names are the netlist's, truncated to its 14 characters.
- With X2 nets and an IPC-D-356 netlist, the netlist is compared with the pad nets.
  The fixed-width fields of IPC-D-356 truncate designators to 6 and net names to 14
  characters; the comparison accounts for this, including KiCad's `#<n>` suffix on
  shortened names.

### Pads

Pads are the flashes and regions with a `.P` attribute (component, pin, pin function). KiCad
marks the inner-layer pads of through-hole pins with `.C` only; they take the pin name of
the outer-layer pad of the same component at the same position. Without X2 attributes no
pads are known.

### Drills

Excellon files are preferred when both Excellon and Gerber drill files are given. Plating
and layer span come from the `.FileFunction` (`Plated,1,4,PTH`), or from `PTH`/`NPTH` in
the file name; hole functions (via, component, mechanical) come from the tool's
`.AperFunction` comment. Slots are read in both the `G85` and the routed (`M15`/`M16`) form.

### Board outline

The outline is built from the centre lines of the profile file's strokes: segments and arcs
are joined end to end (ends closer than $2\ \mu\mathrm{m}$ meet) into closed loops, and
nested loops alternate between board and cut-out. Open chains are reported. Without a
profile, the outline is the bounding box of the copper.

### Components

Positions come from the placement CSV. KiCad, Altium and similar column names are
recognised (`Ref`/`Designator`, `PosX`/`Center-X(mm)`, `Rot`/`Rotation`, `Side`/`Layer`);
units follow the column name or value suffix (`mm`, `mil`, `in`). Without a placement file,
each component is placed at the mean position of its pads.

### Drill and placement origin

Drill and placement files may be written relative to an auxiliary origin. When none of
their points lies on the board and all do relative to the outline's lower-left corner, the
usual auxiliary origin, they are moved there, with a warning that states the offset.

### Stackup YAML

```yaml
layers:            # top to bottom
  - {name: F.Cu,  type: copper,     thickness_mm: 0.035}
  - {name: core,  type: dielectric, thickness_mm: 1.51, material: FR-4, epsilon_r: 4.5}
  - {name: B.Cu,  type: copper,     thickness_mm: 0.035}
```

Optional fields: `material`, `epsilon_r`, `loss_tangent`, `conductivity_s_per_m`. The YAML
stackup takes precedence over the job file.

## IPC-2581

Input: an IPC-2581 XML file (revisions B and C), or a zip archive containing one. Element
names are matched without namespace prefixes.

| IPC-2581 content | model content |
| --- | --- |
| `CadHeader` `units` | MILLIMETER, INCH or MICRON |
| `Layer`, `Stackup/StackupGroup/StackupLayer`, `Spec` | stackup with material, permittivity, loss tangent, conductivity |
| `Step/Profile` (with `Cutout`) | outline |
| `LayerFeature` of conductor layers: `Set` (`net`, `polarity`), `Features`, `Pad` with `PinRef` | copper per net, pads |
| `LayerFeature` of solder mask layers | mask openings |
| `LayerFeature` of `DRILL` and `ROUT` layers: `Hole`, `SlotCavity`, `Span` | drills with span, plating and net |
| `Component` | components |

Shapes: `Circle`, `RectCenter`, `RectRound`, `RectCham`, `Oval`, `Ellipse`, `Diamond`,
`Octagon`, `Hexagon`, `Triangle`, `Donut` (round), `Contour` with `Cutout`, `Polygon` (filled
or hollow), `Polyline`, `Line`, `Arc`, `UserSpecial`, and references to the standard and user
dictionaries, with `LineDesc`/`FillDesc` and their references. Placement applies `Xform`
(scale, mirror, rotation counter-clockwise, offset) and then `Location`. Unsupported shapes
are reported. Copper layers of negative polarity are not supported and are reported.

The rotation of a component on the bottom side is given in the mirrored frame; it is
converted to the board's frame as $(540 - r) \bmod 360$ degrees. This reading of the
standard is checked against KiCad's exports only (the three importers agree on every
component); other generators may differ. Component rotation is not used by the analyses.

Plated holes without a net take the net of the copper at their centre, as for the
fabrication set. `platingStatus` `VIA` gives via holes, `PLATED` component holes,
`NONPLATED` holes of unknown function.

## ODB++

Input: an ODB++ job as a directory, a zip archive, or a tar archive,
plain or gzip-compressed (`.tgz`, `.tar.gz`); members compressed with Unix compress
(`.z`) are read. The job's root is where `matrix/matrix` lies; the first step of the matrix
is read. Names are matched without regard to case.

| ODB++ content | model content |
| --- | --- |
| `matrix/matrix` rows of context `BOARD` | layer order; `SIGNAL`, `POWER_GROUND` and `MIXED` layers are copper |
| layer `attrlist`: `.layer_dielectric`, `.copper_weight`, `.dielectric_constant`, `.loss_tangent`, `.material` | stackup with thickness, permittivity and loss tangent |
| `steps/<step>/profile` | outline |
| `features` of copper layers | copper |
| `features` of `SOLDER_MASK` layers | mask openings |
| `features` of `DRILL` layers, their `.drill` attribute, the rows' start and end layers | drills with span, plating and function |
| `steps/<step>/eda/data`: `NET`, `SNT`, `FID` | nets of copper features, pads with component and pin |
| `comp_+_top` and `comp_+_bot` `components`: `CMP`, `PRP Value` | components |

Feature records: lines (`L`), pads (`P`), arcs (`A`) and surfaces (`S`, with island and hole
contours of segments and arcs); text and barcodes are not copper the analyses use and are
skipped. Standard symbols: round, square, rectangle (with rounded or chamfered corners),
oval, diamond, octagon, round and square donut, ellipse, horizontal and vertical hexagon,
triangle; user symbols from `symbols/`. Unsupported symbols are reported and drawn as their
nominal circle. Pads apply mirroring and rotation (clockwise in ODB++, counter-clockwise in
the model); orientations 0 to 7 and the arbitrary angles of 8 and 9 are supported.

A copper feature belongs to the net that `eda/data` assigns to it; a feature assigned to a
component pin is a pad. Plated holes take the net of the copper at their centre. Without
dielectric thicknesses the stackup is an even 1.6 mm FR-4 board, with a warning. Copper
layers of negative polarity and rout layers are reported and not read.

ODB++ lists the parts that have pads. IPC-2581 also lists footprints without pads (logos,
graphics), and a placement file only the parts marked for assembly.

## Tracks

All three importers keep the drawn tracks besides the merged copper: Gerber draws (`D01`)
with a circular aperture, IPC-2581 `Line`, `Arc` and `Polyline` features with round line
ends, and ODB++ line and arc records with a round symbol. Each track takes the net of the
copper at its middle. The signal-integrity analysis cuts across them (`docs/solvers/si.md`).

## Agreement between the importers

A board exported by the same design software in the different forms gives the same model
to within the discretisation of curves. `docs/validation/importers.md` and
`docs/validation/importers.md` give the measured agreement on the test boards.
