# Canonical board model

Every importer produces the same board model. It holds what the analyses need from a
finished printed circuit board: the board outline, the stackup, the copper of each layer
merged per net, the pads, the drilled holes, the solder mask openings and the components.
This document specifies its file format, version 2. Any tool that writes this format can
feed the analyses; any tool that reads it can use their inputs.

## Conventions

- Lengths are millimetres, as IEEE 754 double precision numbers. The implementation works
  in integer nanometres; a reader may round coordinates to $1\ \mathrm{nm}$ without loss.
- The plane coordinates are those of the source data: $x$ to the right, $y$ upwards, as
  in Gerber and IPC-2581. Nothing is mirrored for bottom layers.
- The vertical coordinate $z$ is zero at the bottom surface of the lowest stackup layer
  and increases upwards. It is not stored: it follows from the stackup thicknesses.
- Copper layers are numbered from 1 at the top to $N$ at the bottom.
- Net 0 is "no net": copper and pads without an electrical net. Its name is empty.

## Storage

A model is a directory, or a zip archive of the same tree:

```
board.json            metadata, stackup, nets, pads, drills, components
outline.wkb           board outline
copper/L<n>.wkb       copper of layer n, one geometry per net
pads/L<n>.wkb         pad shapes on layer n, one geometry per pad
mask/top.wkb          solder mask openings, top side
mask/bottom.wkb       solder mask openings, bottom side
```

### Geometry encoding

Geometry is OGC Well-Known Binary (WKB, ISO/IEC 13249-3), little endian, as a
`MultiPolygon` (type 6). Each polygon has its outer ring first and its holes after it.
Rings are closed (the first point is repeated at the end). Outer rings are counter-clockwise
and holes clockwise. Any WKB reader decodes the files, for example `shapely.wkb.loads` in
Python or `ST_GeomFromWKB` in PostGIS.

A `.wkb` file that holds several geometries is their concatenation. `board.json` gives the
byte `offset` and `length` of each.

Curves are stored as polygons. Arcs and circles are replaced by chords no further than the
arc tolerance (`source.arc_tolerance_mm`) from the true curve. The interior vertices of arcs
that bound filled areas lie slightly outside the arc, so that each polygon has the area of
the shape it replaces; arcs of stroked tracks follow the centre line and keep the exact
track width.

## board.json

```json
{
  "format": "kairoq-board",
  "version": 2,
  "name": "CM5_MINIMA_3",
  "units": "mm",
  "source": { ... },
  "stackup": { "source": "gbrjob", "layers": [ ... ] },
  "copper_layers": [ ... ],
  "mask_layers": [ ... ],
  "outline": { "file": "outline.wkb" },
  "nets": ["", "GND", "+3V3", ...],
  "net_source": "x2",
  "components": [ ... ],
  "pads": [ ... ],
  "drills": [ ... ],
  "warnings": [ ... ]
}
```

A reader rejects a model whose `format` differs or whose `version` is higher than the
versions it knows. Fields not listed here may be added in later versions; readers ignore
unknown fields.

### source

| field | type | meaning |
| --- | --- | --- |
| `importer` | string | `fabrication` (Gerber set) or `ipc2581` |
| `generator` | string | software that wrote the source files, when they state it |
| `files` | string list | file names read |
| `arc_tolerance_mm` | number | largest distance between an arc and its chords |

### stackup

`source` is where the stackup came from: `gbrjob` (Gerber job file), `ipc2581`, `yaml`
(a stackup file supplied with the inputs) or `default` (assumed; see below).

`layers` lists the physical layers from the top of the board to the bottom:

| field | type | meaning |
| --- | --- | --- |
| `name` | string | name in the source |
| `type` | string | `copper`, `dielectric`, `soldermask` |
| `thickness_mm` | number | thickness |
| `material` | string | material name, when given |
| `epsilon_r` | number | relative permittivity, when given |
| `loss_tangent` | number | dielectric loss tangent, when given |
| `conductivity_s_per_m` | number | electrical conductivity of a copper layer, when given |
| `copper_index` | integer | copper layer number (copper layers only) |

The default stackup, used when the inputs contain none, is $N$ copper layers of
$35\ \mu\mathrm{m}$ separated by FR-4 dielectric of equal thickness, $1.6\ \mathrm{mm}$ in
total (or the board thickness of the job file). The importer adds a warning.

### copper_layers

| field | type | meaning |
| --- | --- | --- |
| `index` | integer | 1 (top) to $N$ (bottom) |
| `name` | string | name in the source |
| `side` | string | `top`, `inner`, `bottom` |
| `file` | string | geometry file |
| `nets` | list | `{"net": id, "offset": bytes, "length": bytes}`, one entry per net present, ascending id |

The geometry of a net on a layer is all its copper on that layer merged: tracks, pads, zones
and graphics. Copper of different nets does not overlap in a correct board; the importers
report overlaps as warnings.

### mask_layers

`{"side": "top" | "bottom", "file": ...}`. The geometry is the area of the mask openings.

### nets and net_source

`nets` is the list of net names; a net's id is its position. `net_source` states where the
nets came from:

| value | meaning |
| --- | --- |
| `x2` | Gerber X2 `.N` object attributes |
| `ipc2581` | IPC-2581 `Set` net attributes |
| `ipc356` | copper connectivity, named from an IPC-D-356 netlist |
| `geometric` | copper connectivity only; nets are named `N$1`, `N$2`, ... |

### pads

| field | type | meaning |
| --- | --- | --- |
| `component` | string | reference designator |
| `pin` | string | pin name; may be empty when the source gives none |
| `function` | string | pad function from the source, for example `SMDPad`, `ComponentPad` |
| `net` | integer | net id |
| `layer` | integer | copper layer |
| `position_mm` | `[x, y]` | pad reference point (flash position, or the centre of the bounding box) |
| `file`, `offset`, `length` | | pad shape in `pads/L<layer>.wkb` |

A through-hole pin has one pad entry per copper layer on which it has copper.

### drills

| field | type | meaning |
| --- | --- | --- |
| `position_mm` | `[x, y]` | hole centre, or the first end of a slot |
| `slot_end_mm` | `[x, y]` | second end of a slot (slots only) |
| `diameter_mm` | number | finished hole diameter, or slot width |
| `plated` | boolean | plated through the copper layers it spans |
| `function` | string | `via`, `component`, `mechanical` or `unknown` |
| `from_layer`, `to_layer` | integer | copper layers spanned, `from_layer` ≤ `to_layer` |
| `net` | integer | net id; 0 for non-plated holes and for plated holes on no net |

### tracks

The drawn tracks of the copper layers: strokes of a round tool, straight or circular, as
the design data gives them (version 2). Their copper is also part of the merged copper of
their layer; the tracks keep the centre line and width that the merged copper loses.

| field | type | meaning |
| --- | --- | --- |
| `layer` | integer | copper layer |
| `net` | integer | net id: the net of the copper at the track's middle |
| `width_mm` | number | width (the tool's diameter) |
| `start_mm`, `end_mm` | `[x, y]` | ends |
| `centre_mm`, `ccw` | `[x, y]`, boolean | arcs only: centre, and the direction from start to end |

A version 1 model has no tracks; a reader of version 2 reads it with none.

### components

| field | type | meaning |
| --- | --- | --- |
| `refdes` | string | reference designator |
| `value`, `package` | string | as given by the source |
| `position_mm` | `[x, y]` | placement origin |
| `rotation_deg` | number | rotation, counter-clockwise as seen from the top, for parts on either side |
| `side` | string | `top` or `bottom` |
| `position_source` | string | `placement` (pick-and-place file), `ipc2581`, or `pads` (mean of the pad positions, rotation unknown) |

### warnings

Messages from the import about data that was missing, ambiguous or inconsistent. A model
with warnings is complete as far as the inputs allow.

## Converting other data

A converter needs to write `board.json` with the fields above and the WKB files. Copper must
be merged per net and per layer; overlapping input shapes are allowed only within one net.
The Python package can be used to check the result:

```python
import kairoq
board = kairoq.open("converted-model/")
print(board.summary())
```
