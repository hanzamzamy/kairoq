# Viewer

`kairoq view <input>` (or `kairoq.serve(board)` in Python) starts a local HTTP server
and opens the viewer in the browser. The server listens on `127.0.0.1:8150` by default
(`--host`, `--port`; port 0 picks a free port). The viewer files are compiled into the
program; no network access is needed.

## Use

- Layers: each copper layer, dielectric layer, solder mask and the plated hole barrels can
  be shown or hidden. Layers appear in stackup order.
- Colour: by layer, or by net.
- Net: type a net name, or click copper or a hole barrel to select its net. The selected net
  is highlighted on every layer and the others are dimmed; the panel gives its pad and hole
  counts.
- Pointer: the line at the bottom gives the layer, net and board coordinates under the
  pointer.
- Section: a cutting plane normal to $x$ or $y$, moved with the slider; its position is
  given in board coordinates.
- View: top, bottom and oblique views; the vertical scale enlarges $z$ by 1, 5, 10 or 20 to
  make thin layers visible. The default is true scale.
- Import warnings are listed at the bottom of the panel.
- With a DC solution (`kairoq dc ... --view`): voltage, current density or power density
  on each solved layer, on a linear or logarithmic colour scale. Voltages are scaled per
  net, and the legend gives the range of the selected net. A click gives the values at the
  point; the panel lists the port voltages and currents, the dissipation and the estimated
  error.
- With a thermal solution (`kairoq thermal ... --view`): the temperature on each copper
  layer, with the sources' and heat sinks' results in the panel. For a transient, a time
  control selects one of the fields kept during the run, and a plot shows the highest
  temperature and the sources' temperatures over time with the selected time marked.
- With EMI estimates (`kairoq emi ... --view`): the sheet current of a loop per ampere of
  loop current on each layer of its mesh (a selector chooses the loop when there are
  several), the return-path findings drawn as outlines on their layers (red for crossings
  and splits, orange for partial gaps), and the loops' inductance and largest field and the
  switch nodes' capacitance in the panel.
- With a signal-integrity analysis (`kairoq si ... --view`): the analysed tracks coloured
  by their impedance, or by the differential impedance of pairs, on a range set by the
  bulk of the values (short sections at pads and bends lie outside it); a table of the nets
  with their mean impedances and delays; selecting a net plots its impedance (and its
  partner's) against the distance from its first pad, with its target band.

## Interface

The viewer reads two endpoints, which other clients can use as well.

`GET /api/board` returns `board.json` content without geometry references
(`docs/canonical-model.md`), with a `display` object added:

| field | meaning |
| --- | --- |
| `origin_mm` | `[x, y]` subtracted from all mesh coordinates |
| `stackup_z_mm` | `[bottom, top]` of each stackup layer, in stackup order |
| `meshes` | keys of the available meshes |
| `outline_area_mm2` | board area |

`GET /api/mesh?key=<key>` returns a triangle mesh. Keys: `copper/<n>` (copper layer $n$),
`dielectric/<i>` (stackup layer $i$), `mask/top`, `mask/bottom`, `drills` (barrels of the
plated holes). The body is little endian:

| bytes | content |
| --- | --- |
| 4 | vertex count $V$ (uint32) |
| 4 | index count $I$ (uint32) |
| $12V$ | positions, float32 $x, y, z$ in mm relative to `origin_mm` |
| $4I$ | triangle vertex indices, uint32 |
| $4V$ | vertex ids, uint32: the net id for copper and holes, 0 otherwise |

Curves in display meshes use five times the import arc tolerance.

With a DC solution, `/api/board` also contains `dc` (solved layers, value ranges, port
results, dissipation, estimated error, hole currents), and `GET /api/field?layer=<n>`
returns the solution on copper layer $n$, little endian:

| bytes | content |
| --- | --- |
| 4 | vertex count $V$ (uint32) |
| 4 | triangle count $T$ (uint32) |
| $12V$ | positions, float32, mm relative to `origin_mm`, $z$ at the top of the copper |
| $12T$ | triangle vertex indices, uint32 |
| $4V$ | voltage, float32, V |
| $4T$ | current density, float32, A/mm² |
| $4T$ | power density, float32, W/mm² |
| $4T$ | net id, uint32 |

With a thermal solution, `/api/board` contains `thermal` (layers, temperature range, source
and heat sink results and, for a transient, the time series and the times of the kept
fields), `GET /api/thermal?layer=<n>` returns the final temperature on copper layer $n$
(positions and indices as above, then the temperature per vertex as float32 in degrees
Celsius and a copper flag per triangle as uint32), and
`GET /api/thermal_values?layer=<n>&snapshot=<i>` the temperature per vertex of kept field
$i$ on the same mesh (float32).

With EMI estimates, `/api/board` contains `emi` (loops with their inductance, moment,
layers and largest sheet current; return-path findings with their outlines; switch-node
capacitances), and `GET /api/loop?loop=<i>&layer=<n>` returns the sheet current of loop $i$
on layer $n$: vertex and triangle counts, positions and indices as above, then the
magnitude of the sheet current per triangle per ampere of loop current (float32, A/mm).
