# Test data

## KiCad demonstration boards (`kicad/`)

Three demonstration boards distributed with the KiCad source, exported by
`tools/export_kicad_demos.sh` with kicad-cli 10.0.6. Each directory holds:

| path | content |
| --- | --- |
| `gerber/` | Gerber X2 copper, solder mask and profile; Excellon drill files; Gerber job file; IPC-D-356 netlist; placement CSV |
| `gerber-drill/` | the drill data as Gerber X2 drill files |
| `<name>.xml` or `<name>.ipc2581.zip` | IPC-2581 revision C, plain or zip-compressed |
| `<name>.odb.zip` | ODB++ job, zip-compressed, millimetres with six decimals |
| `SOURCE.txt` | origin, export and licence |

| board | copper layers | size | features |
| --- | --- | --- | --- |
| `ecc83` | 2 | 52 × 46 mm | through-hole parts, one-sided routing |
| `pic_programmer` | 2 | 160 × 99 mm | zones, vias, macro apertures, custom pads |
| `cm5_minima` | 6 | 61 × 64 mm | inner planes, 444 vias, plated slots, fine-pitch parts |

All exports of a board use the zone fill stored in the board file, so that they describe
the same copper.

Licence: the KiCad demonstration files are licensed under the Creative Commons
Attribution-ShareAlike 4.0 International licence (KiCad source, `LICENSE.README`); these
exports are adaptations distributed under the same licence. `cm5_minima` also carries the
licence file of its demo directory (CERN-OHL-S-2.0), copied as `LICENSE.txt`.
Attribution: KiCad project and the authors of the demonstration boards,
<https://gitlab.com/kicad/code/kicad>, `demos/`; `cm5_minima` (CM5 MINIMA REV3) is by
Pierluigi Colangeli. The complete source of each board is its directory under `demos/` at
the tag 10.0.6.

## Synthetic boards

`tests/python/synthetic.py` writes fabrication sets with analytically known geometry
(tracks, via chains, planes with clear polarity) when the tests run; they are not stored.

## Measured boards (`gerber2ems/`)

Two slices of the Antmicro Signal Integrity Test Board with one-port VNA measurements, from
antmicro/gerber2ems (Apache License 2.0, `gerber2ems/LICENSE`; origin and the one file added
here in `gerber2ems/SOURCE.txt`): a 50 Ω microstrip with a short (`stub_short`) and a long
(`stub_long`) open stub, `vna.csv` holding $S_{11}$ and the input impedance at port SP1 from
0.2 to 6 GHz with port SP2 loaded. The drill and placement files are written relative to the
board's lower-left corner.
