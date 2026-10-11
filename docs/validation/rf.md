# Validation: RF structures

## Scope

The RF analysis (`docs/solvers/rf.md`) is checked against:

1. closed forms of lines and fixtures for the de-embedding;
2. a worked example of a textbook and a limiting case for the patch estimate;
3. the native channel solution for the full-wave S-parameters of a microstrip (openEMS);
4. the patch estimate against its full-wave solution;
5. dielectric and conductor losses, and microstrip ports, in the full-wave export.

The tests are in `tests/python/test_deembed.py` and `tests/python/test_rf.py`; the
full-wave ones run when `KAIROQ_OPENEMS_PYTHON` names a Python with openEMS (here
openEMS from the Nix environment of `docs/solvers/openems.md`).

## Results

| case | quantity | reference | result | deviation |
| --- | --- | --- | --- | --- |
| lossy 60 Ω line of 30 mm, 10 mm and 5 mm removed at its ports | S | the 15 mm line | | $7 \times 10^{-15}$ |
| lossless 60 Ω line referred to 60 Ω | $\lvert S_{11} \rvert$ | 0 | $6 \times 10^{-16}$ | |
| fixture of 2 × 4 mm of 40 Ω and 0.2 pF (several wavelengths at 40 GHz), 2x-thru removed at both ports | S | the DUT, 12 mm of 70 Ω with loss | | $3 \times 10^{-14}$ |
| patch, $\varepsilon_r$ 2.2, $h$ 1.588 mm, 10 GHz | $W$, $\varepsilon_\text{eff}$, $\Delta L$, $L$ | Balanis: 1.186 cm, 1.972, 0.081 cm, 0.906 cm | 1.185 cm, 1.9715, 0.0811 cm, 0.9053 cm | below 0.1 % |
| the same | edge resistance, 50 Ω inset | Balanis: 228.3 Ω, 0.3126 cm | 228.4 Ω, 0.3124 cm | 0.05 %, −0.06 % |
| slot of $10^{-4}$ wavelengths over a ground plane | directivity | 3 | 3.0000001 | |
| microstrip 50 Ω, 20 mm on 0.2 mm, ports at its ends | group delay, 1–10 GHz | native channel: 120.7 ps | 122.5 ps | +1.5 % |
| the same | $\lvert S_{21} \rvert$ at 1 and 10 GHz | native: −0.10, −0.71 dB | −0.60, −0.65 dB | see below |
| the same | $\lvert S_{11} \rvert$ to 10 GHz | matched line | 0.06 at most | the ports |
| edge-fed patch 37 × 28.5 mm, 1.5 mm, $\varepsilon_r$ 4.5, lossless | resonant frequency (largest input resistance, feed de-embedded) | openEMS: 2.370 GHz | estimate 2.451 GHz | +3.4 % |
| the same | edge resistance | openEMS: 287 Ω | estimate 324 Ω | +13 % |
| the same | directivity | openEMS: 7.20 dBi | estimate 6.06 dBi | −1.1 dB |
| the same | radiation efficiency | 1 (lossless) | 1.017 | +1.7 % |

- **De-embedding** agrees with the closed forms to rounding. The 2x-thru fixture is several
  wavelengths long at 40 GHz: the square root follows its branch through every multiple of
  $\pi$ of electrical length.
- **Microstrip.** The full-wave and native delays agree within 1.5 %. The full-wave loss
  is higher at low frequency and lower at high frequency: its dielectric loss is a
  conductivity set at 10 GHz (the loss tangent grows as $1/f$ below), and its copper is
  lossless. The reflection of up to 0.06 (0.13 at 20 GHz) is the lumped ports' vertical
  path, which the native channel does not have.
- **Patch.** The transmission-line model places the resonance 3.4 % high and the edge
  resistance 13 % high against the full-wave solution, within the accuracy usually stated
  for it on a substrate of 0.026 wavelengths in the dielectric. Its directivity is that of
  two slots over an infinite ground plane; the full-wave model has a finite ground plane,
  20 mm beyond the copper, and the substrate's fields. The efficiency of the lossless model
  is the full-wave solution's energy balance.
- **Mesh.** With the minimum cell halved, the microstrip's delay changes by −0.09 %
  (122.54 to 122.43 ps), its $\lvert S_{21} \rvert$ at 10 GHz by 0.015 dB and its largest
  $\lvert S_{11} \rvert$ from 0.060 to 0.059 (0.34 to 0.98 million cells, 7735 to 16 641
  time steps); the patch's resonance by one 10 MHz step of the sweep (2.370 to 2.380 GHz),
  its edge resistance by 1.3 % (286.6 to 290.2 Ω), its directivity by 0.02 dB and its energy
  balance from +1.7 % to −0.8 % (0.36 to 1.19 million cells, 24 920 to 67 932 time steps).
  The deviations of the patch estimate above are not the full-wave mesh's.
- **Vias.** The circles of a via's pads and antipads had no mesh lines of their own, so
  the dielectric ring between pad and plane spanned one or two cells of uncertain width.
  Lines at both circles and halfway between leave the stub notches of two via boards
  (a 2.8 mm stub in one dielectric, a 2 mm stub through five planes) at 10.0 and 10.4 GHz;
  tracing the circles at four angles on cells of 0.025 mm moves them to 10.4 and 11.5 GHz
  (`docs/validation/networks.md`, "Via model"): the staircased corners of pads and antipads add capacitance,
  and via structures need the finer mesh.

## A real board

The PCIe transmit pair of `cm5_minima` (`examples/cm5_minima_si.yaml`, channel `PCIE_TX`:
CM5 connector on the bottom to M.2 connector on the top, through vias), exported to 10 GHz
with cells of 0.05 mm at the signal copper and 0.25 mm elsewhere: 6.6 million cells,
port 1 excited, 118 988 time steps to −40 dB (6.6 hours on 16 threads, 1.4 GB).

| quantity | openEMS | native |
| --- | --- | --- |
| $\lvert S_{31} \rvert$ (through) at 1, 2.5, 5, 7.5, 10 GHz | −1.02, −1.17, −1.60, −2.41, −2.57 dB | −0.40, −0.89, −1.66, −2.85, −4.35 dB |
| $\lvert S_{11} \rvert$ at the same | −17.7, −14.1, −18.0, −11.6, −16.2 dB | −28.6, −24.8, −27.1, −17.8, −26.6 dB |
| group delay of $S_{31}$, 1–5 GHz | 362 ps | 314 ps |

- **Delay.** Half of each member's 51 mm runs at 45° (22.4 and 23.6 mm) with arcs at the
  bends. On the rectangular grid such a track is a staircase, whose current path is up to
  $\sqrt{2}$ times as long: at most 9.3 and 9.8 mm more, 59 and 61 ps. The full-wave delay
  exceeds the native one by 48 ps, within that bound; the straight microstrip above agrees
  within 1.5 %. The full-wave model of a board routed at 45° is a reference for loss and
  resonances, not for delay.
- **Loss.** The full-wave loss is higher below 5 GHz and lower above, as on the microstrip:
  its dielectric loss is a conductivity fixed at 5 GHz and its copper is lossless.
- **Reflection.** The full-wave reflection is 7 to 11 dB larger: its lumped ports add the
  path from the pad to the plane, and the staircase edges change the impedance of the
  diagonal sections.

## Losses in the full-wave export

The microstrip above (20 mm of 0.35 mm on 0.2 mm, $\varepsilon_r$ 4.5, loss
tangent 0.02, cells of 0.05 mm), $\lvert S_{21} \rvert$ in dB:

| model | 1 GHz | 2 GHz | 5 GHz | 10 GHz | 15 GHz | 20 GHz |
| --- | --- | --- | --- | --- | --- | --- |
| native channel | −0.10 | −0.18 | −0.39 | −0.71 | −1.02 | −1.31 |
| openEMS, Debye dielectric and copper sheets | −0.12 | −0.19 | −0.44 | −0.82 | −1.21 | −1.64 |
| openEMS, conductivity at 10 GHz and perfect copper | −0.60 | −0.59 | −0.62 | −0.65 | −0.71 | −0.83 |

- Six Debye poles reproduce the Djordjevic-Sarkar permittivity within $2 \times 10^{-4}$
  from 0.1 to 20 GHz; with the copper's conductivity as conducting sheets the full-wave
  loss follows the native one, where the single conductivity gave six times the loss at
  1 GHz and two thirds of it at 20 GHz.
- The full-wave model loses 0.1 dB more at 10 GHz and 0.3 dB more at 20 GHz: radiation of
  the open line and the lumped ports, which the native channel does not have, and the
  conducting sheet's model of the copper's thickness.

### Ports

The same microstrip with microstrip ports at its pads (the reference plane moved back to the
pad along the line, attenuation included) instead of lumped boxes:

| ports | largest $\lvert S_{11} \rvert$ to 10 / 20 GHz | $\lvert S_{21} \rvert$ at 10 / 20 GHz | group delay 1–10 GHz |
| --- | --- | --- | --- |
| native channel (the line's own mismatch) | 0.020 / 0.020 | −0.71 / −1.31 dB | 121.06 ps |
| lumped boxes | 0.058 / 0.113 | −0.82 / −1.64 dB | 122.57 ps |
| microstrip ports | 0.042 / 0.064 | −0.79 / −1.46 dB | 121.89 ps (+0.7 %) |

The vertical path of the lumped port from the pad to the plane goes, and with it a third of
the reflection and half the delay error. The rest of the reflection is the port's feed and
the pad. On the patch above a microstrip port on the feed line moves the
directivity from 7.20 to 6.08 dBi and its energy balance to +6 %; the patch keeps its
lumped port until that is understood.

## Defects found

Three defects of the full-wave export were found by the microstrip, each giving an open
port:

- A lumped port of no extent across the current (a sheet, or a line) has a current probe
  whose loop is degenerate: the port read 370 to 660 Ω. The ports are now boxes as wide as
  the track.
- The smoothing of the mesh moved lines by about $10^{-16}$ mm; a sheet off its line is
  left out of the grid, and the strip was missing. Lines near a sheet's plane or a port's
  face now take it exactly. The loop export uses the same mesh code and has the same
  correction.
- Exported polygons kept collinear vertices; they are now removed.

## Not validated

No measured S-parameters or antenna patterns; no comparison with another full-wave
solver.
