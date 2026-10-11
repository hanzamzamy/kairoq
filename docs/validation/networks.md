# Validation: networks

## Scope

The network analysis (`docs/solvers/networks.md`) is checked against:

1. closed forms of uniform lines, single and coupled, built from cross-sections of the
   same geometry, and the input reflection of a line into a load;
2. the independence of the result from the cut spacing, and the TDR's impedance and
   timing;
3. the analytic mode conversion of a skew, and the quarter-wave resonance of a via stub;
4. reciprocity and passivity on a real channel;
5. Touchstone output read back;
6. the cross-sections of channels: dispersion and the reach of a cut;
7. bodies of revolution in the cross-section solver, and the via model against full-wave
   solutions;
8. measured S-parameters of two stub lines.

The full-wave check of the channel's delay is in `docs/validation/rf.md`. The tests are in
`tests/python/test_networks.py`, `tests/python/test_si.py` and
`tests/cpp/test_via_model.cpp`.

## Lines and networks

| case | quantity | reference | deviation |
| --- | --- | --- | --- |
| microstrip 20 mm, $w = 1$, $h = 1.5$ mm, ports at its pads, to 10 GHz | $\lvert S_{21} \rvert$, phase of $S_{21}$ | uniform line with the per-length parameters of a cross-section of the same geometry | 0.0002 dB, 0.002° at most |
| the same | $S_{11}$ | the same | below 0.01 |
| the same, one port, 50 Ω at the far pad | $S_{11}$ | input impedance of the line into the resistor | below 0.01 |
| the same, cuts every 0.5 and every 3 mm | S-parameters | each other | $10^{-9}$ |
| the same, TDR (45 ps) | impedance on the line | the cross-section: 81.83 Ω | $10^{-4}$ |
| the same | time of the far end's reflection | twice the line's delay: 236 ps | $-2\ \%$ |
| pair $w = s = h = 0.2$ mm, 20 mm, to 10 GHz | $S_{dd21}$ | uniform line of the odd mode of the coupled cross-section | $2.4 \times 10^{-4}$ at most |
| the same | $\lvert S_{cd21} \rvert$ | symmetric pair: 0 | below $10^{-3}$ |
| the pair with 1 mm more on one member | $S_{cd21}$ at 2.5 and 5 GHz | $\lvert \sin(\omega \Delta\tau / 2) \rvert$: −26.5, −20.5 dB | −27.2, −21.2 dB |
| through vias with a 2.835 mm stub, $\varepsilon_r = 4.5$ | frequency of the $S_{21}$ notch | quarter-wave stub: 12.46 GHz | 12.45 GHz |

- The uniform-line comparisons check the assembly: sections, nodes, the conversion to
  S-parameters, the mixed-mode transformation and the coupled two-conductor sections
  (matrix exponential). The cut through the board and the hand-built cut are meshed
  independently; their difference is the solver's.
- The skew's mode conversion follows the closed form for matched lines within 0.7 dB; the
  remainder is the mismatch of the lines to the 50 Ω ports.
- The via notch confirms the via topology (the stub below the inner layer, open at the
  bottom). The via model and its comparison with full-wave solutions are in section "Via
  model".
- The TDR's far-end reflection crosses halfway 2 % before twice the delay: the Gaussian
  step of 45 ps rises over a tenth of the round trip.

### A real channel

The PCIe transmit pair of `cm5_minima`, from the CM5 connector (`Module302`, bottom) to
the M.2 connector (`J701`, top), as a four-port: 149 line sections (92 coupled), 10 via
sections, 238 nodes; 3 s for 400 frequencies on 20 threads.

- Reciprocal to $6 \times 10^{-13}$; passive (largest singular value of $S$ 0.998).
- $S_{dd21}$ at 1 GHz is −0.35 dB, against the 0.34 dB loss of the pair's lines in the
  signal-integrity analysis; −1.1 dB at 5 GHz, −2.4 dB at 10 GHz.
- The differential TDR lies between 85 and 102 Ω for a pair of the 90 Ω class.

## Cross-sections of channels

- **Dispersion.** A cross-section keeps the field energy of each dielectric region at its
  solved frequency and weights it with the region's permittivity at any other frequency.
  On a masked microstrip from 1 MHz to 20 GHz this agrees with new solutions at two
  frequencies per decade (an option) within $1.9 \times 10^{-4}$ in capacitance, 0.35 % in
  conductance (where the interpolation of the anchors is itself approximate) and 0.01 % in
  impedance.
- **Reach of a cut.** A cut reaches twenty times the distance to the shield, at most 25 mm,
  measured beyond the track and its partner; its dielectric ends at the board's outline.
  The field of a microstrip falls as the inverse square of the lateral distance: a
  hand-built cut of the microstrip above, on 1.5 mm, gives 81.50, 82.07 and 82.24 Ω at
  half widths of 8, 15 and 25 mm, and a reach of ten times the distance, at most 2.5 mm,
  gives 77.4 Ω. Measured from the track alone, the reach cut off a partner beyond it: on
  `cm5_minima`, where the HDMI pairs run with a ground strip between their members, one
  cut gave 131 Ω with the partner halved against 107 Ω with it whole.

The pairs of `cm5_minima` give:

| pairs | mean differential impedance | class |
| --- | --- | --- |
| Ethernet | 97.8 Ω | 100 Ω |
| HDMI | 93.2 Ω | 100 Ω |
| MIPI D-PHY | 99.6 Ω | 100 Ω |
| USB | 87.4 Ω | 90 Ω |
| PCIe | 91.3 Ω | 90 Ω |

The classes' geometry itself (tracks of 0.13 mm 0.25 mm apart for 100 Ω, 0.14 mm 0.154 mm
apart for 90 Ω, on the board's stackup) gives 102.2 and 92.5 Ω: the classes agree with the
stackup within 2.7 %; the routed pairs depart further where the routing departs from the
class geometry (ground between the members, other spacings, fan-outs). The test checks the
class geometry within 4 % and the routed means within 10 %. The analysis of all 36 pair
nets takes 86 s (1352 distinct cross-sections).

## Bodies of revolution

The cross-section solver takes the cut's horizontal axis as the radius of a body of
revolution, weights each element by $2\pi r$ and leaves the axis free.

| case | quantity | reference | deviation |
| --- | --- | --- | --- |
| sphere of 1 mm in a shell of 2 mm | capacitance | $4\pi\varepsilon_0 ab/(b-a)$ = 0.2225 pF | $3 \times 10^{-5}$ |
| coaxial line, 0.15 and 0.6 mm, 4 and 8 mm long | capacitance of the 4 mm difference | $2\pi\varepsilon_0 \ell / \ln(b/a)$ = 0.1605 pF | $8 \times 10^{-5}$ |

## Via model

Two boards, each with a track on the top layer to a via, on an inner layer to a second via
10 mm on, and back to the top; the vias are drilled through, so each leaves a stub:

- **8 layers**: dielectrics of 0.2, 0.2, 0.4, 0.4, 0.4, 0.4 and 0.2 mm, ground planes on
  all but L1 and L3, a 2 mm stub from L3 through four pairs of planes, non-functional pads
  of 0.55 mm in antipads of 0.9 mm, four return vias at 1 mm;
- **thick**: four layers, ground on L2 and L4, a 2.8 mm stub in one dielectric from L3 to
  L4, no return vias, 40 × 10 mm.

Frequency of the $S_{21}$ notch:

| case | openEMS, 0.05 mm cells | openEMS, 0.025 mm, circles traced | via model | coaxial model |
| --- | --- | --- | --- | --- |
| 8 layers | 10.4 GHz | 11.5 GHz | 12.3 GHz (+7 %) | 17.9 GHz (+56 %) |
| thick | 10.0 GHz | 10.4 GHz | 8.4 GHz (−19 %) | 12.5 GHz (+20 %) |
| 8 layers, vias blind to L3 (no stub) | no notch; −1.26 dB at 10 GHz | | −1.42 dB | −1.45 dB |

- **The reference.** On cells of 0.05 mm the circles of a pad of 0.55 mm and an antipad of
  0.9 mm are staircased into squares, whose corners narrow the ring between them; tracing
  the circles with mesh lines at four angles on cells of 0.025 mm raises the 8-layer
  notch by 11 %. A larger board (45 × 15 mm) keeps the notch where the board's modes move:
  it is the via's. The blind vias have none: it is the stub's.
- **The model.** The barrel is a ladder: node capacitances from an axisymmetric field
  solution, the pairs' fundamental mode between the barrels, the higher modes'
  inductance. On thin pairs with return vias it places the notch within 7 %; in a single
  thick dielectric on a small board, whose open edges radiate, 19 % low. The quarter-wave
  rule of the coaxial model misses the loading of the pads and planes (8 layers) and holds
  only by chance on the thick board.
- **Field solutions.** The barrel ladder against a coaxial line (a barrel through a thick
  plane, the capacitance of one millimetre more): within 1 %. Complex Hankel functions
  against scipy: $10^{-10}$.
- **Cost.** On `cm5_minima`'s PCIe pair the via model adds 7 s (two barrels' field
  solutions) to the channel's 13 s.

Defects found on the way, corrected: the axisymmetric solver's per-region field energies
lacked the $2\pi r$ weight (the capacitances with a dielectric were wrong; a test now
compares them with the vacuum's in a uniform fill); and the full-wave export gave the
circles of via pads and antipads no mesh lines.

## Measured S-parameters

The open-source tool gerber2ems ships slices of the Antmicro Signal Integrity Test Board
with one-port VNA measurements (`testdata/gerber2ems`, Apache License 2.0): a 50 Ω
microstrip 0.185 mm wide over 0.12 mm of FR-4, with an open stub of 8.4 mm (`stub_short`)
and a longer one (`stub_long`), $S_{11}$ at port SP1 from 0.2 to 6 GHz with SP2 loaded. The
networks model of the channel, with the ports at the placement positions of the simulation
ports:

| case | rms $\lvert S_{11} - S_{11,\text{meas}} \rvert$ | rms $\lvert S_{11,\text{meas}} \rvert$ | largest $\lvert Z \rvert$ | smallest $\lvert Z \rvert$ |
| --- | --- | --- | --- | --- |
| `stub_short` | 0.087 (largest 0.14) | 0.58 | 4.76 GHz (measured 4.84, −1.6 %) | 2.15 GHz (measured 2.28, −6 %) |
| `stub_long` | 0.110 (largest 0.18) | 0.50 | | |

The measurement carries a ripple of about 0.25 GHz from the set-up, averaged out over one
period before the comparison. A 0 Ω part at the stub's end, which would short it, is taken
as not fitted: the measured resonance is that of an open stub, a quarter wave at the
junction, and with the part the model's resonance falls at 4.1 GHz.

The comparison found three defects of the import and the networks model, corrected:

- the drill and placement files of these boards are relative to the board's lower-left
  corner: the via of the stub's ground pad joined nothing, and the simulation ports lay
  off the board. Both are now moved when none of their points lies on the board and all do
  relative to that corner;
- the stub ends on the inside of the line's single track: the model had no junction there,
  and the stub was missing. A track ending on the inside of another is now a junction;
- the arcs and lines of the meanders end within nanometres of each other, but on either
  side of a micrometre: their ends were distinct nodes, and most of a meander was cut off.
  Ends within 1 µm are now one node.

Two further cases were not compared. The measurements of the two meanders exceed
$\lvert S_{11} \rvert = 1$ at 39 % of their points, up to 1.097, which a passive one-port
cannot, and their phase gives a third of the drawn line's delay. The filter is drawn as
copper areas, which the networks model does not cut into lines.

## Performance

Networks: 0.1 s for a synthetic line of 18 sections, 3 s for the PCIe pair of
`cm5_minima` at 400 frequencies; the cuts and their cross-sections (the
signal-integrity analysis of the channel's nets) dominate. All Python tests run in
5 minutes.

## Not validated

Measured data enter for the two stub lines only. Channels of real boards are compared with
full-wave solutions of their structure, not with measurements.
