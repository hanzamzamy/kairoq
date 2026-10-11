# Time-domain response

The waveforms at the ports of a channel (`docs/solvers/networks.md`) with drivers and
receivers: linear sources and loads, or IBIS buffer models. Reflections, ringing,
overshoot, crosstalk with real drivers, and received bit patterns. Validation:
`docs/validation/timedomain.md`.

## Method

**The channel** is its scattering matrix $S(f)$ on the linear sweep. Its impulse
responses $s_{ij}(t)$ follow from the inverse Fourier transform of the sweep, tapered to
zero over its last quarter (a raised cosine from $0.75 f_\text{max}$ to $f_\text{max}$), the
lowest frequency's value extended to DC, at a time step of $1/(2 n f_\text{max})$ with an
oversampling $n$ (4). Band-limited, an instantaneous part of a response (the reflection at
a port's own junction) spreads to either side of $t = 0$; the part before it, within
$3/f_\text{max}$, is added at $t = 0$, which keeps the response causal and its area (its DC
value) exact. The responses repeat with the sweep's period $1/\Delta f$; a simulation lasts
at most half of it (10 ns for 400 points to 20 GHz).

**At each time step** the waves entering the ports $a$ and leaving them $b$ satisfy the
convolution

$$b(t_n) = s(0)\, a(t_n) + \sum_{m \ge 1} s(t_m)\, a(t_{n-m}),$$

with the port voltages and currents (into the channel)
$V = \sqrt{Z_0}\,(a + b)$, $I = (a - b)/\sqrt{Z_0}$ for the channel's reference $Z_0$. The
waves entering are solved together with the terminations by Newton's method:

- **open** ports: no current;
- **linear drivers**: a source between the low and high levels, linear transitions over
  1.25 times the 10–90 % time, behind a resistance;
- **linear receivers**: a resistance to a termination voltage and a capacitance;
- **IBIS buffers**: the pin, the package (series $R_\text{pkg}$, $L_\text{pkg}$; $C_\text{pkg}$ at the
  pin), the die with $C_\text{comp}$, and the buffer's current
  $k_u(t)\, I_\text{pu}(V_\text{pu} - V) + k_d(t)\, I_\text{pd}(V - V_\text{pd}) + I_\text{gc}(V) + I_\text{pc}(V_\text{pc} - V)$
  from its I-V tables; an input has only its clamps and $C_\text{comp}$.

Capacitors and inductors are integrated by backward Euler. The simulation starts from rest
(all voltages zero); the drivers start low.

**Switching coefficients.** An IBIS output switches between its pulldown and its pullup
by $k_u(t)$ and $k_d(t)$. With two V-t waveforms of an edge (into fixtures to different
voltages), at each time the two equations

$$k_u I_\text{pu} + k_d I_\text{pd} + I_\text{clamps} + C_\text{comp}\frac{dV}{dt} = -\frac{V - V_\text{fixture}}{R_\text{fixture}} - C_\text{fixture}\frac{dV}{dt},$$

one per waveform, give both; with one, $k_d = 1 - k_u$; with none, a linear ramp at the
`[Ramp]` slope. A bit pattern applies the rising or falling coefficients from each
transition.

**Bit patterns.** `edge` (low, then high from the delay), `pulse` (one bit high), `prbs7`
to `prbs31` (maximal-length sequences of $x^n + x^m + 1$ from all ones), or a string of
`0` and `1`, at a bit rate after a delay.

## Project file

```yaml
si:
  channels: [{name: CLK, ports: [U1.4, U7.2]}]
  network: {f_max_hz: 20e9, points: 800}       # 800 points: 20 ns, of which 10 ns simulated
  simulations:
    - name: CLK_EDGE
      channel: CLK
      drivers:
        - {port: 1, pattern: edge, ibis: {file: u1.ibs, model: LVCMOS33_8MA, corner: typ}}
      receivers:
        - {port: 2, ibis: {file: u7.ibs, model: IN_33}}
    - name: DATA
      channel: CLK
      drivers: [{port: 1, pattern: prbs7, bit_rate: 1e9, bits: 8, v_high: 3.3, r_ohm: 30, rise_ps: 400}]
      receivers: [{port: 2, r_ohm: 0, c_pf: 3}]
```

| key | default | meaning |
| --- | --- | --- |
| `simulations[]` | | `name`, `channel` (of `si.channels`), `drivers`, `receivers`, `duration_ns` (the pattern and 2 ns), `oversample` (4) |
| `drivers[]` | | `port` (from 1), `pattern` (`edge`), `bit_rate` ($10^9$), `bits` (64, of a PRBS), `delay_ps` (100); linear: `v_low` (0), `v_high` (1), `r_ohm` (50), `rise_ps` (50), `fall_ps` (= `rise_ps`); or `ibis`: `file`, `model`, `corner` (`typ`, `min`, `max`), `package` (true) |
| `receivers[]` | | `port`; linear: `r_ohm` (0: none), `v_term` (0), `c_pf` (0); or `ibis` as for drivers. Other ports are open |

## Use

```
kairoq waveforms project.yaml [<input>...] [-o dir]
```

prints per simulation and port the extremes and the final value, and for single edges the
overshoot and undershoot; with `-o`, `NAME.csv` holds the waveforms (time and one column
per port).

```python
r = kairoq.simulate(board, {"channels": [...], "network": {...}, "simulations": [...]})[0]
r["time_s"], r["voltage"], r["metrics"]
kairoq.ibis_switching("u1.ibs", "LVCMOS33_8MA", rising=True)
```

## Limits

- The channel is linear and its sweep band-limited: edges faster than $0.9/f_\text{max}$ are
  reported, and details of the response finer than about $1/f_\text{max}$ are lost.
- The responses repeat with the sweep's period: the simulation lasts at most half of it.
  Long patterns need finer sweeps; eyes over long patterns follow from the response to a
  single bit (`docs/solvers/eye.md`).
- IBIS: version 3 to 7 keywords of a model's behaviour (`[Pullup]`, `[Pulldown]`,
  `[GND Clamp]`, `[POWER Clamp]`, the references, `[Ramp]`, `[Rising Waveform]`,
  `[Falling Waveform]`, `C_comp`, `[Package]`); not read: IBIS-AMI, `[Driver Schedule]`,
  `[Submodel]`, `[ISSO PU]`/`[ISSO PD]` and supply effects (the rails are ideal); not used:
  the fixtures' inductance and DUT parasitics, the per-pin package values. Unread keywords
  are reported.
- Backward Euler damps resonances of the package slightly at coarse steps.
