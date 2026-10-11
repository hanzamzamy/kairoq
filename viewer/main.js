// kairoq viewer: layered 3D board geometry served by the core (/api/board, /api/mesh).
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

const $ = (id) => document.getElementById(id);

const LAYER_COLOURS = [0xd4a017, 0x3fa7d6, 0xe05a47, 0x59c36a, 0xb57edc, 0xf2a65a, 0x4dd0c4, 0xe86fa8];
const DIELECTRIC_COLOUR = 0x8a7a3d;
const MASK_COLOUR = 0x1f6f3a;
const DRILL_COLOUR = 0xc0c0c0;
const DIM = 0.18;  // brightness of copper outside the selected net

// Viridis, 9 stops.
const VIRIDIS = [[68, 1, 84], [71, 44, 122], [59, 81, 139], [44, 113, 142], [33, 144, 141],
  [39, 173, 129], [92, 200, 99], [170, 220, 50], [253, 231, 37]];

const state = {
  board: null,
  fields: { dc: new Map(), thermal: new Map(), loop: new Map() },  // layer -> { mesh, data }
  loop: 0,  // shown loop of the EMI estimates
  snapshot: -1,  // shown snapshot of a transient
  field: '',
  groups: new Map(),  // mesh key -> { mesh, ids, base colour }
  colourMode: 'layer',
  selectedNet: -1,
  zScale: 1,
  siMeshes: [],  // impedance-coloured tracks
};

// ------------------------------------------------------------------ scene

const canvas = $('view');
const renderer = new THREE.WebGLRenderer({ canvas, antialias: true });
renderer.setPixelRatio(window.devicePixelRatio);
renderer.localClippingEnabled = true;
const scene = new THREE.Scene();
scene.background = new THREE.Color(0x15171a);
const camera = new THREE.PerspectiveCamera(35, 1, 0.01, 10000);
camera.up.set(0, 0, 1);
const controls = new OrbitControls(camera, canvas);
controls.enableDamping = false;
scene.add(new THREE.HemisphereLight(0xffffff, 0x404040, 2.0));
const sun = new THREE.DirectionalLight(0xffffff, 1.6);
sun.position.set(0.4, -0.6, 1.0);
scene.add(sun);
const board = new THREE.Group();  // scaled in z by the vertical scale
scene.add(board);

const sectionPlane = new THREE.Plane(new THREE.Vector3(-1, 0, 0), 0);

function resize() {
  const w = canvas.clientWidth, h = canvas.clientHeight;
  renderer.setSize(w, h, false);
  camera.aspect = w / h;
  camera.updateProjectionMatrix();
  render();
}
window.addEventListener('resize', resize);

let pending = false;
function render() {
  if (pending) return;
  pending = true;
  requestAnimationFrame(() => {
    pending = false;
    renderer.render(scene, camera);
  });
}
controls.addEventListener('change', render);

// ------------------------------------------------------------------ data

async function fetchMesh(key) {
  const res = await fetch(`api/mesh?key=${encodeURIComponent(key)}`);
  if (!res.ok) throw new Error(`${key}: ${await res.text()}`);
  const buf = await res.arrayBuffer();
  const head = new Uint32Array(buf, 0, 2);
  const nv = head[0], ni = head[1];
  let off = 8;
  const positions = new Float32Array(buf, off, nv * 3); off += nv * 12;
  const indices = new Uint32Array(buf, off, ni); off += ni * 4;
  const ids = new Uint32Array(buf, off, nv);
  return { positions, indices, ids };
}

function netColour(id) {
  if (id === 0) return new THREE.Color(0x808080);
  const c = new THREE.Color();
  c.setHSL((id * 0.61803398875) % 1, 0.65, 0.55);
  return c;
}

function describeKey(key) {
  const b = state.board;
  const [kind, arg] = key.split('/');
  if (kind === 'copper') {
    const l = b.copper_layers.find((c) => c.index === Number(arg));
    return { kind, label: `L${arg} ${l ? l.name : ''}`, index: Number(arg) };
  }
  if (kind === 'dielectric') return { kind, label: b.stackup.layers[Number(arg)].name || 'dielectric', index: Number(arg) };
  if (kind === 'mask') return { kind, label: `mask ${arg}` };
  return { kind, label: 'plated holes' };
}

function material(kind, index) {
  const common = { side: THREE.DoubleSide, clippingPlanes: [], flatShading: true };
  if (kind === 'copper') return new THREE.MeshStandardMaterial({ ...common, vertexColors: true, metalness: 0.3, roughness: 0.5 });
  if (kind === 'drills') return new THREE.MeshStandardMaterial({ ...common, vertexColors: true, metalness: 0.5, roughness: 0.4 });
  if (kind === 'dielectric')
    return new THREE.MeshStandardMaterial({ ...common, color: DIELECTRIC_COLOUR, transparent: true, opacity: 0.35, depthWrite: false });
  return new THREE.MeshStandardMaterial({ ...common, color: MASK_COLOUR, transparent: true, opacity: 0.45, depthWrite: false });
}

function baseColour(key, info) {
  if (info.kind === 'copper') return new THREE.Color(LAYER_COLOURS[(info.index - 1) % LAYER_COLOURS.length]);
  return new THREE.Color(DRILL_COLOUR);
}

function recolour() {
  for (const [key, g] of state.groups) {
    if (!g.colours) continue;
    const arr = g.colours.array;
    const cache = new Map();
    for (let v = 0; v < g.ids.length; v++) {
      const id = g.ids[v];
      let c = cache.get(id);
      if (!c) {
        c = state.colourMode === 'net' ? netColour(id) : g.base.clone();
        if (state.selectedNet >= 0 && id !== state.selectedNet) c.multiplyScalar(DIM);
        if (state.selectedNet >= 0 && id === state.selectedNet) c.lerp(new THREE.Color(0xffffff), 0.35);
        cache.set(id, c);
      }
      arr[3 * v] = c.r; arr[3 * v + 1] = c.g; arr[3 * v + 2] = c.b;
    }
    g.colours.needsUpdate = true;
  }
  render();
}

async function addMesh(key) {
  const info = describeKey(key);
  const data = await fetchMesh(key);
  const geom = new THREE.BufferGeometry();
  geom.setAttribute('position', new THREE.BufferAttribute(data.positions, 3));
  geom.setIndex(new THREE.BufferAttribute(data.indices, 1));
  let colours = null;
  if (info.kind === 'copper' || info.kind === 'drills') {
    colours = new THREE.BufferAttribute(new Float32Array(data.ids.length * 3), 3);
    geom.setAttribute('color', colours);
  }
  geom.computeBoundingBox();
  geom.computeBoundingSphere();
  const mesh = new THREE.Mesh(geom, material(info.kind, info.index));
  mesh.userData.key = key;
  mesh.renderOrder = info.kind === 'dielectric' || info.kind === 'mask' ? 1 : 0;
  board.add(mesh);
  state.groups.set(key, { mesh, ids: data.ids, colours, base: baseColour(key, info), info });
}

// ------------------------------------------------------------------ DC fields

// Fields of a solution, from /api/field (DC) or /api/thermal.
async function fetchField(kind, layer) {
  const url = kind === 'loop' ? `api/loop?loop=${state.loop}&layer=${layer}`
    : `api/${kind === 'thermal' ? 'thermal' : 'field'}?layer=${layer}`;
  const res = await fetch(url);
  if (!res.ok) throw new Error(await res.text());
  const buf = await res.arrayBuffer();
  const [nv, nt] = new Uint32Array(buf, 0, 2);
  let off = 8;
  const take = (Type, n) => { const a = new Type(buf, off, n); off += n * 4; return a; };
  if (kind === 'thermal')
    return { positions: take(Float32Array, nv * 3), indices: take(Uint32Array, nt * 3),
      temperature: take(Float32Array, nv), copper: take(Uint32Array, nt) };
  if (kind === 'loop')
    return { positions: take(Float32Array, nv * 3), indices: take(Uint32Array, nt * 3), current: take(Float32Array, nt) };
  return {
    positions: take(Float32Array, nv * 3), indices: take(Uint32Array, nt * 3),
    voltage: take(Float32Array, nv), current_density: take(Float32Array, nt),
    power_density: take(Float32Array, nt), net: take(Uint32Array, nt),
  };
}

function colourMap(t) {
  t = Math.min(1, Math.max(0, t)) * (VIRIDIS.length - 1);
  const i = Math.min(Math.floor(t), VIRIDIS.length - 2), f = t - i;
  const a = VIRIDIS[i], b = VIRIDIS[i + 1];
  return [0, 1, 2].map((k) => (a[k] + f * (b[k] - a[k])) / 255);
}

// Value range of each net in the loaded DC fields: voltages are coloured per net, since
// the drops along one net are small against the differences between nets.
function netRanges() {
  const ranges = new Map();
  for (const f of state.fields.dc.values()) {
    const d = f.data, nt = d.indices.length / 3;
    for (let t = 0; t < nt; t++) {
      let r = ranges.get(d.net[t]);
      if (!r) { r = [Infinity, -Infinity]; ranges.set(d.net[t], r); }
      for (let c = 0; c < 3; c++) {
        const v = d.voltage[d.indices[3 * t + c]];
        if (v < r[0]) r[0] = v;
        if (v > r[1]) r[1] = v;
      }
    }
  }
  return ranges;
}

function fieldScale() {
  if (state.field === 'voltage') {
    const ranges = netRanges();
    const map = (v, net) => {
      const r = ranges.get(net);
      return r && r[1] > r[0] ? (v - r[0]) / (r[1] - r[0]) : 0.5;
    };
    const sel = ranges.get(state.selectedNet);
    return { lo: sel ? sel[0] : NaN, hi: sel ? sel[1] : NaN, map, perNet: true };
  }
  if (state.field === 'temperature') {
    const [lo, hi] = state.board.thermal.range;
    return { lo, hi, map: (v) => (hi > lo ? (v - lo) / (hi - lo) : 0.5), perNet: false };
  }
  const log = $('field-log').checked;
  const hi = state.field === 'loop_current' ? state.board.emi.loops[state.loop].max_current_a_per_mm
    : state.board.dc.range[state.field][1];
  const floor = hi * 1e-4;
  const map = log
    ? (v) => (Math.log10(Math.max(v, floor)) - Math.log10(floor)) / (Math.log10(hi) - Math.log10(floor))
    : (v) => v / hi;
  return { lo: log ? floor : 0, hi, map, perNet: false };
}

function fmt(v) {
  const a = Math.abs(v);
  return a !== 0 && (a < 1e-3 || a >= 1e4) ? v.toExponential(3) : v.toPrecision(4);
}

const fieldKind = (name) => (name === 'temperature' ? 'thermal' : name === 'loop_current' ? 'loop' : 'dc');
const fieldLayers = (kind) => (kind === 'loop' ? state.board.emi.loops[state.loop].layers : state.board[kind].layers);

function value(kind, d, t, c) {
  if (kind === 'thermal') return d.temperature[d.indices[3 * t + c]];
  if (kind === 'loop') return d.current[t];
  if (state.field === 'voltage') return d.voltage[d.indices[3 * t + c]];
  return d[state.field][t];
}

async function showField() {
  state.field = $('field').value;
  for (const set of Object.values(state.fields)) for (const f of set.values()) f.mesh.visible = false;
  if (!state.field) { $('legend-min').textContent = ''; $('legend-max').textContent = ''; render(); return; }
  const kind = fieldKind(state.field);
  const set = state.fields[kind];
  for (const layer of fieldLayers(kind)) {
    if (set.has(layer)) continue;
    const data = await fetchField(kind, layer);
    const nt = data.indices.length / 3;
    const pos = new Float32Array(nt * 9);
    for (let t = 0; t < nt; t++)
      for (let c = 0; c < 3; c++) {
        const v = data.indices[3 * t + c];
        pos[9 * t + 3 * c] = data.positions[3 * v];
        pos[9 * t + 3 * c + 1] = data.positions[3 * v + 1];
        pos[9 * t + 3 * c + 2] = data.positions[3 * v + 2] + 0.001;
      }
    const geom = new THREE.BufferGeometry();
    geom.setAttribute('position', new THREE.BufferAttribute(pos, 3));
    geom.setAttribute('color', new THREE.BufferAttribute(new Float32Array(nt * 9), 3));
    geom.computeBoundingSphere();
    const material = new THREE.MeshBasicMaterial({
      vertexColors: true, side: THREE.DoubleSide, clippingPlanes: [],
      polygonOffset: true, polygonOffsetFactor: -2, polygonOffsetUnits: -4,
    });
    const mesh = new THREE.Mesh(geom, material);
    mesh.userData.fieldLayer = layer;
    mesh.userData.fieldKind = kind;
    board.add(mesh);
    set.set(layer, { mesh, data });
  }
  if (kind === 'thermal' && state.board.thermal.transient && state.snapshot >= 0) await loadSnapshot(state.snapshot);
  const sc = fieldScale();
  for (const [layer, f] of set) {
    const col = f.mesh.geometry.getAttribute('color').array, d = f.data, nt = d.indices.length / 3;
    for (let t = 0; t < nt; t++)
      for (let c = 0; c < 3; c++) col.set(colourMap(sc.map(value(kind, d, t, c), d.net ? d.net[t] : 0)), 9 * t + 3 * c);
    f.mesh.geometry.getAttribute('color').needsUpdate = true;
    const copper = state.groups.get(`copper/${layer}`);
    f.mesh.visible = !copper || copper.mesh.visible;
    f.mesh.material.clippingPlanes = $('section-on').checked ? [sectionPlane] : [];
  }
  const ctx = $('legend').getContext('2d');
  for (let x = 0; x < 240; x++) {
    const [r, g, b] = colourMap(x / 239);
    ctx.fillStyle = `rgb(${255 * r},${255 * g},${255 * b})`;
    ctx.fillRect(x, 0, 1, 14);
  }
  if (sc.perNet && Number.isNaN(sc.lo)) {
    $('legend-min').textContent = 'minimum of each net';
    $('legend-max').textContent = 'maximum';
  } else {
    $('legend-min').textContent = fmt(sc.lo);
    $('legend-max').textContent = fmt(sc.hi);
  }
  render();
}

function probeField(event) {
  if (!state.field) return false;
  const rect = canvas.getBoundingClientRect();
  pointer.set(((event.clientX - rect.left) / rect.width) * 2 - 1, -((event.clientY - rect.top) / rect.height) * 2 + 1);
  raycaster.setFromCamera(pointer, camera);
  const kind = fieldKind(state.field);
  const targets = [...state.fields[kind].values()].filter((f) => f.mesh.visible).map((f) => f.mesh);
  const hit = raycaster.intersectObjects(targets, false)[0];
  if (!hit) return false;
  const layer = hit.object.userData.fieldLayer;
  const d = state.fields[kind].get(layer).data, t = hit.faceIndex;
  const bc = hit.barycoord || new THREE.Vector3(1 / 3, 1 / 3, 1 / 3);
  const interp = (arr) => bc.x * arr[d.indices[3 * t]] + bc.y * arr[d.indices[3 * t + 1]] + bc.z * arr[d.indices[3 * t + 2]];
  const o = state.board.display.origin_mm;
  const at = `L${layer} at (${(hit.point.x + o[0]).toFixed(3)}, ${(hit.point.y + o[1]).toFixed(3)}) mm`;
  if (kind === 'loop') {
    $('probe').textContent = `${at}: ${fmt(d.current[t])} A/mm per ampere of loop current`;
  } else if (kind === 'thermal') {
    $('probe').textContent = `${at}: ${fmt(interp(d.temperature))} °C, ${d.copper[t] ? 'copper' : 'dielectric'}`;
  } else {
    $('probe').textContent = `${at}, ${state.board.nets[d.net[t]]}: ${fmt(interp(d.voltage))} V, ` +
      `${fmt(d.current_density[t])} A/mm², ${fmt(d.power_density[t])} W/mm²`;
  }
  return true;
}

// ------------------------------------------------------------------ transients

// Temperatures of snapshot i on the loaded thermal layers (the last one is the final state,
// as loaded with the mesh).
async function loadSnapshot(i) {
  for (const [layer, f] of state.fields.thermal) {
    f.values = f.values || new Map();
    if (!f.values.has(i)) {
      const res = await fetch(`api/thermal_values?layer=${layer}&snapshot=${i}`);
      if (!res.ok) throw new Error(await res.text());
      f.values.set(i, new Float32Array(await res.arrayBuffer()));
    }
    f.data.temperature = f.values.get(i);
  }
}

async function showSnapshot(i) {
  state.snapshot = i;
  $('time-value').textContent = `${fmt(state.board.thermal.transient.snapshot_time_s[i])} s`;
  drawHistory();
  if (state.field === 'temperature') await showField();
}

// Highest temperature and the sources over time, with the shown time marked.
function drawHistory() {
  const tr = state.board.thermal.transient, cv = $('history'), ctx = cv.getContext('2d');
  const W = cv.width, H = cv.height, L = 34, B = 16;
  ctx.clearRect(0, 0, W, H);
  const t = tr.time_s, tmax = t[t.length - 1];
  const all = [tr.max_c, ...tr.series.map((s) => s.values)];
  let lo = Infinity, hi = -Infinity;
  for (const v of all) for (const x of v) if (Number.isFinite(x)) { lo = Math.min(lo, x); hi = Math.max(hi, x); }
  if (!(hi > lo)) hi = lo + 1;
  const X = (x) => L + (W - L - 4) * x / tmax, Y = (y) => H - B - (H - B - 4) * (y - lo) / (hi - lo);
  ctx.strokeStyle = '#888'; ctx.fillStyle = '#aaa'; ctx.font = '10px sans-serif';
  ctx.strokeRect(L, 4, W - L - 4, H - B - 4);
  ctx.fillText(fmt(hi), 2, 12); ctx.fillText(fmt(lo), 2, H - B);
  ctx.fillText('0', L, H - 4); ctx.fillText(`${fmt(tmax)} s`, W - 46, H - 4);
  const colours = ['#ffffff', '#ffb000', '#4fc3f7', '#e57373', '#81c784', '#ba68c8', '#fff176'];
  all.forEach((v, k) => {
    ctx.strokeStyle = colours[k % colours.length];
    ctx.beginPath();
    v.forEach((y, i) => (i ? ctx.lineTo(X(t[i]), Y(y)) : ctx.moveTo(X(t[i]), Y(y))));
    ctx.stroke();
  });
  if (state.snapshot >= 0) {
    ctx.strokeStyle = '#f44';
    const x = X(tr.snapshot_time_s[state.snapshot]);
    ctx.beginPath(); ctx.moveTo(x, 4); ctx.lineTo(x, H - B); ctx.stroke();
  }
  cv.title = ['highest', ...tr.series.map((s) => s.name)].map((n, k) => `${['white', 'amber', 'blue', 'red', 'green', 'violet', 'yellow'][k % 7]}: ${n}`).join('\n');
}

// Return-path findings: their outlines in red on their layers.
function drawFindings() {
  const emi = state.board.emi;
  if (!emi || !emi.findings.length) return;
  const crossing = new THREE.LineBasicMaterial({ color: 0xff3030 });
  const partial = new THREE.LineBasicMaterial({ color: 0xffa040 });
  for (const f of emi.findings) {
    if (f.outline.length < 3) continue;
    const pts = f.outline.map(([x, y]) => new THREE.Vector3(x, y, f.z_mm + 0.01));
    const line = new THREE.LineLoop(new THREE.BufferGeometry().setFromPoints(pts), f.kind === 'partial' ? partial : crossing);
    line.userData.fieldLayer = f.layer;
    board.add(line);
  }
}

function buildSolutionPanel() {
  const dc = state.board.dc, th = state.board.thermal, emi = state.board.emi;
  if (!dc && !th && !emi) return;
  $('solution').hidden = false;
  $('solution').querySelector('h2').textContent = th || emi ? 'Solution' : 'DC solution';
  for (const opt of $('field').options) {
    if (opt.value === 'temperature') opt.hidden = !th;
    else if (opt.value === 'loop_current') opt.hidden = !emi || !emi.loops.length;
    else if (opt.value) opt.hidden = !dc;
  }
  if (emi && emi.loops.length > 1) {
    $('loop-row').hidden = false;
    emi.loops.forEach((l, i) => {
      const o = document.createElement('option');
      o.value = String(i);
      o.textContent = l.name;
      $('loop').append(o);
    });
    $('loop').addEventListener('change', () => {
      state.loop = Number($('loop').value);
      for (const f of state.fields.loop.values()) board.remove(f.mesh);
      state.fields.loop.clear();
      showField().catch((e) => console.error(e));
    });
  }
  const table = $('ports');
  const row = (cells) => { const tr = document.createElement('tr'); tr.innerHTML = cells; table.append(tr); };
  if (dc) for (const p of dc.ports)
    row(`<td>${p.name}</td><td class="num">${fmt(p.voltage_v)} V</td><td class="num">${fmt(p.current_a)} A</td>`);
  if (th) for (const s of th.sources)
    row(`<td>${s.name}</td><td class="num">${fmt(s.mean_c)} °C</td><td class="num">${fmt(s.power_w)} W</td>`);
  const lines = [];
  if (th) lines.push(`maximum ${fmt(th.max_c)} °C on L${th.max_layer}; heat ${fmt(th.heat_in_w)} W; ` +
    `estimated error ${(100 * th.estimated_error).toFixed(1)} %`);
  if (dc) lines.push(`dissipation ${fmt(dc.power_w)} W; estimated error ${(100 * dc.estimated_error).toFixed(1)} % ` +
    `(energy norm); ${dc.triangles} triangles`);
  $('dc-meta').textContent = lines.join('; ');
  if (emi) {
    for (const l of emi.loops)
      row(`<td>${l.name}</td><td class="num">${fmt(l.inductance_nh)} nH</td><td class="num">${l.max_field_dbuv_m > -900 ? fmt(l.max_field_dbuv_m) + ' dBµV/m' : ''}</td>`);
    for (const c of emi.coupling)
      row(`<td>${c.net}</td><td class="num">${fmt(c.capacitance_pf)} pF</td><td class="num">to ${c.plate}</td>`);
    const crossings = emi.findings.filter((f) => f.kind !== 'partial').length;
    if (emi.findings.length)
      row(`<td>return path</td><td class="num">${crossings} crossings</td><td class="num">${emi.findings.length - crossings} partial</td>`);
    drawFindings();
  }
  if (th && th.heatsinks) for (const h of th.heatsinks)
    row(`<td>${h.name}</td><td class="num">${fmt(h.temperature_c)} °C</td><td class="num">${fmt(h.heat_w)} W</td>`);
  if (th && th.transient) {
    const n = th.transient.snapshot_time_s.length;
    $('time-row').hidden = false;
    $('history').hidden = false;
    $('time').max = String(n - 1);
    $('time').value = String(n - 1);
    state.snapshot = n - 1;
    $('time-value').textContent = `${fmt(th.transient.snapshot_time_s[n - 1])} s`;
    $('time').addEventListener('input', () => showSnapshot(Number($('time').value)).catch((e) => console.error(e)));
    drawHistory();
  }
  $('field').addEventListener('change', () => showField().catch((e) => console.error(e)));
  $('field-log').addEventListener('change', () => showField());
}

// ------------------------------------------------------------------ signal integrity

// Range of values for colours and plots: the 10th to 90th percentile (short sections at
// pads and bends would otherwise set it), at least 5 % either side of the median.
function robustRange(v, p = 0.1) {
  v = v.filter((x) => x > 0).sort((a, b) => a - b);
  if (!v.length) return [0, 1];
  const q = (f) => v[Math.min(v.length - 1, Math.max(0, Math.round(f * (v.length - 1))))];
  const m = q(0.5);
  return [Math.min(q(p), 0.95 * m), Math.max(q(1 - p), 1.05 * m)];
}

function siRange(q) {
  const v = [];
  for (const t of state.board.si.tracks) for (const s of t.samples) v.push(q === 'zd' ? s[2] : s[1]);
  return robustRange(v);
}

// Each analysed track as a ribbon of its width on top of its copper, coloured by the
// sample nearest along it.
function drawSiTracks() {
  for (const m of state.siMeshes) board.remove(m);
  state.siMeshes = [];
  const q = $('si-quantity').value;
  const ctx = $('si-legend').getContext('2d');
  if (!q) { ctx.clearRect(0, 0, 240, 14); $('si-min').textContent = ''; $('si-max').textContent = ''; render(); return; }
  const [lo, hi] = siRange(q);
  const o = state.board.display.origin_mm;
  const byLayer = new Map();
  for (const t of state.board.si.tracks) {
    if (!byLayer.has(t.layer)) byLayer.set(t.layer, { pos: [], col: [] });
    const acc = byLayer.get(t.layer);
    const pts = t.points_mm;
    // Cumulative length at each point, then pieces cut where the nearest sample changes.
    const cum = [0];
    for (let i = 1; i < pts.length; i++) cum.push(cum[i - 1] + Math.hypot(pts[i][0] - pts[i - 1][0], pts[i][1] - pts[i - 1][1]));
    const samples = t.samples;
    const cuts = [];
    for (let k = 1; k < samples.length; k++) cuts.push(0.5 * (samples[k - 1][0] + samples[k][0]));
    const at = (u) => {
      let i = 1;
      while (i < cum.length - 1 && cum[i] < u) i++;
      const f = cum[i] > cum[i - 1] ? (u - cum[i - 1]) / (cum[i] - cum[i - 1]) : 0;
      return [pts[i - 1][0] + f * (pts[i][0] - pts[i - 1][0]), pts[i - 1][1] + f * (pts[i][1] - pts[i - 1][1])];
    };
    const marks = [...new Set([0, ...cum, ...cuts.filter((c) => c > 0 && c < cum[cum.length - 1]), cum[cum.length - 1]])].sort((a, b) => a - b);
    const w = t.width_mm / 2, z = t.z_mm + 0.003;
    for (let k = 1; k < marks.length; k++) {
      const u0 = marks[k - 1], u1 = marks[k];
      if (!(u1 > u0)) continue;
      const mid = 0.5 * (u0 + u1);
      let best = samples[0];
      for (const s of samples) if (Math.abs(s[0] - mid) < Math.abs(best[0] - mid)) best = s;
      const v = q === 'zd' ? best[2] : best[1];
      const c = v > 0 ? colourMap((v - lo) / (hi - lo)) : [0.4, 0.4, 0.4];
      const a = at(u0), b = at(u1);
      const dx = b[0] - a[0], dy = b[1] - a[1], len = Math.hypot(dx, dy) || 1;
      const nx = -dy / len * w, ny = dx / len * w;
      const quad = [[a[0] + nx, a[1] + ny], [a[0] - nx, a[1] - ny], [b[0] - nx, b[1] - ny], [b[0] + nx, b[1] + ny]];
      for (const i of [0, 1, 2, 0, 2, 3]) {
        acc.pos.push(quad[i][0] - o[0], quad[i][1] - o[1], z);
        acc.col.push(...c);
      }
    }
  }
  for (const [layer, acc] of byLayer) {
    const geom = new THREE.BufferGeometry();
    geom.setAttribute('position', new THREE.BufferAttribute(new Float32Array(acc.pos), 3));
    geom.setAttribute('color', new THREE.BufferAttribute(new Float32Array(acc.col), 3));
    geom.computeBoundingSphere();
    const mesh = new THREE.Mesh(geom, new THREE.MeshBasicMaterial({
      vertexColors: true, side: THREE.DoubleSide, polygonOffset: true, polygonOffsetFactor: -3, polygonOffsetUnits: -6,
      clippingPlanes: $('section-on').checked ? [sectionPlane] : [],
    }));
    mesh.userData.fieldLayer = layer;
    const copper = state.groups.get(`copper/${layer}`);
    mesh.visible = !copper || copper.mesh.visible;
    board.add(mesh);
    state.siMeshes.push(mesh);
  }
  for (let x = 0; x < 240; x++) {
    const [r, g, b] = colourMap(x / 239);
    ctx.fillStyle = `rgb(${255 * r},${255 * g},${255 * b})`;
    ctx.fillRect(x, 0, 1, 14);
  }
  $('si-min').textContent = `${fmt(lo)} Ω`;
  $('si-max').textContent = `${fmt(hi)} Ω`;
  render();
}

// Impedance against distance from the net's first pad, for the net and its partner.
function showSiNet(name) {
  const si = state.board.si, cv = $('si-profile');
  for (const tr of $('si-nets').rows) tr.classList.toggle('selected', tr.dataset.net === name);
  const net = si.nets.find((n) => n.name === name);
  if (!net) { cv.hidden = true; $('si-net').textContent = ''; return; }
  const q = $('si-quantity').value === 'zd' && net.partner ? 'z_differential' : 'z0';
  const series = [net.name, net.partner].filter(Boolean).map((n) =>
    si.samples.filter((s) => s.net === n && s.distance_mm >= 0 && s[q] > 0).map((s) => [s.distance_mm, s[q]]).sort((a, b) => a[0] - b[0]));
  const all = series.flat();
  if (!all.length) { cv.hidden = true; return; }
  cv.hidden = false;
  const ctx = cv.getContext('2d'), W = cv.width, H = cv.height, L = 34, B = 16;
  ctx.clearRect(0, 0, W, H);
  let [lo, hi] = robustRange(all.map((p) => p[1]), 0.05);
  const target = net.target && (net.target.differential === (q === 'z_differential')) ? net.target : null;
  if (target) { lo = Math.min(lo, target.ohm * (1 - target.tolerance)); hi = Math.max(hi, target.ohm * (1 + target.tolerance)); }
  const pad = 0.05 * (hi - lo);
  lo -= pad; hi += pad;
  const dmax = Math.max(...all.map((p) => p[0])) || 1;
  // Values outside the range are drawn on its edge.
  const X = (x) => L + (W - L - 4) * x / dmax;
  const Y = (y) => H - B - (H - B - 4) * (Math.min(hi, Math.max(lo, y)) - lo) / (hi - lo);
  ctx.strokeStyle = '#888'; ctx.fillStyle = '#aaa'; ctx.font = '10px sans-serif';
  ctx.strokeRect(L, 4, W - L - 4, H - B - 4);
  ctx.fillText(fmt(hi), 2, 12); ctx.fillText(fmt(lo), 2, H - B);
  ctx.fillText('0', L, H - 4); ctx.fillText(`${fmt(dmax)} mm`, W - 52, H - 4);
  if (target) {
    ctx.strokeStyle = '#4a6';
    for (const f of [1 - target.tolerance, 1, 1 + target.tolerance]) {
      ctx.beginPath(); ctx.moveTo(L, Y(target.ohm * f)); ctx.lineTo(W - 4, Y(target.ohm * f)); ctx.stroke();
    }
  }
  ['#ffffff', '#ffb000'].forEach((colour, k) => {
    if (!series[k]) return;
    ctx.fillStyle = colour;
    for (const [x, y] of series[k]) ctx.fillRect(X(x) - 1, Y(y) - 1, 2, 2);
  });
  const quantity = q === 'z0' ? 'impedance' : 'differential impedance';
  $('si-net').textContent = `${net.name}: ${quantity} along the net (white${net.partner ? `; ${net.partner} amber` : ''}); ` +
    `${fmt(net.length_mm)} mm, ${fmt(net.delay_ps)} ps, ${fmt(net.loss_db)} dB`;
}

function buildSiPanel() {
  const si = state.board.si;
  if (!si) return;
  $('si').hidden = false;
  const table = $('si-nets');
  const head = document.createElement('tr');
  head.innerHTML = '<td>net</td><td class="num">Z Ω</td><td class="num">Zdiff Ω</td><td class="num">ps</td>';
  table.append(head);
  for (const n of si.nets) {
    const tr = document.createElement('tr');
    tr.dataset.net = n.name;
    tr.innerHTML = `<td>${n.name}</td><td class="num">${n.z_mean > 0 ? fmt(n.z_mean) : ''}</td>` +
      `<td class="num">${n.zd_mean > 0 ? fmt(n.zd_mean) : ''}</td><td class="num">${fmt(n.delay_ps)}</td>`;
    tr.title = `Z ${fmt(n.z_min)} to ${fmt(n.z_max)} Ω` + (n.zd_mean > 0 ? `; Zdiff ${fmt(n.zd_min)} to ${fmt(n.zd_max)} Ω` : '') +
      `; ${fmt(n.length_mm)} mm; layers ${n.layers.join(', ')}`;
    tr.addEventListener('click', () => selectNet(state.board.nets.indexOf(n.name)));
    table.append(tr);
  }
  const counts = {};
  for (const f of si.findings) counts[f.kind] = (counts[f.kind] || 0) + 1;
  $('si-meta').textContent = `${si.nets.length} nets, ${si.cuts} cuts, ${si.cross_sections.length} cross-sections; ` +
    (Object.keys(counts).length ? Object.entries(counts).map(([k, v]) => `${v} ${k}`).join(', ') : 'no findings');
  $('si-quantity').addEventListener('change', () => {
    drawSiTracks();
    showSiNet(state.selectedNet > 0 ? state.board.nets[state.selectedNet] : '');
  });
  drawSiTracks();
}

// ------------------------------------------------------------------ panel

function buildPanel() {
  const b = state.board;
  $('title').textContent = b.name || 'board';
  const thick = b.stackup.layers.reduce((s, l) => s + l.thickness_mm, 0);
  $('meta').textContent = `${b.copper_layers.length} copper layers, ${thick.toFixed(3)} mm, ` +
    `${b.nets.length - 1} nets (${b.net_source}), ${b.drills.length} holes`;

  const layers = $('layers');
  for (const [key, g] of state.groups) {
    const row = document.createElement('label');
    const box = document.createElement('input');
    box.type = 'checkbox';
    box.checked = g.info.kind !== 'mask';
    g.mesh.visible = box.checked;
    box.addEventListener('change', () => {
      g.mesh.visible = box.checked;
      const f = g.info.kind === 'copper' && state.field && state.fields[fieldKind(state.field)].get(g.info.index);
      if (f) f.mesh.visible = box.checked;
      if (g.info.kind === 'copper')
        for (const m of state.siMeshes) if (m.userData.fieldLayer === g.info.index) m.visible = box.checked;
      render();
    });
    const swatch = document.createElement('span');
    swatch.className = 'swatch';
    const colour = g.info.kind === 'dielectric' ? DIELECTRIC_COLOUR : g.info.kind === 'mask' ? MASK_COLOUR : g.base.getHex();
    swatch.style.background = `#${colour.toString(16).padStart(6, '0')}`;
    row.append(box, swatch, document.createTextNode(g.info.label));
    layers.append(row);
  }

  const list = $('net-list');
  b.nets.forEach((n, i) => { if (i > 0) { const o = document.createElement('option'); o.value = n; list.append(o); } });

  if (b.warnings.length) {
    $('warnings-section').hidden = false;
    for (const w of b.warnings) { const li = document.createElement('li'); li.textContent = w; $('warnings').append(li); }
  }
}

function selectNet(id) {
  state.selectedNet = id;
  const b = state.board;
  if (id < 0) {
    $('net-info').textContent = '';
    $('net-search').value = '';
  } else {
    const holes = b.drills.filter((d) => d.net === id).length;
    const pads = b.pads.filter((p) => p.net === id).length;
    $('net-info').textContent = `${pads} pads, ${holes} holes`;
    $('net-search').value = b.nets[id] || '(no net)';
  }
  recolour();
  if (state.field === 'voltage') showField();
  if (state.board.si) showSiNet(id > 0 ? b.nets[id] : '');
}

$('net-search').addEventListener('change', (e) => {
  const id = state.board.nets.indexOf(e.target.value);
  selectNet(id > 0 ? id : -1);
});
$('net-clear').addEventListener('click', () => selectNet(-1));
for (const r of document.querySelectorAll('input[name=colour]'))
  r.addEventListener('change', (e) => { state.colourMode = e.target.value; recolour(); });

// ------------------------------------------------------------------ section and view

function boardBox() {
  const box = new THREE.Box3();
  for (const g of state.groups.values()) box.union(g.mesh.geometry.boundingBox);
  return box;
}

function updateSection() {
  const on = $('section-on').checked;
  const axis = $('section-axis').value;
  const box = boardBox();
  const t = Number($('section-pos').value) / 1000;
  const lo = axis === 'x' ? box.min.x : box.min.y, hi = axis === 'x' ? box.max.x : box.max.y;
  const pos = lo + t * (hi - lo);
  sectionPlane.normal.set(axis === 'x' ? -1 : 0, axis === 'y' ? -1 : 0, 0);
  sectionPlane.constant = pos;
  for (const g of state.groups.values()) g.mesh.material.clippingPlanes = on ? [sectionPlane] : [];
  for (const set of Object.values(state.fields))
    for (const f of set.values()) f.mesh.material.clippingPlanes = on ? [sectionPlane] : [];
  for (const m of state.siMeshes) m.material.clippingPlanes = on ? [sectionPlane] : [];
  const origin = state.board.display.origin_mm;
  $('section-value').textContent = on ? `${axis} = ${(pos + (axis === 'x' ? origin[0] : origin[1])).toFixed(3)} mm` : '';
  render();
}
for (const id of ['section-on', 'section-axis', 'section-pos']) $(id).addEventListener('input', updateSection);

function setView(kind) {
  const box = boardBox();
  const size = box.getSize(new THREE.Vector3());
  const centre = box.getCenter(new THREE.Vector3());
  centre.z *= state.zScale;
  const d = Math.max(size.x, size.y) / (2 * Math.tan((camera.fov * Math.PI) / 360)) * 1.15;
  camera.near = d / 1000;
  camera.far = d * 20;
  camera.updateProjectionMatrix();
  if (kind === 'top') camera.position.set(centre.x, centre.y - 1e-3 * d, centre.z + d);
  else if (kind === 'bottom') camera.position.set(centre.x, centre.y + 1e-3 * d, centre.z - d);
  else camera.position.set(centre.x + 0.6 * d, centre.y - 0.8 * d, centre.z + 0.7 * d);
  controls.target.copy(centre);
  controls.update();
  render();
}
$('view-top').addEventListener('click', () => setView('top'));
$('view-bottom').addEventListener('click', () => setView('bottom'));
$('view-iso').addEventListener('click', () => setView('iso'));
$('z-scale').addEventListener('change', (e) => {
  state.zScale = Number(e.target.value);
  board.scale.set(1, 1, state.zScale);
  render();
});

// ------------------------------------------------------------------ picking

const raycaster = new THREE.Raycaster();
const pointer = new THREE.Vector2();

function pick(event) {
  const rect = canvas.getBoundingClientRect();
  pointer.set(((event.clientX - rect.left) / rect.width) * 2 - 1, -((event.clientY - rect.top) / rect.height) * 2 + 1);
  raycaster.setFromCamera(pointer, camera);
  const targets = [...state.groups.values()]
    .filter((g) => g.mesh.visible && (g.info.kind === 'copper' || g.info.kind === 'drills'))
    .map((g) => g.mesh);
  const planes = $('section-on').checked ? [sectionPlane] : [];
  for (const hit of raycaster.intersectObjects(targets, false)) {
    if (planes.length && planes[0].distanceToPoint(hit.point) < 0) continue;
    const g = state.groups.get(hit.object.userData.key);
    return { g, id: g.ids[hit.face.a], point: hit.point };
  }
  return null;
}

let hoverTimer = 0;
canvas.addEventListener('pointermove', (event) => {
  clearTimeout(hoverTimer);
  hoverTimer = setTimeout(() => {
    const hit = pick(event);
    const o = state.board.display.origin_mm;
    if (!hit) { $('status').textContent = ''; return; }
    const name = hit.id ? state.board.nets[hit.id] : '(no net)';
    $('status').textContent = `${hit.g.info.label}  net ${name}  x ${(hit.point.x + o[0]).toFixed(3)}  y ${(hit.point.y + o[1]).toFixed(3)} mm`;
  }, 40);
});

let downAt = null;
canvas.addEventListener('pointerdown', (e) => { downAt = [e.clientX, e.clientY]; });
canvas.addEventListener('pointerup', (e) => {
  if (!downAt || Math.hypot(e.clientX - downAt[0], e.clientY - downAt[1]) > 4) return;
  if (probeField(e)) return;
  const hit = pick(e);
  selectNet(hit && hit.id > 0 ? hit.id : -1);
});

// ------------------------------------------------------------------ start

async function start() {
  resize();
  const res = await fetch('api/board');
  state.board = await res.json();
  const keys = state.board.display.meshes;
  let done = 0;
  await Promise.all(keys.map(async (key) => {
    await addMesh(key);
    $('loading').textContent = `loading ${++done}/${keys.length}`;
  }));
  // Panel order follows the stackup: copper and dielectrics interleaved top to bottom.
  const order = (key) => {
    const g = state.groups.get(key).info;
    if (g.kind === 'mask') return key.endsWith('top') ? -1 : 1e6;
    if (g.kind === 'drills') return 1e6 + 1;
    if (g.kind === 'dielectric') return g.index;
    return state.board.stackup.layers.findIndex((l) => l.copper_index === g.index);
  };
  state.groups = new Map([...state.groups.entries()].sort((a, b) => order(a[0]) - order(b[0])));
  buildPanel();
  buildSolutionPanel();
  buildSiPanel();
  recolour();
  setView('iso');
  $('loading').hidden = true;
}

start().catch((e) => { $('loading').textContent = `error: ${e.message}`; console.error(e); });
