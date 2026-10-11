# Third-party software and data

KairoQ is distributed under the GNU General Public License, version 3 or later
(`LICENSE`). It builds on the following works, whose licences are compatible with it.

## Libraries compiled into the program

Fetched at build time at the pinned versions in `cmake/Dependencies.cmake`.

| library | version | licence | use |
| --- | --- | --- | --- |
| [Clipper2](https://github.com/AngusJohnson/Clipper2) | 2.0.1 | BSL-1.0 | polygon clipping and offsetting |
| [nlohmann/json](https://github.com/nlohmann/json) | 3.12.0 | MIT | JSON |
| [pugixml](https://github.com/zeux/pugixml) | 1.16 | MIT | XML (IPC-2581) |
| [yaml-cpp](https://github.com/jbeder/yaml-cpp) | 0.9.0 | MIT | YAML |
| [miniz](https://github.com/richgel999/miniz) | 3.1.2 | MIT | zip archives |
| [earcut.hpp](https://github.com/mapbox/earcut.hpp) | 3.2.4 | ISC | triangulation of display meshes |
| [cpp-httplib](https://github.com/yhirose/cpp-httplib) | 0.60.0 | MIT | viewer server |
| [CDT](https://github.com/artem-ogre/CDT) | 2.0.0 | MPL-2.0 | constrained Delaunay triangulation and refinement |
| [AMGCL](https://github.com/ddemidov/amgcl) | 1.5.0 | MIT | algebraic multigrid and conjugate gradients |
| [nanobind](https://github.com/wjakob/nanobind) | ≥ 2.4 | BSD-3-Clause | Python bindings |
| [doctest](https://github.com/doctest/doctest) | 2.5.3 | MIT | C++ tests (not distributed) |

## Files in this repository

| path | origin | licence |
| --- | --- | --- |
| `viewer/vendor/three.module.js`, `three.core.js`, `OrbitControls.js` | [three.js](https://threejs.org) r186 (npm `three` 0.186.1), unmodified | MIT (`viewer/vendor/LICENSE-three.txt`) |
| `testdata/kicad/*` | exports of KiCad 10.0.6 demonstration boards (`testdata/README.md`) | CC BY-SA 4.0; `cm5_minima` also CERN-OHL-S-2.0 |
| `testdata/gerber2ems/*` | slices of the Antmicro Signal Integrity Test Board with VNA measurements, from antmicro/gerber2ems (`testdata/gerber2ems/SOURCE.txt`) | Apache-2.0 |
