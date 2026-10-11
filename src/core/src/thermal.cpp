#include "kairoq/thermal.hpp"

#include "fem_util.hpp"
#include "segment_grid.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <unordered_map>

namespace kairoq {

namespace {

using namespace detail;

constexpr double kStefanBoltzmann = 5.670374419e-8;  // W/(m^2 K^4)
constexpr double kKelvin = 273.15;

// Point-in-polygon tests against a polygon set, on a segment grid (coordinates in mm).
class RegionTester {
public:
    explicit RegionTester(const PolygonSet& polys) {
        std::vector<SegmentGrid::Seg> segs;
        double x1 = -INFINITY, y1 = -INFINITY;
        x0_ = y0_ = INFINITY;
        auto bound = [&](const Path64& r) {
            for (const auto& a : r) {
                x0_ = std::min(x0_, a.x * 1e-6);
                y0_ = std::min(y0_, a.y * 1e-6);
                x1 = std::max(x1, a.x * 1e-6);
                y1 = std::max(y1, a.y * 1e-6);
            }
        };
        for (const auto& p : polys) bound(p.outer);
        auto ring = [&](const Path64& r) {
            for (std::size_t k = 0; k < r.size(); ++k) {
                const auto& a = r[k];
                const auto& c = r[(k + 1) % r.size()];
                segs.push_back({a.x * 1e-6 - x0_, a.y * 1e-6 - y0_, c.x * 1e-6 - x0_, c.y * 1e-6 - y0_, true, 0.0});
            }
        };
        for (const auto& p : polys) {
            ring(p.outer);
            for (const auto& h : p.holes) ring(h);
        }
        if (!segs.empty()) grid_ = std::make_unique<SegmentGrid>(std::move(segs), 0.0, 0.0, x1 - x0_, y1 - y0_, 0.2);
    }
    bool inside(double x_nm, double y_nm) const {
        if (!grid_) return false;
        const double x = x_nm * 1e-6 - x0_, y = y_nm * 1e-6 - y0_;
        if (x < 0 || y < 0) return false;
        return grid_->inside(x, y);
    }

private:
    std::unique_ptr<SegmentGrid> grid_;
    double x0_ = 0, y0_ = 0;
};

// Board data that does not change between passes.
struct Model {
    int n = 0;  // copper layers
    // Conduction sheets from the top to the bottom (1-based): one per copper layer and,
    // in thick dielectrics, more between them.
    struct Sheet {
        int copper = 0;       // copper layer, or 0 for a dielectric sheet
        double t_cu = 0.0;    // copper thickness (m)
        double share = 0.0;   // in-plane dielectric conductance (W/K)
        double fill = 0.0;    // in-plane conductance of the resin between the copper (W/K)
        double r_down = 0.0;  // resistance per area to the next sheet (m^2 K/W)
        // Heat capacity per area (J/(m^2 K)): of the copper, of the resin between it, and
        // the share of the dielectric (and, on the outer sheets, of the mask).
        double c_cu = 0.0, c_fill = 0.0, c_share = 0.0;
    };
    std::vector<Sheet> sheets;
    std::vector<int> sheet_of;       // copper layer -> sheet
    double r_top = 0, r_bottom = 0;  // resistance per area from the outer sheets to the surfaces
    PolygonSet domain;               // board outline without non-plated holes
    std::vector<PolygonSet> copper;  // per layer
    std::vector<std::unique_ptr<RegionTester>> copper_test;
    struct Hole {
        std::size_t drill;
        std::vector<int> layers;  // layers whose copper touches the hole wall
        PolygonSet outline;
        double conductance_per_m = 0.0;  // k A of the barrel (W m / K)
    };
    std::vector<Hole> holes;
    std::vector<std::map<int, PolygonSet>> areas;  // per heat spec, per layer
    std::vector<double> zmid;                      // mid-plane height per layer (m)
    struct SinkRegion {
        std::size_t sink;
        bool bottom;
        PolygonSet area;                      // empty: the whole face
        std::unique_ptr<RegionTester> test;   // null for the whole face
    };
    std::vector<SinkRegion> sink_regions;
};

// One mesh: a domain part on a sheet.
struct Piece {
    int sheet = 0;
    int layer = 0;  // copper layer of the sheet, or 0
    std::size_t part = 0;
    TriMesh mesh;
    std::vector<std::uint8_t> copper;
    std::vector<double> conductance;  // in-plane, per triangle (W/K)
    std::vector<double> capacity;     // per area, per triangle (J/(m^2 K))
    std::vector<double> node_area;    // lumped, m^2
    std::size_t offset = 0;           // first global node
    std::unique_ptr<TriangleLocator> locator;
};

using Pieces = std::vector<std::unique_ptr<Piece>>;

// Slices of a dielectric of thickness t (m): graded from `first` at both faces, doubling up
// to `largest`, which fills the middle.
std::vector<double> slice_thicknesses(double t, double first, double largest) {
    first = std::min(first, largest);
    std::vector<double> side;
    double rest = t, step = first;
    while (step < largest && rest > 2.0 * step * (1.0 + 1e-9)) {
        side.push_back(step);
        rest -= 2.0 * step;
        step = std::min(2.0 * step, largest);
    }
    const int middle = std::max(1, static_cast<int>(std::ceil(rest / largest - 1e-9)));
    std::vector<double> out(side);
    for (int k = 0; k < middle; ++k) out.push_back(rest / middle);
    out.insert(out.end(), side.rbegin(), side.rend());
    return out;
}

Model build_model(const Board& b, const ThermalSettings& s, double tol) {
    Model m;
    m.n = static_cast<int>(b.copper.size());
    const int n = m.n;
    m.zmid.assign(n + 1, 0.0);
    const auto zs = stackup_z(b);
    // Vertical elements from the top: copper layers and dielectric slices. Dielectric is
    // cut into slices that start at dielectric_first_slice at both faces and double
    // towards the middle up to dielectric_slice; a sheet sits on every copper layer, at
    // every boundary between slices and on an outer face of dielectric. Each slice gives
    // half its in-plane conduction (and heat capacity) to the sheet above and half to the
    // sheet below, and its through-thickness resistance between them. Solder mask outside
    // the outermost element is a series resistance to the surface.
    struct Element {
        bool copper;
        int index;
        double t, kxy, kz, c;
    };
    std::vector<Element> el;
    double c_mask_above = 0, c_mask_below = 0;
    bool seen = false;
    for (std::size_t i = 0; i < b.stackup.size(); ++i) {
        const auto& l = b.stackup[i];
        const double t = l.thickness_nm * 1e-9;
        if (l.type == LayerType::Copper) {
            m.zmid[l.copper_index] = 0.5 * (zs[i].top_nm + zs[i].bottom_nm) * 1e-9;
            el.push_back({true, l.copper_index, t, 0, 0, s.copper_heat_capacity});
            seen = true;
            continue;
        }
        if (t <= 0) continue;
        const bool mask = l.type == LayerType::SolderMask;
        double kxy = mask ? s.mask_k : s.dielectric_k_xy, kz = mask ? s.mask_k : s.dielectric_k_z;
        double c = mask ? s.mask_heat_capacity : s.dielectric_heat_capacity;
        if (auto it = s.layer_k.find(l.name); it != s.layer_k.end()) {
            kxy = it->second[0];
            kz = it->second[1];
        }
        if (auto it = s.layer_heat_capacity.find(l.name); it != s.layer_heat_capacity.end()) c = it->second;
        if (mask) {
            (seen ? m.r_bottom : m.r_top) += t / kz;
            (seen ? c_mask_below : c_mask_above) += t * c;
            continue;
        }
        for (double tau : slice_thicknesses(t, s.dielectric_first_slice_nm * 1e-9, s.dielectric_slice_nm * 1e-9))
            el.push_back({false, 0, tau, kxy, kz, c});
        seen = true;
    }
    m.sheets.assign(1, Model::Sheet{});
    m.sheet_of.assign(n + 1, 0);
    std::optional<Element> pending;  // slice waiting for the sheet below it
    int cur = 0;
    auto new_sheet = [&](int copper, double t_cu, double c_cu) {
        Model::Sheet sh;
        sh.copper = copper;
        sh.t_cu = t_cu;
        sh.c_cu = c_cu * t_cu;
        m.sheets.push_back(sh);
        const int k = static_cast<int>(m.sheets.size()) - 1;
        if (pending) {
            for (int j : {cur, k}) {
                m.sheets[j].share += 0.5 * pending->t * pending->kxy;
                m.sheets[j].c_share += 0.5 * pending->t * pending->c;
            }
            m.sheets[cur].r_down += pending->t / pending->kz;
            pending.reset();
        }
        cur = k;
        if (copper) m.sheet_of[copper] = k;
    };
    for (const auto& e : el) {
        if (e.copper) {
            new_sheet(e.index, e.t, e.c);
            continue;
        }
        if (cur == 0 || pending) new_sheet(0, 0.0, 0.0);  // top face, or between slices
        pending = e;
    }
    if (pending) new_sheet(0, 0.0, 0.0);  // bottom face
    // Inner copper layers are embedded in resin, taken to be the dielectric below them;
    // on outer layers the space between the copper is air, or mask filling it.
    for (std::size_t i = 0; i < el.size(); ++i) {
        if (!el[i].copper) continue;
        auto& sh = m.sheets[m.sheet_of[el[i].index]];
        const bool outer = i == 0 || i + 1 >= el.size();
        if (outer) {
            if (s.outer_gaps_mask) {
                sh.fill = s.mask_k * el[i].t;
                sh.c_fill = s.mask_heat_capacity * el[i].t;
            }
            continue;
        }
        if (el[i - 1].copper || el[i + 1].copper) continue;
        sh.fill = el[i + 1].kxy * el[i].t;
        sh.c_fill = el[i + 1].c * el[i].t;
    }
    m.sheets[1].c_share += c_mask_above;
    m.sheets.back().c_share += c_mask_below;
    for (int l = 1; l <= n; ++l)
        if (m.sheet_of[l] == 0 || m.sheets[m.sheet_of[l]].t_cu <= 0)
            throw std::runtime_error("copper layer " + std::to_string(l) + " has no thickness in the stackup");
    for (std::size_t k = 1; k + 1 < m.sheets.size(); ++k)
        if (m.sheets[k].r_down <= 0) throw std::runtime_error("no dielectric between two copper layers");

    // Domain: the outline without non-plated holes.
    Paths64 npth;
    for (const auto& d : b.drills)
        if (!d.plated) npth.push_back(hole_outline(d, tol));
    m.domain = npth.empty() ? b.outline : difference(b.outline, union_polygons(npth));

    m.copper.resize(n + 1);
    m.copper_test.resize(n + 1);
    for (const auto& cl : b.copper) {
        Paths64 all;
        for (const auto& nc : cl.nets) {
            const auto p = to_paths(nc.polygons);
            all.insert(all.end(), p.begin(), p.end());
        }
        m.copper[cl.index] = intersection(union_polygons(all), m.domain);
        m.copper_test[cl.index] = std::make_unique<RegionTester>(m.copper[cl.index]);
    }

    // Plated holes: the layers whose copper touches them, and the barrel.
    const double tp = s.plating_thickness_nm;
    for (std::size_t d = 0; d < b.drills.size(); ++d) {
        const auto& dr = b.drills[d];
        if (!dr.plated) continue;
        Model::Hole h;
        h.drill = d;
        for (int l = std::clamp(dr.from_layer, 1, n); l <= std::clamp(dr.to_layer, 1, n); ++l)
            if (m.copper_test[l]->inside(double(dr.at.x), double(dr.at.y))) h.layers.push_back(l);
        if (h.layers.empty()) continue;
        h.outline = union_polygons({hole_outline(dr, tol)});
        const double slot = dr.slot_end ? std::hypot(double(dr.slot_end->x - dr.at.x), double(dr.slot_end->y - dr.at.y)) : 0.0;
        const double area = tp * (kPi * (dr.diameter_nm + tp) + 2.0 * slot) * 1e-18;  // m^2
        h.conductance_per_m = s.copper_k * area;
        m.holes.push_back(std::move(h));
    }

    // Heat areas.
    for (const auto& hs : s.sources) {
        const PortShape shape = resolve_port(b, hs.area, tol);
        std::map<int, PolygonSet> clipped;
        for (const auto& [layer, poly] : shape.by_layer) {
            auto c = intersection(poly, m.domain);
            if (!c.empty()) clipped[layer] = std::move(c);
        }
        if (clipped.empty()) throw std::runtime_error(hs.area.name + ": the area lies outside the board");
        m.areas.push_back(std::move(clipped));
    }

    // Heat sink regions.
    for (std::size_t k = 0; k < s.heatsinks.size(); ++k)
        for (const auto& r : s.heatsinks[k].regions) {
            Model::SinkRegion sr;
            sr.sink = k;
            sr.bottom = r.bottom;
            if (!r.area.empty()) {
                sr.area = intersection(r.area, m.domain);
                if (sr.area.empty()) throw std::runtime_error("heat sink " + s.heatsinks[k].name + ": a region lies outside the board");
                sr.test = std::make_unique<RegionTester>(sr.area);
            }
            m.sink_regions.push_back(std::move(sr));
        }
    return m;
}

// Heat flux (W/m^2) from a sheet at temperature t through the series resistance r (mask,
// outer dielectric) to a surface cooled by convection and radiation into ambient ta. The
// surface temperature follows by bisection.
double surface_flux(const SurfaceSpec& sp, double r, double t, double ta) {
    auto out = [&](double ts) {
        const double tk = ts + kKelvin, tak = ta + kKelvin;
        return sp.convection_w_m2k * (ts - ta) + sp.emissivity * kStefanBoltzmann * (tk * tk * tk * tk - tak * tak * tak * tak);
    };
    if (r <= 0) return out(t);
    double lo = std::min(t, ta), hi = std::max(t, ta);
    for (int i = 0; i < 60; ++i) {
        const double mid = 0.5 * (lo + hi);
        // Flux through r equals the flux from the surface.
        if ((t - mid) / r > out(mid)) lo = mid;
        else hi = mid;
    }
    return (t - 0.5 * (lo + hi)) / r;
}

double surface_slope(const SurfaceSpec& sp, double r, double t, double ta) {
    return (surface_flux(sp, r, t + 0.01, ta) - surface_flux(sp, r, t - 0.01, ta)) / 0.02;
}

// Interpolation of a nodal field over the pieces of each sheet; owns its meshes.
struct Field {
    std::shared_ptr<const Pieces> pieces;
    std::vector<double> nodal;  // per global node
    double fallback = 0.0;
    std::vector<int> sheet_of;  // copper layer -> sheet
    double at(int layer, double x, double y) const { return at_sheet(sheet_of[layer], x, y); }
    double at_sheet(int sheet, double x, double y) const {
        for (const auto& p : *pieces) {
            if (p->sheet != sheet) continue;
            std::array<double, 3> w;
            const long t = p->locator->locate(x, y, w, 50000.0);
            if (t < 0) continue;
            const auto& tr = p->mesh.triangles[static_cast<std::size_t>(t)];
            return w[0] * nodal[p->offset + tr[0]] + w[1] * nodal[p->offset + tr[1]] + w[2] * nodal[p->offset + tr[2]];
        }
        return fallback;
    }
};

double power_at(const HeatSpec& h, double t) {
    const auto& p = h.power_schedule;
    if (p.empty()) return h.power_w.value_or(0.0);
    if (t <= p.front()[0]) return p.front()[1];
    for (std::size_t k = 1; k < p.size(); ++k)
        if (t <= p[k][0]) {
            const double span = p[k][0] - p[k - 1][0];
            return span > 0 ? p[k - 1][1] + (p[k][1] - p[k - 1][1]) * (t - p[k - 1][0]) / span : p[k][1];
        }
    return p.back()[1];
}

double steady_power(const HeatSpec& h) { return h.power_w ? *h.power_w : power_at(h, 0.0); }

double peak_power(const HeatSpec& h) {
    double p = h.power_w.value_or(0.0);
    for (const auto& q : h.power_schedule) p = std::max(p, q[1]);
    return p;
}

// Unknowns, couplings and loads on a fixed mesh.
struct Layout {
    std::shared_ptr<Pieces> pieces;
    std::size_t nnodes = 0, ntri = 0, ndof = 0;
    std::vector<std::size_t> dof;  // node -> unknown
    std::vector<std::vector<long>> hole_node;
    std::vector<long> area_dof;  // fixed-temperature source -> its unknown
    std::vector<std::vector<std::pair<std::size_t, double>>> area_nodes;  // per source: node, lumped area
    std::vector<long> junction;                                            // per source: unknown, or -1
    std::vector<std::size_t> sink_dof;                                     // per heat sink
    std::vector<std::vector<std::pair<std::size_t, double>>> source_weight;  // per source: unknown, share
    // Faces: per node the area cooled by convection and radiation (0 inside the board), and
    // the contacts with the heat sinks (node, sink, conductance W/K).
    std::vector<double> face_area;
    std::vector<std::uint8_t> face;  // 0 none, 1 top, 2 bottom
    struct Contact {
        std::size_t node, sink;
        double g;
    };
    std::vector<Contact> contacts;
    std::vector<double> capacity;  // per unknown, J/K
    std::vector<double> fixed;     // per unknown: temperature, or NaN
    std::size_t unlocated = 0;
};

// The matrix without the faces, and its restriction to the free unknowns.
struct System {
    Csr k;                           // static conductance matrix, every diagonal present
    std::vector<double> rhs_static;  // loads that do not change (sink bodies to ambient)
    std::vector<long> free_index;    // unknown -> free row, or -1
    std::vector<std::size_t> free_dofs;
    std::vector<std::ptrdiff_t> ptr, col;  // free part
    std::vector<std::size_t> pos;          // index of each free entry in k
    std::vector<std::size_t> diag;         // index in the free part of each row's diagonal
    std::vector<double> b_fixed;           // per free row: coupling to the fixed temperatures
    // Per node and side (0 the sheet below, 1 the sheet above): the unknowns and weights of
    // the point it couples to, and the conductance per area (0: none).
    struct Counterpart {
        std::array<std::uint32_t, 3> dof{};
        std::array<float, 3> w{};
        float gz = 0;
    };
    std::vector<std::array<Counterpart, 2>> counterpart;
};

}  // namespace

namespace {

class Run {
public:
    Run(const Board& b, const ThermalSettings& s) : b_(b), s_(s) {}
    ThermalResult run();

private:
    const Board& b_;
    const ThermalSettings& s_;
    double tol_ = 1000.0;
    std::unique_ptr<Model> m_;
    int ns_ = 0;
    double ta_ = 25.0;
    bool coupled_ = false, radiation_ = false;
    std::vector<PolygonSet> constraints_;
    std::size_t nparts_ = 0;
    std::vector<std::vector<std::array<double, 3>>> seeds_;
    Layout lay_;
    System sys_;
    ThermalResult res_;
    std::vector<std::vector<std::array<double, 3>>> dc_seeds_;
    std::vector<double> peak_t_;  // nodal temperature at the peak powers (transients)

    std::size_t owner(int sheet, std::size_t part) const { return static_cast<std::size_t>(sheet) * nparts_ + part; }
    const SurfaceSpec& face_spec(std::uint8_t f) const { return f == 1 ? s_.top : s_.bottom; }
    double face_r(std::uint8_t f) const { return f == 1 ? m_->r_top : m_->r_bottom; }

    std::shared_ptr<Pieces> mesh_pieces(bool& any_unrefined);
    void make_layout(std::shared_ptr<Pieces> pieces);
    void make_system();
    std::vector<double> nodal(const std::vector<double>& x) const;  // unknowns -> nodes
    std::vector<double> joule_load(const DcResult& dc, double& placed, std::size_t& lost) const;
    std::vector<double> source_load(const std::vector<double>& power) const;
    // Face terms linearised about the nodal temperatures t0 with the slopes g (per node,
    // W/(m^2 K)): the diagonal per unknown and the load.
    void face_terms(const std::vector<double>& t0, const std::vector<double>& g, std::vector<double>& diag,
                    std::vector<double>& load) const;
    std::vector<double> face_slopes(const std::vector<double>& t0) const;
    std::vector<double> free_values(const std::vector<double>& diag) const;
    // Heat into each fixed unknown (out of the board), from the full system.
    std::vector<double> fixed_heat(const std::vector<double>& x, const std::vector<double>& diag,
                                   const std::vector<double>& load) const;
    // Adds the estimate of the vertical coupling's error per triangle of a piece.
    void add_vertical_estimate(const Piece& pc, const std::vector<double>& x, const std::vector<char>& skip, std::vector<double>& eta2,
                               double& total) const;
    void steady(const std::vector<double>& power, std::vector<double>& x, bool adapt);
    void report(const std::vector<double>& x, const std::vector<double>& power, double joule_w, bool final_fields);
    void transient(std::vector<double> x);
    std::vector<double> initial_state() const;
    Field field(const std::vector<double>& x) const { return Field{lay_.pieces, nodal(x), ta_, m_->sheet_of}; }
};

std::shared_ptr<Pieces> Run::mesh_pieces(bool& any_unrefined) {
    auto pieces = std::make_shared<Pieces>();
    std::vector<std::pair<int, std::size_t>> jobs;
    for (int k = 1; k <= ns_; ++k)
        for (std::size_t p = 0; p < nparts_; ++p) jobs.emplace_back(k, p);
    pieces->resize(jobs.size());
    std::atomic<std::size_t> unrefined{0};
    parallel_for(jobs.size(), [&](std::size_t j) {
        auto pc = std::make_unique<Piece>();
        pc->sheet = jobs[j].first;
        pc->layer = m_->sheets[pc->sheet].copper;
        pc->part = jobs[j].second;
        pc->mesh = mesh_region(m_->domain[pc->part], constraints_[pc->sheet], s_.mesh, seeds_[owner(pc->sheet, pc->part)]);
        unrefined += pc->mesh.unrefined;
        const TriMesh& tm = pc->mesh;
        const auto& sh = m_->sheets[pc->sheet];
        pc->copper.resize(tm.triangles.size());
        pc->conductance.resize(tm.triangles.size());
        pc->capacity.resize(tm.triangles.size());
        pc->node_area.assign(tm.x.size(), 0.0);
        for (std::size_t t = 0; t < tm.triangles.size(); ++t) {
            const auto& tr = tm.triangles[t];
            const double cx = (tm.x[tr[0]] + tm.x[tr[1]] + tm.x[tr[2]]) / 3, cy = (tm.y[tr[0]] + tm.y[tr[1]] + tm.y[tr[2]]) / 3;
            pc->copper[t] = pc->layer > 0 && m_->copper_test[pc->layer]->inside(cx, cy);
            pc->conductance[t] = (pc->copper[t] ? s_.copper_k * sh.t_cu : sh.fill) + sh.share;
            pc->capacity[t] = (pc->copper[t] ? sh.c_cu : sh.c_fill) + sh.c_share;
            const double a = triangle_area(tm, t) * 1e-18;
            for (auto v : tr) pc->node_area[v] += a / 3.0;
        }
        pc->locator = std::make_unique<TriangleLocator>(pc->mesh);
        (*pieces)[j] = std::move(pc);
    });
    any_unrefined = unrefined > 0;
    return pieces;
}

void Run::make_layout(std::shared_ptr<Pieces> pieces_sp) {
    Layout L;
    L.pieces = std::move(pieces_sp);
    const Pieces& pieces = *L.pieces;
    const int n = m_->n;
    for (auto& pc : pieces) {
        pc->offset = L.nnodes;
        L.nnodes += pc->mesh.x.size();
        L.ntri += pc->mesh.triangles.size();
    }
    const std::size_t nsrc = s_.sources.size();
    // Tied nodes: hole walls per layer, fixed-temperature areas.
    UnionFind uf(L.nnodes);
    L.hole_node.assign(m_->holes.size(), std::vector<long>(n + 1, -1));
    std::vector<long> area_node(nsrc, -1);
    L.area_nodes.assign(nsrc, {});
    for (auto& pc : pieces) {
        NodeGrid grid(pc->mesh, 2.0e5);
        for (std::size_t h = 0; h < m_->holes.size(); ++h) {
            const auto& hl = m_->holes[h];
            if (std::find(hl.layers.begin(), hl.layers.end(), pc->layer) == hl.layers.end()) continue;
            for (auto v : grid.inside(hl.outline)) {
                const std::size_t g = pc->offset + v;
                if (L.hole_node[h][pc->layer] < 0) L.hole_node[h][pc->layer] = static_cast<long>(g);
                else uf.unite(g, static_cast<std::size_t>(L.hole_node[h][pc->layer]));
            }
        }
        for (std::size_t a = 0; a < nsrc; ++a) {
            auto it = m_->areas[a].find(pc->layer);
            if (pc->layer == 0 || it == m_->areas[a].end()) continue;
            for (auto v : grid.inside(it->second)) {
                const std::size_t g = pc->offset + v;
                L.area_nodes[a].push_back({g, pc->node_area[v]});
                if (s_.sources[a].temperature_c) {
                    if (area_node[a] < 0) area_node[a] = static_cast<long>(g);
                    else uf.unite(g, static_cast<std::size_t>(area_node[a]));
                }
            }
        }
    }
    for (std::size_t a = 0; a < nsrc; ++a)
        if (L.area_nodes[a].empty()) throw std::runtime_error(s_.sources[a].area.name + " covers no mesh node");
    L.dof.assign(L.nnodes, 0);
    std::unordered_map<std::size_t, std::size_t> root_dof;
    for (std::size_t g = 0; g < L.nnodes; ++g) {
        auto [it, ins] = root_dof.emplace(uf.find(g), root_dof.size());
        L.dof[g] = it->second;
    }
    L.ndof = root_dof.size();
    L.junction.assign(nsrc, -1);
    for (std::size_t a = 0; a < nsrc; ++a)
        if (s_.sources[a].junction_k_per_w && !s_.sources[a].temperature_c) L.junction[a] = static_cast<long>(L.ndof++);
    for (std::size_t k = 0; k < s_.heatsinks.size(); ++k) L.sink_dof.push_back(L.ndof++);
    L.area_dof.assign(nsrc, -1);
    for (std::size_t a = 0; a < nsrc; ++a)
        if (area_node[a] >= 0) L.area_dof[a] = static_cast<long>(L.dof[static_cast<std::size_t>(area_node[a])]);

    // Fixed temperatures.
    L.fixed.assign(L.ndof, NAN);
    for (std::size_t a = 0; a < nsrc; ++a)
        if (s_.sources[a].temperature_c) L.fixed[static_cast<std::size_t>(L.area_dof[a])] = *s_.sources[a].temperature_c;
    for (std::size_t k = 0; k < s_.heatsinks.size(); ++k)
        if (s_.heatsinks[k].temperature_c) L.fixed[L.sink_dof[k]] = *s_.heatsinks[k].temperature_c;

    // Faces: per triangle of the outer sheets, a heat sink region or the open face.
    L.face_area.assign(L.nnodes, 0.0);
    L.face.assign(L.nnodes, 0);
    std::map<std::pair<std::size_t, std::size_t>, double> contact_area;  // (node, sink) -> area
    for (const auto& pc : pieces) {
        const bool top = pc->sheet == 1, bottom = pc->sheet == ns_;
        if (!top && !bottom) continue;
        const TriMesh& tm = pc->mesh;
        for (std::size_t t = 0; t < tm.triangles.size(); ++t) {
            const auto& tr = tm.triangles[t];
            const double cx = (tm.x[tr[0]] + tm.x[tr[1]] + tm.x[tr[2]]) / 3, cy = (tm.y[tr[0]] + tm.y[tr[1]] + tm.y[tr[2]]) / 3;
            const double a3 = triangle_area(tm, t) * 1e-18 / 3.0;
            // A sheet that is both the top and the bottom one (no dielectric) cannot occur.
            const bool is_bottom = bottom && !top;
            long sink = -1;
            for (const auto& r : m_->sink_regions)
                if (r.bottom == is_bottom && (!r.test || r.test->inside(cx, cy))) {
                    sink = static_cast<long>(r.sink);
                    break;
                }
            for (auto v : tr) {
                const std::size_t g = pc->offset + v;
                L.face[g] = is_bottom ? 2 : 1;
                if (sink < 0) L.face_area[g] += a3;
                else contact_area[{g, static_cast<std::size_t>(sink)}] += a3;
            }
        }
    }
    for (const auto& [key, area] : contact_area) {
        const auto& hs = s_.heatsinks[key.second];
        const double r = face_r(L.face[key.first]) + hs.tim_m2k_per_w;
        L.contacts.push_back({key.first, key.second, r > 0 ? area / r : area * 1e12});
    }
    for (std::size_t k = 0; k < s_.heatsinks.size(); ++k) {
        bool any = false;
        for (const auto& c : L.contacts) any = any || c.sink == k;
        if (!any) throw std::runtime_error("heat sink " + s_.heatsinks[k].name + " touches no part of the board");
    }

    // Loads of the sources: by lumped area over their area, or into the junction.
    L.source_weight.assign(nsrc, {});
    for (std::size_t a = 0; a < nsrc; ++a) {
        if (s_.sources[a].temperature_c) continue;
        if (L.junction[a] >= 0) {
            L.source_weight[a].push_back({static_cast<std::size_t>(L.junction[a]), 1.0});
            continue;
        }
        double total = 0;
        for (const auto& [g, ar] : L.area_nodes[a]) total += ar;
        std::map<std::size_t, double> w;
        for (const auto& [g, ar] : L.area_nodes[a]) w[L.dof[g]] += ar / total;
        L.source_weight[a].assign(w.begin(), w.end());
    }

    // Heat capacity per unknown.
    L.capacity.assign(L.ndof, 0.0);
    for (const auto& pc : pieces) {
        const TriMesh& tm = pc->mesh;
        for (std::size_t t = 0; t < tm.triangles.size(); ++t) {
            const double c3 = pc->capacity[t] * triangle_area(tm, t) * 1e-18 / 3.0;
            for (auto v : tm.triangles[t]) L.capacity[L.dof[pc->offset + v]] += c3;
        }
    }
    for (std::size_t a = 0; a < nsrc; ++a) {
        const double c = s_.sources[a].heat_capacity_j_per_k;
        if (c <= 0 || s_.sources[a].temperature_c) continue;
        for (const auto& [d, w] : L.source_weight[a]) L.capacity[d] += c * w;
    }
    for (std::size_t k = 0; k < s_.heatsinks.size(); ++k) L.capacity[L.sink_dof[k]] += s_.heatsinks[k].heat_capacity_j_per_k;
    for (std::size_t d = 0; d < L.ndof; ++d)
        if (!std::isnan(L.fixed[d])) L.capacity[d] = 0.0;
    lay_ = std::move(L);
}

void Run::make_system() {
    const Layout& L = lay_;
    const Pieces& pieces = *L.pieces;
    System S;
    std::vector<Triplet> trip;
    trip.reserve(9 * L.ntri + 34 * L.nnodes);
    for (std::size_t d = 0; d < L.ndof; ++d) trip.push_back({d, d, 0.0});
    // In-plane conduction.
    for (const auto& pc : pieces) {
        const TriMesh& mm = pc->mesh;
        for (std::size_t t = 0; t < mm.triangles.size(); ++t) {
            const auto& tr = mm.triangles[t];
            const double x1 = mm.x[tr[0]], y1 = mm.y[tr[0]], x2 = mm.x[tr[1]], y2 = mm.y[tr[1]], x3 = mm.x[tr[2]], y3 = mm.y[tr[2]];
            const double bb[3] = {y2 - y3, y3 - y1, y1 - y2};
            const double cc[3] = {x3 - x2, x1 - x3, x2 - x1};
            const double a2 = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
            const double g = pc->conductance[t];
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    trip.push_back({L.dof[pc->offset + tr[i]], L.dof[pc->offset + tr[j]], g * (bb[i] * bb[j] + cc[i] * cc[j]) / (2.0 * a2)});
        }
    }
    // Vertical conduction between neighbouring sheets on non-matching meshes: each node
    // couples, over half its lumped area, to the point below (or above) it.
    std::size_t unlocated = 0;
    S.counterpart.assign(L.nnodes, {});
    for (int l = 1; l < ns_; ++l) {
        const double gz = 1.0 / m_->sheets[l].r_down;
        for (int dir = 0; dir < 2; ++dir) {
            const int from = dir == 0 ? l : l + 1, to = dir == 0 ? l + 1 : l;
            for (const auto& pa : pieces) {
                if (pa->sheet != from) continue;
                for (std::size_t v = 0; v < pa->mesh.x.size(); ++v) {
                    const double x = pa->mesh.x[v], y = pa->mesh.y[v];
                    std::array<double, 3> w;
                    long tri = -1;
                    const Piece* pb = nullptr;
                    for (const auto& cand : pieces) {
                        if (cand->sheet != to) continue;
                        tri = cand->locator->locate(x, y, w, 4.0 * s_.mesh.boundary_tolerance_nm + 1000.0);
                        if (tri >= 0) {
                            pb = cand.get();
                            break;
                        }
                    }
                    if (tri < 0) {
                        ++unlocated;
                        continue;
                    }
                    const double g = 0.5 * gz * pa->node_area[v];
                    const auto& tr = pb->mesh.triangles[static_cast<std::size_t>(tri)];
                    std::array<std::size_t, 4> ids{L.dof[pa->offset + v], L.dof[pb->offset + tr[0]], L.dof[pb->offset + tr[1]],
                                                   L.dof[pb->offset + tr[2]]};
                    const std::array<double, 4> u{1.0, -w[0], -w[1], -w[2]};
                    for (int i = 0; i < 4; ++i)
                        for (int j = 0; j < 4; ++j) trip.push_back({ids[i], ids[j], g * u[i] * u[j]});
                    auto& cp = S.counterpart[pa->offset + v][dir == 0 ? 0 : 1];
                    for (int i = 0; i < 3; ++i) {
                        cp.dof[static_cast<std::size_t>(i)] = static_cast<std::uint32_t>(ids[static_cast<std::size_t>(i) + 1]);
                        cp.w[static_cast<std::size_t>(i)] = static_cast<float>(w[static_cast<std::size_t>(i)]);
                    }
                    cp.gz = static_cast<float>(gz);
                }
            }
        }
    }
    lay_.unlocated = unlocated;
    auto couple = [&](std::size_t a, std::size_t c, double g) {
        trip.push_back({a, a, g});
        trip.push_back({c, c, g});
        trip.push_back({a, c, -g});
        trip.push_back({c, a, -g});
    };
    // Barrels of plated holes.
    for (std::size_t h = 0; h < m_->holes.size(); ++h) {
        const auto& hl = m_->holes[h];
        for (std::size_t k = 1; k < hl.layers.size(); ++k) {
            const int la = hl.layers[k - 1], lb = hl.layers[k];
            if (L.hole_node[h][la] < 0 || L.hole_node[h][lb] < 0) continue;
            couple(L.dof[static_cast<std::size_t>(L.hole_node[h][la])], L.dof[static_cast<std::size_t>(L.hole_node[h][lb])],
                   hl.conductance_per_m / std::abs(m_->zmid[la] - m_->zmid[lb]));
        }
    }
    // Junctions: the junction-to-case resistance spread over the area by lumped area.
    for (std::size_t a = 0; a < s_.sources.size(); ++a) {
        if (L.junction[a] < 0) continue;
        double total = 0;
        for (const auto& [g, ar] : L.area_nodes[a]) total += ar;
        const double gj = 1.0 / *s_.sources[a].junction_k_per_w;
        for (const auto& [g, ar] : L.area_nodes[a]) couple(static_cast<std::size_t>(L.junction[a]), L.dof[g], gj * ar / total);
    }
    // Heat sinks: the interface to the sink, and the sink body to ambient.
    S.rhs_static.assign(L.ndof, 0.0);
    for (const auto& c : L.contacts) couple(L.dof[c.node], L.sink_dof[c.sink], c.g);
    for (std::size_t k = 0; k < s_.heatsinks.size(); ++k) {
        const auto& hs = s_.heatsinks[k];
        if (hs.temperature_c) continue;
        const double g = 1.0 / *hs.resistance_k_per_w;
        trip.push_back({L.sink_dof[k], L.sink_dof[k], g});
        S.rhs_static[L.sink_dof[k]] += g * ta_;
    }
    S.k = to_csr(L.ndof, trip);
    trip.clear();
    trip.shrink_to_fit();

    // Restriction to the free unknowns.
    S.free_index.assign(L.ndof, -1);
    for (std::size_t d = 0; d < L.ndof; ++d)
        if (std::isnan(L.fixed[d])) {
            S.free_index[d] = static_cast<long>(S.free_dofs.size());
            S.free_dofs.push_back(d);
        }
    const std::size_t nf = S.free_dofs.size();
    S.ptr.assign(nf + 1, 0);
    S.diag.assign(nf, 0);
    S.b_fixed.assign(nf, 0.0);
    for (std::size_t r = 0; r < nf; ++r) {
        const std::size_t d = S.free_dofs[r];
        for (auto k = S.k.ptr[d]; k < S.k.ptr[d + 1]; ++k) {
            const auto c = static_cast<std::size_t>(S.k.col[k]);
            if (S.free_index[c] >= 0) {
                if (c == d) S.diag[r] = S.col.size();
                S.col.push_back(S.free_index[c]);
                S.pos.push_back(static_cast<std::size_t>(k));
            } else {
                S.b_fixed[r] += S.k.val[static_cast<std::size_t>(k)] * L.fixed[c];
            }
        }
        S.ptr[r + 1] = static_cast<std::ptrdiff_t>(S.col.size());
    }
    sys_ = std::move(S);
}

std::vector<double> Run::nodal(const std::vector<double>& x) const {
    std::vector<double> t(lay_.nnodes);
    for (std::size_t g = 0; g < lay_.nnodes; ++g) t[g] = x[lay_.dof[g]];
    return t;
}

std::vector<double> Run::face_slopes(const std::vector<double>& t0) const {
    std::vector<double> g(lay_.nnodes, 0.0);
    parallel_for((lay_.nnodes + 4095) / 4096, [&](std::size_t chunk) {
        for (std::size_t v = chunk * 4096; v < std::min(lay_.nnodes, (chunk + 1) * 4096); ++v)
            if (lay_.face_area[v] > 0) g[v] = surface_slope(face_spec(lay_.face[v]), face_r(lay_.face[v]), t0[v], ta_);
    });
    return g;
}

void Run::face_terms(const std::vector<double>& t0, const std::vector<double>& g, std::vector<double>& diag,
                     std::vector<double>& load) const {
    // q(T) ~ q(t0) + g (T - t0) per node: g a on the diagonal, a (g t0 - q(t0)) in the load.
    std::vector<double> q(lay_.nnodes, 0.0);
    parallel_for((lay_.nnodes + 4095) / 4096, [&](std::size_t chunk) {
        for (std::size_t v = chunk * 4096; v < std::min(lay_.nnodes, (chunk + 1) * 4096); ++v)
            if (lay_.face_area[v] > 0) q[v] = surface_flux(face_spec(lay_.face[v]), face_r(lay_.face[v]), t0[v], ta_);
    });
    for (std::size_t v = 0; v < lay_.nnodes; ++v) {
        const double a = lay_.face_area[v];
        if (a <= 0) continue;
        diag[lay_.dof[v]] += a * g[v];
        load[lay_.dof[v]] += a * (g[v] * t0[v] - q[v]);
    }
}

std::vector<double> Run::free_values(const std::vector<double>& diag) const {
    std::vector<double> val(sys_.col.size());
    for (std::size_t i = 0; i < val.size(); ++i) val[i] = sys_.k.val[sys_.pos[i]];
    for (std::size_t r = 0; r < sys_.free_dofs.size(); ++r) val[sys_.diag[r]] += diag[sys_.free_dofs[r]];
    return val;
}

std::vector<double> Run::fixed_heat(const std::vector<double>& x, const std::vector<double>& diag,
                                    const std::vector<double>& load) const {
    std::vector<double> q(lay_.ndof, 0.0);
    for (std::size_t d = 0; d < lay_.ndof; ++d) {
        if (std::isnan(lay_.fixed[d])) continue;
        double kt = diag[d] * x[d];
        for (auto k = sys_.k.ptr[d]; k < sys_.k.ptr[d + 1]; ++k) kt += sys_.k.val[static_cast<std::size_t>(k)] * x[static_cast<std::size_t>(sys_.k.col[k])];
        q[d] = kt - load[d];
    }
    return q;
}

std::vector<double> Run::source_load(const std::vector<double>& power) const {
    std::vector<double> f(lay_.ndof, 0.0);
    for (std::size_t a = 0; a < s_.sources.size(); ++a)
        for (const auto& [d, w] : lay_.source_weight[a]) f[d] += power[a] * w;
    return f;
}

std::vector<double> Run::joule_load(const DcResult& dc, double& placed, std::size_t& lost) const {
    // Joule heat of the DC solution, deposited at the triangle centroids; barrel heat half at
    // each end.
    std::vector<double> f(lay_.ndof, 0.0);
    placed = 0;
    lost = 0;
    const Pieces& pieces = *lay_.pieces;
    auto put = [&](int layer, double x, double y, double p, double slack) {
        for (const auto& pc : pieces) {
            if (pc->layer != layer) continue;
            std::array<double, 3> w;
            const long k = pc->locator->locate(x, y, w, slack);
            if (k < 0) continue;
            const auto& tt = pc->mesh.triangles[static_cast<std::size_t>(k)];
            for (int c = 0; c < 3; ++c) f[lay_.dof[pc->offset + tt[c]]] += p * w[c];
            placed += p;
            return;
        }
        ++lost;
    };
    for (const auto& lf : dc.layers)
        for (std::size_t t = 0; t < lf.triangles.size(); ++t) {
            const auto& tr = lf.triangles[t];
            put(lf.layer, (lf.x[tr[0]] + lf.x[tr[1]] + lf.x[tr[2]]) / 3, (lf.y[tr[0]] + lf.y[tr[1]] + lf.y[tr[2]]) / 3, lf.power_w[t],
                4.0 * s_.mesh.boundary_tolerance_nm + 1000.0);
        }
    for (const auto& br : dc.barrels) {
        const double p = br.current_a * br.current_a * br.resistance_ohm;
        const auto& dr = b_.drills[br.drill];
        for (int l : {br.layer_a, br.layer_b}) put(l, double(dr.at.x), double(dr.at.y), 0.5 * p, 1.0e5);
    }
    return f;
}

std::vector<double> Run::initial_state() const {
    std::vector<double> x(lay_.ndof, ta_);
    for (std::size_t d = 0; d < lay_.ndof; ++d)
        if (!std::isnan(lay_.fixed[d])) x[d] = lay_.fixed[d];
    return x;
}

// The vertical coupling lumps the energy G (u - u')^2 of a triangle at its nodes. For a
// linear difference w the lumped value exceeds the exact integral by A/12 times the sum over
// the edges of (w_i - w_j)^2; that excess, per side at the half conductance per area with
// which each side is assembled, estimates the error of the term. It is large where thin
// dielectric slices, stiffly coupled and nearly without in-plane conduction, change
// steeply under coarse elements, which the gradient recovery of the sheets does not see.
void Run::add_vertical_estimate(const Piece& pc, const std::vector<double>& x, const std::vector<char>& skip, std::vector<double>& eta2,
                                double& total) const {
    const TriMesh& m = pc.mesh;
    for (std::size_t t = 0; t < m.triangles.size(); ++t) {
        if (skip[t]) continue;
        const auto& tr = m.triangles[t];
        const double area = triangle_area(m, t) * 1e-18;
        for (std::size_t side = 0; side < 2; ++side) {
            std::array<double, 3> w{};
            double gz = INFINITY;
            for (std::size_t i = 0; i < 3; ++i) {
                const auto& cp = sys_.counterpart[pc.offset + tr[i]][side];
                gz = std::min(gz, double(cp.gz));
                w[i] = x[lay_.dof[pc.offset + tr[i]]];
                for (std::size_t j = 0; j < 3; ++j) w[i] -= double(cp.w[j]) * x[cp.dof[j]];
            }
            if (!(gz > 0)) continue;
            const double e = 0.5 * gz * area / 12.0 *
                             ((w[0] - w[1]) * (w[0] - w[1]) + (w[1] - w[2]) * (w[1] - w[2]) + (w[2] - w[0]) * (w[2] - w[0]));
            eta2[t] += e;
            total += e;
        }
    }
}

// Steady solution: adaptive passes (when `adapt`) and the iteration of the electrothermal
// coupling and of the faces' nonlinearity. x is the initial guess on the current layout
// and receives the solution.
void Run::steady(const std::vector<double>& power, std::vector<double>& x, bool adapt) {
    std::optional<Field> previous;  // temperature of the previous coupling iteration
    std::optional<Field> guess;     // temperature of the previous adaptive pass
    if (!adapt) previous = field(x);
    bool mesh_final = !adapt;
    std::vector<double> jl(lay_.ndof, 0.0);
    // Anderson acceleration of the fixed point T -> G(T) on the final mesh: the next state
    // combines the last few so that their residuals G(T) - T cancel as far as they can.
    std::vector<double> x_in;                       // the state the last iteration started from
    std::vector<std::vector<double>> hist_g, hist_f;  // G(T) and G(T) - T, the latest last
    const std::size_t depth = 3;
    std::vector<double> changes;
    for (int iteration = 0;; ++iteration) {
        // DC solution with the conductivity at the current temperature.
        if (coupled_) {
            const auto td = Clock::now();
            DcSettings dcs = *s_.dc;
            if (!previous) {
                dcs.temperature_c = ta_;
            } else {
                const Field* f = &*previous;
                dcs.temperature_field = [f](int layer, double xx, double yy) { return f->at(layer, xx, yy); };
                if (!dc_seeds_.empty()) {
                    dcs.seeds = dc_seeds_;
                    dcs.max_passes = 1;
                }
            }
            res_.dc = solve_dc(b_, dcs);
            if (dc_seeds_.empty()) dc_seeds_ = res_.dc->mesh_seeds;
            res_.seconds_dc += since(td);
        }

        for (int pass = 0;; ++pass) {
            if (!mesh_final) {
                const auto tm = Clock::now();
                bool unrefined = false;
                auto pieces = mesh_pieces(unrefined);
                if (unrefined && pass == 0 && iteration == 0) res_.warnings.push_back("some mesh refinements could not be performed");
                res_.seconds_mesh += since(tm);
                const auto ta0 = Clock::now();
                make_layout(std::move(pieces));
                make_system();
                res_.seconds_assemble += since(ta0);
                // Initial guess (and the state the faces are linearised about): the previous
                // pass, or the previous iteration, or ambient.
                x = initial_state();
                const Field* from = guess ? &*guess : previous ? &*previous : nullptr;
                if (from) {
                    for (const auto& pc : *lay_.pieces)
                        for (std::size_t v = 0; v < pc->mesh.x.size(); ++v) {
                            const std::size_t d = lay_.dof[pc->offset + v];
                            if (std::isnan(lay_.fixed[d])) x[d] = from->at_sheet(pc->sheet, pc->mesh.x[v], pc->mesh.y[v]);
                        }
                }
            }
            const auto ta0 = Clock::now();
            const std::vector<double> t0 = nodal(x);
            std::vector<double> diag(lay_.ndof, 0.0), load = source_load(power);
            for (std::size_t d = 0; d < lay_.ndof; ++d) load[d] += sys_.rhs_static[d];
            face_terms(t0, face_slopes(t0), diag, load);
            double placed = 0;
            if (coupled_) {
                std::size_t lost = 0;
                jl = joule_load(*res_.dc, placed, lost);
                if (lost && iteration == 0 && pass == 0) res_.warnings.push_back(std::to_string(lost) + " Joule heat elements outside the thermal mesh");
                for (std::size_t d = 0; d < lay_.ndof; ++d) load[d] += jl[d];
            }
            res_.heat_joule_w = placed;
            const std::size_t nf = sys_.free_dofs.size();
            std::vector<std::vector<double>> rhs(1, std::vector<double>(nf));
            std::vector<std::vector<double>> sol(1, std::vector<double>(nf));
            for (std::size_t r = 0; r < nf; ++r) {
                rhs[0][r] = load[sys_.free_dofs[r]] - sys_.b_fixed[r];
                sol[0][r] = x[sys_.free_dofs[r]];
            }
            const auto val = free_values(diag);
            res_.unknowns = nf;
            res_.triangles = lay_.ntri;
            res_.seconds_assemble += since(ta0);
            const auto ts = Clock::now();
            const SpdSolver solver(nf, sys_.ptr, sys_.col, val, s_.tolerance, s_.max_iterations);
            const SolveStats st = solver.solve(rhs[0], sol[0]);
            res_.iterations = st.iterations;
            res_.residual = st.residual;
            if (nf > 0 && !(st.residual <= s_.tolerance * 10))
                res_.warnings.push_back("linear solver did not reach the tolerance: residual " + std::to_string(st.residual));
            res_.seconds_solve += since(ts);
            for (std::size_t r = 0; r < nf; ++r) x[sys_.free_dofs[r]] = sol[0][r];

            if (mesh_final) break;
            // Error estimate in the energy norm of the whole problem (T - ambient).
            const Pieces& pieces = *lay_.pieces;
            const std::vector<double> T = nodal(x);
            std::vector<std::vector<double>> eta2(pieces.size());
            double eta2_total = 0, energy_plane = 0, area_total = 0;
            for (std::size_t k = 0; k < pieces.size(); ++k) {
                const Piece& pc = *pieces[k];
                std::vector<double> u(pc.mesh.x.size());
                for (std::size_t v = 0; v < u.size(); ++v) u[v] = T[pc.offset + v];
                std::vector<char> skip(pc.mesh.triangles.size(), 0);
                for (std::size_t t = 0; t < skip.size(); ++t) {
                    const auto& tr = pc.mesh.triangles[t];
                    const auto d0 = lay_.dof[pc.offset + tr[0]];
                    skip[t] = d0 == lay_.dof[pc.offset + tr[1]] && d0 == lay_.dof[pc.offset + tr[2]];
                }
                eta2[k] = zz_estimate(pc.mesh, u, pc.conductance, skip, eta2_total, energy_plane, area_total, &pc.copper);
                add_vertical_estimate(pc, x, skip, eta2[k], eta2_total);
            }
            double energy = 0;
            for (std::size_t d = 0; d < lay_.ndof; ++d) {
                double kt = diag[d] * (x[d] - ta_);
                for (auto k = sys_.k.ptr[d]; k < sys_.k.ptr[d + 1]; ++k)
                    kt += sys_.k.val[static_cast<std::size_t>(k)] * (x[static_cast<std::size_t>(sys_.k.col[k])] - ta_);
                energy += (x[d] - ta_) * kt;
            }
            energy = std::max(energy, energy_plane);
            res_.estimated_error = energy > 0 ? std::sqrt(eta2_total / energy) : 0.0;
            res_.passes.push_back({lay_.ntri, res_.estimated_error, energy});
            bool stop = res_.estimated_error <= s_.target_error || pass + 1 >= s_.max_passes;
            if (!stop) {
                std::vector<RefinePart> rp;
                for (std::size_t k = 0; k < pieces.size(); ++k) rp.push_back({&pieces[k]->mesh, &eta2[k], owner(pieces[k]->sheet, pieces[k]->part)});
                bool budget = false;
                const std::size_t marked =
                    mark_and_seed(rp, eta2_total, energy, s_.target_error, lay_.ntri, s_.max_triangles, s_.mesh.min_edge_nm, seeds_, budget);
                if (budget) res_.warnings.push_back("refinement limited by the triangle budget (" + std::to_string(s_.max_triangles) + ")");
                stop = marked == 0;
            }
            if (stop) {
                if (res_.estimated_error > s_.target_error)
                    res_.warnings.push_back("estimated error " + std::to_string(res_.estimated_error) + " above the target after " +
                                            std::to_string(pass + 1) + " passes");
                mesh_final = true;
                break;
            }
            guess = field(x);
        }

        // Change since the previous iteration; stop when settled.
        double change = INFINITY;
        if (previous) {
            change = 0;
            for (const auto& pc : *lay_.pieces)
                for (std::size_t v = 0; v < pc->mesh.x.size(); ++v)
                    change = std::max(change, std::abs(x[lay_.dof[pc->offset + v]] - previous->at_sheet(pc->sheet, pc->mesh.x[v], pc->mesh.y[v])));
        }
        res_.coupling_iterations = iteration + 1;
        res_.coupling_change_k = change;
        // Below 20 C - 1/alpha the copper's resistivity would be negative: a state the
        // coupling can settle in only when the heating has no physical steady state.
        if (coupled_) {
            const double floor_c = 20.0 - 1.0 / s_.dc->temperature_coefficient;
            const double tmin = *std::min_element(x.begin(), x.end());
            if (tmin <= floor_c) {
                res_.warnings.push_back("thermal runaway: the coupled iteration went below " + std::to_string(floor_c) +
                                        " C, where the copper's resistivity vanishes; the heating has no steady state at this current");
                break;
            }
        }
        if (!(coupled_ || radiation_) || change <= s_.coupling_tolerance_k) break;
        changes.push_back(change);
        const std::size_t nc = changes.size();
        if (nc >= 4 && changes[nc - 1] > changes[nc - 2] && changes[nc - 2] > changes[nc - 3] && changes[nc - 3] > changes[nc - 4]) {
            res_.warnings.push_back("thermal runaway: the temperature grows with each iteration (change " + std::to_string(change) +
                                    " K); the heating has no steady state at this current");
            break;
        }
        if (iteration + 1 >= s_.max_coupling_iterations) {
            res_.warnings.push_back("temperature did not settle: change " + std::to_string(change) + " K after " +
                                    std::to_string(iteration + 1) + " iterations");
            break;
        }
        if (s_.coupling_acceleration && previous && x_in.size() == x.size()) {
            const std::size_t n = x.size();
            std::vector<double> f(n);
            for (std::size_t d = 0; d < n; ++d) f[d] = x[d] - x_in[d];
            hist_g.push_back(x);
            hist_f.push_back(f);
            if (hist_g.size() > depth + 1) {
                hist_g.erase(hist_g.begin());
                hist_f.erase(hist_f.begin());
            }
            const std::size_t m = hist_g.size() - 1;
            if (m > 0) {
                // min |f_k - dF gamma| over the differences of consecutive residuals.
                std::vector<std::vector<double>> dF(m, std::vector<double>(n)), dG(m, std::vector<double>(n));
                for (std::size_t j = 0; j < m; ++j)
                    for (std::size_t d = 0; d < n; ++d) {
                        dF[j][d] = hist_f[j + 1][d] - hist_f[j][d];
                        dG[j][d] = hist_g[j + 1][d] - hist_g[j][d];
                    }
                std::vector<double> A(m * m, 0.0), rhs(m, 0.0);
                for (std::size_t a = 0; a < m; ++a) {
                    for (std::size_t c = 0; c < m; ++c)
                        for (std::size_t d = 0; d < n; ++d) A[a * m + c] += dF[a][d] * dF[c][d];
                    for (std::size_t d = 0; d < n; ++d) rhs[a] += dF[a][d] * f[d];
                }
                double tr = 0;
                for (std::size_t a = 0; a < m; ++a) tr += A[a * m + a];
                for (std::size_t a = 0; a < m; ++a) A[a * m + a] += 1e-10 * tr + 1e-300;
                for (std::size_t c = 0; c < m; ++c)
                    for (std::size_t q = c + 1; q < m; ++q) {
                        const double r = A[q * m + c] / A[c * m + c];
                        for (std::size_t k = c; k < m; ++k) A[q * m + k] -= r * A[c * m + k];
                        rhs[q] -= r * rhs[c];
                    }
                std::vector<double> gamma(m);
                for (std::size_t c = m; c-- > 0;) {
                    double v = rhs[c];
                    for (std::size_t k = c + 1; k < m; ++k) v -= A[c * m + k] * gamma[k];
                    gamma[c] = v / A[c * m + c];
                }
                std::vector<double> next = x;
                bool finite = true;
                for (std::size_t d = 0; d < n; ++d) {
                    for (std::size_t j = 0; j < m; ++j) next[d] -= gamma[j] * dG[j][d];
                    finite = finite && std::isfinite(next[d]);
                }
                if (finite) x = std::move(next);
            }
        }
        x_in = x;
        previous = field(x);
    }
    if (lay_.unlocated) res_.warnings.push_back(std::to_string(lay_.unlocated) + " nodes without a counterpart on the neighbouring sheet");
}

// Results of a state x: sources, heat sinks, the heat balance and, with final_fields, the
// fields per copper layer.
void Run::report(const std::vector<double>& x, const std::vector<double>& power, double joule_w, bool final_fields) {
    const Pieces& pieces = *lay_.pieces;
    const std::vector<double> T = nodal(x);
    res_.heat_sources_w = 0;
    for (std::size_t a = 0; a < s_.sources.size(); ++a)
        if (!s_.sources[a].temperature_c) res_.heat_sources_w += power[a];
    res_.heat_joule_w = joule_w;
    res_.heat_top_w = res_.heat_bottom_w = 0;
    for (std::size_t v = 0; v < lay_.nnodes; ++v) {
        if (lay_.face_area[v] <= 0) continue;
        const double q = lay_.face_area[v] * surface_flux(face_spec(lay_.face[v]), face_r(lay_.face[v]), T[v], ta_);
        (lay_.face[v] == 1 ? res_.heat_top_w : res_.heat_bottom_w) += q;
    }
    // Heat into the sinks through the interface, and into fixed areas.
    res_.heatsinks.clear();
    res_.heat_sinks_w = 0;
    for (std::size_t k = 0; k < s_.heatsinks.size(); ++k) {
        ThermalHeatsinkResult h;
        h.name = s_.heatsinks[k].name;
        h.temperature_c = x[lay_.sink_dof[k]];
        h.temperature_set = s_.heatsinks[k].temperature_c.has_value();
        for (const auto& c : lay_.contacts)
            if (c.sink == k) h.heat_w += c.g * (T[c.node] - h.temperature_c);
        res_.heat_sinks_w += h.heat_w;
        res_.heatsinks.push_back(h);
    }
    res_.heat_fixed_w = 0;
    res_.sources.clear();
    for (std::size_t a = 0; a < s_.sources.size(); ++a) {
        ThermalSourceResult r;
        r.name = s_.sources[a].area.name;
        r.temperature_set = s_.sources[a].temperature_c.has_value();
        double sum = 0, area = 0, mx = -INFINITY;
        for (const auto& [g, ar] : lay_.area_nodes[a]) {
            sum += T[g] * ar;
            area += ar;
            mx = std::max(mx, T[g]);
        }
        r.mean_c = area > 0 ? sum / area : NAN;
        r.max_c = mx;
        if (lay_.junction[a] >= 0) r.junction_c = x[static_cast<std::size_t>(lay_.junction[a])];
        r.power_w = power[a];
        res_.sources.push_back(r);
    }
    // Heat taken out by fixed areas: what the faces and sinks of the board do not take,
    // per area by its share of the conduction into it.
    {
        std::vector<double> diag(lay_.ndof, 0.0), load = source_load(power);
        for (std::size_t d = 0; d < lay_.ndof; ++d) load[d] += sys_.rhs_static[d];
        face_terms(T, face_slopes(T), diag, load);
        if (coupled_ && res_.dc) {
            double placed = 0;
            std::size_t lost = 0;
            const auto jl = joule_load(*res_.dc, placed, lost);
            for (std::size_t d = 0; d < lay_.ndof; ++d) load[d] += jl[d];
        }
        const auto q = fixed_heat(x, diag, load);
        for (std::size_t a = 0; a < s_.sources.size(); ++a)
            if (s_.sources[a].temperature_c) {
                res_.sources[a].power_w = q[static_cast<std::size_t>(lay_.area_dof[a])];
                res_.heat_fixed_w += res_.sources[a].power_w;
            }
    }
    if (!final_fields) return;
    res_.layers.clear();
    res_.max_c = -INFINITY;
    std::map<int, std::size_t> slot;
    for (const auto& pc : pieces) {
        if (pc->layer == 0) continue;
        if (!slot.count(pc->layer)) {
            slot[pc->layer] = res_.layers.size();
            ThermalLayerField f;
            f.layer = pc->layer;
            res_.layers.push_back(std::move(f));
        }
        ThermalLayerField& f = res_.layers[slot[pc->layer]];
        const auto base = static_cast<std::uint32_t>(f.x.size());
        for (std::size_t v = 0; v < pc->mesh.x.size(); ++v) {
            f.x.push_back(pc->mesh.x[v]);
            f.y.push_back(pc->mesh.y[v]);
            const double t = T[pc->offset + v];
            f.temperature_c.push_back(t);
            if (t > res_.max_c) {
                res_.max_c = t;
                res_.max_layer = pc->layer;
                res_.max_x_nm = pc->mesh.x[v];
                res_.max_y_nm = pc->mesh.y[v];
            }
        }
        for (std::size_t t = 0; t < pc->mesh.triangles.size(); ++t) {
            const auto& tr = pc->mesh.triangles[t];
            f.triangles.push_back({base + tr[0], base + tr[1], base + tr[2]});
            f.copper.push_back(pc->copper[t]);
        }
    }
}

// Backward Euler from the state x at t = 0 on the current mesh. The faces are linearised
// with fixed slopes (those of the steady solution at the peak power), their deviation from
// the exact flux carried in the load, so that the matrix depends on the step alone.
void Run::transient(std::vector<double> x) {
    const TransientSettings& tr = *s_.transient;
    const double duration = tr.duration_s;
    const double max_step = tr.max_step_s > 0 ? tr.max_step_s : duration / 50.0;
    const double first = tr.step_s > 0 ? std::min(tr.step_s, max_step) : max_step / std::pow(2.0, 6);
    ThermalTransientResult out;
    const std::size_t nsrc = s_.sources.size(), nsink = s_.heatsinks.size();
    out.source_mean_c.assign(nsrc, {});
    out.source_max_c.assign(nsrc, {});
    out.source_junction_c.assign(nsrc, {});
    out.heatsink_c.assign(nsink, {});
    out.heatsink_w.assign(nsink, {});

    const std::vector<double> g = face_slopes(peak_t_.empty() ? nodal(x) : peak_t_);
    // Breakpoints: the points of the schedules and the end.
    std::set<double> breaks{duration};
    for (const auto& h : s_.sources)
        for (const auto& p : h.power_schedule)
            if (p[0] > 0 && p[0] < duration) breaks.insert(p[0]);
    std::vector<double> snaps = tr.snapshot_times_s;
    if (snaps.empty()) {
        const int n = 20;
        for (int k = 1; k <= n; ++k) snaps.push_back(first * std::pow(duration / first, double(k) / n));
    }
    std::sort(snaps.begin(), snaps.end());
    std::size_t next_snap = 0;
    // Copper layers in the order of the result fields, for the snapshots.
    std::vector<std::vector<std::size_t>> layer_nodes;
    {
        std::map<int, std::size_t> slot;
        for (const auto& pc : *lay_.pieces) {
            if (pc->layer == 0) continue;
            if (!slot.count(pc->layer)) {
                slot[pc->layer] = layer_nodes.size();
                layer_nodes.emplace_back();
            }
            auto& ln = layer_nodes[slot[pc->layer]];
            for (std::size_t v = 0; v < pc->mesh.x.size(); ++v) ln.push_back(pc->offset + v);
        }
    }
    auto snapshot = [&](double t, const std::vector<double>& T) {
        out.snapshot_time_s.push_back(t);
        std::vector<std::vector<double>> s(layer_nodes.size());
        for (std::size_t l = 0; l < layer_nodes.size(); ++l) {
            s[l].reserve(layer_nodes[l].size());
            for (auto g2 : layer_nodes[l]) s[l].push_back(T[g2]);
        }
        out.snapshot_c.push_back(std::move(s));
    };

    std::vector<double> jl(lay_.ndof, 0.0);
    double joule_w = 0;
    std::vector<double> t_dc;  // nodal temperature of the last DC solution
    auto update_dc = [&](const std::vector<double>& xs) {
        const auto td = Clock::now();
        DcSettings dcs = *s_.dc;
        const Field f = field(xs);
        dcs.temperature_field = [&f](int layer, double xx, double yy) { return f.at(layer, xx, yy); };
        if (!dc_seeds_.empty()) {
            dcs.seeds = dc_seeds_;
            dcs.max_passes = 1;
        }
        res_.dc = solve_dc(b_, dcs);
        std::size_t lost = 0;
        jl = joule_load(*res_.dc, joule_w, lost);
        t_dc = nodal(xs);
        ++out.dc_solutions;
        res_.seconds_dc += since(td);
    };
    if (coupled_) update_dc(x);

    auto record = [&](double t, const std::vector<double>& xs, double in_w, double out_w, double stored_w) {
        const std::vector<double> T = nodal(xs);
        out.time_s.push_back(t);
        double mx = -INFINITY;
        for (const auto& ln : layer_nodes)
            for (auto g2 : ln) mx = std::max(mx, T[g2]);
        out.max_c.push_back(mx);
        for (std::size_t a = 0; a < nsrc; ++a) {
            double sum = 0, area = 0, m2 = -INFINITY;
            for (const auto& [g2, ar] : lay_.area_nodes[a]) {
                sum += T[g2] * ar;
                area += ar;
                m2 = std::max(m2, T[g2]);
            }
            out.source_mean_c[a].push_back(sum / area);
            out.source_max_c[a].push_back(m2);
            out.source_junction_c[a].push_back(lay_.junction[a] >= 0 ? xs[static_cast<std::size_t>(lay_.junction[a])] : NAN);
        }
        for (std::size_t k = 0; k < nsink; ++k) {
            const double ts = xs[lay_.sink_dof[k]];
            double q = 0;
            for (const auto& c : lay_.contacts)
                if (c.sink == k) q += c.g * (T[c.node] - ts);
            out.heatsink_c[k].push_back(ts);
            out.heatsink_w[k].push_back(q);
        }
        out.heat_in_w.push_back(in_w);
        out.heat_out_w.push_back(out_w);
        out.heat_stored_w.push_back(stored_w);
    };
    std::vector<double> power(nsrc);
    for (std::size_t a = 0; a < nsrc; ++a) power[a] = power_at(s_.sources[a], 0.0);
    {
        double in0 = joule_w;
        for (std::size_t a = 0; a < nsrc; ++a)
            if (!s_.sources[a].temperature_c) in0 += power[a];
        record(0.0, x, in0, NAN, NAN);
        snapshot(0.0, nodal(x));
    }

    // One solver per step length.
    struct StepSolver {
        std::unique_ptr<SpdSolver> solver;
        std::vector<double> val;  // the matrix, for the residual of the previous state
    };
    std::map<double, StepSolver> solvers;
    const bool increments = tr.step_tolerance > 0;
    std::vector<double> diag_faces(lay_.ndof, 0.0);
    for (std::size_t v = 0; v < lay_.nnodes; ++v) diag_faces[lay_.dof[v]] += lay_.face_area[v] * g[v];
    auto solver_for = [&](double dt) -> const StepSolver& {
        auto it = solvers.find(dt);
        if (it != solvers.end()) return it->second;
        if (solvers.size() >= 3) solvers.erase(solvers.begin());
        std::vector<double> d(diag_faces);
        for (std::size_t k = 0; k < lay_.ndof; ++k) d[k] += lay_.capacity[k] / dt;
        StepSolver ss;
        ss.val = free_values(d);
        ++out.matrix_setups;
        const auto ts = Clock::now();
        ss.solver = std::make_unique<SpdSolver>(sys_.free_dofs.size(), sys_.ptr, sys_.col, ss.val, increments ? tr.step_tolerance : s_.tolerance,
                                                s_.max_iterations);
        res_.seconds_solve += since(ts);
        out.seconds_setup += since(ts);
        return solvers.emplace(dt, std::move(ss)).first->second;
    };

    double t = 0, dt_base = first;
    int since_reset = 0;
    const std::size_t nf = sys_.free_dofs.size();
    std::vector<double> rhs(nf), sol(nf);
    while (t < duration * (1 - 1e-12)) {
        const double next_break = *breaks.upper_bound(t * (1 + 1e-12) + 1e-15);
        double dt = std::min(dt_base, next_break - t);
        if (next_break - (t + dt) < 0.5 * dt) dt = next_break - t;  // no sliver before a breakpoint
        const double t1 = t + dt;
        for (std::size_t a = 0; a < nsrc; ++a) power[a] = power_at(s_.sources[a], t1);
        if (coupled_) {
            double change = 0;
            const auto T = nodal(x);
            for (std::size_t v = 0; v < T.size(); ++v) change = std::max(change, std::abs(T[v] - t_dc[v]));
            if (change > tr.dc_update_k) update_dc(x);
        }
        const auto ta0 = Clock::now();
        const std::vector<double> t0 = nodal(x);
        std::vector<double> diag(lay_.ndof, 0.0), load = source_load(power);
        for (std::size_t d = 0; d < lay_.ndof; ++d) load[d] += sys_.rhs_static[d] + jl[d] + lay_.capacity[d] / dt * x[d];
        face_terms(t0, g, diag, load);
        for (std::size_t r = 0; r < nf; ++r) {
            rhs[r] = load[sys_.free_dofs[r]] - sys_.b_fixed[r];
            sol[r] = x[sys_.free_dofs[r]];
        }
        res_.seconds_assemble += since(ta0);
        const StepSolver& sv = solver_for(dt);
        const auto ts = Clock::now();
        SolveStats st;
        if (increments) {
            // The change from the previous state: A d = rhs - A x_prev, to a tolerance relative
            // to that residual, so that the accuracy follows the step's change.
            std::vector<double> r0(rhs), dlt(nf, 0.0);
            for (std::size_t r = 0; r < nf; ++r)
                for (auto k = sys_.ptr[r]; k < sys_.ptr[r + 1]; ++k)
                    r0[r] -= sv.val[static_cast<std::size_t>(k)] * sol[static_cast<std::size_t>(sys_.col[static_cast<std::size_t>(k)])];
            double norm = 0;
            for (double v : r0) norm += v * v;
            if (norm > 0) st = sv.solver->solve(r0, dlt);
            for (std::size_t r = 0; r < nf; ++r) sol[r] += dlt[r];
        } else {
            st = sv.solver->solve(rhs, sol);
        }
        res_.seconds_solve += since(ts);
        out.solver_iterations += st.iterations;
        if (nf > 0 && !(st.residual <= 10 * (increments ? tr.step_tolerance : s_.tolerance)))
            res_.warnings.push_back("linear solver did not reach the tolerance at t = " + std::to_string(t1) + " s");
        std::vector<double> x1(x);
        for (std::size_t r = 0; r < nf; ++r) x1[sys_.free_dofs[r]] = sol[r];

        // Balance of the step: heat in, heat out (faces as linearised, sinks to ambient,
        // fixed temperatures) and the rate of change of the stored heat.
        double in_w = joule_w, stored = 0, out_w = 0;
        for (std::size_t a = 0; a < nsrc; ++a)
            if (!s_.sources[a].temperature_c) in_w += power[a];
        for (std::size_t d = 0; d < lay_.ndof; ++d) stored += lay_.capacity[d] * (x1[d] - x[d]) / dt;
        {
            const auto T1 = nodal(x1);
            for (std::size_t v = 0; v < lay_.nnodes; ++v)
                if (lay_.face_area[v] > 0)
                    out_w += lay_.face_area[v] * (surface_flux(face_spec(lay_.face[v]), face_r(lay_.face[v]), t0[v], ta_) + g[v] * (T1[v] - t0[v]));
            for (std::size_t k = 0; k < nsink; ++k) {
                const auto& hs = s_.heatsinks[k];
                if (hs.temperature_c) {
                    for (const auto& c : lay_.contacts)
                        if (c.sink == k) out_w += c.g * (T1[c.node] - x1[lay_.sink_dof[k]]);
                } else {
                    out_w += (x1[lay_.sink_dof[k]] - ta_) / *hs.resistance_k_per_w;
                }
            }
            // Fixed areas: the conduction into them, less the Joule heat deposited there.
            for (std::size_t a = 0; a < nsrc; ++a) {
                if (!s_.sources[a].temperature_c) continue;
                const auto dfix = static_cast<std::size_t>(lay_.area_dof[a]);
                double kt = 0;
                for (auto k = sys_.k.ptr[dfix]; k < sys_.k.ptr[dfix + 1]; ++k)
                    kt += sys_.k.val[static_cast<std::size_t>(k)] * x1[static_cast<std::size_t>(sys_.k.col[k])];
                // Faces of the fixed area itself are counted with the faces above.
                out_w += kt - jl[dfix];
            }
        }
        x = std::move(x1);
        t = t1;
        record(t, x, in_w, out_w, stored);
        ++out.steps;
        while (next_snap < snaps.size() && snaps[next_snap] <= t * (1 + 1e-9)) {
            if (out.snapshot_time_s.empty() || out.snapshot_time_s.back() != t) snapshot(t, nodal(x));
            ++next_snap;
        }
        // Step length: back to the first step at a breakpoint of a schedule, doubled every
        // steps_per_doubling steps otherwise.
        if (std::abs(t - next_break) <= 1e-12 * std::max(1.0, t) && next_break < duration) {
            dt_base = first;
            since_reset = 0;
        } else if (++since_reset % std::max(1, tr.steps_per_doubling) == 0) {
            dt_base = std::min(2.0 * dt_base, max_step);
        }
    }
    if (out.snapshot_time_s.back() != t) snapshot(t, nodal(x));
    report(x, power, joule_w, true);
    res_.transient = std::move(out);
}

ThermalResult Run::run() {
    tol_ = b_.source.arc_tolerance_nm > 0 ? b_.source.arc_tolerance_nm : 1000.0;
    ta_ = s_.ambient_c;
    for (const auto& hs : s_.sources) {
        const bool powered = hs.power_w.has_value() || !hs.power_schedule.empty();
        if (powered == hs.temperature_c.has_value())
            throw std::runtime_error(hs.area.name + ": give either a power (power_w or power_schedule) or temperature_c");
        for (std::size_t k = 1; k < hs.power_schedule.size(); ++k)
            if (hs.power_schedule[k][0] < hs.power_schedule[k - 1][0])
                throw std::runtime_error(hs.area.name + ": the times of power_schedule must not decrease");
    }
    std::set<std::string> sink_names;
    for (const auto& h : s_.heatsinks) {
        if (!sink_names.insert(h.name).second) throw std::runtime_error("heat sink " + h.name + " is defined twice");
        if (h.temperature_c.has_value() == h.resistance_k_per_w.has_value())
            throw std::runtime_error("heat sink " + h.name + ": give either temperature_c or resistance_k_per_w");
        if (h.resistance_k_per_w && !(*h.resistance_k_per_w > 0))
            throw std::runtime_error("heat sink " + h.name + ": resistance_k_per_w must be positive");
        if (h.regions.empty()) throw std::runtime_error("heat sink " + h.name + " has no region");
        if (h.tim_m2k_per_w < 0) throw std::runtime_error("heat sink " + h.name + ": negative interface resistance");
    }
    bool any_fixed = !s_.heatsinks.empty();
    for (const auto& hs : s_.sources) any_fixed = any_fixed || hs.temperature_c.has_value();
    const bool any_surface = s_.top.convection_w_m2k > 0 || s_.top.emissivity > 0 || s_.bottom.convection_w_m2k > 0 ||
                             s_.bottom.emissivity > 0;
    if (!any_fixed && !any_surface) throw std::runtime_error("no path for the heat: no surface cooling, heat sink or fixed temperature");
    if (s_.transient) {
        const auto& tr = *s_.transient;
        if (!(tr.duration_s > 0)) throw std::runtime_error("transient: duration_s must be positive");
        if (tr.step_s < 0 || tr.max_step_s < 0) throw std::runtime_error("transient: step lengths must not be negative");
    }
    radiation_ = s_.top.emissivity > 0 || s_.bottom.emissivity > 0;
    coupled_ = s_.dc && !s_.dc->ports.empty();

    auto t0 = Clock::now();
    m_ = std::make_unique<Model>(build_model(b_, s_, tol_));
    const int n = m_->n;
    ns_ = static_cast<int>(m_->sheets.size()) - 1;
    // Mesh constraints per sheet: copper outlines, hole walls and heat areas on the copper
    // sheets, heat sink regions on the outer sheets.
    constraints_.assign(ns_ + 1, {});
    for (int l = 1; l <= n; ++l) {
        PolygonSet& c = constraints_[m_->sheet_of[l]];
        c = m_->copper[l];
        for (const auto& h : m_->holes)
            if (std::find(h.layers.begin(), h.layers.end(), l) != h.layers.end()) c.insert(c.end(), h.outline.begin(), h.outline.end());
        for (const auto& a : m_->areas)
            if (auto it = a.find(l); it != a.end()) c.insert(c.end(), it->second.begin(), it->second.end());
    }
    for (const auto& r : m_->sink_regions) {
        PolygonSet& c = constraints_[r.bottom ? ns_ : 1];
        c.insert(c.end(), r.area.begin(), r.area.end());
    }
    nparts_ = m_->domain.size();
    seeds_.assign(static_cast<std::size_t>(ns_ + 1) * nparts_, {});
    res_.seconds_mesh += since(t0);

    std::vector<double> power(s_.sources.size());
    std::vector<double> x;
    if (!s_.transient) {
        for (std::size_t a = 0; a < power.size(); ++a) power[a] = steady_power(s_.sources[a]);
        steady(power, x, true);
        report(x, power, res_.heat_joule_w, true);
        return std::move(res_);
    }
    // Transient: the mesh from the steady solution at the peak powers, whose temperature
    // also gives the slopes of the faces' linearisation; then the start.
    const auto tsteady = Clock::now();
    for (std::size_t a = 0; a < power.size(); ++a) power[a] = peak_power(s_.sources[a]);
    steady(power, x, true);
    std::vector<double> start;
    for (std::size_t a = 0; a < power.size(); ++a) power[a] = power_at(s_.sources[a], 0.0);
    if (s_.transient->start_steady) {
        start = x;
        steady(power, start, false);
    } else {
        start = initial_state();
    }
    // The faces are linearised about the peak state; the transient keeps those slopes.
    peak_t_ = nodal(x);
    const double steady_s = since(tsteady);
    transient(std::move(start));
    res_.transient->seconds_steady = steady_s;
    return std::move(res_);
}

}  // namespace

ThermalResult solve_thermal(const Board& b, const ThermalSettings& s) {
    Run r(b, s);
    return r.run();
}

}  // namespace kairoq
