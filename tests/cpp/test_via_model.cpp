#include "via_model.hpp"

#include <doctest/doctest.h>

#include <cmath>

using namespace kairoq;

TEST_CASE("Hankel functions of the second kind against reference values") {
    // scipy.special.hankel2
    struct Ref {
        cplx z, h0, h1;
    };
    const Ref refs[] = {
        {{1e-4, 0}, {0.9999999975000001, 5.937289069709336}, {5.0000000327317196e-05, 6366.198036455762}},
        {{0.05, -0.001}, {0.986587188421374, 1.979209109565159}, {-0.22872379741865195, 12.784267921907036}},
        {{1, 0}, {0.7651976865579664, -0.088256964215677}, {0.44005058574493355, 0.7812128213002889}},
        {{3.7, -0.02}, {-0.3909789614984973, -0.10499640944587393}, {0.053952532661849746, -0.40847522411297904}},
        {{9.5, -0.1}, {-0.17463006068978656, -0.15581092817985454}, {0.14694825013458906, -0.18313718976725105}},
        {{10.5, -0.1}, {-0.2143871741023453, 0.060078095675383816}, {-0.07035331995536519, -0.21187044144875813}},
        {{25, -0.3}, {0.07074283845839906, 0.09468558472957102}, {-0.09331340832256146, 0.0726661355113718}},
        {{80, 0}, {-0.06974216551221005, 0.055620339089770016}, {-0.05605729667571261, -0.06939591378458808}},
    };
    for (const auto& r : refs) {
        cplx h0, h1;
        hankel2(r.z, h0, h1);
        CHECK(std::abs(h0 - r.h0) < 1e-10 * std::abs(r.h0));
        CHECK(std::abs(h1 - r.h1) < 1e-10 * std::abs(r.h1));
    }
}

TEST_CASE("Barrel ladder through a thick plane: the coaxial capacitance per length") {
    // A barrel of 0.3 mm through a plane with a 0.6 mm hole, between two pad-less layers
    // 0.5 mm away, in vacuum. One more millimetre of plane adds 2 pi eps0 / ln(b/a) to the
    // nodes' total capacitance.
    auto cap = [](double t_mm) {
        ViaCut c;
        c.drill_nm = 0.3e6;
        const double t = t_mm * 1e6;
        c.copper = {{t + 1.0e6, t + 1.035e6, 0, 0}, {0.5e6, 0.5e6 + t, 0, 0.6e6}, {-0.035e6, 0, 0, 0}};
        c.from = 0;
        c.to = 2;
        const ViaLadder lad = via_ladder(c, 1e9, 0.01);
        double total = 0;
        for (double v : lad.capacitance) total += v;
        return total;
    };
    const double d = cap(3.0) - cap(2.0);
    const double ref = 2 * 3.14159265358979323846 * 8.8541878128e-12 * 1e-3 / std::log(2.0);
    CHECK(d == doctest::Approx(ref).epsilon(0.01));
}
