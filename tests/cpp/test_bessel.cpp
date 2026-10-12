#include "bessel.hpp"

#include <doctest/doctest.h>

#include <cmath>

using namespace kairoq::detail;

// Reference values from the Cephes routines (scipy.special j0, k0, i0), each branch of the
// implementation and its ends.
TEST_CASE("J0 against reference values") {
    const double ref[][2] = {{0, 1}, {0.5, 0.93846980724081297}, {2.4048255576957729, -9.5868825549168066e-17}, {5, -0.17759677131433829},
                             {7.99, 0.17399001312793252}, {8, 0.1716508071375539}, {12, 0.047689310796833549}, {18.3, 0.04233583847140502},
                             {24.99, 0.095008236967548135}, {25, 0.096266783275958126}, {40, 0.0073668905842372923},
                             {173.2, -0.056442243340822935}};
    for (const auto& r : ref) {
        CHECK(std::abs(bessel_j0(r[0]) - r[1]) < 2e-14);
        CHECK(bessel_j0(-r[0]) == bessel_j0(r[0]));
    }
}

TEST_CASE("K0 and I0 against reference values") {
    const double ref[][3] = {{0.001, 7.0236888005623825, 1.0000002500000158}, {0.05, 3.1142340294719917, 1.0006250976630324},
                             {0.7, 0.66051985991510165, 1.1263030183068092},  {2, 0.1138938727495334, 2.2795853023360668},
                             {7.5, 0.00024917761635611437, 268.16131151518931}, {15, 9.8195364823964346e-08, 339649.37329791381},
                             {20, 5.7412378153365238e-10, 43558282.559553534}};
    for (const auto& r : ref) {
        CHECK(std::abs(bessel_k0(r[0]) / r[1] - 1) < 1e-14);
        CHECK(std::abs(bessel_i0(r[0]) / r[2] - 1) < 1e-14);
    }
}
