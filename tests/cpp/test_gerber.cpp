#include "kairoq/gerber.hpp"

#include <doctest/doctest.h>

#include <cmath>
#include <string>

using namespace kairoq;

namespace {

const char* kHeader = "%FSLAX46Y46*%\n%MOMM*%\n%LPD*%\nG01*\n";

double mm2(const PolygonSet& p) { return area_nm2(p) / 1e12; }

double total_mm2(const GerberImage& img) {
    double a = 0;
    for (const auto& [net, polys] : render_by_net(img.features)) a += mm2(polys);
    return a;
}

}  // namespace

TEST_CASE("macro expressions follow arithmetic precedence") {
    const std::vector<double> v{0, 1.5, 2.0, 4.0};
    CHECK(evaluate_macro_expression("$1+$1", v) == doctest::Approx(3.0));
    CHECK(evaluate_macro_expression("$2x$3+1", v) == doctest::Approx(9.0));
    CHECK(evaluate_macro_expression("1+$2x$3", v) == doctest::Approx(9.0));
    CHECK(evaluate_macro_expression("($1+0.5)x$3", v) == doctest::Approx(8.0));
    CHECK(evaluate_macro_expression("-$2/4", v) == doctest::Approx(-0.5));
    CHECK(evaluate_macro_expression("$3-$2-1", v) == doctest::Approx(1.0));
    CHECK(evaluate_macro_expression("$9", v) == doctest::Approx(0.0));
}

TEST_CASE("attribute fields decode unicode escapes") {
    CHECK(decode_attribute_field("A\\u002CB") == "A,B");
    CHECK(decode_attribute_field("\\u00B5") == "\xC2\xB5");
    CHECK(decode_attribute_field("plain\\x") == "plain\\x");
}

TEST_CASE("a stroked track has the stadium area and carries its net") {
    const std::string g = std::string(kHeader) +
                          "%TF.FileFunction,Copper,L1,Top*%\n"
                          "%TA.AperFunction,Conductor*%\n%ADD10C,0.500000*%\n%TD*%\n"
                          "D10*\n%TO.N,VCC*%\nX0Y0D02*\nX10000000Y0D01*\n%TD*%\nM02*\n";
    const GerberImage img = parse_gerber(g);
    CHECK(img.warnings.empty());
    CHECK(img.has_net_attributes);
    CHECK(img.file_function() == std::vector<std::string>{"Copper", "L1", "Top"});
    REQUIRE(img.features.size() == 1);
    CHECK(img.features[0].net == "VCC");
    CHECK(img.features[0].function == "Conductor");
    CHECK(img.features[0].stroke_width_nm == doctest::Approx(500000));
    const auto by_net = render_by_net(img.features);
    REQUIRE(by_net.count("VCC") == 1);
    CHECK(mm2(by_net.at("VCC")) == doctest::Approx(10 * 0.5 + kPi * 0.0625).epsilon(1e-6));
}

TEST_CASE("regions with arcs and clear polarity") {
    // Disc of radius 2 mm as a region of two half arcs, then a cleared 1 mm square.
    const std::string g = std::string(kHeader) +
                          "G75*\nG36*\nX2000000Y0D02*\nG03*\nX-2000000Y0I-2000000J0D01*\nX2000000Y0I2000000J0D01*\nG01*\nG37*\n"
                          "%LPC*%\n%ADD10R,1.0X1.0*%\nD10*\nX0Y0D03*\n%LPD*%\nM02*\n";
    const GerberImage img = parse_gerber(g);
    CHECK(img.warnings.empty());
    CHECK(total_mm2(img) == doctest::Approx(kPi * 4.0 - 1.0).epsilon(1e-5));
}

TEST_CASE("standard apertures with holes, obrounds and polygons") {
    const std::string g = std::string(kHeader) +
                          "%ADD10C,2.0X1.0*%\n%ADD11O,3.0X1.0*%\n%ADD12P,2.0X4X45*%\n%ADD13R,2.0X1.0X0.5*%\n"
                          "D10*\nX0Y0D03*\nD11*\nX10000000Y0D03*\nD12*\nX20000000Y0D03*\nD13*\nX30000000Y0D03*\nM02*\n";
    const GerberImage img = parse_gerber(g);
    REQUIRE(img.features.size() == 4);
    auto area = [&](int i) { return mm2(union_polygons(img.features[i].shape)); };
    CHECK(area(0) == doctest::Approx(kPi * (1.0 - 0.25)).epsilon(1e-5));
    CHECK(area(1) == doctest::Approx(2.0 * 1.0 + kPi * 0.25).epsilon(1e-5));
    CHECK(area(2) == doctest::Approx(2.0).epsilon(1e-6));  // square with 2 mm diagonal
    CHECK(area(3) == doctest::Approx(2.0 - kPi * 0.0625).epsilon(1e-5));
}

TEST_CASE("the KiCad RoundRect macro gives a rounded rectangle") {
    const std::string g = std::string(kHeader) +
                          "%AMRoundRect*\n0 Rectangle with rounded corners*\n0 $1 Rounding radius*\n"
                          "4,1,4,$2,$3,$4,$5,$6,$7,$8,$9,$2,$3,0*\n1,1,$1+$1,$2,$3*\n1,1,$1+$1,$4,$5*\n"
                          "1,1,$1+$1,$6,$7*\n1,1,$1+$1,$8,$9*\n20,1,$1+$1,$2,$3,$4,$5,0*\n20,1,$1+$1,$4,$5,$6,$7,0*\n"
                          "20,1,$1+$1,$6,$7,$8,$9,0*\n20,1,$1+$1,$8,$9,$2,$3,0*%\n"
                          // 2 x 1 mm with 0.25 mm radius: corner centres at +-0.75, +-0.25.
                          "%ADD20RoundRect,0.250000X-0.750000X-0.250000X0.750000X-0.250000X0.750000X0.250000X-0.750000X0.250000X0*%\n"
                          "D20*\nX0Y0D03*\nM02*\n";
    const GerberImage img = parse_gerber(g);
    CHECK(img.warnings.empty());
    REQUIRE(img.features.size() == 1);
    const double r = 0.25;
    CHECK(mm2(union_polygons(img.features[0].shape)) == doctest::Approx(2.0 * 1.0 - (4 - kPi) * r * r).epsilon(1e-4));
}

TEST_CASE("macro primitives: centre line rotation, outline, thermal, exposure off") {
    const std::string g = std::string(kHeader) +
                          "%AMBOX*\n21,1,2,1,0,0,90*\n1,0,0.5,0,0*%\n"  // 2x1 rotated, minus a 0.5 hole
                          "%AMTRI*\n4,1,3,0,0,1,0,0,1,0,0,0*%\n"
                          "%AMTH*\n7,0,0,2,1,0.2,0*%\n"
                          "%ADD10BOX*%\n%ADD11TRI*%\n%ADD12TH*%\n"
                          "D10*\nX0Y0D03*\nD11*\nX0Y0D03*\nD12*\nX0Y0D03*\nM02*\n";
    const GerberImage img = parse_gerber(g);
    CHECK(img.warnings.empty());
    REQUIRE(img.features.size() == 3);
    const auto box = union_polygons(img.features[0].shape);
    CHECK(mm2(box) == doctest::Approx(2.0 - kPi * 0.0625).epsilon(1e-5));
    const Rect64 b = bounds(box);
    CHECK((b.bottom - b.top) == doctest::Approx(2e6).epsilon(1e-6));  // rotated: tall
    CHECK(mm2(union_polygons(img.features[1].shape)) == doctest::Approx(0.5).epsilon(1e-9));
    // Thermal: annulus 2/1 mm minus two 0.2 mm gaps crossing it.
    const double annulus = kPi * (1.0 - 0.25);
    const double gaps = 4 * 0.2 * 0.5 - 0.0;  // four ring segments 0.5 mm long, 0.2 mm wide (approximately)
    CHECK(mm2(union_polygons(img.features[2].shape)) == doctest::Approx(annulus - gaps).epsilon(0.01));
}

TEST_CASE("inch units with trailing zeros omitted") {
    const std::string g = "%FSTAX24Y24*%\n%MOIN*%\n%ADD10R,0.1X0.1*%\nD10*\nX1Y1D03*\nM02*\n";
    const GerberImage img = parse_gerber(g);
    REQUIRE(img.features.size() == 1);
    // X1 with 2.4 format and trailing zeros omitted is 10.0000 in.
    CHECK(img.features[0].at.x == std::llround(10.0 * 25.4e6));
    CHECK(mm2(union_polygons(img.features[0].shape)) == doctest::Approx(2.54 * 2.54).epsilon(1e-6));
}

TEST_CASE("aperture blocks and step and repeat") {
    const std::string g = std::string(kHeader) +
                          "%ADD10R,1.0X1.0*%\n"
                          "%ABD20*%\nD10*\nX0Y0D03*\nX2000000Y0D03*\n%AB*%\n"
                          "%SRX3Y2I10.0J5.0*%\nD20*\nX0Y0D03*\n%SR*%\nM02*\n";
    const GerberImage img = parse_gerber(g);
    CHECK(img.warnings.empty());
    CHECK(img.features.size() == 2 * 3 * 2);
    CHECK(total_mm2(img) == doctest::Approx(12.0).epsilon(1e-9));
}

TEST_CASE("X2 attributes in KiCad compatibility comments") {
    const std::string g = std::string(kHeader) + "G04 #@! TF.FileFunction,Copper,L2,Bot*\n%ADD10C,1.0*%\nD10*\n"
                                                 "G04 #@! TO.N,GND*\nG04 #@! TO.P,R1,2*\nX0Y0D03*\nM02*\n";
    const GerberImage img = parse_gerber(g);
    CHECK(img.file_function() == std::vector<std::string>{"Copper", "L2", "Bot"});
    REQUIRE(img.features.size() == 1);
    CHECK(img.features[0].net == "GND");
    CHECK(img.features[0].component == "R1");
    CHECK(img.features[0].pin == "2");
    CHECK(img.features[0].is_pad());
}

TEST_CASE("single-quadrant arcs choose the short centre") {
    const std::string g = std::string(kHeader) + "%ADD10C,0.1*%\nD10*\nG74*\nX1000000Y0D02*\nG03*\nX0Y1000000I1000000J0D01*\nM02*\n";
    const GerberImage img = parse_gerber(g);
    REQUIRE(img.features.size() == 1);
    REQUIRE(img.features[0].arc_centre);
    CHECK(*img.features[0].arc_centre == Point64(0, 0));
}
