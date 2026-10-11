#include "kairoq/drill.hpp"
#include "kairoq/fabrication.hpp"

#include <doctest/doctest.h>

using namespace kairoq;

TEST_CASE("Excellon with decimal coordinates and X2 comments (KiCad)") {
    const char* text =
        "M48\n; DRILL file KiCad\n; FORMAT={-:-/ absolute / metric / decimal}\n"
        "; #@! TF.FileFunction,Plated,1,4,PTH\nFMAT,2\nMETRIC\n"
        "; #@! TA.AperFunction,Plated,PTH,ViaDrill\nT1C0.300\n"
        "; #@! TA.AperFunction,Plated,PTH,ComponentDrill\nT2C1.000\n%\nG90\nG05\n"
        "T1\nX126.365Y-125.095\nX1.0Y2.0\nT2\nX10.0Y-5.0G85X12.0Y-5.0\nM30\n";
    const DrillFile d = parse_excellon(text);
    CHECK(d.warnings.empty());
    REQUIRE(d.plated);
    CHECK(*d.plated);
    CHECK(d.from_layer == 1);
    CHECK(d.to_layer == 4);
    REQUIRE(d.tools.size() == 2);
    CHECK(d.tools[0].function == "ViaDrill");
    CHECK(d.tools[0].diameter_nm == doctest::Approx(300000));
    CHECK(d.tools[1].function == "ComponentDrill");
    REQUIRE(d.hits.size() == 3);
    CHECK(d.hits[0].at == Point64(126365000, -125095000));
    CHECK(d.hits[2].at == Point64(10000000, -5000000));
    REQUIRE(d.hits[2].slot_end);
    CHECK(*d.hits[2].slot_end == Point64(12000000, -5000000));
}

TEST_CASE("Excellon integer formats: leading and trailing zero suppression") {
    // INCH,LZ: leading zeros kept, trailing suppressed, 2.4 format.
    const DrillFile lz = parse_excellon("M48\nINCH,LZ\nT1C0.0400\n%\nT1\nX01Y-0015\nM30\n");
    REQUIRE(lz.hits.size() == 1);
    CHECK(lz.hits[0].at.x == std::llround(1.0 * 25.4e6));
    CHECK(lz.hits[0].at.y == std::llround(-0.15 * 25.4e6));
    CHECK(lz.tools[0].diameter_nm == doctest::Approx(0.04 * 25.4e6));
    // METRIC,TZ: trailing zeros kept, leading suppressed, 3.3 format.
    const DrillFile tz = parse_excellon("M48\nMETRIC,TZ\nT1C0.8\n%\nT1\nX12500Y-300\nM30\n");
    REQUIRE(tz.hits.size() == 1);
    CHECK(tz.hits[0].at == Point64(12500000, -300000));
}

TEST_CASE("Excellon routed slots") {
    const DrillFile d = parse_excellon("M48\nMETRIC\nT1C1.0\n%\nT1\nG00X1.0Y1.0\nM15\nG01X3.0Y1.0\nM16\nG05\nX5.0Y5.0\nM30\n");
    REQUIRE(d.hits.size() == 2);
    REQUIRE(d.hits[0].slot_end);
    CHECK(d.hits[0].at == Point64(1000000, 1000000));
    CHECK(*d.hits[0].slot_end == Point64(3000000, 1000000));
    CHECK_FALSE(d.hits[1].slot_end);
}

TEST_CASE("IPC-D-356 records (KiCad, inch units)") {
    const char* text =
        "P  CODE 00\nP  UNITS CUST 0\n"
        "317NET-(P3-P1)      C1    -1    D0394PA00X+055750Y-039250X0787Y0787R270S0\n"
        "327GND              U1    -14         A01X+001000Y-000500X0100Y0200R000S0\n999\n";
    const Netlist n = parse_ipc356(text);
    CHECK(n.warnings.empty());
    REQUIRE(n.points.size() == 2);
    CHECK(n.points[0].net == "NET-(P3-P1)");
    CHECK(n.points[0].refdes == "C1");
    CHECK(n.points[0].pin == "1");
    CHECK(n.points[0].plated);
    CHECK(n.points[0].access == 0);
    CHECK(n.points[0].hole_nm == doctest::Approx(0.0394 * 25.4e6));
    CHECK(n.points[0].at == Point64(141605000, -99695000));
    CHECK(n.points[1].net == "GND");
    CHECK(n.points[1].pin == "14");
    CHECK(n.points[1].access == 1);
    CHECK(n.points[1].hole_nm == 0.0);
}

TEST_CASE("Gerber job file stackup") {
    const char* text = R"({"Header":{"GenerationSoftware":{"Vendor":"KiCad","Application":"Pcbnew","Version":"10"}},
        "GeneralSpecs":{"ProjectId":{"Name":"demo"},"LayerNumber":2,"BoardThickness":1.6},
        "FilesAttributes":[{"Path":"a.gtl","FileFunction":"Copper,L1,Top","FilePolarity":"Positive"}],
        "MaterialStackup":[{"Type":"Legend","Name":"silk"},{"Type":"SolderMask","Thickness":0.01,"Name":"m"},
          {"Type":"Copper","Thickness":0.035,"Name":"F.Cu"},
          {"Type":"Dielectric","Thickness":1.51,"Material":"FR4","DielectricConstant":"4.5","LossTangent":0.02,"Name":"core"},
          {"Type":"Copper","Thickness":0.035,"Name":"B.Cu"}]})";
    const JobFile j = parse_job_file(text);
    CHECK(j.generator == "KiCad Pcbnew 10");
    CHECK(j.project == "demo");
    CHECK(j.layer_count == 2);
    REQUIRE(j.files.size() == 1);
    CHECK(j.files[0].function == std::vector<std::string>{"Copper", "L1", "Top"});
    REQUIRE(j.stackup.size() == 4);
    CHECK(j.stackup[0].type == LayerType::SolderMask);
    CHECK(j.stackup[1].copper_index == 1);
    CHECK(j.stackup[2].thickness_nm == doctest::Approx(1.51e6));
    CHECK(*j.stackup[2].epsilon_r == doctest::Approx(4.5));
    CHECK(j.stackup[3].copper_index == 2);
}

TEST_CASE("placement CSV: KiCad and Altium headers") {
    const auto k = parse_placement_csv("Ref,Val,Package,PosX,PosY,Rot,Side\n\"C1\",\"10uF\",\"CP\",141.605,-99.695,90,top\n"
                                       "\"R1\",\"1k\",\"R\",1,2,0,bottom\n");
    REQUIRE(k.size() == 2);
    CHECK(k[0].refdes == "C1");
    CHECK(k[0].at == Point64(141605000, -99695000));
    CHECK(k[0].rotation_deg == 90);
    CHECK(k[1].side == Side::Bottom);
    const auto a = parse_placement_csv("Designator,Comment,Layer,Footprint,Center-X(mil),Center-Y(mil),Rotation\n"
                                       "U1,MCU,TopLayer,QFN,100,200,45\n");
    REQUIRE(a.size() == 1);
    CHECK(a[0].at == Point64(2540000, 5080000));
    CHECK(a[0].side == Side::Top);
}

TEST_CASE("stackup YAML") {
    const auto s = parse_stackup_yaml(
        "layers:\n"
        "  - {name: top, type: copper, thickness_mm: 0.035}\n"
        "  - {name: core, type: dielectric, thickness_mm: 1.5, material: FR-4, epsilon_r: 4.3}\n"
        "  - {name: bottom, type: copper, thickness_mm: 0.07}\n");
    REQUIRE(s.size() == 3);
    CHECK(s[0].copper_index == 1);
    CHECK(s[2].copper_index == 2);
    CHECK(s[2].thickness_nm == doctest::Approx(70000));
    CHECK(s[1].material == "FR-4");
    CHECK_THROWS(parse_stackup_yaml("layers:\n  - {name: x, type: copper}\n"));
}
