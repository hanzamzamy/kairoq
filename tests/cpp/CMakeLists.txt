add_executable(kairoq_tests
    main.cpp
    test_geometry.cpp
    test_gerber.cpp
    test_readers.cpp
    test_model.cpp
    test_mesh2d.cpp
    test_peec.cpp
    test_cavity.cpp
    test_via_model.cpp
    test_layered_green.cpp
    test_bessel.cpp
    test_sha256.cpp)
target_link_libraries(kairoq_tests PRIVATE kairoq::core kairoq_doctest)
target_include_directories(kairoq_tests PRIVATE ${PROJECT_SOURCE_DIR}/src/core/src)
target_compile_definitions(kairoq_tests PRIVATE KAIROQ_TESTDATA="${PROJECT_SOURCE_DIR}/testdata")
add_test(NAME cpp COMMAND kairoq_tests)
