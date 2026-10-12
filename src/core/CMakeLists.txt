set(VIEWER_DIR ${PROJECT_SOURCE_DIR}/viewer)
set(VIEWER_FILES
    ${VIEWER_DIR}/index.html
    ${VIEWER_DIR}/main.js
    ${VIEWER_DIR}/style.css
    ${VIEWER_DIR}/vendor/three.module.js
    ${VIEWER_DIR}/vendor/three.core.js
    ${VIEWER_DIR}/vendor/OrbitControls.js)
set(VIEWER_SOURCE ${CMAKE_CURRENT_BINARY_DIR}/viewer_files.cpp)
kairoq_embed_files(OUTPUT ${VIEWER_SOURCE} NAME viewer_files BASE ${VIEWER_DIR} FILES ${VIEWER_FILES})
set(OPENEMS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/resources/openems)
set(OPENEMS_SOURCE ${CMAKE_CURRENT_BINARY_DIR}/openems_files.cpp)
kairoq_embed_files(OUTPUT ${OPENEMS_SOURCE} NAME openems_files BASE ${OPENEMS_DIR} FILES ${OPENEMS_DIR}/simulate.py ${OPENEMS_DIR}/simulate_ports.py)

set(BUILD_INFO_SOURCE ${CMAKE_CURRENT_BINARY_DIR}/build_info.cpp)
add_custom_target(kairoq_build_info
    COMMAND ${CMAKE_COMMAND} -DOUTPUT=${BUILD_INFO_SOURCE} -DSOURCE_DIR=${PROJECT_SOURCE_DIR} -DVERSION=${PROJECT_VERSION}
            "-DCOMPILER=${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}" "-DPLATFORM=${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR}"
            -P ${PROJECT_SOURCE_DIR}/cmake/BuildInfo.cmake
    BYPRODUCTS ${BUILD_INFO_SOURCE}
    COMMENT "Recording the build information")

add_library(kairoq_core STATIC
    src/archive.cpp
    src/board.cpp
    src/build.cpp
    src/cross_section.cpp
    src/dc.cpp
    src/deembed.cpp
    src/drill.cpp
    src/emi.cpp
    src/elmer.cpp
    src/eye.cpp
    src/emi_loop.cpp
    src/fabrication.cpp
    src/feature.cpp
    src/fem_util.cpp
    src/geometry.cpp
    src/ibis.cpp
    src/gerber.cpp
    src/ipc2581.cpp
    src/mesh.cpp
    src/mesh2d.cpp
    src/network.cpp
    src/model_io.cpp
    src/odbpp.cpp
    src/openems.cpp
    src/open.cpp
    src/peec.cpp
    src/project.cpp
    src/readers.cpp
    src/rf.cpp
    src/server.cpp
    src/si.cpp
    src/thermal.cpp
    src/cavity.cpp
    src/pdn.cpp
    src/via_model.cpp
    src/layered_green.cpp
    src/emi_capacitance.cpp
    src/timedomain.cpp
    src/provenance.cpp
    src/sha256.cpp
    ${BUILD_INFO_SOURCE}
    ${VIEWER_SOURCE}
    ${OPENEMS_SOURCE})
add_library(kairoq::core ALIAS kairoq_core)
add_dependencies(kairoq_core kairoq_build_info)

target_include_directories(kairoq_core PUBLIC include)
target_link_libraries(kairoq_core
    PUBLIC Clipper2 nlohmann_json::nlohmann_json
    PRIVATE pugixml::pugixml yaml-cpp::yaml-cpp miniz kairoq_earcut kairoq_httplib kairoq_cdt kairoq_amgcl)
find_package(Threads REQUIRED)
target_link_libraries(kairoq_core PUBLIC Threads::Threads)
if(WIN32)
    target_link_libraries(kairoq_core PRIVATE ws2_32)
endif()
target_compile_options(kairoq_core PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wno-missing-field-initializers>)
