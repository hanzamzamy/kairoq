# Third-party dependencies, pinned by version and SHA-256, built from source and linked
# statically so that a build does not depend on the libraries installed on the host.
# For an offline build, point FETCHCONTENT_SOURCE_DIR_<NAME> at an unpacked copy of
# the archive.

include(FetchContent)

set(FETCHCONTENT_QUIET ON)

# Clipper2: polygon clipping and offsetting (BSL-1.0).
set(CLIPPER2_UTILS OFF CACHE BOOL "" FORCE)
set(CLIPPER2_EXAMPLES OFF CACHE BOOL "" FORCE)
set(CLIPPER2_TESTS OFF CACHE BOOL "" FORCE)
set(CLIPPER2_USINGZ OFF CACHE STRING "" FORCE)
FetchContent_Declare(Clipper2
    URL https://github.com/AngusJohnson/Clipper2/archive/refs/tags/Clipper2_2.0.1.tar.gz
    URL_HASH SHA256=2a3693aceab4aed3e39b743e038d87701acc53cf05ed7b2013aab3e0aec5287e
    SOURCE_SUBDIR CPP
    DOWNLOAD_EXTRACT_TIMESTAMP ON)

# nlohmann/json: JSON reading and writing (MIT).
set(JSON_BuildTests OFF CACHE BOOL "" FORCE)
FetchContent_Declare(nlohmann_json
    URL https://github.com/nlohmann/json/releases/download/v3.12.0/json.tar.xz
    URL_HASH SHA256=42f6e95cad6ec532fd372391373363b62a14af6d771056dbfc86160e6dfff7aa
    DOWNLOAD_EXTRACT_TIMESTAMP ON)

# pugixml: XML parsing for IPC-2581 (MIT).
set(PUGIXML_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
FetchContent_Declare(pugixml
    URL https://github.com/zeux/pugixml/archive/refs/tags/v1.16.tar.gz
    URL_HASH SHA256=357bcab8877dc9943f355d3a72daba1b053238ba955f50fa81586afb65090219
    DOWNLOAD_EXTRACT_TIMESTAMP ON)

# yaml-cpp: stackup and project files (MIT).
set(YAML_CPP_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(YAML_CPP_BUILD_TOOLS OFF CACHE BOOL "" FORCE)
set(YAML_CPP_BUILD_CONTRIB OFF CACHE BOOL "" FORCE)
set(YAML_CPP_INSTALL OFF CACHE BOOL "" FORCE)
FetchContent_Declare(yaml-cpp
    URL https://github.com/jbeder/yaml-cpp/archive/refs/tags/yaml-cpp-0.9.0.tar.gz
    URL_HASH SHA256=25cb043240f828a8c51beb830569634bc7ac603978e0f69d6b63558dadefd49a
    DOWNLOAD_EXTRACT_TIMESTAMP ON)

# miniz: zip archives for the canonical model and compressed IPC-2581 (MIT).
FetchContent_Declare(miniz
    URL https://github.com/richgel999/miniz/archive/refs/tags/3.1.2.tar.gz
    URL_HASH SHA256=98468f8924934b723276680f85238b6c78bf1f8b49b4459cc9b7214a20e2e9fb
    DOWNLOAD_EXTRACT_TIMESTAMP ON)

# earcut.hpp: polygon triangulation for display meshes (ISC). Header only.
FetchContent_Declare(earcut
    URL https://github.com/mapbox/earcut.hpp/archive/refs/tags/v3.2.4.tar.gz
    URL_HASH SHA256=70c00e90811950b2c3d973ecfd46010680b946f7952306ffa2785c0368fe2db0
    DOWNLOAD_EXTRACT_TIMESTAMP ON
    SOURCE_SUBDIR do-not-configure)

# cpp-httplib: the local viewer server (MIT). Header only.
FetchContent_Declare(httplib
    URL https://github.com/yhirose/cpp-httplib/archive/refs/tags/v0.60.0.tar.gz
    URL_HASH SHA256=f8d6bc6e2a958f9f21bde11ca408bd987cdf91b99ed9fc63c3e4ea3439376e71
    DOWNLOAD_EXTRACT_TIMESTAMP ON
    SOURCE_SUBDIR do-not-configure)

# CDT: constrained Delaunay triangulation and refinement for solver meshes (MPL-2.0).
# Header only.
FetchContent_Declare(cdt
    URL https://github.com/artem-ogre/CDT/archive/refs/tags/2.0.0.tar.gz
    URL_HASH SHA256=81ed0c7a8cbdc059b26948a160c3ce68c7a147976334093fe3175c8fa161a447
    DOWNLOAD_EXTRACT_TIMESTAMP ON
    SOURCE_SUBDIR do-not-configure)

# AMGCL: algebraic multigrid preconditioned solvers (MIT). Header only; used without
# Boost.
FetchContent_Declare(amgcl
    URL https://github.com/ddemidov/amgcl/archive/refs/tags/1.5.0.tar.gz
    URL_HASH SHA256=6794cab752f0f15bc40db9c8a08a5ad761adf7d74b3d6652b9e2d4bf7f19ab31
    DOWNLOAD_EXTRACT_TIMESTAMP ON
    SOURCE_SUBDIR do-not-configure)

FetchContent_MakeAvailable(Clipper2 nlohmann_json pugixml yaml-cpp miniz earcut httplib cdt amgcl)

add_library(kairoq_cdt INTERFACE)
target_include_directories(kairoq_cdt SYSTEM INTERFACE ${cdt_SOURCE_DIR}/CDT/include)

add_library(kairoq_amgcl INTERFACE)
target_include_directories(kairoq_amgcl SYSTEM INTERFACE ${amgcl_SOURCE_DIR})
target_compile_definitions(kairoq_amgcl INTERFACE AMGCL_NO_BOOST)
find_package(OpenMP COMPONENTS CXX)
if(OpenMP_CXX_FOUND)
    target_link_libraries(kairoq_amgcl INTERFACE OpenMP::OpenMP_CXX)
endif()

if(NOT TARGET pugixml::pugixml AND TARGET pugixml)
    add_library(pugixml::pugixml ALIAS pugixml)
endif()
if(NOT TARGET yaml-cpp::yaml-cpp AND TARGET yaml-cpp)
    add_library(yaml-cpp::yaml-cpp ALIAS yaml-cpp)
endif()

add_library(kairoq_earcut INTERFACE)
target_include_directories(kairoq_earcut SYSTEM INTERFACE ${earcut_SOURCE_DIR}/include)

add_library(kairoq_httplib INTERFACE)
target_include_directories(kairoq_httplib SYSTEM INTERFACE ${httplib_SOURCE_DIR})

if(KAIROQ_TESTS)
    # doctest: C++ unit tests (MIT). Header only.
    FetchContent_Declare(doctest
        URL https://github.com/doctest/doctest/archive/refs/tags/v2.5.3.tar.gz
        URL_HASH SHA256=174ebc4e769928959614789c5b4e9c3d0a0f81a62bb608756b127bfebfb21331
        DOWNLOAD_EXTRACT_TIMESTAMP ON
        SOURCE_SUBDIR do-not-configure)
    FetchContent_MakeAvailable(doctest)
    add_library(kairoq_doctest INTERFACE)
    target_include_directories(kairoq_doctest SYSTEM INTERFACE ${doctest_SOURCE_DIR})
endif()
