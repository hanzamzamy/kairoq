cmake_minimum_required(VERSION 3.24)

project(kairoq VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(KAIROQ_CLI "Build the command-line tool" ON)
option(KAIROQ_PYTHON "Build the Python module" OFF)
option(KAIROQ_TESTS "Build the C++ tests" ON)

include(cmake/Dependencies.cmake)
include(cmake/EmbedFiles.cmake)

add_subdirectory(src/core)

if(KAIROQ_CLI)
    add_subdirectory(src/cli)
endif()

if(KAIROQ_PYTHON)
    add_subdirectory(src/python)
endif()

if(KAIROQ_TESTS)
    enable_testing()
    add_subdirectory(tests/cpp)
endif()
