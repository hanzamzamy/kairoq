find_package(Python 3.10 COMPONENTS Interpreter Development.Module REQUIRED)
execute_process(
    COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_ROOT)
find_package(nanobind CONFIG REQUIRED)

nanobind_add_module(_core NB_STATIC bindings.cpp)
target_link_libraries(_core PRIVATE kairoq::core)

# Development builds: place the module next to the package sources' build copy, so
# that PYTHONPATH=<build>/src/python imports it together with the pure Python files.
set_target_properties(_core PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/kairoq)
file(GLOB PACKAGE_FILES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/kairoq/*.py)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/kairoq)
foreach(f IN LISTS PACKAGE_FILES)
    get_filename_component(name ${f} NAME)
    configure_file(${f} ${CMAKE_CURRENT_BINARY_DIR}/kairoq/${name} COPYONLY)
endforeach()

install(TARGETS _core LIBRARY DESTINATION kairoq)
if(DEFINED SKBUILD_SCRIPTS_DIR AND KAIROQ_CLI)
    install(TARGETS kairoq RUNTIME DESTINATION ${SKBUILD_SCRIPTS_DIR})
endif()
