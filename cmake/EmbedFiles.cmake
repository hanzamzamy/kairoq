# Embed files into a C++ source as byte arrays.
#
# kairoq_embed_files(OUTPUT <file.cpp> NAME <name> BASE <dir> FILES <file>...)
# generates <file.cpp> defining the arrays kairoq::<name> and kairoq::<name>_count
# (declared in include/kairoq/embedded.hpp), each file keyed by its
# path relative to BASE. The source is regenerated when an input changes.

add_executable(kairoq_embed ${CMAKE_CURRENT_LIST_DIR}/embed.cpp)

function(kairoq_embed_files)
    cmake_parse_arguments(ARG "" "OUTPUT;NAME;BASE" "FILES" ${ARGN})
    add_custom_command(
        OUTPUT ${ARG_OUTPUT}
        COMMAND kairoq_embed ${ARG_OUTPUT} ${ARG_NAME} ${ARG_BASE} ${ARG_FILES}
        DEPENDS kairoq_embed ${ARG_FILES}
        COMMENT "Embedding ${ARG_BASE}"
        VERBATIM)
endfunction()
