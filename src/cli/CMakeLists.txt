add_executable(kairoq main.cpp)
target_link_libraries(kairoq PRIVATE kairoq::core)
# Python wheels install the command into the scripts directory (src/python).
if(NOT SKBUILD)
    install(TARGETS kairoq RUNTIME DESTINATION bin)
endif()
