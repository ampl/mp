# CMake build script for the SMPL writer.

add_library(amplsmpswriter smpswriter.cc)
target_link_libraries(amplsmpswriter util)
add_dependencies(amplsmpswriter arith_h)
add_executable(smpswriter main.cc)
target_link_libraries(smpswriter amplsmpswriter)
