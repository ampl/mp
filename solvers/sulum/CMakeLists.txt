# CMake build script for the AMPL/Sulum solver.

add_definitions(-DSLM_NO_BOOST)
include_directories(amplsulum-static ${SULUM_INCLUDE_DIRS})
add_library(amplsulum-static sulum.cc)
add_dependencies(amplsulum-static arith_h)
target_link_libraries(amplsulum-static ${SULUM_LIBRARIES} util)

add_library(amplsulum SHARED ../util/solver-c.cc)
target_link_libraries(amplsulum amplsulum-static)

add_executable(sulum main.cc)
target_link_libraries(sulum amplsulum-static)
