cmake_minimum_required(VERSION 3.10)

project(NLWriter2)


set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)


set(NLW2_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(NLW2_LIB_FILES
        ${NLW2_DIR}/src/nl-writer2.cc
        ${NLW2_DIR}/src/nl-utils2.cc
        ${NLW2_DIR}/src/dtoa.c)
set(NLW2_INC_FILES
              ${NLW2_DIR}/include/mp/nlsol.h
              ${NLW2_DIR}/include/mp/nl-header.h
              ${NLW2_DIR}/include/mp/nl-feeder2.h
              ${NLW2_DIR}/include/mp/nl-writer2.h
              ${NLW2_DIR}/include/mp/nl-writer2.hpp
              ${NLW2_DIR}/include/mp/nl-writer2-misc.h
              ${NLW2_DIR}/include/mp/nl-utils2.h
              ${NLW2_DIR}/include/mp/sol-reader2.h
              ${NLW2_DIR}/include/mp/sol-reader2.hpp
              ${NLW2_DIR}/include/mp/sol-handler2.h
              )

# Create library
set(NLW2_LIB_NAME "nlw2")
add_library(${NLW2_LIB_NAME} STATIC
  ${NLW2_LIB_FILES} ${NLW2_INC_FILES})
target_include_directories(
    ${NLW2_LIB_NAME} PUBLIC ${NLW2_DIR}/include)

if (MSVC)               ## Set the same as for your parent project
##  target_compile_options(
##      ${NLW2_LIB_NAME} PRIVATE ${AMPL_MSVC_COMPILE_OPTIONS})
endif()

# NLW2 example
set(NLSOL_EX_FILES
  ${NLW2_DIR}/examples/nlsol_ex.cc
  ${NLW2_DIR}/examples/nlsol_ex_mdl.h
  ${NLW2_DIR}/examples/nlsol_ex_nl.h
  ${NLW2_DIR}/examples/nlsol_ex_sol.h)

option(BUILD_EXAMPLES_NLW2
  "Build NLW2 examples (library not affected)" ON)
if (BUILD_EXAMPLES_NLW2)         # Actually build the example
  add_executable(nl-writer-example ${NLSOL_EX_FILES})
  target_include_directories(nl-writer-example PUBLIC
    ${NLW2_DIR}/include)
  target_link_libraries(nl-writer-example ${NLW2_LIB_NAME})
endif()

# Export library include path
set(NLW2_INCLUDE_PATH ${NLW2_DIR}/include PARENT_SCOPE)
# Export source files
set(NLW2_SRC
  ${NLW2_LIB_FILES} ${NLW2_INC_FILES} ${NLSOL_EX_FILES}
  PARENT_SCOPE)
