cmake_minimum_required(VERSION 3.10)

project(NLWriter2)


set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)


set(NLW2_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(NLW2_LIB_FILES
        ${NLW2_DIR}/src/nl-writer2.cc
        ${NLW2_DIR}/src/nl-utils2.cc
        ${NLW2_DIR}/src/dtoa.c)

# Create library
set(NLW2_LIB_NAME "nlw2")
add_library(${NLW2_LIB_NAME} STATIC ${NLW2_LIB_FILES})
target_include_directories(
    ${NLW2_LIB_NAME} PUBLIC ${NLW2_DIR}/include)

if (MSVC)               ## Set the same as for your parent project
##  target_compile_options(
##      ${NLW2_LIB_NAME} PRIVATE ${AMPL_MSVC_COMPILE_OPTIONS})
endif()

# NLW2 example
add_executable(nl-writer-example
            ${NLW2_DIR}/examples/nlsol_ex.cc)
target_include_directories(nl-writer-example PUBLIC
            ${NLW2_DIR}/include)
target_link_libraries(nl-writer-example ${NLW2_LIB_NAME})
