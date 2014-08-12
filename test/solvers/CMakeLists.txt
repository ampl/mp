add_library(solver-test ../args.cc solver-test.cc)
target_link_libraries(solver-test asl)
target_compile_definitions(solver-test PRIVATE
  MP_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../data")

include(CheckCXXSourceCompiles)
check_cxx_source_compiles(
  "#include <memory>
  int main() { std::unique_ptr<int> p; }" HAVE_UNIQUE_PTR)
if (HAVE_UNIQUE_PTR)
  target_compile_definitions(solver-test PUBLIC MP_USE_UNIQUE_PTR)
endif ()

if (TARGET gecode)
  add_mp_test(gecode-test PREFIX gecode-test.cc
    LIBS amplgecode-static solver-test)
endif ()

if (TARGET ilogcp)
  include_directories(${PROJECT_SOURCE_DIR}/thirdparty/gtest)
  # Include the gtest source in test instead of linking with the library to
  # avoid a problem when gtest is built with a different runtime than ilogcp.
  add_prefix(ILOGCP_TEST_SRC ../ ${GTEST_SRC})
  add_executable(ilogcp-test ilogcp-test.cc ${ILOGCP_TEST_SRC})
  target_link_libraries(ilogcp-test amplilogcp-static solver-test)
  add_test(NAME ilogcp-test COMMAND $<TARGET_FILE:ilogcp-test>)
endif ()

if (TARGET jacop)
  add_mp_test(jacop-test jacop-test.cc LIBS ampljacop-static solver-test)
  if (MSVC)
    set_target_properties(jacop-test PROPERTIES LINK_FLAGS "/DELAYLOAD:jvm.dll")
    target_link_libraries(jacop-test delayimp)
  endif ()
endif ()

if (TARGET amplsmpswriter)
  add_mp_test(smpswriter-test ../args.cc smpswriter-test.cc LIBS amplsmpswriter)
endif ()
