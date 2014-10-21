check_cxx_source_compiles("
  #include <thread>
  int main() {}" MP_THREAD)
if (MP_THREAD)
  add_definitions(-DMP_THREAD=1)
endif ()

if (TARGET gecode)
  add_mp_test(gecode-test PREFIX gecode-test.cc
    LIBS amplgecode-static solver-impl-test)
endif ()

if (TARGET ilogcp)
  include_directories(${PROJECT_SOURCE_DIR}/thirdparty/gtest)
  # Include the gtest source in test instead of linking with the library to
  # avoid a problem when gtest is built with a different runtime than ilogcp.
  add_prefix(ILOGCP_TEST_SRC ../ ${GTEST_SRC})
  add_executable(ilogcp-test ilogcp-test.cc feature.h solver-impl-test.h
                 ${ILOGCP_TEST_SRC})
  target_link_libraries(ilogcp-test amplilogcp-static)
  add_test(NAME ilogcp-test COMMAND $<TARGET_FILE:ilogcp-test>)
endif ()

if (TARGET jacop)
  add_mp_test(jacop-test jacop-test.cc solver-impl-test.h LIBS ampljacop-static)
  if (MSVC)
    set_target_properties(jacop-test PROPERTIES LINK_FLAGS "/DELAYLOAD:jvm.dll")
    target_link_libraries(jacop-test delayimp)
  endif ()
endif ()

if (TARGET smpswriter)
  add_mp_test(smpswriter-test smpswriter-test.cc LIBS amplsmpswriter)
endif ()
