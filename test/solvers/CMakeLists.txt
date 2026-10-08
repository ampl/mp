check_cxx_source_compiles("
  #include <thread>
  int main() {}" MP_THREAD)
if (MP_THREAD)
  add_definitions(-DMP_THREAD=1)
endif ()

if (TARGET gecode)
  add_mp_test(gecode-mp-test gecode-test.cc feature.h nl-solver-test.h
    LIBS amplgecode-static)
endif ()

if (TARGET ilogcp)
  # Include the gtest source in test instead of linking with the library to
  # avoid a problem when gtest is built with a different runtime than ilogcp.
  set(ILOGCP_TEST_SRC
      ${MP_GTEST_DIR}/googletest/src/gtest-all.cc
      ${MP_GTEST_DIR}/googlemock/src/gmock-all.cc
      ${PROJECT_SOURCE_DIR}/test/gtest-extra.cc
      ${PROJECT_SOURCE_DIR}/test/util.cc
      ${PROJECT_SOURCE_DIR}/test/test-main.cc)
  add_executable(ilogcp-test ilogcp-test.cc feature.h nl-solver-test.h
                 ${ILOGCP_TEST_SRC})
  target_include_directories(ilogcp-test PRIVATE
      ${MP_GTEST_INCLUDE_DIRS} ${MP_GTEST_SOURCE_DIRS})
  target_compile_definitions(ilogcp-test PRIVATE FMT_USE_FILE_DESCRIPTORS)
  target_link_libraries(ilogcp-test amplilogcp-static)
  add_test(NAME ilogcp-test COMMAND $<TARGET_FILE:ilogcp-test>)
  
  add_executable(numberofmap-speed-test numberofmap-speed-test.cc)
  target_link_libraries(numberofmap-speed-test amplilogcp-static)
endif ()

if (TARGET jacop)
  add_mp_test(jacop-test jacop-test.cc feature.h nl-solver-test.h
    LIBS ampljacop-static)
  if (MSVC)
    set_target_properties(jacop-test PROPERTIES LINK_FLAGS "/DELAYLOAD:jvm.dll")
    target_link_libraries(jacop-test delayimp)
  endif ()
endif ()

if (TARGET localsolver)
  add_mp_test(localsolver-test localsolver-test.cc feature.h nl-solver-test.h
    LIBS ampllocalsolver-static)
endif ()
