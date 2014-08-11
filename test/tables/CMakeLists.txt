find_package(ODBC)
if (ODBC_FOUND)
  include_directories(${ODBC_INCLUDE_DIRS})

  set(CMAKE_REQUIRED_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/..)
  set(CMAKE_REQUIRED_LIBRARIES ${ODBC_LIBRARIES})

  add_ampl_test(odbc-test odbc-test.cc odbc.cc LIBS function ${ODBC_LIBRARIES})

  include(CheckCXXSourceRuns)
  check_cxx_source_runs(
    "#include \"tests/odbc.cc\"
     int main() {
       return odbc::Env().FindDriver(\"access\").empty() ? 1 : 0;
     }" HAVE_ACCESS_DRIVER)

  if (HAVE_ACCESS_DRIVER)
    cxx_test_with_flags(access-test ""
      "gtest;function;util;${ODBC_LIBRARIES}"
      ../test-main.cc access-test.cc odbc.cc)
  endif ()

  check_cxx_source_runs(
    "#include \"tests/odbc.cc\"
     int main() {
       return odbc::Env().FindDriver(\"excel\").empty() ? 1 : 0;
     }" HAVE_EXCEL_DRIVER)

  if (HAVE_EXCEL_DRIVER)
    add_ampl_test(excel-test excel-test.cc odbc.cc
      LIBS function ${ODBC_LIBRARIES})
  endif ()

  check_cxx_source_runs(
    "#include \"tests/odbc.cc\"
     int main() {
       return odbc::Env().FindDriver(\"*.csv\").empty() ? 1 : 0;
     }" HAVE_ODBC_TEXT_DRIVER)

  if (HAVE_ODBC_TEXT_DRIVER)
    add_ampl_test(csv-test csv-test.cc odbc.cc LIBS function ${ODBC_LIBRARIES})
  endif ()

  if (ENABLE_MYSQL_TESTS)
    check_cxx_source_runs(
      "#include \"tests/odbc.cc\"
      int main() {
        return odbc::Env().FindDriver(\"mysql\").empty() ? 1 : 0;
      }" HAVE_MYSQL_DRIVER)

    if (HAVE_MYSQL_DRIVER)
      if (WIN32)
        set(MYSQL_TEST_LIBRARIES Ws2_32.lib)
      endif ()
      add_ampl_test(mysql-test mysql-test.cc odbc.cc
        LIBS function ${ODBC_LIBRARIES} ${MYSQL_TEST_LIBRARIES})
    endif ()
  endif ()
endif ()

add_ampl_test(tableproxy-test tableproxy-test.cc LIBS function)
add_dependencies(tableproxy-test tableproxy)
