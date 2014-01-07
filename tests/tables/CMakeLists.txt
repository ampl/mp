find_package(ODBC)
if (ODBC_FOUND)
  include_directories(${ODBC_INCLUDE_DIRS})

  set(CMAKE_REQUIRED_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/..)
  set(CMAKE_REQUIRED_LIBRARIES ${ODBC_LIBRARIES})

  cxx_test_with_flags(odbc-test ""
    "gtest;function;${ODBC_LIBRARIES}" ../test-main.cc odbc-test.cc odbc.cc)

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
    cxx_test_with_flags(excel-test ""
      "gtest;function;util;${ODBC_LIBRARIES}"
      ../test-main.cc excel-test.cc odbc.cc)
  endif ()

  check_cxx_source_runs(
    "#include \"tests/odbc.cc\"
     int main() {
       return odbc::Env().FindDriver(\"*.csv\").empty() ? 1 : 0;
     }" HAVE_ODBC_TEXT_DRIVER)

  if (HAVE_ODBC_TEXT_DRIVER)
    cxx_test_with_flags(csv-test ""
      "gtest;function;util;${ODBC_LIBRARIES}"
      ../test-main.cc csv-test.cc odbc.cc)
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
      cxx_test_with_flags(mysql-test ""
        "gtest;function;util;${ODBC_LIBRARIES};${MYSQL_TEST_LIBRARIES}"
        ../test-main.cc mysql-test.cc odbc.cc)
    endif ()
  endif ()
endif ()

cxx_test_with_flags(tableproxy-test ""
  "gtest;function;util" ../test-main.cc tableproxy-test.cc ../util.cc)
add_dependencies(tableproxy-test tableproxy)
