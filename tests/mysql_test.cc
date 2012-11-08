/*
 Tests of the MySQL ODBC connection.

 Copyright (C) 2012 AMPL Optimization LLC

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization LLC disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Victor Zverovich
 */

#include "gtest/gtest.h"
#include "tests/config.h"
#include "tests/function.h"
#include "tests/odbc.h"

#ifdef _WIN32
# include <winsock2.h>
# include <process.h>
#define getpid _getpid
#else
# include <sys/types.h>
# include <unistd.h>
#endif

#include "solvers/funcadd.h"

#undef snprintf

using fun::Table;

#define SERVER "callisto.local"

namespace {

class MySQLTest : public ::testing::Test {
 protected:
  static fun::Library lib_;
  odbc::Env env_;
  std::string connection_;

  static void SetUpTestCase() {
    lib_.Load();
  }

  void SetUp() {
    connection_ = "DRIVER={" + env_.FindDriver("mysql") +
        "}; SERVER=" SERVER "; DATABASE=test;";
  }
};

fun::Library MySQLTest::lib_("../tables/ampltabl.dll");

TEST_F(MySQLTest, Read) {
  Table t("", "ODBC", connection_.c_str(), "SQL=SELECT VERSION();");
  t.AddCol("VERSION()");
  EXPECT_EQ(DB_Done, lib_.GetHandler("odbc")->Read(&t));
  EXPECT_EQ(nullptr, t.error_message());
  EXPECT_EQ(1, t.num_rows());
  EXPECT_TRUE(t.GetString(0) != nullptr);
}

TEST_F(MySQLTest, Write) {
  // Create a unique table name from the hostname and pid. This is necessary
  // to avoid clashes between tests running in parallel on different machines
  // and accessing the same database server.
  const int BUFFER_SIZE = 256;
  char hostname[BUFFER_SIZE] = "";
  ASSERT_EQ(0, gethostname(hostname, BUFFER_SIZE));
  int pid = getpid();
  char table_name[BUFFER_SIZE] = "";
  // The table name contains space to check quotation.
  snprintf(table_name, BUFFER_SIZE, "%s %d", hostname, pid);

  Table t(table_name, "ODBC", connection_.c_str(), "verbose");
  t.AddCol("Test");
  EXPECT_EQ(DB_Done, lib_.GetHandler("odbc")->Write(&t));
  EXPECT_STREQ(nullptr, t.error_message());

  // Drop the table.
  odbc::Connection con(env_);
  con.Connect(connection_.c_str());
  char sql[BUFFER_SIZE];
  snprintf(sql, BUFFER_SIZE, "DROP TABLE `%s`", table_name);
  odbc::Statement stmt(con);
  stmt.Execute(sql);
}

// TODO(viz): more tests
}
