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
#include "solvers/funcadd.h"
#undef snprintf

#ifdef _WIN32
# include <process.h>
# define getpid _getpid
# define snprintf _snprintf
#else
# include <sys/types.h>
# include <unistd.h>
#endif

using fun::Table;

#define SERVER "callisto.local"

namespace {

class MySQLTest : public ::testing::Test {
 protected:
  static fun::Library lib_;
  odbc::Env env_;
  std::string connection_;
  std::string table_name_;
  enum {BUFFER_SIZE = 256};

  static void SetUpTestCase() {
    lib_.Load();
  }

  void SetUp();
  void TearDown();
};

void MySQLTest::SetUp() {
  connection_ = "DRIVER={" + env_.FindDriver("mysql") +
      "}; SERVER=" SERVER "; DATABASE=test;";

  // Create a unique table name from the hostname and pid. This is necessary
  // to avoid clashes between tests running in parallel on different machines
  // and accessing the same database server.
  char hostname[BUFFER_SIZE] = "";
  ASSERT_EQ(0, gethostname(hostname, BUFFER_SIZE));
  int pid = getpid();
  char table_name[BUFFER_SIZE] = "";
  // The table name contains space to check quotation.
  snprintf(table_name, BUFFER_SIZE, "%s %d", hostname, pid);
  table_name_ = table_name;
}

void MySQLTest::TearDown() {
  // Drop the table.
  odbc::Connection con(env_);
  con.Connect(connection_.c_str());
  char sql[BUFFER_SIZE];
  snprintf(sql, BUFFER_SIZE, "DROP TABLE `%s`", table_name_.c_str());
  odbc::Statement stmt(con);
  try {
    stmt.Execute(sql);
  } catch (const std::exception &) {}  // Ignore errors.
}

fun::Library MySQLTest::lib_("../tables/ampltabl.dll");

TEST_F(MySQLTest, Read) {
  Table t("", 1, 1, "ODBC", connection_.c_str(), "SQL=SELECT VERSION();");
  t.SetColName(0, "VERSION()");
  EXPECT_EQ(DB_Done, lib_.GetHandler("odbc")->Read(&t));
  EXPECT_EQ(nullptr, t.error_message());
  EXPECT_EQ(1, t.num_rows());
  EXPECT_TRUE(t.GetString(0, 0) != nullptr);
}

TEST_F(MySQLTest, Write) {
  {
    Table t(table_name_.c_str(), 2, 1, "ODBC", connection_.c_str());
    t.SetColName(0, "Character Name");
    t.SetString(0, "Arthur Dent");
    t.SetString(1, "Ford Prefect");
    EXPECT_EQ(DB_Done, lib_.GetHandler("odbc")->Write(&t));
    EXPECT_STREQ(nullptr, t.error_message());
  }
  {
    Table t(table_name_.c_str(), 2, 1, "ODBC", connection_.c_str());
    t.SetColName(0, "Character Name");
    EXPECT_EQ(DB_Done, lib_.GetHandler("odbc")->Read(&t));
    EXPECT_STREQ(nullptr, t.error_message());
    ASSERT_EQ(2, t.num_rows());
    EXPECT_STREQ("Arthur Dent", t.GetString(0, 0));
    EXPECT_STREQ("Ford Prefect", t.GetString(1, 0));
  }
}

TEST_F(MySQLTest, Rewrite) {
  {
    // The first write creates a table.
    Table t(table_name_.c_str(), 1, 1, "ODBC", connection_.c_str());
    t.SetColName(0, "Test");
    t.SetString(0, "foo");
    EXPECT_EQ(DB_Done, lib_.GetHandler("odbc")->Write(&t));
    EXPECT_STREQ(nullptr, t.error_message());
  }
  {
    Table t(table_name_.c_str(), 1, 1, "ODBC", connection_.c_str());
    t.SetColName(0, "Test");
    EXPECT_EQ(DB_Done, lib_.GetHandler("odbc")->Read(&t));
    EXPECT_STREQ(nullptr, t.error_message());
    ASSERT_EQ(1, t.num_rows());
    EXPECT_STREQ("foo", t.GetString(0, 0));
  }
  {
    // The second write should drop the table and create a new one.
    Table t(table_name_.c_str(), 1, 1, "ODBC", connection_.c_str());
    t.SetColName(0, "Character");
    t.SetString(0, "Zaphod");
    EXPECT_EQ(DB_Done, lib_.GetHandler("odbc")->Write(&t));
    EXPECT_STREQ(nullptr, t.error_message());
  }
  {
    Table t(table_name_.c_str(), 1, 1, "ODBC", connection_.c_str());
    t.SetColName(0, "Character");
    EXPECT_EQ(DB_Done, lib_.GetHandler("odbc")->Read(&t));
    EXPECT_STREQ(nullptr, t.error_message());
    ASSERT_EQ(1, t.num_rows());
    EXPECT_STREQ("Zaphod", t.GetString(0, 0));
  }
}

// TODO(viz): more tests
}
