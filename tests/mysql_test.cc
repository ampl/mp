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

// A socket environment.
class SocketEnv {
 public:
  SocketEnv() {
#ifdef WIN32
    WSADATA data = {};
    if (WSAStartup(MAKEWORD(2, 2), &data))
      throw std::runtime_error("WSAStartup failed");
#endif
  }

  ~SocketEnv() {
#ifdef WIN32
    WSACleanup();
#endif
  }

  void GetHostname(char *name, int len) const {
    if (gethostname(name, len))
      throw std::runtime_error("gethostname failed");
  }
};

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
  SocketEnv().GetHostname(hostname, BUFFER_SIZE);
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
  Table t("", 1);
  t = "VERSION()";
  lib_.GetHandler("odbc")->Read(connection_, &t, "SQL=SELECT VERSION();");
  EXPECT_EQ(1, t.num_rows());
  EXPECT_TRUE(t.GetString(0, 0) != nullptr);
}

TEST_F(MySQLTest, Write) {
  Table t1(table_name_, 1);
  t1 = "Character Name",
       "Arthur Dent",
       "Ford Prefect";
  lib_.GetHandler("odbc")->Write(connection_, t1);
  Table t2(table_name_, 1);
  t2 = "Character Name";
  lib_.GetHandler("odbc")->Read(connection_, &t2);
  ASSERT_EQ(2, t2.num_rows());
  EXPECT_STREQ("Arthur Dent", t2.GetString(0, 0));
  EXPECT_STREQ("Ford Prefect", t2.GetString(1, 0));
}

TEST_F(MySQLTest, Rewrite) {
  Table t1(table_name_, 1);
  t1 = "Test",
       "foo";
  // The first write creates a table.
  lib_.GetHandler("odbc")->Write(connection_, t1);
  {
    Table t(table_name_, 1);
    t = "Test";
    lib_.GetHandler("odbc")->Read(connection_, &t);
    ASSERT_EQ(1, t.num_rows());
    EXPECT_STREQ("foo", t.GetString(0, 0));
  }
  // The second write should drop the table and create a new one.
  Table t2(table_name_, 1);
  t2 = "Character",
       "Zaphod";
  lib_.GetHandler("odbc")->Write(connection_, t2);
  {
    Table t(table_name_, 1);
    t = "Character";
    lib_.GetHandler("odbc")->Read(connection_, &t);
    ASSERT_EQ(1, t.num_rows());
    EXPECT_STREQ("Zaphod", t.GetString(0, 0));
  }
}

// TODO(viz): more tests
}
