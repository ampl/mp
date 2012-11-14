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

#include <climits>

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
  const fun::Handler *handler_;
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
  handler_ = lib_.GetHandler("odbc");
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
  handler_->Read(connection_, &t, "SQL=SELECT VERSION();");
  EXPECT_EQ(1u, t.num_rows());
  EXPECT_TRUE(t(0, 0).string() != nullptr);
}

TEST_F(MySQLTest, Write) {
  Table t1(table_name_, 1);
  t1 = "Character Name",
       "Arthur Dent",
       "Ford Prefect";
  handler_->Write(connection_, t1);
  Table t2(table_name_, 1);
  t2 = "Character Name";
  handler_->Read(connection_, &t2);
  EXPECT_EQ(t2, t1);
}

TEST_F(MySQLTest, Rewrite) {
  Table t1(table_name_, 1);
  t1 = "Test",
       "foo";
  // The first write creates a table.
  handler_->Write(connection_, t1);
  Table t2(table_name_, 1);
  t2 = "Test";
  handler_->Read(connection_, &t2);
  ASSERT_EQ(t1, t2);
  // The second write should drop the table and create a new one.
  Table t3(table_name_, 1);
  t3 = "Character",
       "Zaphod";
  handler_->Write(connection_, t3);
  Table t4(table_name_, 1);
  t4 = "Character";
  handler_->Read(connection_, &t4);
  ASSERT_EQ(t3, t4);
}

TEST_F(MySQLTest, Append) {
  Table t1(table_name_, 1);
  t1 = "Test",
       "Beeblebrox";
  // The first write creates a table.
  handler_->Write(connection_, t1);
  Table t2(table_name_, 1);
  t2 = "Test";
  handler_->Read(connection_, &t2);
  ASSERT_EQ(t2, t1);
  // The second write appends to the table.
  Table t3(table_name_, 1);
  t3 = "Test",
       "Zaphod";
  handler_->Write(connection_, t3, true);
  Table t4(table_name_, 1);
  t4 = "Test";
  handler_->Read(connection_, &t4);
  Table t5(table_name_, 1);
  t5= "Test",
      "Zaphod",
      "Beeblebrox";
  ASSERT_EQ(t5, t4);
}

TEST_F(MySQLTest, AdjustColNames) {
  Table t(table_name_, 3);
  t = "Time:a",         "Strcol:b", "Mixed:c",
      20121112143000.0, "e",        "f",
      20121112150000.0, "111",      "222";
  handler_->Write(connection_, t);
  Table t2(table_name_, 3);
  t2 = "a", "b", "c";
  handler_->Read(connection_, &t2);
  Table t3(table_name_, 3);
  t3 = "a",             "b", "c",
      20121112143000.0, "e", "f",
      20121112150000.0, 111, 222;
  EXPECT_EQ(t3, t2);
}

#define EXPECT_ERROR(statement, message) \
do { \
  bool catched = false; \
  try { \
    statement; \
  } catch (const std::runtime_error &e) { \
    EXPECT_STREQ(message, e.what()); \
    catched = true; \
  } \
  EXPECT_TRUE(catched); \
} while (false)

TEST_F(MySQLTest, EmptyColName) {
  Table t(table_name_, 2);
  t = "a", "";
  EXPECT_ERROR(handler_->Write(connection_, t),
      "Column 2's name is the empty string.");
}

TEST_F(MySQLTest, QuoteInTableName) {
  table_name_ += '`';
  Table t(table_name_, 1);
  t = "c", "v";
  handler_->Write(connection_, t);
}

TEST_F(MySQLTest, QuoteInColumnName) {
  Table t(table_name_, 1);
  t = "c`", "v";
  handler_->Write(connection_, t);
}

TEST_F(MySQLTest, SpecialCharInColumnName) {
  for (unsigned c = 1; c <= UCHAR_MAX; ++c) {
    if (c == ' ')
      continue; // MySQL doesn't allow trailing spaces in column names.
    char col_name[2] = {static_cast<char>(c)};
    Table t(table_name_, 1);
    t = col_name, "v";
    if (!std::isprint(c)) {
      char buffer[BUFFER_SIZE] = "";
      snprintf(buffer, BUFFER_SIZE, "Name \"\\x%02x\" contains invalid "
          "character with code %d ('\\x%02x')", c, c, c);
      EXPECT_ERROR(handler_->Write(connection_, t), buffer);
      continue;
    }
    try {
      handler_->Write(connection_, t);
    } catch (...) {
      std::cout << "Failed on character " << c << std::endl;
      throw;
    }
    Table in(table_name_, 1);
    in = col_name;
    handler_->Read(connection_, &in);
    EXPECT_EQ(t, in);
  }
}

// TODO(viz): more tests
}
