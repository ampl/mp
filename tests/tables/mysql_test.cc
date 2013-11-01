/*
 Tests of the MySQL ODBC connection.

 Copyright (C) 2012 AMPL Optimization Inc

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization Inc disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Victor Zverovich
 */

#include <cctype>
#include <climits>

#include "gtest/gtest.h"
#include "tests/config.h"
#include "tests/function.h"
#include "tests/tables/odbc.h"
#include "solvers/util/format.h"
#include "solvers/funcadd.h"

#ifdef _WIN32
# include <process.h>
# define getpid _getpid
#else
# include <sys/types.h>
# include <unistd.h>
#endif

using fun::Handler;
using fun::Table;

#define SERVER "callisto.local"

namespace {

// A socket environment.
class SocketEnv {
 public:
  SocketEnv() {
#ifdef _WIN32
    WSADATA data = {};
    if (WSAStartup(MAKEWORD(2, 2), &data))
      throw std::runtime_error("WSAStartup failed");
#endif
  }

  ~SocketEnv() {
#ifdef _WIN32
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
  const Handler *handler_;
  odbc::Env env_;
  std::vector<std::string> strings_;
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
  strings_.push_back("ODBC");
  strings_.push_back("DRIVER={" + env_.FindDriver("mysql") +
      "}; SERVER=" SERVER "; DATABASE=test; USER=test;");

  // Create a unique table name from the hostname and pid. This is necessary
  // to avoid clashes between tests running in parallel on different machines
  // and accessing the same database server.
  char hostname[BUFFER_SIZE] = "";
  SocketEnv().GetHostname(hostname, BUFFER_SIZE);
  // The table name contains space to check quotation.
  table_name_ = str(fmt::Format("{} {}") << hostname << getpid());
}

void MySQLTest::TearDown() {
  // Drop the table.
  odbc::Connection con(env_);
  con.Connect(strings_[1].c_str());
  odbc::Statement stmt(con);
  try {
    stmt.Execute(c_str(fmt::Format("DROP TABLE `{}`") << table_name_));
  } catch (const std::exception &) {}  // Ignore errors.
}

fun::Library MySQLTest::lib_("../../tables/ampltabl.dll");

TEST_F(MySQLTest, Read) {
  Table t("", 0, 1, strings_);
  t.AddString("SQL=SELECT VERSION();");
  t = "VERSION()";
  handler_->Read(&t);
  EXPECT_EQ(1u, t.num_rows());
  EXPECT_TRUE(t(0, 0).string() != nullptr);
}

TEST_F(MySQLTest, Write) {
  Table t1(table_name_, 1, 0, strings_);
  t1 = "Character Name",
       "Arthur Dent",
       "Ford Prefect";
  handler_->Write(t1);
  Table t2(table_name_, 1, 0, strings_);
  t2 = "Character Name";
  handler_->Read(&t2);
  EXPECT_EQ(t1, t2);
}

TEST_F(MySQLTest, Rewrite) {
  Table t1(table_name_, 1, 0, strings_);
  t1 = "Test",
       "foo";
  // The first write creates a table.
  handler_->Write(t1);
  Table t2(table_name_, 1, 0, strings_);
  t2 = "Test";
  handler_->Read(&t2);
  ASSERT_EQ(t1, t2);
  // The second write should drop the table and create a new one.
  Table t3(table_name_, 1, 0, strings_);
  t3 = "Character",
       "Zaphod";
  handler_->Write(t3);
  Table t4(table_name_, 1, 0, strings_);
  t4 = "Character";
  handler_->Read(&t4);
  ASSERT_EQ(t3, t4);
}

TEST_F(MySQLTest, WriteInOut) {
  Table t1(table_name_, 1, 0, strings_);
  t1 = "Name",
       "Beeblebrox";
  // The first write creates a table.
  handler_->Write(t1);
  Table t2(table_name_, 1, 0, strings_);
  t2 = "Name";
  handler_->Read(&t2);
  ASSERT_EQ(t1, t2);
  // The second write appends data to the table.
  Table t3(table_name_, 1, 0, strings_);
  t3 = "Name",
       "Zaphod";
  handler_->Write(t3, Handler::INOUT);
  Table t4(table_name_, 1, 0, strings_);
  t4 = "Name";
  handler_->Read(&t4);
  Table t5(table_name_, 1, 0, strings_);
  t5= "Name",
      "Zaphod",
      "Beeblebrox";
  ASSERT_EQ(t5, t4);
}

TEST_F(MySQLTest, Append) {
  Table t1(table_name_, 1, 0, strings_);
  t1 = "Name",
       "Zaphod";
  // The first write creates a table.
  handler_->Write(t1);
  Table t2(table_name_, 1, 0, strings_);
  t2 = "Name";
  handler_->Read(&t2);
  ASSERT_EQ(t1, t2);
  // The second write appends data to the table.
  Table t3(table_name_, 1, 0, strings_);
  t3.AddString("write=append");
  t3 = "Name",
       "Beeblebrox";
  handler_->Write(t3, Handler::INOUT);
  Table t4(table_name_, 1, 0, strings_);
  t4 = "Name";
  handler_->Read(&t4);
  Table t5(table_name_, 1, 0, strings_);
  t5= "Name",
      "Zaphod",
      "Beeblebrox";
  ASSERT_EQ(t5, t4);
}

TEST_F(MySQLTest, AdjustColNames) {
  Table t(table_name_, 1, 2, strings_);
  t = "Time:a",         "Strcol:b", "Mixed:c",
      20121112143000.0, "e",        "f",
      20121112150000.0, "111",      "222";
  handler_->Write(t);
  Table t2(table_name_, 1, 2, strings_);
  t2 = "a", "b", "c";
  handler_->Read(&t2);
  Table t3(table_name_, 1, 2, strings_);
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
  Table t(table_name_, 1, 1, strings_);
  t = "a", "";
  EXPECT_ERROR(handler_->Write(t),
      "Column 2's name is the empty string.");
  EXPECT_ERROR(handler_->Read(&t),
      "Column 2's name is the empty string.");
}

TEST_F(MySQLTest, QuoteInTableName) {
  table_name_ += '`';
  Table t(table_name_, 1, 0, strings_);
  t = "c", "v";
  handler_->Write(t);
}

TEST_F(MySQLTest, QuoteInColumnName) {
  Table t(table_name_, 1, 0, strings_);
  t = "c`", "v";
  handler_->Write(t);
}

TEST_F(MySQLTest, InvalidCharInTableName) {
  table_name_ += '\t';
  Table t(table_name_, 1, 0, strings_);
  t = "c", "v";
  std::string error = str(fmt::Format(
      "Table name contains invalid character with code {}")
      << static_cast<int>('\t'));
  EXPECT_ERROR(handler_->Write(t), error.c_str());
  EXPECT_ERROR(handler_->Read(&t), error.c_str());
}

TEST_F(MySQLTest, LowerCaseLettersInTableName) {
  for (int c = 1; c <= UCHAR_MAX; ++c) {
    if (std::isalpha(c) && std::tolower(c) == c)
      table_name_ += static_cast<char>(c);
  }
  Table t(table_name_, 1, 0, strings_);
  t = "c", "v";
  handler_->Write(t);
  Table in(table_name_, 1, 0, strings_);
  in = "c";
  handler_->Read(&in);
  EXPECT_EQ(t, in);
}

TEST_F(MySQLTest, UpperCaseLettersInTableName) {
  for (int c = 1; c <= UCHAR_MAX; ++c) {
    if (std::isalpha(c) && std::tolower(c) != c)
      table_name_ += static_cast<char>(c);
  }
  Table t(table_name_, 1, 0, strings_);
  t = "c", "v";
  handler_->Write(t);
  Table in(table_name_, 1, 0, strings_);
  in = "c";
  handler_->Read(&in);
  EXPECT_EQ(t, in);
}

TEST_F(MySQLTest, DigitsInTableName) {
  for (int c = 1; c <= UCHAR_MAX; ++c) {
    if (std::isdigit(c))
      table_name_ += static_cast<char>(c);
  }
  Table t(table_name_, 1, 0, strings_);
  t = "c", "v";
  handler_->Write(t);
  Table in(table_name_, 1, 0, strings_);
  in = "c";
  handler_->Read(&in);
  EXPECT_EQ(t, in);
}

TEST_F(MySQLTest, SpecialCharsInTableName) {
  for (int c = 1; c <= UCHAR_MAX; ++c) {
    if (std::isprint(c) && !std::isalnum(c))
      table_name_ += static_cast<char>(c);
  }
  Table t(table_name_, 1, 0, strings_);
  t = "c", "v";
  handler_->Write(t);
  Table in(table_name_, 1, 0, strings_);
  in = "c";
  handler_->Read(&in);
  EXPECT_EQ(t, in);
}

TEST_F(MySQLTest, InvalidCharsInColumnName) {
  for (int c = 1; c <= UCHAR_MAX; ++c) {
    if (std::isprint(c))
      continue;
    char col_name[2] = {static_cast<char>(c)};
    Table t(table_name_, 1, 0, strings_);
    t = col_name, "v";
    std::string error = str(fmt::Format(
          "Column 1's name contains invalid character with code {}")
          << static_cast<int>(c));
    EXPECT_ERROR(handler_->Write(t), error.c_str());
    EXPECT_ERROR(handler_->Read(&t), error.c_str());
  }
}

TEST_F(MySQLTest, AlphaNumericColumnName) {
  std::string col_name;
  for (int c = 1; c <= UCHAR_MAX; ++c) {
    if (std::isalnum(c))
      col_name += static_cast<char>(c);
  }
  Table t(table_name_, 1, 0, strings_);
  t = col_name.c_str(), "v";
  handler_->Write(t);
  Table in(table_name_, 1, 0, strings_);
  in = col_name.c_str();
  handler_->Read(&in);
  EXPECT_EQ(t, in);
}

TEST_F(MySQLTest, SpecialCharsInColumnName) {
  std::string col_name;
  for (int c = 1; c <= UCHAR_MAX; ++c) {
    if (std::isprint(c) && !std::isalnum(c))
      col_name += static_cast<char>(c);
  }
  Table t(table_name_, 1, 0, strings_);
  t = col_name.c_str(), "v";
  handler_->Write(t);
  Table in(table_name_, 1, 0, strings_);
  in = col_name.c_str();
  handler_->Read(&in);
  EXPECT_EQ(t, in);
}
}
