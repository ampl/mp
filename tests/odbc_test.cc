/*
 Tests of the AMPL ODBC table handler.

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

#include <cstring>
#include <algorithm>
#include <deque>
#include <vector>

#include <sql.h>
#include <sqlext.h>

#include "gtest/gtest.h"
#include "tests/function.h"
#include "solvers/asl.h"
#include "tests/config.h"

using fun::Table;

namespace {

// An ODBC environment.
class Env {
 private:
  SQLHENV env_;

  bool GetDiag(SQLSMALLINT rec_number, SQLCHAR *sql_state,
      SQLINTEGER &native_error, std::vector<SQLCHAR> &message,
      std::ostream &os) const;

  bool Check(const char *func_name, SQLRETURN ret) const;

 public:
  Env();
  ~Env();

  bool FindDriver(const char *name) const;
};

bool Env::GetDiag(SQLSMALLINT rec_number, SQLCHAR *sql_state,
    SQLINTEGER &native_error, std::vector<SQLCHAR> &message,
    std::ostream &os) const {
  SQLSMALLINT length = 0;
  SQLRETURN ret = SQLGetDiagRec(SQL_HANDLE_ENV,
          env_, rec_number, sql_state, &native_error, 0, 0, &length);
  if (ret == SQL_NO_DATA)
    return false;
  if (ret != SQL_SUCCESS && ret != SQL_SUCCESS_WITH_INFO) {
    os << "SQLGetDiagRec returned error code " << ret << "\n";
    return false;
  }
  message.resize(length + 1);
  ret = SQLGetDiagRec(SQL_HANDLE_ENV, env_, rec_number,
          sql_state, &native_error, &message[0], length + 1, &length);
  if (ret != SQL_SUCCESS && ret != SQL_SUCCESS_WITH_INFO) {
    os << "SQLGetDiagRec returned error code " << ret << "\n";
    return false;
  }
  return true;
}

bool Env::Check(const char *func_name, SQLRETURN ret) const {
  if (ret == SQL_SUCCESS)
    return true;
  std::ostringstream os;
  os << func_name << " returned error code " << ret << "\n";
  SQLCHAR sql_state[6] = "";
  SQLINTEGER native_error = 0;
  std::vector<SQLCHAR> message;
  for (SQLSMALLINT i = 1;
      GetDiag(i, sql_state, native_error, message, os); ++i) {
    os << "SQLState: " << sql_state << "\n";
    os << "Native Error: " << native_error << "\n";
    os << "Message: " << &message[0] << "\n";
  }
  if (ret != SQL_SUCCESS_WITH_INFO)
    throw std::runtime_error(os.str());
  std::cout << os.str();
  return true;
}

Env::Env() : env_(SQL_NULL_HANDLE) {
  Check("SQLAllocHandle", SQLAllocHandle(
      SQL_HANDLE_ENV, SQL_NULL_HANDLE, &env_));
  Check("SQLSetEnvAttr", SQLSetEnvAttr(env_, SQL_ATTR_ODBC_VERSION,
      reinterpret_cast<SQLPOINTER>(SQL_OV_ODBC3), 0));
}

Env::~Env() {
  try {
    Check("SQLFreeHandle", SQLFreeHandle(SQL_HANDLE_ENV, env_));
  } catch (const std::exception &e) {
    std::cout << e.what();
  }
}

bool Env::FindDriver(const char *name) const {
  for (SQLUSMALLINT direction = SQL_FETCH_FIRST;; direction = SQL_FETCH_NEXT) {
    const int BUFFER_SIZE = 256;
    SQLCHAR desc[BUFFER_SIZE] = "", attr[BUFFER_SIZE] = "";
    SQLSMALLINT desc_length = 0, attr_length = 0;
    SQLRETURN ret = SQLDrivers(env_, direction,
        desc, BUFFER_SIZE, &desc_length, attr, BUFFER_SIZE, &attr_length);
    if (ret == SQL_NO_DATA) break;
    Check("SQLDrivers", ret);
    SQLCHAR *desc_end = desc + std::min<int>(BUFFER_SIZE, desc_length);
    std::transform(desc, desc_end, desc, std::ptr_fun<int, int>(std::tolower));
    return std::search(desc, desc_end, name, name + std::strlen(name)) !=
        desc_end;
  }
  return true;
}

class ODBCTest : public ::testing::Test {
 protected:
  static fun::Library lib_;

  static void SetUpTestCase() {
    lib_.Load();
  }

  void Read(Table &t) {
    lib_.GetHandler("odbc")->Read(t);
  }
};

fun::Library ODBCTest::lib_("../tables/ampltabl.dll");

TEST_F(ODBCTest, Loaded) {
  EXPECT_EQ("", lib_.error());
  EXPECT_TRUE(lib_.GetHandler("odbc"));
}

TEST_F(ODBCTest, ReadMySQL) {
  if (!Env().FindDriver("mysql")) {
    std::cerr << "Skipping MySQL test." << std::endl;
    return;
  }
  Table t("", "ODBC",
      "DRIVER={MySQL ODBC 5.2 Driver}; SOCKET={/var/run/mysqld/mysqld.sock};",
      "SQL=SELECT VERSION();");
  t.AddCol("VERSION()");
  Read(t);
  EXPECT_EQ(1, t.num_rows());
  EXPECT_TRUE(t.GetString(0) != nullptr);
}

// TODO: more tests
}
