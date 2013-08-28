/*
 A simple C++ interface to ODBC.

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

#include "tests/odbc.h"

#include <sql.h>
#include <sqlext.h>

#include <cctype>
#include <cstring>

#include <algorithm>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace {

bool GetDiag(SQLSMALLINT handle_type, SQLHANDLE handle,
    SQLSMALLINT rec_number, SQLCHAR *sql_state, SQLINTEGER &native_error,
    std::vector<SQLCHAR> &message, std::ostream &os) {
  SQLSMALLINT length = 0;
  SQLRETURN ret = SQLGetDiagRec(handle_type,
          handle, rec_number, sql_state, &native_error, 0, 0, &length);
  if (ret == SQL_NO_DATA)
    return false;
  if (ret != SQL_SUCCESS && ret != SQL_SUCCESS_WITH_INFO) {
    os << "SQLGetDiagRec returned error code " << ret << "\n";
    return false;
  }
  message.resize(length + 1);
  ret = SQLGetDiagRec(handle_type, handle, rec_number,
          sql_state, &native_error, &message[0], length + 1, &length);
  if (ret != SQL_SUCCESS && ret != SQL_SUCCESS_WITH_INFO) {
    os << "SQLGetDiagRec returned error code " << ret << "\n";
    return false;
  }
  return true;
}

void Check(SQLSMALLINT handle_type, SQLHANDLE handle,
    const char *func_name, SQLRETURN ret, bool nothrow = false) {
  if (ret == SQL_SUCCESS)
    return;
  std::ostringstream os;
  os << func_name << " returned error code " << ret << "\n";
  SQLCHAR sql_state[6] = "";
  SQLINTEGER native_error = 0;
  std::vector<SQLCHAR> message;
  for (SQLSMALLINT i = 1;
      GetDiag(handle_type, handle, i, sql_state, native_error, message, os);
      ++i) {
    os << "SQLState: " << sql_state << "\n";
    os << "Native Error: " << native_error << "\n";
    os << "Message: " << &message[0] << "\n";
  }
  if (ret != SQL_SUCCESS_WITH_INFO && !nothrow)
    throw std::runtime_error(os.str());
  std::cout << os.str();
}

void FreeHandle(SQLSMALLINT handle_type, SQLHANDLE handle) {  // throw()
  Check(handle_type, handle, "SQLFreeHandle",
      SQLFreeHandle(handle_type, handle), true);
}

SQLCHAR *ConvertString(const char *s) {
  return reinterpret_cast<SQLCHAR*>(const_cast<char*>(s));
}
}

namespace odbc {

void Env::Check(const char *func_name, SQLRETURN ret) const {
  ::Check(SQL_HANDLE_ENV, env_, func_name, ret);
}

Env::Env() : env_(SQL_NULL_HANDLE) {
  Check("SQLAllocHandle", SQLAllocHandle(
      SQL_HANDLE_ENV, SQL_NULL_HANDLE, &env_));
  Check("SQLSetEnvAttr", SQLSetEnvAttr(env_, SQL_ATTR_ODBC_VERSION,
      reinterpret_cast<SQLPOINTER>(SQL_OV_ODBC3), 0));
}

Env::~Env() {
  FreeHandle(SQL_HANDLE_ENV, env_);
}

std::string Env::FindDriver(const char *name) const {
  for (SQLUSMALLINT direction = SQL_FETCH_FIRST;; direction = SQL_FETCH_NEXT) {
    const int BUFFER_SIZE = 256;
    SQLCHAR desc[BUFFER_SIZE] = "", attr[BUFFER_SIZE] = "";
    SQLSMALLINT desc_length = 0, attr_length = 0;
    SQLRETURN ret = SQLDrivers(env_, direction,
        desc, BUFFER_SIZE, &desc_length, attr, BUFFER_SIZE, &attr_length);
    if (ret == SQL_NO_DATA) break;
    Check("SQLDrivers", ret);
    std::string driver_name(reinterpret_cast<char*>(desc));
    SQLCHAR *desc_end = desc + std::min<int>(BUFFER_SIZE, desc_length);
    std::transform(desc, desc_end, desc, std::ptr_fun<int, int>(std::tolower));
    if (std::search(desc, desc_end, name, name + std::strlen(name)) != desc_end)
        return driver_name;
  }
  return std::string();
}

void Connection::Check(
    const char *func_name, SQLRETURN ret, bool nothrow) const {
  ::Check(SQL_HANDLE_DBC, dbc_, func_name, ret, nothrow);
}

Connection::Connection(const Env &env) : dbc_(SQL_NULL_HANDLE) {
  env.Check("SQLAllocHandle", SQLAllocHandle(
      SQL_HANDLE_DBC, env.env_, &dbc_));
}

Connection::~Connection() {
  Check("SQLDisconnect", SQLDisconnect(dbc_), true);
  FreeHandle(SQL_HANDLE_DBC, dbc_);
}

void Connection::Connect(const char *connection_string) {
  SQLSMALLINT length = 0;
  Check("SQLDriverConnect", SQLDriverConnect(dbc_, 0,
      ConvertString(connection_string), SQL_NTS,
      0, 0, &length, SQL_DRIVER_NOPROMPT));
}

void Statement::Check(
    const char *func_name, SQLRETURN ret, bool nothrow) const {
  ::Check(SQL_HANDLE_STMT, stmt_, func_name, ret, nothrow);
}

Statement::Statement(const Connection &c) : stmt_(SQL_NULL_HANDLE) {
  c.Check("SQLAllocHandle", SQLAllocHandle(
      SQL_HANDLE_STMT, c.dbc_, &stmt_));
}

Statement::~Statement() {
  FreeHandle(SQL_HANDLE_STMT, stmt_);
}

void Statement::Execute(const char *sql_statement) {
  Check("SQLExecDirect", SQLExecDirect(
      stmt_, ConvertString(sql_statement), SQL_NTS));
}
}
