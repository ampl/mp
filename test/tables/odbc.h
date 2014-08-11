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

#ifndef TESTS_ODBC_H_
#define TESTS_ODBC_H_

#ifdef _WIN32
// sqltypes.h needs HWND on Windows.
# include <windows.h>
# undef max
#endif

#include <sqltypes.h>

#include <algorithm>
#include <string>
#include <vector>

namespace odbc {

// An ODBC environment.
class Env {
 private:
  SQLHENV env_;

  // Do not implement.
  Env(const Env &);
  Env &operator=(const Env &);

  friend class Connection;

  void Check(const char *func_name, SQLRETURN ret) const;

 public:
  Env();
  ~Env();

  std::string FindDriver(const char *name) const;
};

class Connection {
 private:
  SQLHDBC dbc_;

  // Do not implement.
  Connection(const Connection &);
  Connection &operator=(const Connection &);

  friend class Statement;

  void Check(const char *func_name, SQLRETURN ret, bool nothrow = false) const;

 public:
  Connection(const Env &env);
  ~Connection();

  void Connect(const char *connection_string);
};

class Statement {
 private:
  SQLHSTMT stmt_;

  // Do not implement.
  Statement(const Statement &);
  Statement &operator=(const Statement &);

  void Check(const char *func_name, SQLRETURN ret, bool nothrow = false) const;

 public:
  Statement(const Connection &c);
  ~Statement();

  void Execute(const char *sql_statement);
};
}

#endif  // TESTS_ODBC_H_
