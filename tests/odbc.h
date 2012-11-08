/*
 A simple C++ interface to ODBC.

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

#ifndef TESTS_ODBC_H_
#define TESTS_ODBC_H_

#include <sqltypes.h>

#include <ostream>
#include <string>
#include <vector>

namespace odbc {

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

  std::string FindDriver(const char *name) const;
};
}

#endif  // TESTS_ODBC_H_
