/*
 AMPL solver utilities.

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

#ifndef SOLVERS_UTIL_UTIL_H_
#define SOLVERS_UTIL_UTIL_H_

#include <stdexcept>
#include <string>

struct expr;

namespace ampl {

/// Returns the name of an operation with the specified code.
const char *GetOpName(int opcode);

// Returns true if the expressions e1 and e2 are structurally equal.
bool Equal(const expr *e1, const expr *e2);

/// A general solver error.
class Error : public std::runtime_error {
 public:
  explicit Error(const std::string &message) : std::runtime_error(message) {}
};

/// An exception that is thrown when an ASL expression not supported
/// by the solver is encountered.
class UnsupportedExprError : public Error {
 public:
  explicit UnsupportedExprError(const char *expr) :
    Error(std::string("unsupported expression: ") + expr) {}
};

/// An exception that is thrown when an incomplete constraint expression
/// is encountered.
class IncompleteConstraintExprError : public Error {
 public:
  explicit IncompleteConstraintExprError(const char *expr) :
    Error(std::string("incomplete constraint expression using ") + expr) {}
};
}

#endif  // SOLVERS_UTIL_UTIL_H_
