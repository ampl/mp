/*-------------------------------------------------------------------------*/
/* AMPL solver utilities                                  Victor Zverovich */
/*                                                                         */
/* Name           : util.h                                                 */
/* Title          : AMPL solver utilities                                  */
/* By             : Victor Zverovich                                       */
/* Date           : May 2012                                               */
/*                                                                         */
/* Various utilities for AMPL solver drivers.                              */
/* The same_expr function is based on the one by Robert Fourer.            */
/*-------------------------------------------------------------------------*/

#ifndef AMPL_SOLVERS_UTIL_H
#define AMPL_SOLVERS_UTIL_H

#include <stdexcept>

// Operation types
enum {
  // Unary operation
  OPTYPE_UNARY = 1,

  // Binary operation
  OPTYPE_BINARY = 2,

  // Variable-argument function such as min or max
  OPTYPE_VARARG = 3,

  // Piecewise-linear term
  OPTYPE_PLTERM = 4,

  // The if-then-else expression
  OPTYPE_IF = 5,

  // The sum expression
  OPTYPE_SUM = 6,

  // Function call
  OPTYPE_FUNCALL = 7,

  // String
  OPTYPE_STRING = 8,

  // Number
  OPTYPE_NUMBER = 9,

  // Variable
  OPTYPE_VARIABLE = 10,

  // The count expression
  OPTYPE_COUNT = 11
};

struct expr;

const char *get_opname(int opcode);
bool same_expr(const expr *e1, const expr *e2);

// A general solver error.
class Error : public std::runtime_error {
 public:
  Error(const std::string &message) : std::runtime_error(message) {}
};

// An exception that is thrown when an ASL expression not supported
// by the solver is encountered.
class UnsupportedExprError : public Error {
 public:
  UnsupportedExprError(const char *expr) :
    Error(std::string("unsupported expression: ") + expr) {}
};

// An exception that is thrown when an incomplete constraint expression
// is encountered.
class IncompleteConstraintExprError : public Error {
 public:
  IncompleteConstraintExprError(const char *expr) :
    Error(std::string("incomplete constraint expression using ") + expr) {}
};

#endif  // AMPL_SOLVERS_UTIL_H
