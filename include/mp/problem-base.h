/*
 Base types and functions related to the optimization problem.

 This header is used to decouple the .nl reader from a specific problem
 representation as the reader can be used to construct different types of
 problems. So instead of including problem.h and expr.h from nl.h, this
 header containing common definitions is included from problem.h, expr.h
 and nl.h.

 Copyright (C) 2014 AMPL Optimization Inc

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

#ifndef MP_PROBLEM_BASE_H_
#define MP_PROBLEM_BASE_H_

#include <cassert>
#include <cstddef>  // for std::size_t

namespace mp {

class Expr;

namespace expr {

// Expression kinds.
enum Kind {
  // An unknown expression.
  UNKNOWN = 0,

  FIRST_EXPR,

  // To simplify checks, numeric expression kinds are in a range
  // [FIRST_NUMERIC, LAST_NUMERIC].
  FIRST_NUMERIC = FIRST_EXPR,
  VARIABLE = FIRST_NUMERIC,

  // Unary expressions.
  FIRST_UNARY,
  FLOOR = FIRST_UNARY,
  CEIL,
  ABS,
  MINUS,
  TANH,
  TAN,
  SQRT,
  SINH,
  SIN,
  LOG10,
  LOG,
  EXP,
  COSH,
  COS,
  ATANH,
  ATAN,
  ASINH,
  ASIN,
  ACOSH,
  ACOS,
  POW2,
  LAST_UNARY = POW2,

  // Binary expressions.
  FIRST_BINARY,
  ADD = FIRST_BINARY,
  SUB,
  MUL,
  DIV,
  INT_DIV,
  MOD,
  POW,
  POW_CONST_BASE,
  POW_CONST_EXP,
  LESS,
  ATAN2,
  PRECISION,
  ROUND,
  TRUNC,
  LAST_BINARY = TRUNC,

  IF,
  IFSYM,
  PLTERM,
  CALL,

  // Variable argument expressions.
  FIRST_VARARG,
  MIN = FIRST_VARARG,
  MAX,
  LAST_VARARG = MAX,

  SUM,
  COUNT,
  NUMBEROF,
  NUMBEROF_SYM,
  LAST_NUMERIC,

  // CONSTANT belongs both to numeric and logical expressions therefore
  // the [FIRST_NUMERIC, LAST_NUMERIC] and [FIRST_LOGICAL, LAST_LOGICAL]
  // ranges overlap at CONSTANT = LAST_NUMERIC = FIRST_LOGICAL.
  CONSTANT = LAST_NUMERIC,

  // To simplify checks, logical expression kinds are in a range
  // [FIRST_LOGICAL, LAST_LOGICAL].
  FIRST_LOGICAL = CONSTANT,
  NOT,

  // Binary logical expressions.
  FIRST_BINARY_LOGICAL,
  OR = FIRST_BINARY_LOGICAL,
  AND,
  IFF,
  LAST_BINARY_LOGICAL = IFF,

  // Relational expressions.
  FIRST_RELATIONAL,
  LT = FIRST_RELATIONAL,  // <
  LE,                     // <=
  EQ,                     // =
  GE,                     // >=
  GT,                     // >
  NE,                     // !=
  LAST_RELATIONAL = NE,

  FIRST_LOGICAL_COUNT,
  ATLEAST = FIRST_LOGICAL_COUNT,
  ATMOST,
  EXACTLY,
  NOT_ATLEAST,
  NOT_ATMOST,
  NOT_EXACTLY,
  LAST_LOGICAL_COUNT = NOT_EXACTLY,

  IMPLICATION,

  FIRST_ITERATED_LOGICAL,
  FORALL = FIRST_ITERATED_LOGICAL,
  EXISTS,
  LAST_ITERATED_LOGICAL = EXISTS,

  ALLDIFF,
  LAST_LOGICAL = ALLDIFF,

  STRING,
  LAST_EXPR = STRING
};

// Maximum opcode index.
enum { MAX_OPCODE = 81 };

class OpCodeInfo {
 private:
  static const OpCodeInfo INFO[MAX_OPCODE + 1];

 public:
  expr::Kind kind;
  expr::Kind first_kind;  // First member of a kind.

  friend const OpCodeInfo &GetOpCodeInfo(int opcode);
};

inline const OpCodeInfo &GetOpCodeInfo(int opcode) {
  assert(opcode >= 0 && opcode <= MAX_OPCODE);
  return OpCodeInfo::INFO[opcode];
}

int opcode(expr::Kind kind);
}  // namespace expr

class Expr;

namespace internal {
// Expression information.
class ExprInfo {
 private:
  static const ExprInfo INFO[];

  friend class mp::Expr;
  friend int expr::opcode(expr::Kind kind);

 public:
  int opcode;
  int precedence;
  const char *str;
};
}

inline int expr::opcode(expr::Kind kind) {
  assert(kind >= expr::UNKNOWN && kind <= expr::LAST_EXPR);
  return internal::ExprInfo::INFO[kind].opcode;
}

namespace func {
// Function type.
enum Type {
  NUMERIC  = 0,
  SYMBOLIC = 1  // Accepts symbolic arguments.
};
}

namespace var {
// Variable type.
enum Type { CONTINUOUS, INTEGER };
}

namespace obj {
// Objective type.
enum Type { MIN = 0, MAX = 1 };
}

// Complementarity namespace. It would make more sense to call it compl,
// but the latter is a reserved word in C++.
namespace comp {
// Flags for complementarity constraints.
enum { INF_LB = 1, INF_UB = 2 };
}

namespace suf {
// Suffix kinds.
enum {
  VAR     =  0,  // Applies to variables.
  CON     =  1,  // Applies to constraints.
  OBJ     =  2,  // Applies to objectives.
  PROBLEM =  3,  // Applies to problems.
  MASK    =  3,  // Mask for the above.
  FLOAT   =  4,  // Suffix values are floating-point numbers.
  IODECL  =  8,  // Tell AMPL to make this an INOUT suffix.
  OUTPUT  = 16,  // Output suffix: return values to AMPL.
  INPUT   = 32,  // Input suffix: values were received from AMPL.
  OUTONLY = 64   // Output only: reject as an input value.
};
}
}  // namespace mp

#endif  // MP_PROBLEM_BASE_H_
