/*
 Base types and functions related to the optimization problem.

 This header is used to decouple the .nl reader from a specific problem
 representation as the reader can be used to construct different types of
 problems. So instead of including problem.h and expr.h from nl.h, this
 header containing common stuff is included from problem.h, expr.h and nl.h.

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
#include <cstddef>

#include "mp/opcode.hd"

namespace mp {

class Expr;

namespace expr {

// Expression kinds.
enum Kind {
  // An unknown expression.
  UNKNOWN = 0,

  EXPR_START,

  // To simplify checks, numeric expression kinds are in a range
  // [NUMERIC_START, NUMERIC_END].
  NUMERIC_START = EXPR_START,
  VARIABLE = NUMERIC_START,
  UNARY,
  BINARY,
  IF,
  PLTERM,
  CALL,
  VARARG,
  SUM,
  COUNT,
  NUMBEROF,
  NUMERIC_END,

  // CONSTANT belongs both to numeric and logical expressions therefore
  // the [NUMERIC_START, NUMERIC_END] and [LOGICAL_START, LOGICAL_END]
  // ranges overlap at CONSTANT = NUMERIC_END = LOGICAL_START.
  CONSTANT = NUMERIC_END,

  // To simplify checks, logical expression kinds are in a range
  // [LOGICAL_START, LOGICAL_END].
  LOGICAL_START = CONSTANT,
  NOT,
  BINARY_LOGICAL,
  RELATIONAL,
  LOGICAL_COUNT,
  IMPLICATION,
  ITERATED_LOGICAL,
  ALLDIFF,
  LOGICAL_END = ALLDIFF,

  STRING,
  EXPR_END = STRING
};

Kind kind(int opcode);

// Expression information.
class Info {
 private:
  static const Info INFO[N_OPS];

  friend class mp::Expr;

 public:
  expr::Kind kind;
  int precedence;
  const char *str;

  // Returns the expression kind for the opcode which should be in the
  // range [0, N_OPS).
  friend Kind kind(int opcode) {
    assert(opcode >= 0 && opcode < N_OPS);
    return INFO[opcode].kind;
  }
};
}  // namespace expr

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

// A reference to an immutable array.
template <typename T>
class ArrayRef {
 private:
  const T *data_;
  std::size_t size_;

 public:
  ArrayRef(const T *data, std::size_t size) : data_(data), size_(size) {}

  template <typename U>
  ArrayRef(ArrayRef<U> other) : data_(other.data()), size_(other.size()) {}

  template <std::size_t SIZE>
  ArrayRef(const T (&data)[SIZE]) : data_(data), size_(SIZE) {}

  const T *data() const { return data_; }
  std::size_t size() const { return size_; }

  const T &operator[](std::size_t i) const { return data_[i]; }
};

template <typename T>
ArrayRef<T> MakeArrayRef(const T *data, std::size_t size) {
  return ArrayRef<T>(data, size);
}
}  // namespace mp

#endif  // MP_PROBLEM_BASE_H_
