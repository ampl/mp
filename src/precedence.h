/*
 Operator precedence.

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

#ifndef MP_PRECEDENCE_H_
#define MP_PRECEDENCE_H_

#include "mp/common.h"

namespace mp {
namespace prec {
enum Precedence {
  UNKNOWN,
  CONDITIONAL,       // if-then-else
  IFF,               // <==>
  IMPLICATION,       // ==> else
  LOGICAL_OR,        // or ||
  LOGICAL_AND,       // and &&
  NOT,               // not !
  RELATIONAL,        // < <= = == >= > != <>
  PIECEWISE_LINEAR,  // a piecewise-linear expression
  ADDITIVE,          // + - less
  ITERATIVE,         // sum prod min max
  MULTIPLICATIVE,    // * / div mod
  EXPONENTIATION,    // ^
  UNARY,             // + - (unary)
  CALL,              // a function call including functional forms of
                     // min and max
  PRIMARY            // variable, string or constant
};
}

prec::Precedence precedence(expr::Kind kind);

namespace internal {
class PrecInfo {
 private:
  static const prec::Precedence INFO[expr::LAST_EXPR + 1];
  friend prec::Precedence mp::precedence(expr::Kind kind);
};
}

inline prec::Precedence precedence(expr::Kind kind) {
  assert(kind >= expr::UNKNOWN && kind <= expr::LAST_EXPR);
  return internal::PrecInfo::INFO[kind];
}
}  // namespace mp

#endif  // MP_PRECEDENCE_H_
