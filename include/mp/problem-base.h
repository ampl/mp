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

#define MP_DISPATCH(call) static_cast<Impl*>(this)->call

// Suppresses warnings about unused variables.
#define MP_UNUSED(x) (void)(x)
#define MP_UNUSED2(x, y) MP_UNUSED(x); MP_UNUSED(y)
#define MP_UNUSED3(x, y, z) MP_UNUSED2(x, y); MP_UNUSED(z)

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
  NUM_KINDS,     // The number of suffix kinds.
  MASK    =  3,  // Mask for the above.
  FLOAT   =  4,  // Suffix values are floating-point numbers.
  IODECL  =  8,  // Tell AMPL to make this an INOUT suffix.
  OUTPUT  = 16,  // Output suffix: return values to AMPL.
  INPUT   = 32,  // Input suffix: values were received from AMPL.
  OUTONLY = 64   // Output only: reject as an input value.
};
}

// Information about an optimization problem.
struct ProblemInfo {
  // Total number of variables.
  int num_vars;

  // Number of algebraic constraints including ranges and equality constraints.
  // It doesn't include logical constraints.
  int num_algebraic_cons;

  // Total number of objectives.
  int num_objs;

  // Number of ranges (constraints with -Infinity < LHS < RHS < Infinity).
  int num_ranges;

  // Number of equality constraints or -1 if unknown (AMPL prior to 19970627).
  int num_eqns;

  // Number of logical constraints.
  int num_logical_cons;

  int num_integer_vars() const {
    return num_linear_binary_vars + num_linear_integer_vars +
        num_nl_integer_vars_in_both + num_nl_integer_vars_in_cons +
        num_nl_integer_vars_in_objs;
  }

  int num_continuous_vars() const { return num_vars - num_integer_vars(); }

  // Nonlinear and complementarity information
  // -----------------------------------------

  // Total number of nonlinear constraints.
  int num_nl_cons;

  // Total number of nonlinear objectives.
  int num_nl_objs;

  // Total number of complementarity conditions.
  int num_compl_conds;

  // Number of nonlinear complementarity conditions.
  int num_nl_compl_conds;

  // Number of complementarities involving double inequalities
  // (for ASL_cc_simplify).
  int num_compl_dbl_ineqs;

  // Number of complemented variables with a nonzero lower bound
  // (for ASL_cc_simplify).
  int num_compl_vars_with_nz_lb;

  // Information about network constraints
  // -------------------------------------

  // Number of nonlinear network constraints.
  int num_nl_net_cons;

  // Number of linear network constraints.
  int num_linear_net_cons;

  // Information about nonlinear variables
  // -------------------------------------

  // Number of nonlinear variables in constraints including nonlinear
  // variables in both constraints and objectives.
  int num_nl_vars_in_cons;

  // Number of nonlinear variables in objectives including nonlinear
  // variables in both constraints and objectives.
  int num_nl_vars_in_objs;

  // Number of nonlinear variables in both constraints and objectives.
  int num_nl_vars_in_both;

  // Miscellaneous
  // -------------

  // Number of linear network variables (arcs).
  int num_linear_net_vars;

  // Number of functions.
  int num_funcs;

  // Information about discrete variables
  // ------------------------------------

  // Number of linear binary variables.
  int num_linear_binary_vars;

  // Number of linear non-binary integer variables.
  int num_linear_integer_vars;

  // Number of integer nonlinear variables in both constraints and objectives.
  int num_nl_integer_vars_in_both;

  // Number of integer nonlinear variables just in constraints.
  int num_nl_integer_vars_in_cons;

  // Number of integer nonlinear variables just in objectives.
  int num_nl_integer_vars_in_objs;

  // Information about nonzeros
  // --------------------------

  // Number of nonzeros in constraints' Jacobian.
  std::size_t num_con_nonzeros;

  // Number of nonzeros in all objective gradients.
  std::size_t num_obj_nonzeros;

  // Information about names
  // -----------------------

  // Length of longest constraint name (if stub.row exists).
  int max_con_name_len;

  // Length of longest variable name (if stub.col exists).
  int max_var_name_len;

  // Information about common expressions
  // ------------------------------------

  int num_common_exprs_in_both;
  int num_common_exprs_in_cons;
  int num_common_exprs_in_objs;

  // Number of common expressions that only appear in a single constraint
  // and don't appear in objectives.
  int num_common_exprs_in_single_cons;

  // Number of common expressions that only appear in a single objective
  // and don't appear in constraints.
  int num_common_exprs_in_single_objs;

  int num_common_exprs() const {
    return num_common_exprs_in_both + num_common_exprs_in_cons +
        num_common_exprs_in_objs + num_common_exprs_in_single_cons +
        num_common_exprs_in_single_objs;
  }
};
}  // namespace mp

#endif  // MP_PROBLEM_BASE_H_
