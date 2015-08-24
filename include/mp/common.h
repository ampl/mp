/*
 Common declarations

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

#ifndef MP_COMMON_H_
#define MP_COMMON_H_

#include <cassert>
#include <cstddef>  // for std::size_t

#define MP_DISPATCH(call) static_cast<Impl*>(this)->call

/** The mp namespace. */
namespace mp {

namespace internal {
// Suppresses warnings about unused variables.
inline void Unused(...) {}
}

/** Expression information. */
namespace expr {

/** Expression kind. */
enum Kind {
  /** An unknown expression. */
  UNKNOWN = 0,

  /** The first expression kind other than the unknown expression kind. */
  FIRST_EXPR,

  /**
    \rst
    .. _numeric-expr:

    The first numeric expression kind. Numeric expression kinds are in
    the range ``[FIRST_NUMERIC, LAST_NUMERIC]``.
    \endrst
   */
  FIRST_NUMERIC = FIRST_EXPR,

  /**
    \rst
    .. _number:

    A number such as ``42`` or ``-1.23e-4``.
    \endrst
   */
  NUMBER = FIRST_NUMERIC,

  /**
    \rst
    The first reference expression kind. Reference expression kinds are in
    the range ``[FIRST_REFERENCE, LAST_REFERENCE]``.
    \endrst
   */
  FIRST_REFERENCE,

  /**
    \rst
    .. _variable:

    A reference to a variable.
    \endrst
   */
  VARIABLE = FIRST_REFERENCE,

  /**
    \rst
    .. _common-expr:

    A reference to a common expression.
    \endrst
   */
  COMMON_EXPR,

  /** The last reference expression kind. */
  LAST_REFERENCE = COMMON_EXPR,

  /**
    \rst
    .. _unary-expr:

    The first unary numeric expression kind.
    Examples: ``-x``, ``sin(x)``, where ``x`` is a variable.
    \endrst
   */
  FIRST_UNARY,

  /** A unary minus expression. */
  MINUS = FIRST_UNARY,

  /** The `abs <http://ampl.com/BOOK/CHAPTERS/24-refman.pdf#page=6>`_ function. */
  ABS,

  /** The ``floor`` function. */
  FLOOR,

  /** The ``ceil`` function. */
  CEIL,

  /** The ``sqrt`` function. */
  SQRT,

  POW2,
  EXP,
  LOG,
  LOG10,
  SIN,
  SINH,
  COS,
  COSH,
  TAN,
  TANH,
  ASIN,
  ASINH,
  ACOS,
  ACOSH,
  ATAN,
  ATANH,
  
  /** The last unary numeric expression kind. */
  LAST_UNARY = ATANH,

  // Binary expressions.
  /**
    \rst
    .. _binary-expr:

    A binary numeric expression.
    Examples: ``x / y``, ``atan2(x, y)``, where ``x`` and ``y`` are variables.
    \endrst
   */
  FIRST_BINARY,
  ADD = FIRST_BINARY,
  SUB,
  LESS,
  MUL,
  DIV,
  INT_DIV,
  MOD,
  POW,
  POW_CONST_BASE,
  POW_CONST_EXP,
  ATAN2,
  PRECISION,
  ROUND,
  TRUNC,
  LAST_BINARY = TRUNC,

  /**
    \rst
    .. _if-expr:

    An if-then-else expression.
    Example: ``if x != 0 then y else z``, where ``x``, ``y`` and ``z`` are
    variables.
    \endrst
   */
  IF,

  /**
    \rst
    .. _plterm:

    A piecewise-linear term.
    Example: ``<<0; -1, 1>> x``, where ``x`` is a variable.
    \endrst
   */
  PLTERM,

  /**
    \rst
    .. _call-expr:

    A function call expression.
    Example: ``f(x)``, where ``f`` is a function and ``x`` is a variable.
    \endrst
   */
  CALL,

  // Iterated expressions.
  // The term "iterated" in the context of operators and expressions comes
  // from the article "AMPL: A Mathematical Programming Language" and is
  // used to denote operators indexed over sets.
  FIRST_ITERATED,

  /**
    \rst
    .. _vararg-expr:

    A varag expression (``min`` or ``max``). Example: ``min{i in I} x[i]``,
    where ``I`` is a set and ``x`` is a variable.
    \endrst
   */
  FIRST_VARARG = FIRST_ITERATED,
  MIN = FIRST_VARARG,
  MAX,
  LAST_VARARG = MAX,

  /**
    \rst
    .. _sum-expr:

    A sum expression. Example: ``sum{i in I} x[i]``, where ``I`` is a set and
    ``x`` is a variable.
    \endrst
   */
  SUM,

  /**
    \rst
    .. _numberof-expr:

    A numberof expression. Example: ``numberof 42 in ({i in I} x[i])``,
    where ``I`` is a set and ``x`` is a variable.
    \endrst
   */
  NUMBEROF,
  LAST_ITERATED = NUMBEROF,

  /**
    \rst
    .. _numberof-sym-expr:

    A symbolic numberof expression.
    Example: ``numberof (if x != 0 then 'a' else 'b') in ('a', 'b', 'c')``,
    where ``x`` is a variable.
    \endrst
   */
  NUMBEROF_SYM,

  /**
    \rst
    .. _count-expr:

    A count expression.
    Example: ``count{i in I} (x[i] >= 0)``, where ``I`` is a set and ``x``
    is a variable.
    \endrst
   */
  COUNT,
  LAST_NUMERIC = COUNT,

  // To simplify checks, logical expression kinds are in the range
  // [FIRST_LOGICAL, LAST_LOGICAL].
  /**
    \rst
    .. _logical-expr:

    A logical expression.
    \endrst
   */
  FIRST_LOGICAL,

  /**
    \rst
    .. _bool:

    A Boolean (logical) constant such as ``0`` or ``1``.
    \endrst
   */
  BOOL = FIRST_LOGICAL,

  /**
    \rst
    .. _not-expr:

    A logical NOT expression.
    Example: ``not a``, where ``a`` is a logical expression.
    \endrst
   */
  NOT,

  // Binary logical expressions.
  /**
    \rst
    .. _binary-logical-expr:

    A binary logical expression.
    Examples: ``a || b``, ``a && b``, where ``a`` and ``b`` are logical
    expressions.
    \endrst
   */
  FIRST_BINARY_LOGICAL,
  OR = FIRST_BINARY_LOGICAL,
  AND,
  IFF,
  LAST_BINARY_LOGICAL = IFF,

  // Relational expressions.
  /**
    \rst
    .. _relational-expr:

    A relational expression.
    Examples: ``a < b``, ``a != b``, where ``a`` and ``b`` are numeric
    expressions.
    \endrst
   */
  FIRST_RELATIONAL,
  LT = FIRST_RELATIONAL,  // <
  LE,                     // <=
  EQ,                     // =
  GE,                     // >=
  GT,                     // >
  NE,                     // !=
  LAST_RELATIONAL = NE,

  /**
    \rst
    .. _logical-count-expr:

    A logical count expression.
    Examples: ``atleast 1 (a < b, a != b)``, where ``a`` and ``b`` are
    numeric expressions.
    \endrst
   */
  FIRST_LOGICAL_COUNT,
  ATLEAST = FIRST_LOGICAL_COUNT,
  ATMOST,
  EXACTLY,
  NOT_ATLEAST,
  NOT_ATMOST,
  NOT_EXACTLY,
  LAST_LOGICAL_COUNT = NOT_EXACTLY,

  /**
    \rst
    .. _implication-expr:

    An implication expression.
    Example: ``a ==> b else c``, where ``a``, ``b`` and ``c`` are logical
    expressions.
    \endrst
   */
  IMPLICATION,

  // Iterated logical expressions.
  /**
    \rst
    .. _iterated-logical-expr:

    An iterated logical expression.
    Example: ``exists{i in I} x[i] >= 0``, where ``I`` is a set and ``x`` is a
    variable.
    \endrst
   */
  FIRST_ITERATED_LOGICAL,
  EXISTS = FIRST_ITERATED_LOGICAL,
  FORALL,
  LAST_ITERATED_LOGICAL = FORALL,

  // Pairwise expressions.
  /**
    \rst
    .. _pairwise-expr:

    A pairwise expression (``alldiff`` or ``!alldiff``).
    Example: ``alldiff{i in I} x[i]``, where ``I`` is a set and ``x`` is a
    variable.
    \endrst
   */
  FIRST_PAIRWISE,
  ALLDIFF = FIRST_PAIRWISE,
  NOT_ALLDIFF,
  LAST_PAIRWISE = NOT_ALLDIFF,
  LAST_LOGICAL = LAST_PAIRWISE,

  // String expressions.
  /**
    \rst
    .. _string:

    A string such as ``"abc"``.
    \endrst
   */
  STRING,

  /**
    \rst
    .. _ifsym:

    A symbolic if-then-else expression.
    Example: ``if x != 0 then 'a' else 0``, where ``x`` is a variable.
    \endrst
   */
  IFSYM,

  LAST_EXPR = IFSYM
};

// Maximum opcode.
enum { MAX_OPCODE = 82 };

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

// Returns the string representation of this expression kind.
// Expressions of different kinds can have identical strings.
// For example, POW, POW_CONST_BASE and POW_CONST_EXP all have
// the same representation "^".
const char *str(expr::Kind kind);
}  // namespace expr

namespace internal {

template<bool B, typename T = void>
struct enable_if {};

template<class T>
struct enable_if<true, T> { typedef T type; };

// Returns true if ExprType is of kind k.
template <typename ExprType>
inline bool Is(expr::Kind k) {
  int kind = k;
  // If FIRST_KIND == LAST_KIND, then a decent optimizing compiler simplifies
  // this to kind == ExprType::FIRST_KIND (checked with GCC 4.8.2).
  // No need to do it ourselves.
  return ExprType::FIRST_KIND <= kind && kind <= ExprType::LAST_KIND;
}

int precedence(expr::Kind kind);

// Expression information.
class ExprInfo {
 private:
  static const ExprInfo INFO[];

  friend int expr::opcode(expr::Kind kind);
  friend const char *expr::str(expr::Kind kind);

 public:
  int opcode;
  int precedence;
  const char *str;

  friend int precedence(expr::Kind kind) {
    assert(kind >= expr::UNKNOWN && kind <= expr::LAST_EXPR);
    return internal::ExprInfo::INFO[kind].precedence;
  }
};
}

inline int expr::opcode(expr::Kind kind) {
  assert(kind >= expr::UNKNOWN && kind <= expr::LAST_EXPR);
  return internal::ExprInfo::INFO[kind].opcode;
}

inline const char *expr::str(expr::Kind kind) {
  assert(kind >= expr::UNKNOWN && kind <= expr::LAST_EXPR);
  return internal::ExprInfo::INFO[kind].str;
}

/** Function information. */
namespace func {
/** Function type. */
enum Type {
  /** A numeric function. */
  NUMERIC  = 0,
  /** A symbolic function - accepts numeric and string arguments. */
  SYMBOLIC = 1
};
}

namespace var {
// Variable type.
enum Type { CONTINUOUS, INTEGER };
}

/** Objective information. */
namespace obj {
/** Objective type. */
enum Type {
  MIN = 0, /**< A minimization objective. */
  MAX = 1  /**< A maximization objective. */
};
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
  MASK    =  3,  // Mask for suffix kind.
  FLOAT   =  4,  // Suffix values are floating-point numbers.
  IODECL  =  8,  // Tell AMPL to make this an INOUT suffix.
  OUTPUT  = 16,  // Output suffix: return values to AMPL.
  INPUT   = 32,  // Input suffix: values were received from AMPL.
  OUTONLY = 64   // Output only: reject as an input value.
};
}

namespace sol {
// Solution status.
enum Status {
  UNKNOWN     =  -1,

  // An optimal solution found for an optimization problem or a feasible
  // solution found for a satisfaction problem.
  SOLVED      =   0,

  // Solution returned but it can be non-optimal or even infeasible.
  UNSOLVED    = 100,

  // Problem is infeasible.
  INFEASIBLE  = 200,

  // Problem is unbounded.
  UNBOUNDED   = 300,

  // Stopped by a limit, e.g. on iterations or time.
  LIMIT       = 400,

  FAILURE     = 500,

  // Interrupted by the user.
  INTERRUPTED = 600
};
}

/** Information about an optimization problem. */
struct ProblemInfo {
  /** Total number of variables. */
  int num_vars;

  /**
    Number of algebraic constraints including ranges and equality constraints.
    It doesn't include logical constraints.
   */
  int num_algebraic_cons;

  /** Total number of objectives. */
  int num_objs;

  /** Number of ranges (constraints with -Infinity < LHS < RHS < Infinity). */
  int num_ranges;

  /**
    Number of equality constraints or -1 if unknown (AMPL prior to 19970627).
   */
  int num_eqns;

  /** Number of logical constraints. */
  int num_logical_cons;

  /** Returns the number of integer variables (includes binary variable). */
  int num_integer_vars() const {
    return num_linear_binary_vars + num_linear_integer_vars +
        num_nl_integer_vars_in_both + num_nl_integer_vars_in_cons +
        num_nl_integer_vars_in_objs;
  }

  /** Returns the number of continuous variables. */
  int num_continuous_vars() const { return num_vars - num_integer_vars(); }

  // Nonlinear and complementarity information
  // -----------------------------------------

  /** Total number of nonlinear constraints. */
  int num_nl_cons;

  /** Total number of nonlinear objectives. */
  int num_nl_objs;

  /** Total number of complementarity conditions. */
  int num_compl_conds;

  /** Number of nonlinear complementarity conditions. */
  int num_nl_compl_conds;

  /**
    Number of complementarities involving double inequalities
    (for ``ASL_cc_simplify``).
   */
  int num_compl_dbl_ineqs;

  /**
    Number of complemented variables with a nonzero lower bound
    (for ASL_cc_simplify).
   */
  int num_compl_vars_with_nz_lb;

  // Information about network constraints
  // -------------------------------------

  /** Number of nonlinear network constraints. */
  int num_nl_net_cons;

  /** Number of linear network constraints. */
  int num_linear_net_cons;

  // Information about nonlinear variables
  // -------------------------------------

  /**
    Number of nonlinear variables in constraints including nonlinear
    variables in both constraints and objectives.
   */
  int num_nl_vars_in_cons;

  /**
    Number of nonlinear variables in objectives including nonlinear
    variables in both constraints and objectives.
   */
  int num_nl_vars_in_objs;

  /** Number of nonlinear variables in both constraints and objectives. */
  int num_nl_vars_in_both;

  // Miscellaneous
  // -------------

  /** Number of linear network variables (arcs). */
  int num_linear_net_vars;

  /** Number of functions. */
  int num_funcs;

  // Information about discrete variables
  // ------------------------------------

  /** Number of linear binary variables. */
  int num_linear_binary_vars;

  /** Number of linear non-binary integer variables. */
  int num_linear_integer_vars;

  /**
    Number of integer nonlinear variables in both constraints and objectives.
   */
  int num_nl_integer_vars_in_both;

  /** Number of integer nonlinear variables just in constraints. */
  int num_nl_integer_vars_in_cons;

  /** Number of integer nonlinear variables just in objectives. */
  int num_nl_integer_vars_in_objs;

  // Information about nonzeros
  // --------------------------

  /** Number of nonzeros in constraints' Jacobian. */
  std::size_t num_con_nonzeros;

  /** Number of nonzeros in all objective gradients. */
  std::size_t num_obj_nonzeros;

  // Information about names
  // -----------------------

  /** Length of longest constraint name (if ``stub.row`` exists). */
  int max_con_name_len;

  /** Length of longest variable name (if ``stub.col`` exists). */
  int max_var_name_len;

  // Information about common expressions
  // ------------------------------------

  int num_common_exprs_in_both;
  int num_common_exprs_in_cons;
  int num_common_exprs_in_objs;

  /**
    Number of common expressions that only appear in a single constraint
    and don't appear in objectives.
   */
  int num_common_exprs_in_single_cons;

  /**
    Number of common expressions that only appear in a single objective
    and don't appear in constraints.
   */
  int num_common_exprs_in_single_objs;

  /** Returns the total number of common expressions. */
  int num_common_exprs() const {
    return num_common_exprs_in_both + num_common_exprs_in_cons +
        num_common_exprs_in_objs + num_common_exprs_in_single_cons +
        num_common_exprs_in_single_objs;
  }
};
}  // namespace mp

#endif  // MP_COMMON_H_
