/*
 A C++ interface to AMPL expressions.

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

#include "aslexpr.h"
#include "expr-writer.h"
#include "precedence.h"

#include "asl/aslexpr-visitor.h"

#include <cstdio>
#include <cstring>

using std::size_t;

using mp::asl::Cast;
using mp::asl::Expr;
using mp::asl::NumericConstant;
using mp::asl::NumericExpr;
using mp::asl::LogicalExpr;
namespace prec = mp::prec;
namespace asl = mp::asl;

namespace {

// Compares expressions for equality.
class ExprEqual : public mp::asl::ExprVisitor<ExprEqual, bool, bool> {
 private:
  Expr expr_;

 public:
  explicit ExprEqual(Expr e) : expr_(e) {}

  template <typename T>
  bool VisitNumericConstant(T c) { return Cast<T>(expr_).value() == c.value(); }

  bool VisitVariable(asl::Variable v) {
    return Cast<asl::Variable>(expr_).index() == v.index();
  }

  template <typename E>
  bool VisitUnary(E e) {
    return Equal(Cast<E>(expr_).arg(), e.arg());
  }

  template <typename E>
  bool VisitBinary(E e) {
    E binary = Cast<E>(expr_);
    return Equal(binary.lhs(), e.lhs()) && Equal(binary.rhs(), e.rhs());
  }

  template <typename E>
  bool VisitIf(E e) {
    E if_expr = Cast<E>(expr_);
    return Equal(if_expr.condition(), e.condition()) &&
           Equal(if_expr.true_expr(), e.true_expr()) &&
           Equal(if_expr.false_expr(), e.false_expr());
  }

  bool VisitPLTerm(asl::PiecewiseLinearExpr e) {
    asl::PiecewiseLinearExpr pl = Cast<asl::PiecewiseLinearExpr>(expr_);
    int num_breakpoints = pl.num_breakpoints();
    if (num_breakpoints != e.num_breakpoints())
      return false;
    for (int i = 0; i < num_breakpoints; ++i) {
      if (pl.slope(i) != e.slope(i) || pl.breakpoint(i) != e.breakpoint(i))
        return false;
    }
    return pl.slope(num_breakpoints) == e.slope(num_breakpoints) &&
           Equal(pl.arg(), e.arg());
  }

  bool VisitCall(asl::CallExpr e) {
    asl::CallExpr call = Cast<asl::CallExpr>(expr_);
    int num_args = call.num_args();
    if (call.function() != e.function() || num_args != e.num_args())
      return false;
    for (int i = 0; i < num_args; ++i) {
      Expr arg = call[i], other_arg = e[i];
      if (arg.kind() != other_arg.kind())
        return false;
      if (NumericExpr num_arg = Cast<NumericExpr>(arg)) {
        if (!Equal(num_arg, Cast<NumericExpr>(other_arg)))
          return false;
      } else if (std::strcmp(
              Cast<asl::StringLiteral>(arg).value(),
              Cast<asl::StringLiteral>(other_arg).value()) != 0)
        return false;
    }
    return true;
  }

  template <typename E>
  bool VisitVarArg(E e) {
    E vararg = Cast<E>(expr_);
    typename E::iterator i = vararg.begin(), iend = vararg.end();
    typename E::iterator j = e.begin(), jend = e.end();
    for (; i != iend; ++i, ++j) {
      if (j == jend || !Equal(*i, *j))
        return false;
    }
    return j == jend;
  }

  bool VisitSum(asl::SumExpr e) { return VisitVarArg(e); }
  bool VisitCount(asl::CountExpr e) { return VisitVarArg(e); }
  bool VisitNumberOf(asl::NumberOfExpr e) { return VisitVarArg(e); }

  bool VisitLogicalConstant(asl::LogicalConstant c) {
    return VisitNumericConstant(c);
  }

  bool VisitNot(asl::NotExpr e) { return VisitUnary(e); }

  bool VisitBinaryLogical(asl::BinaryLogicalExpr e) { return VisitBinary(e); }
  bool VisitRelational(asl::RelationalExpr e) { return VisitBinary(e); }
  bool VisitLogicalCount(asl::LogicalCountExpr e) { return VisitBinary(e); }

  bool VisitImplication(asl::ImplicationExpr e) { return VisitIf(e); }

  bool VisitIteratedLogical(asl::IteratedLogicalExpr e) {
    return VisitVarArg(e);
  }

  bool VisitAllDiff(asl::PairwiseExpr e) { return VisitVarArg(e); }
};
}  // namespace

#ifdef MP_USE_UNORDERED_MAP

using asl::internal::HashCombine;

namespace {
// Computes a hash value for an expression.
class ExprHasher : public mp::asl::ExprVisitor<ExprHasher, size_t, size_t> {
 private:
  static size_t Hash(Expr e) {
    return HashCombine<int>(0, e.kind());
  }

  template <typename T>
  static size_t Hash(Expr e, const T &value) {
    return HashCombine(Hash(e), value);
  }

 public:
  size_t VisitNumericConstant(NumericConstant c) { return Hash(c, c.value()); }
  size_t VisitVariable(asl::Variable v) { return Hash(v, v.index()); }

  size_t VisitUnary(asl::UnaryExpr e) { return Hash(e, e.arg()); }

  template <typename E>
  size_t VisitBinary(E e) { return HashCombine(Hash(e, e.lhs()), e.rhs()); }

  template <typename E>
  size_t VisitIf(E e) {
    size_t hash = HashCombine(Hash(e), e.condition());
    return HashCombine(HashCombine(hash, e.true_expr()), e.false_expr());
  }

  size_t VisitPLTerm(asl::PiecewiseLinearExpr e) {
    size_t hash = Hash(e);
    int num_breakpoints = e.num_breakpoints();
    for (int i = 0; i < num_breakpoints; ++i) {
      hash = HashCombine(hash, e.slope(i));
      hash = HashCombine(hash, e.breakpoint(i));
    }
    hash = HashCombine(hash, e.slope(num_breakpoints));
    return HashCombine(hash, e.arg());
  }

  size_t VisitCall(asl::CallExpr e) {
    // Function name is hashed as a pointer. This works because the function
    // object is the same for all calls to the same function.
    size_t hash = Hash(e, e.function().name());
    for (int i = 0, n = e.num_args(); i < n; ++i)
      hash = HashCombine(hash, e[i]);
    return hash;
  }

  template <typename E>
  size_t VisitVarArg(E e) {
    size_t hash = Hash(e);
    for (typename E::iterator i = e.begin(), end = e.end(); i != end; ++i)
      hash = HashCombine(hash, *i);
    return hash;
  }

  size_t VisitSum(asl::SumExpr e) { return VisitVarArg(e); }
  size_t VisitCount(asl::CountExpr e) { return VisitVarArg(e); }
  size_t VisitNumberOf(asl::NumberOfExpr e) { return VisitVarArg(e); }

  size_t VisitLogicalConstant(asl::LogicalConstant c) {
    return Hash(c, c.value());
  }

  size_t VisitNot(asl::NotExpr e) { return Hash(e, e.arg()); }

  size_t VisitBinaryLogical(asl::BinaryLogicalExpr e) { return VisitBinary(e); }
  size_t VisitRelational(asl::RelationalExpr e) { return VisitBinary(e); }

  size_t VisitLogicalCount(asl::LogicalCountExpr e) {
    NumericExpr rhs = e.rhs();
    return HashCombine(Hash(e, e.lhs()), rhs);
  }

  size_t VisitImplication(asl::ImplicationExpr e) { return VisitIf(e); }

  size_t VisitIteratedLogical(asl::IteratedLogicalExpr e) {
    return VisitVarArg(e);
  }

  size_t VisitAllDiff(asl::PairwiseExpr e) { return VisitVarArg(e); }

  size_t VisitStringLiteral(asl::StringLiteral s) {
    size_t hash = Hash(s);
    for (const char *value = s.value(); *value; ++value)
      hash = HashCombine(hash, *value);
    return hash;
  }
};
}

namespace std {
template <>
struct hash<Expr> {
  std::size_t operator()(Expr e) const;
};
}

size_t std::hash<Expr>::operator()(Expr expr) const {
  ExprHasher hasher;
  NumericExpr n = Cast<NumericExpr>(expr);
  return n ? hasher.Visit(n) :
             hasher.VisitStringLiteral(Cast<asl::StringLiteral>(expr));
}

size_t std::hash<NumericExpr>::operator()(NumericExpr expr) const {
  return ExprHasher().Visit(expr);
}

size_t std::hash<LogicalExpr>::operator()(LogicalExpr expr) const {
  return ExprHasher().Visit(expr);
}
#else
# if defined(AMPL_NO_UNORDERED_MAP_WARNING)
  // Do nothing.
# elif defined(_MSC_VER)
#  pragma message("warning: unordered_map not available, numberof may be slow")
# else
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wpedantic"
#  warning "unordered_map not available, numberof may be slow"
#  pragma clang diagnostic pop
# endif
#endif

namespace mp {
namespace asl {

const de VarArgExpr::END = de();

bool Equal(NumericExpr e1, NumericExpr e2) {
  if (e1.kind() != e2.kind())
    return false;
  return ExprEqual(e1).Visit(e2);
}

bool Equal(LogicalExpr e1, LogicalExpr e2) {
  if (e1.kind() != e2.kind())
    return false;
  return ExprEqual(e1).Visit(e2);
}

#ifdef MP_USE_UNORDERED_MAP
size_t internal::HashNumberOfArgs::operator()(NumberOfExpr e) const {
  size_t hash = 0;
  for (int i = 1, n = e.num_args(); i < n; ++i)
    hash = HashCombine(hash, e[i]);
  return hash;
}
#endif

bool internal::EqualNumberOfArgs::operator()(
    NumberOfExpr lhs, NumberOfExpr rhs) const {
  int num_args = lhs.num_args();
  if (num_args != rhs.num_args())
    return false;
  for (int i = 1; i < num_args; ++i) {
    if (!Equal(lhs[i], rhs[i]))
      return false;
  }
  return true;
}
}
}
