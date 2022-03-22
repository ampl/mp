/*
 Expression classes

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

#include "mp/expr.h"
#include "mp/expr-visitor.h"
#include "expr-writer.h"

#include <cstring>

using mp::Cast;

namespace {

/// Compares expressions for equality.
class ExprComparator : public mp::ExprVisitor<ExprComparator, bool> {
 private:
  mp::Expr expr_;

 public:
  explicit ExprComparator(Expr e) : expr_(e) {}

  template <typename T>
  bool VisitNumericConstant(T c) { return Cast<T>(expr_).value() == c.value(); }

  bool VisitVariable(Variable v) {
    return Cast<Variable>(expr_).index() == v.index();
  }

  bool VisitCommonExpr(CommonExpr e) {
    return Cast<Variable>(expr_).index() == e.index();
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
           Equal(if_expr.then_expr(), e.then_expr()) &&
           Equal(if_expr.else_expr(), e.else_expr());
  }

  bool VisitPLTerm(PLTerm e);

  bool VisitCall(CallExpr e);

  template <typename E>
  bool VisitVarArg(E e);

  bool VisitSum(SumExpr e) { return VisitVarArg(e); }
  bool VisitNumberOf(NumberOfExpr e) { return VisitVarArg(e); }
  bool VisitCount(CountExpr e) { return VisitVarArg(e); }

  bool VisitLogicalConstant(LogicalConstant c) {
    return VisitNumericConstant(c);
  }

  bool VisitNot(NotExpr e) { return VisitUnary(e); }
  bool VisitBinaryLogical(BinaryLogicalExpr e) { return VisitBinary(e); }
  bool VisitRelational(RelationalExpr e) { return VisitBinary(e); }
  bool VisitLogicalCount(LogicalCountExpr e) { return VisitBinary(e); }
  bool VisitImplication(ImplicationExpr e) { return VisitIf(e); }
  bool VisitIteratedLogical(IteratedLogicalExpr e) { return VisitVarArg(e); }
  bool VisitAllDiff(PairwiseExpr e) { return VisitVarArg(e); }
};

bool ExprComparator::VisitPLTerm(PLTerm e) {
  PLTerm pl = Cast<PLTerm>(expr_);
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

bool ExprComparator::VisitCall(CallExpr e) {
  CallExpr call = Cast<CallExpr>(expr_);
  int num_args = call.num_args();
  if (call.function() != e.function() || num_args != e.num_args())
    return false;
  for (int i = 0; i < num_args; ++i) {
    Expr arg = call.arg(i), other_arg = e.arg(i);
    if (arg.kind() != other_arg.kind())
      return false;
    if (NumericExpr num_arg = Cast<NumericExpr>(arg)) {
      if (!Equal(num_arg, Cast<NumericExpr>(other_arg)))
        return false;
    } else if (std::strcmp(
            Cast<StringLiteral>(arg).value(),
            Cast<StringLiteral>(other_arg).value()) != 0)
      return false;
  }
  return true;
}

template <typename E>
bool ExprComparator::VisitVarArg(E e) {
  E vararg = Cast<E>(expr_);
  typename E::iterator i = vararg.begin(), iend = vararg.end();
  typename E::iterator j = e.begin(), jend = e.end();
  for (; i != iend; ++i, ++j) {
    if (j == jend || !Equal(*i, *j))
      return false;
  }
  return j == jend;
}
}  // namespace

void mp::format(fmt::BasicFormatter<char> &f, const char *&, NumericExpr e) {
  fmt::MemoryWriter writer;
  ExprWriter<internal::ExprTypes>(writer).Visit(e);
  f.writer() << fmt::StringRef(writer.data(), writer.size());
}

bool mp::Equal(Expr e1, Expr e2) {
  if (e1.kind() != e2.kind())
    return false;
  return ExprComparator(e1).Visit(e2);
}

#ifdef MP_USE_HASH

using std::size_t;
using mp::internal::HashCombine;

inline size_t HashCombine(size_t seed, mp::Reference r) {
  return HashCombine<mp::Expr>(seed, r);
}

template <mp::expr::Kind FIRST, mp::expr::Kind LAST>
inline size_t HashCombine(size_t seed, mp::BasicExpr<FIRST, LAST> e) {
  return HashCombine<mp::Expr>(seed, e);
}

namespace {

/// Computes a hash value for an expression.
class ExprHasher : public mp::ExprVisitor<ExprHasher, size_t> {
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
  size_t VisitVariable(Variable v) { return Hash(v, v.index()); }
  size_t VisitCommonExpr(CommonExpr e) { return Hash(e, e.index()); }

  size_t VisitUnary(UnaryExpr e) { return Hash(e, e.arg()); }

  template <typename E>
  size_t VisitBinary(E e) { return HashCombine(Hash(e, e.lhs()), e.rhs()); }

  template <typename E>
  size_t VisitIf(E e) {
    size_t hash = HashCombine(Hash(e), e.condition());
    return HashCombine(HashCombine(hash, e.then_expr()), e.else_expr());
  }

  size_t VisitPLTerm(PLTerm e) {
    size_t hash = Hash(e);
    int num_breakpoints = e.num_breakpoints();
    for (int i = 0; i < num_breakpoints; ++i) {
      hash = HashCombine(hash, e.slope(i));
      hash = HashCombine(hash, e.breakpoint(i));
    }
    hash = HashCombine(hash, e.slope(num_breakpoints));
    return HashCombine(hash, e.arg());
  }

  size_t VisitCall(CallExpr e) {
    // Function name is hashed as a pointer. This works because the function
    // object is the same for all calls to the same function.
    size_t hash = Hash(e, e.function().name());
    for (int i = 0, n = e.num_args(); i < n; ++i)
      hash = HashCombine(hash, e.arg(i));
    return hash;
  }

  template <typename E>
  size_t VisitVarArg(E e) {
    size_t hash = Hash(e);
    for (typename E::iterator i = e.begin(), end = e.end(); i != end; ++i)
      hash = HashCombine(hash, *i);
    return hash;
  }

  size_t VisitSum(SumExpr e) { return VisitVarArg(e); }
  size_t VisitCount(CountExpr e) { return VisitVarArg(e); }
  size_t VisitNumberOf(NumberOfExpr e) { return VisitVarArg(e); }
  size_t VisitLogicalConstant(LogicalConstant c) { return Hash(c, c.value()); }
  size_t VisitNot(NotExpr e) { return Hash(e, e.arg()); }
  size_t VisitBinaryLogical(BinaryLogicalExpr e) { return VisitBinary(e); }
  size_t VisitRelational(RelationalExpr e) { return VisitBinary(e); }

  size_t VisitLogicalCount(LogicalCountExpr e) {
    return HashCombine(Hash(e, e.lhs()), NumericExpr(e.rhs()));
  }

  size_t VisitImplication(ImplicationExpr e) { return VisitIf(e); }
  size_t VisitIteratedLogical(IteratedLogicalExpr e) { return VisitVarArg(e); }
  size_t VisitAllDiff(PairwiseExpr e) { return VisitVarArg(e); }

  size_t VisitStringLiteral(StringLiteral s) {
    size_t hash = Hash(s);
    for (const char *value = s.value(); *value; ++value)
      hash = HashCombine(hash, *value);
    return hash;
  }
};
}  // namespace

size_t std::hash<mp::Expr>::operator()(mp::Expr expr) const {
  return ExprHasher().Visit(expr);
}

#endif  // MP_USE_HASH
