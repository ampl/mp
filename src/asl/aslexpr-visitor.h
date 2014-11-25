/*
 ASL expression visitor

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

#ifndef MP_ASL_ASLEXPR_VISITOR_H_
#define MP_ASL_ASLEXPR_VISITOR_H_

#include "asl/aslexpr.h"
#include "mp/expr-visitor.h"

namespace mp {
namespace asl {
// ASL expression types.
struct ExprTypes {
  typedef mp::asl::NumericExpr NumericExpr;
  typedef mp::asl::LogicalExpr LogicalExpr;
  typedef mp::asl::NumericConstant NumericConstant;
  typedef mp::asl::Variable Variable;
  typedef mp::asl::UnaryExpr UnaryExpr;
  typedef mp::asl::BinaryExpr BinaryExpr;
  typedef mp::asl::IfExpr IfExpr;
  typedef mp::asl::PiecewiseLinearExpr PLTerm;
  typedef mp::asl::CallExpr CallExpr;
  typedef mp::asl::VarArgExpr VarArgExpr;
  typedef mp::asl::SumExpr SumExpr;
  typedef mp::asl::NumberOfExpr NumberOfExpr;
  typedef mp::asl::CountExpr CountExpr;
  typedef mp::asl::LogicalConstant LogicalConstant;
  typedef mp::asl::NotExpr NotExpr;
  typedef mp::asl::BinaryLogicalExpr BinaryLogicalExpr;
  typedef mp::asl::RelationalExpr RelationalExpr;
  typedef mp::asl::LogicalCountExpr LogicalCountExpr;
  typedef mp::asl::ImplicationExpr ImplicationExpr;
  typedef mp::asl::IteratedLogicalExpr IteratedLogicalExpr;
  typedef mp::asl::PairwiseExpr PairwiseExpr;
};

// ASL expression visitor.
template <typename Impl, typename Result, typename LResult = Result>
class ExprVisitor :
    public mp::ExprVisitor<Impl, Result, LResult, mp::asl::ExprTypes> {};

// Returns true iff e is a zero constant.
inline bool IsZero(NumericExpr e) {
  NumericConstant c = Cast<NumericConstant>(e);
  return c && c.value() == 0;
}

// Expression converter.
// Converts logical count expressions to corresponding relational expressions.
// For example "atleast" is converted to "<=".
template <typename Impl, typename Result, typename LResult = Result>
class ExprConverter : public ExprVisitor<Impl, Result, LResult> {
 private:
  std::vector< ::expr> exprs_;

  RelationalExpr Convert(LogicalCountExpr e, expr::Kind kind) {
    exprs_.push_back(*e.impl_);
    ::expr *result = &exprs_.back();
    result->op = reinterpret_cast<efunc*>(opcode(kind));
    return Expr::Create<RelationalExpr>(result);
  }

 public:
  LResult VisitAtLeast(LogicalCountExpr e) {
    return MP_DISPATCH(VisitLE(Convert(e, expr::LE)));
  }
  LResult VisitAtMost(LogicalCountExpr e) {
    return MP_DISPATCH(VisitGE(Convert(e, expr::GE)));
  }
  LResult VisitExactly(LogicalCountExpr e) {
    return MP_DISPATCH(VisitEQ(Convert(e, expr::EQ)));
  }
  LResult VisitNotAtLeast(LogicalCountExpr e) {
    return MP_DISPATCH(VisitGT(Convert(e, expr::GT)));
  }
  LResult VisitNotAtMost(LogicalCountExpr e) {
    return MP_DISPATCH(VisitLT(Convert(e, expr::LT)));
  }
  LResult VisitNotExactly(LogicalCountExpr e) {
    return MP_DISPATCH(VisitNE(Convert(e, expr::NE)));
  }
};
}  // namespace asl
}  // namespace mp

#endif  // MP_ASL_ASLEXPR_VISITOR_H_
