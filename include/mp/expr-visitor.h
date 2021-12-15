/*
 Expression visitor

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

#ifndef MP_EXPR_VISITOR_H_
#define MP_EXPR_VISITOR_H_

#include "mp/basic-expr-visitor.h"
#include "mp/expr.h"

namespace mp {

/// An expression visitor.
///
/// To use ExprVisitor define a subclass that implements some or all of the
/// Visit* methods with the same signatures as the methods in ExprVisitor,
/// for example, VisitDiv(BinaryExpr).
/// Specify the subclass name as the Impl template parameter. Then calling
/// BasicExprVisitor::Visit for some expression will dispatch to a Visit*
/// method specific to the expression type. For example, if the expression is
/// a division then VisitDiv(BinaryExpr) method of a subclass will be called.
/// If the subclass doesn't contain a method with this signature, then
/// a corresponding method of BasicExprVisitor will be called.
///
/// Example:
///  class MyExprVisitor : public ExprVisitor<MyExprVisitor, double> {
///   public:
///    double VisitAdd(BinaryExpr e) { return Visit(e.lhs()) + Visit(e.rhs()); }
///    double VisitConstant(NumericConstant n) { return n.value(); }
///  };
template <typename Impl, typename Result>
class ExprVisitor :
    public BasicExprVisitor<Impl, Result, internal::ExprTypes> {};

/// Expression converter.
/// Converts logical count expressions to corresponding relational expressions.
/// For example "atleast" is converted to "<=".
template <typename Impl, typename Result>
class ExprConverter : public ExprVisitor<Impl, Result> {
 private:
  ExprFactory factory_;

  RelationalExpr Convert(LogicalCountExpr e, expr::Kind kind) {
    return factory_.MakeRelational(kind, e.lhs(), e.rhs());
  }

 public:
  Result VisitAtLeast(LogicalCountExpr e) {
    return MP_DISPATCH(VisitLE(Convert(e, expr::LE)));
  }
  Result VisitAtMost(LogicalCountExpr e) {
    return MP_DISPATCH(VisitGE(Convert(e, expr::GE)));
  }
  Result VisitExactly(LogicalCountExpr e) {
    return MP_DISPATCH(VisitEQ(Convert(e, expr::EQ)));
  }
  Result VisitNotAtLeast(LogicalCountExpr e) {
    return MP_DISPATCH(VisitGT(Convert(e, expr::GT)));
  }
  Result VisitNotAtMost(LogicalCountExpr e) {
    return MP_DISPATCH(VisitLT(Convert(e, expr::LT)));
  }
  Result VisitNotExactly(LogicalCountExpr e) {
    return MP_DISPATCH(VisitNE(Convert(e, expr::NE)));
  }
};
}  // namespace mp

#endif  // MP_EXPR_VISITOR_H_
