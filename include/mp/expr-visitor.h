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

#include "mp/error.h"
#include "mp/problem-base.h"

namespace mp {

class ExprTypes;

// An exception that is thrown when an unsuported expression is encountered.
class UnsupportedExprError : public UnsupportedError {
 private:
  explicit UnsupportedExprError(fmt::StringRef message)
    : UnsupportedError(message) {}

 public:
  static UnsupportedExprError CreateFromExprString(fmt::StringRef expr) {
    return UnsupportedExprError(std::string("unsupported: ") + expr.c_str());
  }
};

// An exception that is thrown when an invalid numeric expression
// is encountered.
class InvalidNumericExprError : public Error {
 public:
  template <typename NumericExpr>
  explicit InvalidNumericExprError(NumericExpr e) :
    Error("invalid numeric expression: {}", e.kind()) {}
};

// An exception that is thrown when an invalid logical or constraint
// expression is encountered.
class InvalidLogicalExprError : public Error {
 public:
  template <typename LogicalExpr>
  explicit InvalidLogicalExprError(LogicalExpr e) :
    Error("invalid logical expression: {}", e.kind()) {}
};

// An expression visitor.
// To use ExprVisitor define a subclass that implements some or all of the
// Visit* methods with the same signatures as the methods in ExprVisitor,
// for example, VisitDiv(BinaryExpr).
// Specify the subclass name as the Impl template parameter. Then calling
// ExprVisitor::Visit for some expression will dispatch to a Visit* method
// specific to the expression type. For example, if the expression is
// a division then VisitDiv(BinaryExpr) method of a subclass will be called.
// If the subclass doesn't contain a method with this signature, then
// a corresponding method of ExprVisitor will be called.
//
// Example:
//  class MyExprVisitor : public ExprVisitor<MyExprVisitor, double, void> {
//   public:
//    double VisitAdd(BinaryExpr e) { return Visit(e.lhs()) + Visit(e.rhs()); }
//    double VisitConstant(NumericConstant n) { return n.value(); }
//  };
//
// ExprVisitor uses the curiously recurring template pattern:
// http://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
template <typename Impl, typename Result, typename LResult = Result,
          typename ExprTypes = mp::ExprTypes>
class ExprVisitor {
 public:
  typedef typename ExprTypes::NumericExpr NumericExpr;
  typedef typename ExprTypes::LogicalExpr LogicalExpr;
  typedef typename ExprTypes::NumericConstant NumericConstant;
  typedef typename ExprTypes::Variable Variable;
  typedef typename ExprTypes::UnaryExpr UnaryExpr;
  typedef typename ExprTypes::BinaryExpr BinaryExpr;
  typedef typename ExprTypes::IfExpr IfExpr;
  typedef typename ExprTypes::PLTerm PLTerm;
  typedef typename ExprTypes::CallExpr CallExpr;
  typedef typename ExprTypes::VarArgExpr VarArgExpr;
  typedef typename ExprTypes::SumExpr SumExpr;
  typedef typename ExprTypes::CountExpr CountExpr;
  typedef typename ExprTypes::NumberOfExpr NumberOfExpr;
  typedef typename ExprTypes::LogicalConstant LogicalConstant;
  typedef typename ExprTypes::NotExpr NotExpr;
  typedef typename ExprTypes::BinaryLogicalExpr BinaryLogicalExpr;
  typedef typename ExprTypes::RelationalExpr RelationalExpr;
  typedef typename ExprTypes::LogicalCountExpr LogicalCountExpr;
  typedef typename ExprTypes::ImplicationExpr ImplicationExpr;
  typedef typename ExprTypes::IteratedLogicalExpr IteratedLogicalExpr;
  typedef typename ExprTypes::PairwiseExpr PairwiseExpr;

  Result Visit(NumericExpr e);
  LResult Visit(LogicalExpr e);

  Result VisitUnhandledNumericExpr(NumericExpr e) {
    throw UnsupportedExprError::CreateFromExprString(e.opstr());
  }

  LResult VisitUnhandledLogicalExpr(LogicalExpr e) {
    throw UnsupportedExprError::CreateFromExprString(e.opstr());
  }

  Result VisitInvalidNumericExpr(NumericExpr e) {
    throw InvalidNumericExprError(e);
  }

  LResult VisitInvalidLogicalExpr(LogicalExpr e) {
    throw InvalidLogicalExprError(e);
  }

  Result VisitNumericConstant(NumericConstant c) {
    return MP_DISPATCH(VisitUnhandledNumericExpr(c));
  }

  Result VisitVariable(Variable v) {
    return MP_DISPATCH(VisitUnhandledNumericExpr(v));
  }

  // Visits a unary expression or a function taking one argument.
  Result VisitUnary(UnaryExpr e) {
    return MP_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitMinus(UnaryExpr e) {
    return MP_DISPATCH(VisitUnary(e));
  }

  Result VisitAbs(UnaryExpr e) {
    return MP_DISPATCH(VisitUnary(e));
  }

  Result VisitFloor(UnaryExpr e) {
    return MP_DISPATCH(VisitUnary(e));
  }

  Result VisitCeil(UnaryExpr e) {
    return MP_DISPATCH(VisitUnary(e));
  }

  Result VisitSqrt(UnaryExpr e) {
    return MP_DISPATCH(VisitUnary(e));
  }

  Result VisitPow2(UnaryExpr e) {
    return MP_DISPATCH(VisitUnary(e));
  }

  Result VisitExp(UnaryExpr e) {
    return MP_DISPATCH(VisitUnary(e));
  }

  Result VisitLog(UnaryExpr e) {
    return MP_DISPATCH(VisitUnary(e));
  }

  Result VisitLog10(UnaryExpr e) {
    return MP_DISPATCH(VisitUnary(e));
  }

  Result VisitSin(UnaryExpr e) {
    return MP_DISPATCH(VisitUnary(e));
  }

  Result VisitSinh(UnaryExpr e) {
    return MP_DISPATCH(VisitUnary(e));
  }

  Result VisitCos(UnaryExpr e) {
    return MP_DISPATCH(VisitUnary(e));
  }

  Result VisitCosh(UnaryExpr e) {
    return MP_DISPATCH(VisitUnary(e));
  }

  Result VisitTan(UnaryExpr e) {
    return MP_DISPATCH(VisitUnary(e));
  }

  Result VisitTanh(UnaryExpr e) {
    return MP_DISPATCH(VisitUnary(e));
  }

  Result VisitAsin(UnaryExpr e) {
    return MP_DISPATCH(VisitUnary(e));
  }

  Result VisitAsinh(UnaryExpr e) {
    return MP_DISPATCH(VisitUnary(e));
  }

  Result VisitAcos(UnaryExpr e) {
    return MP_DISPATCH(VisitUnary(e));
  }

  Result VisitAcosh(UnaryExpr e) {
    return MP_DISPATCH(VisitUnary(e));
  }

  Result VisitAtan(UnaryExpr e) {
    return MP_DISPATCH(VisitUnary(e));
  }

  Result VisitAtanh(UnaryExpr e) {
    return MP_DISPATCH(VisitUnary(e));
  }

  // Visits a binary expression or a function taking two arguments.
  Result VisitBinary(BinaryExpr e) {
    return MP_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitAdd(BinaryExpr e) {
    return MP_DISPATCH(VisitBinary(e));
  }

  Result VisitSub(BinaryExpr e) {
    return MP_DISPATCH(VisitBinary(e));
  }

  Result VisitLess(BinaryExpr e) {
    return MP_DISPATCH(VisitBinary(e));
  }

  Result VisitMul(BinaryExpr e) {
    return MP_DISPATCH(VisitBinary(e));
  }

  Result VisitDiv(BinaryExpr e) {
    return MP_DISPATCH(VisitBinary(e));
  }

  Result VisitIntDiv(BinaryExpr e) {
    return MP_DISPATCH(VisitBinary(e));
  }

  Result VisitMod(BinaryExpr e) {
    return MP_DISPATCH(VisitBinary(e));
  }

  Result VisitPow(BinaryExpr e) {
    return MP_DISPATCH(VisitBinary(e));
  }

  Result VisitPowConstBase(BinaryExpr e) {
    return MP_DISPATCH(VisitBinary(e));
  }

  Result VisitPowConstExp(BinaryExpr e) {
    return MP_DISPATCH(VisitBinary(e));
  }

  // Visits a function taking two arguments.
  Result VisitBinaryFunc(BinaryExpr e) {
    return MP_DISPATCH(VisitBinary(e));
  }

  Result VisitAtan2(BinaryExpr e) {
    return MP_DISPATCH(VisitBinaryFunc(e));
  }

  Result VisitPrecision(BinaryExpr e) {
    return MP_DISPATCH(VisitBinaryFunc(e));
  }

  Result VisitRound(BinaryExpr e) {
    return MP_DISPATCH(VisitBinaryFunc(e));
  }

  Result VisitTrunc(BinaryExpr e) {
    return MP_DISPATCH(VisitBinaryFunc(e));
  }

  Result VisitIf(IfExpr e) {
    return MP_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitPiecewiseLinear(PLTerm e) {
    return MP_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitCall(CallExpr e) {
    return MP_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitVarArg(VarArgExpr e) {
    return MP_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitMin(VarArgExpr e) {
    return MP_DISPATCH(VisitVarArg(e));
  }

  Result VisitMax(VarArgExpr e) {
    return MP_DISPATCH(VisitVarArg(e));
  }

  Result VisitSum(SumExpr e) {
    return MP_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitNumberOf(NumberOfExpr e) {
    return MP_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  Result VisitCount(CountExpr e) {
    return MP_DISPATCH(VisitUnhandledNumericExpr(e));
  }

  LResult VisitLogicalConstant(LogicalConstant c) {
    return MP_DISPATCH(VisitUnhandledLogicalExpr(c));
  }

  LResult VisitNot(NotExpr e) {
    return MP_DISPATCH(VisitUnhandledLogicalExpr(e));
  }

  LResult VisitBinaryLogical(BinaryLogicalExpr e) {
    return MP_DISPATCH(VisitUnhandledLogicalExpr(e));
  }

  LResult VisitOr(BinaryLogicalExpr e) {
    return MP_DISPATCH(VisitBinaryLogical(e));
  }

  LResult VisitAnd(BinaryLogicalExpr e) {
    return MP_DISPATCH(VisitBinaryLogical(e));
  }

  LResult VisitIff(BinaryLogicalExpr e) {
    return MP_DISPATCH(VisitBinaryLogical(e));
  }

  LResult VisitRelational(RelationalExpr e) {
    return MP_DISPATCH(VisitUnhandledLogicalExpr(e));
  }

  LResult VisitLT(RelationalExpr e) {
    return MP_DISPATCH(VisitRelational(e));
  }

  LResult VisitLE(RelationalExpr e) {
    return MP_DISPATCH(VisitRelational(e));
  }

  LResult VisitEQ(RelationalExpr e) {
    return MP_DISPATCH(VisitRelational(e));
  }

  LResult VisitGE(RelationalExpr e) {
    return MP_DISPATCH(VisitRelational(e));
  }

  LResult VisitGT(RelationalExpr e) {
    return MP_DISPATCH(VisitRelational(e));
  }

  LResult VisitNE(RelationalExpr e) {
    return MP_DISPATCH(VisitRelational(e));
  }

  LResult VisitLogicalCount(LogicalCountExpr e) {
    return MP_DISPATCH(VisitUnhandledLogicalExpr(e));
  }

  LResult VisitAtLeast(LogicalCountExpr e) {
    return MP_DISPATCH(VisitLogicalCount(e));
  }

  LResult VisitAtMost(LogicalCountExpr e) {
    return MP_DISPATCH(VisitLogicalCount(e));
  }

  LResult VisitExactly(LogicalCountExpr e) {
    return MP_DISPATCH(VisitLogicalCount(e));
  }

  LResult VisitNotAtLeast(LogicalCountExpr e) {
    return MP_DISPATCH(VisitLogicalCount(e));
  }

  LResult VisitNotAtMost(LogicalCountExpr e) {
    return MP_DISPATCH(VisitLogicalCount(e));
  }

  LResult VisitNotExactly(LogicalCountExpr e) {
    return MP_DISPATCH(VisitLogicalCount(e));
  }

  LResult VisitImplication(ImplicationExpr e) {
    return MP_DISPATCH(VisitUnhandledLogicalExpr(e));
  }

  LResult VisitIteratedLogical(IteratedLogicalExpr e) {
    return MP_DISPATCH(VisitUnhandledLogicalExpr(e));
  }

  LResult VisitForAll(IteratedLogicalExpr e) {
    return MP_DISPATCH(VisitIteratedLogical(e));
  }

  LResult VisitExists(IteratedLogicalExpr e) {
    return MP_DISPATCH(VisitIteratedLogical(e));
  }

  LResult VisitAllDiff(PairwiseExpr e) {
    return MP_DISPATCH(VisitUnhandledLogicalExpr(e));
  }

  LResult VisitNotAllDiff(PairwiseExpr e) {
    return MP_DISPATCH(VisitUnhandledLogicalExpr(e));
  }
};

template <typename Impl, typename Result, typename LResult, typename ExprTypes>
Result ExprVisitor<Impl, Result, LResult, ExprTypes>::Visit(NumericExpr e) {
  // All expressions except OPNUMBEROFs and OPIFSYM are supported.
  switch (e.kind()) {
  case expr::CONSTANT:
    return MP_DISPATCH(VisitNumericConstant(
                         internal::Cast<NumericConstant>(e)));
  case expr::VARIABLE:
    return MP_DISPATCH(VisitVariable(internal::Cast<Variable>(e)));

  // Unary expressions.
  case expr::MINUS:
    return MP_DISPATCH(VisitMinus(internal::Cast<UnaryExpr>(e)));
  case expr::ABS:
    return MP_DISPATCH(VisitAbs(internal::Cast<UnaryExpr>(e)));
  case expr::FLOOR:
    return MP_DISPATCH(VisitFloor(internal::Cast<UnaryExpr>(e)));
  case expr::CEIL:
    return MP_DISPATCH(VisitCeil(internal::Cast<UnaryExpr>(e)));
  case expr::SQRT:
    return MP_DISPATCH(VisitSqrt(internal::Cast<UnaryExpr>(e)));
  case expr::POW2:
    return MP_DISPATCH(VisitPow2(internal::Cast<UnaryExpr>(e)));
  case expr::EXP:
    return MP_DISPATCH(VisitExp(internal::Cast<UnaryExpr>(e)));
  case expr::LOG:
    return MP_DISPATCH(VisitLog(internal::Cast<UnaryExpr>(e)));
  case expr::LOG10:
    return MP_DISPATCH(VisitLog10(internal::Cast<UnaryExpr>(e)));
  case expr::SIN:
    return MP_DISPATCH(VisitSin(internal::Cast<UnaryExpr>(e)));
  case expr::SINH:
    return MP_DISPATCH(VisitSinh(internal::Cast<UnaryExpr>(e)));
  case expr::COS:
    return MP_DISPATCH(VisitCos(internal::Cast<UnaryExpr>(e)));
  case expr::COSH:
    return MP_DISPATCH(VisitCosh(internal::Cast<UnaryExpr>(e)));
  case expr::TAN:
    return MP_DISPATCH(VisitTan(internal::Cast<UnaryExpr>(e)));
  case expr::TANH:
    return MP_DISPATCH(VisitTanh(internal::Cast<UnaryExpr>(e)));
  case expr::ASIN:
    return MP_DISPATCH(VisitAsin(internal::Cast<UnaryExpr>(e)));
  case expr::ASINH:
    return MP_DISPATCH(VisitAsinh(internal::Cast<UnaryExpr>(e)));
  case expr::ACOS:
    return MP_DISPATCH(VisitAcos(internal::Cast<UnaryExpr>(e)));
  case expr::ACOSH:
    return MP_DISPATCH(VisitAcosh(internal::Cast<UnaryExpr>(e)));
  case expr::ATAN:
    return MP_DISPATCH(VisitAtan(internal::Cast<UnaryExpr>(e)));
  case expr::ATANH:
    return MP_DISPATCH(VisitAtanh(internal::Cast<UnaryExpr>(e)));

  // Binary expressions.
  case expr::ADD:
    return MP_DISPATCH(VisitAdd(internal::Cast<BinaryExpr>(e)));
  case expr::SUB:
    return MP_DISPATCH(VisitSub(internal::Cast<BinaryExpr>(e)));
  case expr::LESS:
    return MP_DISPATCH(VisitLess(internal::Cast<BinaryExpr>(e)));
  case expr::MUL:
    return MP_DISPATCH(VisitMul(internal::Cast<BinaryExpr>(e)));
  case expr::DIV:
    return MP_DISPATCH(VisitDiv(internal::Cast<BinaryExpr>(e)));
  case expr::INT_DIV:
    return MP_DISPATCH(VisitIntDiv(internal::Cast<BinaryExpr>(e)));
  case expr::MOD:
    return MP_DISPATCH(VisitMod(internal::Cast<BinaryExpr>(e)));
  case expr::POW:
    return MP_DISPATCH(VisitPow(internal::Cast<BinaryExpr>(e)));
  case expr::POW_CONST_BASE:
    return MP_DISPATCH(VisitPowConstBase(internal::Cast<BinaryExpr>(e)));
  case expr::POW_CONST_EXP:
    return MP_DISPATCH(VisitPowConstExp(internal::Cast<BinaryExpr>(e)));
  case expr::ATAN2:
    return MP_DISPATCH(VisitAtan2(internal::Cast<BinaryExpr>(e)));
  case expr::PRECISION:
    return MP_DISPATCH(VisitPrecision(internal::Cast<BinaryExpr>(e)));
  case expr::ROUND:
    return MP_DISPATCH(VisitRound(internal::Cast<BinaryExpr>(e)));
  case expr::TRUNC:
    return MP_DISPATCH(VisitTrunc(internal::Cast<BinaryExpr>(e)));

  case expr::IF:
    return MP_DISPATCH(VisitIf(internal::Cast<IfExpr>(e)));
  case expr::PLTERM:
    return MP_DISPATCH(VisitPiecewiseLinear(internal::Cast<PLTerm>(e)));
  case expr::CALL:
    return MP_DISPATCH(VisitCall(internal::Cast<CallExpr>(e)));
  case expr::MIN:
    return MP_DISPATCH(VisitMin(internal::Cast<VarArgExpr>(e)));
  case expr::MAX:
    return MP_DISPATCH(VisitMax(internal::Cast<VarArgExpr>(e)));
  case expr::SUM:
    return MP_DISPATCH(VisitSum(internal::Cast<SumExpr>(e)));
  case expr::NUMBEROF:
    return MP_DISPATCH(VisitNumberOf(internal::Cast<NumberOfExpr>(e)));
  case expr::COUNT:
    return MP_DISPATCH(VisitCount(internal::Cast<CountExpr>(e)));
  default:
    // Normally this branch shouldn't be executed.
    return MP_DISPATCH(VisitInvalidNumericExpr(e));
  }
}

template <typename Impl, typename Result, typename LResult, typename ExprTypes>
LResult ExprVisitor<Impl, Result, LResult, ExprTypes>::Visit(LogicalExpr e) {
  switch (e.kind()) {
  case expr::CONSTANT:
    return MP_DISPATCH(VisitLogicalConstant(
                         internal::Cast<LogicalConstant>(e)));
  case expr::NOT:
    return MP_DISPATCH(VisitNot(internal::Cast<NotExpr>(e)));
  case expr::OR:
    return MP_DISPATCH(VisitOr(internal::Cast<BinaryLogicalExpr>(e)));
  case expr::AND:
    return MP_DISPATCH(VisitAnd(internal::Cast<BinaryLogicalExpr>(e)));
  case expr::IFF:
    return MP_DISPATCH(VisitIff(internal::Cast<BinaryLogicalExpr>(e)));
  case expr::LT:
    return MP_DISPATCH(VisitLT(internal::Cast<RelationalExpr>(e)));
  case expr::LE:
    return MP_DISPATCH(VisitLE(internal::Cast<RelationalExpr>(e)));
  case expr::EQ:
    return MP_DISPATCH(VisitEQ(internal::Cast<RelationalExpr>(e)));
  case expr::GE:
    return MP_DISPATCH(VisitGE(internal::Cast<RelationalExpr>(e)));
  case expr::GT:
    return MP_DISPATCH(VisitGT(internal::Cast<RelationalExpr>(e)));
  case expr::NE:
    return MP_DISPATCH(VisitNE(internal::Cast<RelationalExpr>(e)));
  case expr::ATLEAST:
    return MP_DISPATCH(VisitAtLeast(internal::Cast<LogicalCountExpr>(e)));
  case expr::ATMOST:
    return MP_DISPATCH(VisitAtMost(internal::Cast<LogicalCountExpr>(e)));
  case expr::EXACTLY:
    return MP_DISPATCH(VisitExactly(internal::Cast<LogicalCountExpr>(e)));
  case expr::NOT_ATLEAST:
    return MP_DISPATCH(VisitNotAtLeast(internal::Cast<LogicalCountExpr>(e)));
  case expr::NOT_ATMOST:
    return MP_DISPATCH(VisitNotAtMost(internal::Cast<LogicalCountExpr>(e)));
  case expr::NOT_EXACTLY:
    return MP_DISPATCH(VisitNotExactly(internal::Cast<LogicalCountExpr>(e)));
  case expr::IMPLICATION:
    return MP_DISPATCH(VisitImplication(internal::Cast<ImplicationExpr>(e)));
  case expr::FORALL:
    return MP_DISPATCH(VisitForAll(internal::Cast<IteratedLogicalExpr>(e)));
  case expr::EXISTS:
    return MP_DISPATCH(VisitExists(internal::Cast<IteratedLogicalExpr>(e)));
  case expr::ALLDIFF:
    return MP_DISPATCH(VisitAllDiff(internal::Cast<PairwiseExpr>(e)));
  case expr::NOT_ALLDIFF:
    return MP_DISPATCH(VisitNotAllDiff(internal::Cast<PairwiseExpr>(e)));
  default:
    // Normally this branch shouldn't be executed.
    return MP_DISPATCH(VisitInvalidLogicalExpr(e));
  }
}

// TODO: test

}  // namespace mp

#endif  // MP_EXPR_VISITOR_H_
