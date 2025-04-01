#ifndef QP2PASSES_HPP
#define QP2PASSES_HPP

#include "mp/error.h"
#include "mp/flat/qp2passes.h"
#include "mp/flat/bucketaccum.h"

#include <vector>

namespace mp {

////////////////////////////////////////////////////////////////////
/////////////////// IMPLEMENTATIONS /////////////////////
//////////////////////////////////////////////

void QP2Passes::Process(Expr expr) {
  assert(expr::Kind::SUM == expr.kind);             //for now
  top_expr_ = Cast<internal::ExprTypes::SumExpr>(expr);
  RunPass1();
  if (Pass2SeemsWorth())
    RunPass2Full();
  else
    RunPass2Buckets();
}

EExpr QP2Passes::GetResult() {
  return buckets_.ExtractSum();
}

void QP2Passes::InitPass1() {
  visitor_.InitPass1();
  is_term_qp_.resize(GetTopExpr().num_args());
}

void QP2Passes::RunPass1() {
  InitPass1();
  auto e0 = GetTopExpr();
  int term_index=0;
  for (auto term_iter=e0.begin(), term_end=e0.end();
       term_end!=term_iter; ++term_iter, ++term_index) {

  }
}

bool QP2Passes::Pass2SeemsWorth() const {
  // @todo
  return false;
}

void QP2Passes::RunPass2Full() {
  InitPass2Full();
  // ...
  ExtractPass2ResultIntoBuckets();
}

void QP2Passes::InitPass2Full() {
}

void QP2Passes::ExtractPass2ResultIntoBuckets() {
}

void QP2Passes::RunPass2Buckets() {
}

void QP2PassVisitor::InitPass1() {
  ++timestamp_;              // To distinguish active factor variables
  var_timestamps_.resize(GetQP2P().GetFlattener().num_vars_orig());
}


///////////////////////////////////////////////////////////////////

QP2PassNodeResult QP2PassVisitor::Visit(Expr e, double f) {
  double f_save = factor_;
  factor_ *= f;
  auto result = Base::Visit(e);
  factor_ = f_save;
  return result;
}


QP2PassNodeResult QP2PassVisitor::VisitUnsupported(Expr )
{ return 1000; }


QP2PassNodeResult QP2PassVisitor::VisitMinus(UnaryExpr e)
{ return Visit(e.arg(), -1.0); }

QP2PassNodeResult QP2PassVisitor::VisitAdd(BinaryExpr e) {
  auto deg1 = Visit(e.lhs());
  if (CheckDegree(deg1)) {
    auto deg2 = Visit(e.rhs());
    CheckAndMax(deg2, deg1);
  }
  return deg1;
}

QP2PassNodeResult QP2PassVisitor::VisitSub(BinaryExpr e) {
  auto deg1 = Visit(e.lhs());
  if (CheckDegree(deg1)) {
    auto deg2 = Visit(e.rhs(), -1.0);
    CheckAndMax(deg2, deg1);
  }
  return deg1;
}

QP2PassNodeResult QP2PassVisitor::VisitDiv(BinaryExpr e) {
  auto isc = IsConst(e.rhs());
  if (isc.first) {
    if (!isc.second)
      MP_RAISE("Division by 0");
    return Visit(e.lhs(), 1.0/isc.second);
  }
  return 1000;
}

QP2PassNodeResult QP2PassVisitor::VisitSum(
    internal::ExprTypes::SumExpr expr) {
  int degree_max = 0;
  for (auto term_iter=expr.begin(), term_end=expr.end();
       term_end!=term_iter; ++term_iter) {
    auto degree1 = Visit(*term_iter);
    if (!CheckAndMax(degree1, degree_max))
      return degree_max;
  }
  return degree_max;
}


QP2PassNodeResult QP2PassVisitor::VisitNumericConstant(
    NumericConstant c)
{ return (int)c.value(); }


bool QP2PassVisitor::CheckDegree(int degree) {
  return (degree <= mode_);
}

bool QP2PassVisitor::CheckAndMax(
    int degree_from_node, int &deg_max) {
  if (deg_max < degree_from_node)
    deg_max = degree_from_node;
  return CheckDegree(deg_max);
}

std::pair<bool, double> QP2PassVisitor::IsConst(Expr e) {
  // @todo Could walk in "mode 0"
  // but we hope AMPL has constant as a constant expr.
  // Or, what if presolve=0 and some variables are fixed?
  // @todo although non-critical
  // - we just process it via buckets.
  if (internal::Is<NumericConstant>(e.kind())) {
    auto ec = internal::UncheckedCast<NumericConstant>(e);
    return {true, ec.value()};
  }
  return {false, {}};
}

std::unique_ptr<BasicQP2Passes>
MakeQP2Passes(BasicProblemFlattener& flt) {
  return std::make_unique<QP2Passes>(flt);
}

}  // namespace mp

#endif // QP2PASSES_HPP
