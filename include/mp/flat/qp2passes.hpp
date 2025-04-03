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
  assert(expr::Kind::SUM == expr.kind());             //for now
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
    if (!(is_term_qp_[term_index]             // degree > 2
          = (visitor_.Visit(*term_iter) <= 2))) {
      GetBuckets().Add(                       // insert in buckets
          GetFlattener().VisitVirtual(*term_iter) );
    }
  }
}

bool QP2Passes::Pass2SeemsWorth() const {
  return              // @todo a parameter?
      0.25*visitor_.NumSourceTermsQP()
         > double(visitor_.NumQPVars())*visitor_.NumQPVars();
}

void QP2Passes::RunPass2Full() {
  InitPass2Full();
  CollectMarkedTerms();
  ExtractPass2ResultIntoBuckets();
}

void QP2Passes::InitPass2Full() {
  visitor_.InitPass2();
}

void QP2Passes::CollectMarkedTerms() {
  auto e0 = GetTopExpr();
  int term_index=0;
  for (auto term_iter=e0.begin(), term_end=e0.end();
       term_end!=term_iter; ++term_iter, ++term_index) {
    if (is_term_qp_[term_index]) {            // degree <= 2
      auto deg = visitor_.Visit(*term_iter);
      assert (deg <= 2);
    }
  }
}

void QP2Passes::ExtractPass2ResultIntoBuckets() {
  // @todo use const_term_
  // ...
}

void QP2Passes::RunPass2Buckets() {
  // ...
}

void QP2PassVisitor::InitPass1() {
  pass_ = 1;
  ++timestamp_;              // To distinguish active factor variables
  ts_lin_.resize(GetQP2P().GetFlattener().num_vars_orig());
  ts_qp_.resize(GetQP2P().GetFlattener().num_vars_orig());
  vars_lin_.clear();
  vars_qp_.clear();
  n_source_terms_qp_ = 0;
}

void QP2PassVisitor::InitPass2() {
  pass_ = 2;
  const_term_ = 0.0;
  assert (1.0 == factor_);
  coefs_lin_.resize(GetQP2P().GetFlattener().num_vars_orig());
  // 0 out necessary elements in coefs_lin_
  for (auto v: vars_lin_) {
    assert(v < coefs_lin_.size());
    coefs_lin_[v] = 0.0;
  }
  std::sort(vars_lin_.begin(), vars_lin_.end());
  coefs_qp_.clear();
  coefs_qp_.resize(NumQPVars());
  vperm_qp_.resize(GetQP2P().GetFlattener().num_vars_orig());
  std::sort(vars_qp_.begin(), vars_qp_.end());
  for (auto i = vars_qp_.size(); i--; ) {
    assert(vars_qp_[i] < vperm_qp_.size());
    vperm_qp_[vars_qp_[i]] = i;
  }
}


///////////////////////////////////////////////////////////////////

QP2PassNodeResult QP2PassVisitor::Visit(Expr e, double f) {
  auto f_save = factor_;
  factor_ *= f;
  auto result = Base::Visit(e);
  factor_ = f_save;
  return result;
}

QP2PassNodeResult QP2PassVisitor::VisitAtmostAffine(
    Expr e, AffineExpr& ae) {
  auto f_save = factor_;
  factor_ = 1.0;
  auto pae_save = p_ae_;
  p_ae_ = &ae;
  ae.clear();
  auto mode_save = mode_;
  if (mode_>1)                  // at most linear mode
    mode_ = 1;
  auto result = Base::Visit(e);
  mode_ = mode_save;
  p_ae_ = pae_save;
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

QP2PassNodeResult QP2PassVisitor::VisitMul(
    BinaryExpr expr) {
  auto iscL = IsConst(expr.lhs());
  if (iscL.first)    // const factor - stay in current mode
    return VisitFactor(expr.rhs(), iscL.second);
  auto iscR = IsConst(expr.rhs());
  if (iscR.first)
    return VisitFactor(expr.lhs(), iscR.second);
  AffineExpr aeL;
  int degL = VisitAtmostAffine(expr.lhs(), aeL);
  if (degL > mode_)
    return 1000;               // abort top-level term
  if (!degL)                   // can happen
    return VisitFactor(expr.rhs(), aeL.constant_term());
  AffineExpr aeR;
  int degR = VisitAtmostAffine(expr.rhs(), aeR);
  if (degL + degR > mode_)
    return 1000;               // abort top-level term
  if (!degR)                   // can happen
    return VisitFactor(expr.lhs(), aeR.constant_term());
  if (!ProcessAffineFactors(aeL, aeR))
    return 1000;
  return degL + degR;
}

QP2PassNodeResult QP2PassVisitor::VisitPowConstExp(
    BinaryExpr expr) {
  auto c = Cast<NumericConstant>(expr.rhs()).value();
  if (2.0==c && GetQP2P().GetFlattener().IfQuadratizePow2()) {
    return DoVisitPow2(expr.lhs());
  }
  return 1000;
}

QP2PassNodeResult QP2PassVisitor::VisitPow2(
    UnaryExpr expr) {
  if (GetQP2P().GetFlattener().IfQuadratizePow2()) {
    return DoVisitPow2(expr);
  }
  return 1000;
}

QP2PassNodeResult QP2PassVisitor::VisitPow(
    BinaryExpr expr) {
  auto iscR = IsConst(expr.rhs());
  if (iscR.first) {
    if (2.0==iscR.second &&
        GetQP2P().GetFlattener().IfQuadratizePow2()) {
      return DoVisitPow2(expr.lhs());
    }
  }
  return 1000;
}

QP2PassNodeResult QP2PassVisitor::DoVisitPow2(
    Expr expr) {
  auto isc = IsConst(expr);
  if (isc.first) {
    DoAddConst(factor_ * isc.second * isc.second);
    return 0;
  }
  AffineExpr aeL;
  int degL = VisitAtmostAffine(expr, aeL);
  if (degL*2 > mode_)
    return 1000;               // abort top-level term
  if (!degL) {                 // can happen
    DoAddConst(factor_
               * aeL.constant_term() * aeL.constant_term());
    return 0;
  }
  if (!ProcessAffineFactors(aeL, aeL))
    return 1000;
  return degL*2;
}

void QP2PassVisitor::DoAddConst(long double c) {
  if (GetPAffineExpr()) {
    GetPAffineExpr()->add_to_constant(c);
  } else {
    if (2==pass_)
      const_term_ += c;
  }
}

QP2PassNodeResult QP2PassVisitor::VisitNumericConstant(
    NumericConstant c) {
  DoAddConst(factor_ * c.value());
  return 0;                   // degree 0
}

QP2PassNodeResult QP2PassVisitor::VisitVariable(
    Reference var) {
  const auto& flt = GetQP2P().GetFlattener();
  if (flt.var_orig_lb(var.index())    // var is constant
      == flt.var_orig_ub(var.index())) {
    DoAddConst(factor_ * flt.var_orig_lb(var.index()));
    return 0;
  }                                   // else: var is var
  if (GetPAffineExpr())               // we are in a factor
    GetPAffineExpr()->add_term(factor_, var.index());
  else {
    if (1==pass_) {
      NoteLinVar(var.index());
    } else {
      assert(2==pass_);
      AddLinTerm(factor_, var.index());
    }
  }
  return 1;                   // degree 1
}

QP2PassNodeResult QP2PassVisitor::VisitCommonExpr(
    Reference ) {
  // ProblemFlattener converts them to EExpr's.
  // To reuse, we need access to auxiliary vars.
  return 1000;        // abort this top-level term
}


bool QP2PassVisitor::CheckDegree(int degree) {
  return (degree <= mode_);
}

bool QP2PassVisitor::CheckAndMax(
    int degree_from_node, int &deg_max) {
  if (deg_max < degree_from_node)
    deg_max = degree_from_node;
  return CheckDegree(deg_max);
}


bool QP2PassVisitor::ProcessAffineFactors(
    AffineExpr& aeL, AffineExpr& aeR) {
  aeL.sort_terms();
  if (&aeL != &aeR)
    aeR.sort_terms();
  if (1==pass_) {
    return
        EstimateOutmultiplication(aeL, aeR);
  }
  assert(2==pass_);
  MultiplyOut(aeL, aeR);
  return true;
}

bool QP2PassVisitor::EstimateOutmultiplication(
    AffineExpr& aeL, AffineExpr& aeR) {
  assert(1==pass_);
  if (aeL.size() && aeR.size()) {      // QP terms: estimate
    assert(2==mode_);
    if (double(aeL.size()) * aeR.size()
        > GetQP2P().GetFlattener().MultOutCard())  // cvt:multoutcard
      return false;
    for (auto v: aeL.GetBody().vars())
      NoteQPVar(v);
    for (auto v: aeR.GetBody().vars())
      NoteQPVar(v);
  }
  if (aeL.constant_term()) {           // linear terms from aeR
    for (auto v: aeR.GetBody().vars())
      NoteLinVar(v);
  }
  if (aeR.constant_term()) {           // lin terms from aeL
    for (auto v: aeR.GetBody().vars())
      NoteLinVar(v);
  }
  DoAddConst(                          // constant term
      factor_ * aeL.constant_term() * aeR.constant_term());
  return true;
}

void QP2PassVisitor::MultiplyOut(
    AffineExpr& aeL, AffineExpr& aeR) {
  assert(2==pass_);
  if (aeL.size() && aeR.size()) {      // QP terms: estimate
    assert(2==mode_);
    for (auto i1 = aeL.size(); i1--; )
      for (auto i2 = aeR.size(); i2--; )
        AddQPTerm(factor_ * aeL.coef(i1) * aeR.coef(i2),
                  aeL.var(i1), aeR.var(i2));
  }
  if (aeL.constant_term()) {           // linear terms from aeR
    for (auto i2 = aeR.size(); i2--; )
      AddLinTerm(factor_ * aeL.constant_term() * aeR.coef(i2),
                 aeR.var(i2));
  }
  if (aeR.constant_term()) {           // lin terms from aeL
    for (auto i1 = aeL.size(); i1--; )
      AddLinTerm(factor_ * aeR.constant_term() * aeL.coef(i1),
                 aeL.var(i1));
  }
  DoAddConst(                          // constant term
      factor_ * aeL.constant_term() * aeR.constant_term());
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

QP2PassNodeResult QP2PassVisitor::VisitFactor(
    Expr e, double f) {
  if (f)                         // factor non-0
    return Visit(e, f);
  return 0;
}

void QP2PassVisitor::NoteLinVar(int v) {
  assert(v < ts_lin_.size());
  if (GetTimeStamp() != ts_lin_[v]) {
    ts_lin_[v] = GetTimeStamp();
    vars_lin_.push_back(v);
  }
}
void QP2PassVisitor::NoteQPVar(int v) {
  assert(v < ts_qp_.size());
  if (GetTimeStamp() != ts_qp_[v]) {
    ts_qp_[v] = GetTimeStamp();
    vars_qp_.push_back(v);
  }
}
void QP2PassVisitor::AddLinTerm(double c, int v) {
  assert(v < coefs_lin_.size());
  coefs_lin_[v] += c;
}
void QP2PassVisitor::AddQPTerm(double c, int v1, int v2) {
  assert(v1 < vperm_qp_.size());
  assert(v2 < vperm_qp_.size());
  coefs_qp_.add_to(vperm_qp_[v1], vperm_qp_[v2], c);
}


std::unique_ptr<BasicQP2Passes>
MakeQP2Passes(BasicProblemFlattener& flt) {
  return std::make_unique<QP2Passes>(flt);
}

}  // namespace mp

#endif // QP2PASSES_HPP
