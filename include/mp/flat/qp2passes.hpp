#ifndef QP2PASSES_HPP
#define QP2PASSES_HPP

#include <vector>

#include "mp/error.h"
#include "mp/expr-linear.h"
#include "mp/flat/qp2passes.h"
#include "mp/utils-matrix.h"
#include "mp/expr-visitor.h"

namespace mp {

/// Node visitor result.
/// Polynomial degree of the node.
using QP2PassNodeResult = int;

/// Linear/QP expression visitor
class QP2PassVisitor
    : public
      ExprVisitor<
          QP2PassVisitor,
          QP2PassNodeResult> {
public:
  /// Typedef base class
  using Base = ExprVisitor<
      QP2PassVisitor,
      QP2PassNodeResult>;

  /// Construct
  QP2PassVisitor(BasicProblemFlattener& flt)
      : flattener_(flt) { }

  /// Init pass 1
  void InitPass1();
  /// Init pass 2
  void InitPass2();
  /// Call to say the 2nd pass is not going ahead
  void CancelPass2();

  /// Entry point for the top-level expression
  /// @return whether \a expr has degree <= 2.
  bool Process(Expr expr);

  /// An estimate on the number of source QP terms
  unsigned long long NumSourceTermsQP() const
  { return n_source_terms_qp_; }

  /// Number of QP vars
  auto NumQPVars() const { return vars_qp_.size(); }

  /// Extract result of Pass 2
  EExpr GetPass2Result();

  /// Call after the model is flattened.
  void Shrink();

  /// @section The below Visit... are public for CRTP

  /// Visit at most linear
  /// @param ae: the affine expression to store the subtree
  /// @note resets factor_=1.0 for the subtree
  QP2PassNodeResult VisitAtmostAffine(Expr e, AffineExpr& ae);

  /// Any unsupported expr - we process this top-level term
  /// via buckets
  QP2PassNodeResult VisitUnsupported(Expr e);

  /// Unary minus
  QP2PassNodeResult VisitMinus(UnaryExpr e);
  /// Add
  QP2PassNodeResult VisitAdd(BinaryExpr e);
  /// Sub
  QP2PassNodeResult VisitSub(BinaryExpr e);
  /// Div
  QP2PassNodeResult VisitDiv(BinaryExpr e);
  /// Visit SumExpr
  QP2PassNodeResult VisitSum(internal::ExprTypes::SumExpr expr);
  /// Mul
  QP2PassNodeResult VisitMul(BinaryExpr e);
  /// PowConstExp
  QP2PassNodeResult VisitPowConstExp(BinaryExpr e);
  /// Pow2
  QP2PassNodeResult VisitPow2(UnaryExpr e);
  /// Pow
  QP2PassNodeResult VisitPow(BinaryExpr e);

  /// Constant
  QP2PassNodeResult VisitNumericConstant(NumericConstant );
  /// Variable
  QP2PassNodeResult VisitVariable(Reference );
  /// Defined variable
  QP2PassNodeResult VisitCommonExpr(Reference );

protected:
  const BasicProblemFlattener& GetFlattener() const
  { return flattener_; }
  BasicProblemFlattener& GetFlattener()
  { return flattener_; }

  /// Reuse dispatcher
  using Base::Visit;
  /// Overlaod the dispatching Visit()
  QP2PassNodeResult Visit(Expr e);
  /// Overload the dispatching Visit()
  /// to multiply the constant factor for its subtree
  QP2PassNodeResult Visit(Expr e, double f);

  /// Variable index
  QP2PassNodeResult VisitVariableIndex(int i);
  /// Overload VisitVariableIndex()
  /// to multiply by the constant factor
  QP2PassNodeResult VisitVariableIndex(int i, double f);

  bool CheckDegree(int degree);
  bool CheckAndMax(int degree_from_node, int& deg_max);

  QP2PassNodeResult DoVisitPow2(Expr e);
  /// Add to the top-level or current ae's constant term
  void DoAddConst(long double c);

  /// @todo consider binary terms (cvt:prod)
  bool ProcessAffineFactors(AffineExpr& aeL, AffineExpr& aeR);
  bool EstimateOutmultiplication(AffineExpr& aeL, AffineExpr& aeR);
  void MultiplyOut(AffineExpr& aeL, AffineExpr& aeR);

  std::pair<bool, double> IsConst(Expr );

  const AffineExpr* GetPAffineExpr() const { return p_ae_; }
  AffineExpr* GetPAffineExpr() { return p_ae_; }

  /// Note top-level var
  void NoteLinVar(int v);
  void NoteQPVar(int v);
  /// Add top-level lin/qp term
  void AddLinTerm(double c, int v);
  void AddQPTerm(double c, int v1, int v2);

  unsigned int GetTimeStamp() const { return timestamp_; }

private:
  BasicProblemFlattener& flattener_;

  int pass_ {};         // 1 or 2
  int mode_ {};         // 2: root (accepts QP), 1: only linear
  bool pass1_5_{};

  long double factor_ {1.0};
  AffineExpr* p_ae_{};  // pointer to chosen ae being filled

  unsigned int timestamp_ {0};
  std::vector<unsigned int> ts_lin_, ts_qp_;
  SmallVec<int, 64> vars_lin_, vars_qp_;
  unsigned long long n_source_terms_qp_ {};

  /// Pass 2
  long double const_term_ {};   // the top-level constant term
  std::vector<double> coefs_lin_dense_;

  TMatrix<double, 16> coefs_qp_;
  std::vector<int> vperm_qp_;   // inverse of vars_qp_
};


////////////////////////////////////////////////////////////////////
/////////////////// IMPLEMENTATIONS /////////////////////
//////////////////////////////////////////////

#ifdef DEBUG_QP2PASSES
/// Write algebraic expression (linear + non-linear.)
template <typename ExprTypes,
         typename LinearExpr, typename NumericExpr,
         typename Namer>
void WriteExpr(fmt::Writer &w, const LinearExpr &linear,
               NumericExpr nonlinear, Namer);
#endif

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
  return std::move(result_);
}

void QP2Passes::InitPass1() {
  visitor_.InitPass1();
  ResizeWithExtraCapacity(is_term_qp_, GetTopExpr().num_args());
  n_qp_terms_ = 0;
}

void QP2Passes::RunPass1() {
  InitPass1();
  auto e0 = GetTopExpr();
  int term_index=0;
  for (auto term_iter=e0.begin(), term_end=e0.end();
       term_end!=term_iter; ++term_iter, ++term_index) {
    if ((is_term_qp_[term_index]             // degree > 2
          = visitor_.Process(*term_iter))) {
      ++n_qp_terms_;
    }
  }
}

bool QP2Passes::Pass2SeemsWorth() const {
  return              // @todo a parameter?
      n_qp_terms_
      && visitor_.NumSourceTermsQP()
             > 0.25*visitor_.NumQPVars()
                   *visitor_.NumQPVars();
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
      auto deg2 = visitor_.Process(*term_iter);
      assert (deg2);
      MP_UNUSED(deg2);
    }
  }
}

void QP2Passes::ExtractPass2ResultIntoBuckets() {
  if ((int)n_qp_terms_ < GetTopExpr().num_args()) { // need buckets
    auto e0 = GetTopExpr();
    BucketAccumulator<EExpr> buckets
        {e0.num_args() - n_qp_terms_ + 1};
    buckets.Add(visitor_.GetPass2Result());     // before VisitVirtual()
    int term_index=0;
    for (auto term_iter=e0.begin(), term_end=e0.end();
         term_end!=term_iter; ++term_iter, ++term_index) {
      if (!is_term_qp_[term_index]) {           // degree > 2
        buckets.Add(
            GetFlattener().VisitVirtual(*term_iter) );
      }
    }
    result_ = buckets.ExtractSum();
  } else
    result_ = visitor_.GetPass2Result();
}

void QP2Passes::RunPass2Buckets() {
  visitor_.CancelPass2();
  auto e0 = GetTopExpr();
  BucketAccumulator<EExpr> buckets(e0.num_args());
  for (auto term_iter=e0.begin(), term_end=e0.end();
       term_end!=term_iter; ++term_iter) {
    buckets.Add(
        GetFlattener().VisitVirtual(*term_iter) );
  }
  result_ = buckets.ExtractSum();
}


////////////////////////////////////////////////////////////////
/// QP2PassVisitor methods
////////////////////////////////////////////////////////////////

bool QP2PassVisitor::Process(Expr expr) {
  assert(1==pass_ || 2==pass_);
  if (1==pass_) {
    pass1_5_ = false;
    auto deg1_0 = Visit(expr);
    if (deg1_0 > 2)
      return false;
    pass1_5_ = true;
    auto deg1_5 = Visit(expr);
    assert(deg1_0 == deg1_5);
    MP_UNUSED(deg1_5);
    return true;
  }
  return Visit(expr) <= 2;
}

void QP2PassVisitor::InitPass1() {
  mode_ = 2;                 // quadratics
  assert(0 == pass_);        // We are not inside another run
  assert(!p_ae_);
  pass_ = 1;
  ++timestamp_;              // To distinguish active factor variables
  ResizeWithExtraCapacity(ts_lin_, GetFlattener().num_vars_flat());
  ResizeWithExtraCapacity(ts_qp_, GetFlattener().num_vars_flat());
  vars_lin_.clear();
  vars_qp_.clear();
  n_source_terms_qp_ = 0;
}

void QP2PassVisitor::InitPass2() {
  pass_ = 2;
  const_term_ = 0.0;
  assert (1.0 == factor_);
  ResizeWithExtraCapacity(coefs_lin_dense_, GetFlattener().num_vars_flat());
  // 0 out necessary elements in coefs_lin_
  for (auto v: vars_lin_) {
    assert(v < (int)coefs_lin_dense_.size());
    coefs_lin_dense_[v] = 0.0;
  }
  std::sort(vars_lin_.begin(), vars_lin_.end());
  coefs_qp_.clear();
  ResizeWithExtraCapacity(coefs_qp_, NumQPVars());
  ResizeWithExtraCapacity(vperm_qp_, GetFlattener().num_vars_flat());
  std::sort(vars_qp_.begin(), vars_qp_.end());
  for (auto i = vars_qp_.size(); i--; ) {
    assert(vars_qp_[i] < (int)vperm_qp_.size());
    vperm_qp_[vars_qp_[i]] = i;
  }
}

void QP2PassVisitor::CancelPass2() {
  assert(1 == pass_);
  pass_ = 0;
}


///////////////////////////////////////////////////////////////////

QP2PassNodeResult QP2PassVisitor::Visit(Expr e) {
  assert(mode_);
  assert(pass_);
#ifdef DEBUG_QP2PASSES
  static int depth=0;
  ++depth;
  {
    fmt::MemoryWriter wrt;
    wrt << "QP/Pass " << pass_ << ": ";
    WriteExpr<typename internal::ExprTypes>(
        wrt, LinearExpr{}, Cast<NumericExpr>(e),
        GetFlattener().GetOrigProblem().GetVarNamer());
    fmt::print("{:{}}{}\n", "", depth*2, wrt.str());
  }
#endif
  auto result = Base::Visit(e);
#ifdef DEBUG_QP2PASSES
  {
    fmt::MemoryWriter wrt;
    wrt << "DONE QP/Pass " << pass_ << " on: ";
    WriteExpr<typename internal::ExprTypes>(
        wrt, LinearExpr{}, Cast<NumericExpr>(e),
        GetFlattener().GetOrigProblem().GetVarNamer());
    fmt::print("{:{}}{};  result = {}\n",
               "", depth*2, wrt.str(), result);
  }
  --depth;
#endif
  return result;
}

QP2PassNodeResult QP2PassVisitor::Visit(Expr e, double f) {
  if (f) {
    auto f_save = factor_;
    factor_ *= f;
    auto result = Visit(e);
    factor_ = f_save;
    return result;
  }
  return 0;
}

QP2PassNodeResult QP2PassVisitor::VisitVariableIndex(int i, double f) {
  if (f) {
    auto f_save = factor_;
    factor_ *= f;
    auto result = VisitVariableIndex(i);
    factor_ = f_save;
    return result;
  }
  return 0;
}

QP2PassNodeResult QP2PassVisitor::VisitVariableIndex(int index) {
  const auto& flt = GetFlattener();
  if (flt.var_lb_flat(index)    // var is constant
      >= flt.var_ub_flat(index)) {
    DoAddConst(factor_ * flt.var_lb_flat(index));
    return 0;
  }                                   // else: var is var
  if (GetPAffineExpr())               // we are in a factor
    GetPAffineExpr()->add_term(factor_, index);
  else {
    if (1==pass_) {
      if (pass1_5_)
        NoteLinVar(index);
    } else {
      assert(2==pass_);
      AddLinTerm(factor_, index);
    }
  }
  return 1;                   // degree 1
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
  auto result = Visit(e);
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
    return Visit(expr.rhs(), iscL.second);
  auto iscR = IsConst(expr.rhs());
  if (iscR.first)
    return Visit(expr.lhs(), iscR.second);
  AffineExpr aeL;
  int degL = VisitAtmostAffine(expr.lhs(), aeL);
  if (degL > mode_)
    return 1000;               // abort top-level term
  if (!degL)                   // can happen
    return Visit(expr.rhs(), aeL.constant_term());
  AffineExpr aeR;
  int degR = VisitAtmostAffine(expr.rhs(), aeR);
  if (degL + degR > mode_)
    return 1000;               // abort top-level term
  if (!degR)                   // can happen
    return Visit(expr.lhs(), aeR.constant_term());
  if (!ProcessAffineFactors(aeL, aeR))
    return 1000;
  return degL + degR;
}

QP2PassNodeResult QP2PassVisitor::VisitPowConstExp(
    BinaryExpr expr) {
  auto c = Cast<NumericConstant>(expr.rhs()).value();
  if (2.0==c && GetFlattener().IfQuadratizePow2()) {  // #276
    return DoVisitPow2(expr.lhs());
  }
  return 1000;
}

QP2PassNodeResult QP2PassVisitor::VisitPow2(
    UnaryExpr expr) {
  if (GetFlattener().IfQuadratizePow2()) {
    return DoVisitPow2(expr.arg());
  }
  return 1000;
}

QP2PassNodeResult QP2PassVisitor::VisitPow(
    BinaryExpr expr) {
  auto iscR = IsConst(expr.rhs());
  if (iscR.first) {
    if (2.0==iscR.second &&
        GetFlattener().IfQuadratizePow2()) {
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
    Reference var)
{ return VisitVariableIndex(var.index()); }

QP2PassNodeResult QP2PassVisitor::VisitCommonExpr(
    Reference r) {
  if (!GetFlattener().CanElimCommonExpr(r)) {
    return 1000;        // abandon
    /// @todo Use the explicit result variable of the common expr.
    /// @todo this requires re-entrable Visitor.
    // const auto ee = GetFlattener().VisitVirtual(r);
    // // this asserts that Flattener converted the common expr
    // // to a single result variable (it is not substituted):
    // return VisitVariableIndex(ee.get_representing_variable());
  }
  // Treat a DV like a sum
  auto ce = GetFlattener().GetCommonExpr(r.index());
  const auto& le = ce.linear_expr();
  int degree_max = 0;
  for (const auto& term: le) {
    auto degree1 = VisitVariableIndex(term.var_index(), term.coef());
    if (!CheckAndMax(degree1, degree_max))
      return degree_max;
  }
  if (auto e = ce.nonlinear_expr()) {
    auto degree1 = Visit(e);
    if (!CheckAndMax(degree1, degree_max))
      return degree_max;
  }
  return degree_max;
}


EExpr QP2PassVisitor::GetPass2Result() {
  assert(2==pass_);
  pass_ = 0;          // to indicate the run is finished
  EExpr result;
  result.constant_term(const_term_);
  result.GetLinTerms().reserve(vars_lin_.size());
  for (auto vl: vars_lin_)
    if (auto c = coefs_lin_dense_[vl])
      result.GetLinTerms().add_term(c, vl);
  assert(result.GetLinTerms().is_sorted());
  for (std::size_t i=0; i<vars_qp_.size(); ++i)
    for (std::size_t j=i; j<vars_qp_.size(); ++j)
      if (auto c = coefs_qp_(i, j))
        result.GetQPTerms().add_term(
            c, vars_qp_[i], vars_qp_[j]);
  assert(result.GetQPTerms().is_sorted());
  return result;
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
    if (!pass1_5_) {
      if (double(aeL.size()) * aeR.size()
          > GetFlattener().MultOutCard())  // cvt:multoutcard
        return false;
      if (GetFlattener().LogicalizeProd2BinVars()
          && GetFlattener().IsBinaryOrNegatedBinary(aeL)
          && GetFlattener().IsBinaryOrNegatedBinary(aeR))
        return false;
    } else {
      for (auto v: aeL.GetBody().vars())
        NoteQPVar(v);
      for (auto v: aeR.GetBody().vars())
        NoteQPVar(v);
      n_source_terms_qp_
          += ((unsigned long long)(aeL.size())) * aeR.size();
    }
  }
  if (pass1_5_) {
    if (aeL.constant_term()) {           // linear terms from aeR
      for (auto v: aeR.GetBody().vars())
        NoteLinVar(v);
    }
    if (aeR.constant_term()) {           // lin terms from aeL
      for (auto v: aeL.GetBody().vars())
        NoteLinVar(v);
    }
  }
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

void QP2PassVisitor::NoteLinVar(int v) {
  assert(1 == pass_);
  assert(pass1_5_);
  assert(v < (int)ts_lin_.size());
  if (GetTimeStamp() != ts_lin_[v]) {
    ts_lin_[v] = GetTimeStamp();
    vars_lin_.push_back(v);
  }
}
void QP2PassVisitor::NoteQPVar(int v) {
  assert(1 == pass_);
  assert(pass1_5_);
  assert(v < (int)ts_qp_.size());
  if (GetTimeStamp() != ts_qp_[v]) {
    ts_qp_[v] = GetTimeStamp();
    vars_qp_.push_back(v);
  }
}
void QP2PassVisitor::AddLinTerm(double c, int v) {
  assert(2 == pass_);
  assert(GetTimeStamp() == ts_lin_[v]);
  assert(v < (int)coefs_lin_dense_.size());
  coefs_lin_dense_[v] += c;
}
void QP2PassVisitor::AddQPTerm(double c, int v1, int v2) {
  assert(2 == pass_);
  assert(GetTimeStamp() == ts_qp_[v1]);
  assert(GetTimeStamp() == ts_qp_[v2]);
  assert(v1 < (int)vperm_qp_.size());
  assert(v2 < (int)vperm_qp_.size());
  coefs_qp_.add_to(vperm_qp_[v1], vperm_qp_[v2], c);
}

void QP2PassVisitor::Shrink() {
  coefs_qp_.clear();
  coefs_qp_.shrink_to_fit();
}


QP2PassVisitor*
MakeQP2PassVisitor(BasicProblemFlattener& flt) {
  return new QP2PassVisitor(flt);
}

void DeleteQP2PassVisitor(QP2PassVisitor* pv)
{ delete pv; }

void Shrink(QP2PassVisitor& v) { v.Shrink(); }

}  // namespace mp

#endif // QP2PASSES_HPP
