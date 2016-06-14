/*
 Stochastic programming support

 Copyright (C) 2016 AMPL Optimization Inc

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

#include "sp.h"
#include "mp/expr-visitor.h"

#include <cstring>    // std::strcmp
#include <algorithm>  // std::max

namespace mp {
namespace {

// Returns the name of a logical constraint.
inline std::string lcon_name(int index) {
  return fmt::format("_slogcon[{}]", index + 1);
}

template <typename Impl>
class RandomConstExprExtractor: public ExprVisitor<Impl, double> {
 private:
  int scenario_;

 protected:
  const SPAdapter &sp_;

  int core_var_index(int var_index) const {
    return sp_.core_var_index(var_index);
  }

 public:
  explicit RandomConstExprExtractor(const SPAdapter &sp, int scenario):
    scenario_(scenario), sp_(sp) {}

  double VisitNumericConstant(NumericConstant n) { return n.value(); }

  double VisitVariable(Reference v) {
    if (auto rv = sp_.random_var(v.index()))
      return rv.realization(scenario_);
    return mp::ExprVisitor<Impl, double>::VisitVariable(v);
  }
};

// Extracts an affine expression for a single scenario from an expression
// containing random variables.
template <typename Impl>
class BasicRandomAffineExprExtractor: public RandomConstExprExtractor<Impl> {
 private:
  LinearExpr linear_;
  double coef_;

  typedef RandomConstExprExtractor<Impl> Base;

  void AddTerm(int var_index, double coef) {
    linear_.AddTerm(var_index, coef);
  }

  double DoAddTerm(Expr coef, Expr var) {
    AddTerm(Cast<Reference>(var).index(), coef_ * Base(*this).Visit(coef));
    return 0;
  }

 public:
  BasicRandomAffineExprExtractor(const SPAdapter &sp, int scenario):
    Base(sp, scenario), coef_(1) {}

  const LinearExpr &linear_expr() const { return linear_; }

  double VisitUnary(UnaryExpr e) {
    if (e.kind() != expr::MINUS)
      return Base::VisitUnary(e);
    double saved_coef = coef_;
    coef_ = -coef_;
    double result = -this->Visit(e.arg());
    coef_ = saved_coef;
    return result;
  }

  double VisitBinary(BinaryExpr e);

  // TODO
};

template <typename Impl>
double BasicRandomAffineExprExtractor<Impl>::VisitBinary(BinaryExpr e) {
  switch (e.kind()) {
  case expr::MUL:
    if (e.rhs().kind() == expr::VARIABLE)
      return DoAddTerm(e.lhs(), e.rhs());
    if (e.lhs().kind() == expr::VARIABLE)
      return DoAddTerm(e.rhs(), e.lhs());
    throw UnsupportedError("nonlinear expression");
    break;
  default:
    return Base::VisitBinary(e);
  }
}

class RandomAffineExprExtractor:
    public BasicRandomAffineExprExtractor<RandomAffineExprExtractor> {
 public:
  RandomAffineExprExtractor(const SPAdapter &sp, int scenario):
    BasicRandomAffineExprExtractor<RandomAffineExprExtractor>(sp, scenario) {}
};

// Collects the list of variables that appear in a nonlinear expression.
class VariableCollector:
    public BasicRandomAffineExprExtractor<VariableCollector> {
 private:
  SparseMatrix<int> &vars_in_nonlinear_;
  std::vector<bool> visited_vars_;
  int con_index_;

 public:
  VariableCollector(const SPAdapter &sp, SparseMatrix<int> &vars_in_nonlinear):
    BasicRandomAffineExprExtractor<VariableCollector>(sp, 0),
    vars_in_nonlinear_(vars_in_nonlinear), visited_vars_(sp.num_vars()),
    con_index_(0) {}

  void Collect();

  double VisitNumericConstant(NumericConstant) { return 0; }

  void AddTerm(int var_index, double) {
    int core_var_index = this->core_var_index(var_index);
    if (core_var_index >= 0 && !visited_vars_[core_var_index]) {
      vars_in_nonlinear_.add_index(core_var_index);
      visited_vars_[core_var_index] = true;
    }
  }

  void VisitMultiplied(double, NumericExpr expr) { this->Visit(expr); }
};

void VariableCollector::Collect() {
  int num_stage1_cons = sp_.stage(0).num_cons();
  int num_stage2_cons = sp_.stage(1).num_cons();
  vars_in_nonlinear_.resize_major(num_stage2_cons);
  if (num_stage2_cons == 0) return;
  for (;;) {
    auto expr = sp_.con(num_stage1_cons + con_index_).nonlinear_expr();
    if (expr)
      Visit(expr);
    ++con_index_;
    vars_in_nonlinear_.start(con_index_) = vars_in_nonlinear_.num_elements();
    if (con_index_ == num_stage2_cons) break;
    for (int i = vars_in_nonlinear_.start(con_index_ - 1),
         n = vars_in_nonlinear_.start(con_index_); i < n; ++i) {
      visited_vars_[vars_in_nonlinear_.index(i)] = false;
    }
  }
}

class NullSuffix {
 public:
  int value(int) const { return 0; }
};
}  // namespace

namespace internal {

// Extracts an affine expression from a nonlinear one.
class AffineExprExtractor: public ExprVisitor<AffineExprExtractor, void> {
 private:
  const SPAdapter &sp_;
  std::vector<double> &coefs_;
  double constant_;
  double coef_;

  void AddTerm(int var_index, double coef) {
    coefs_[sp_.var_orig2core_[var_index]] += coef * coef_;
  }

  void VisitMultiplied(double multiplier, NumericExpr expr) {
    double saved_coef = coef_;
    coef_ *= multiplier;
    Visit(expr);
    coef_ = saved_coef;
  }

 public:
  AffineExprExtractor(const SPAdapter &sp, std::vector<double> &coefs):
    sp_(sp), coefs_(coefs), constant_(0), coef_(1){}

  double constant() const { return constant_; }

  void VisitNumericConstant(NumericConstant n) {
    constant_ += coef_ * n.value();
  }

  void VisitVariable(Reference v) { AddTerm(v.index(), 1); }

  void VisitCommonExpr(Reference e);

  void VisitBinary(BinaryExpr e);

  void VisitSum(IteratedExpr e) {
    for (auto arg: e)
      Visit(arg);
  }
};

void AffineExprExtractor::VisitCommonExpr(Reference e) {
  auto common_expr = sp_.problem_.common_expr(e.index());
  for (auto term: common_expr.linear_expr())
    AddTerm(term.var_index(), term.coef());
  Visit(common_expr.nonlinear_expr());
}

void AffineExprExtractor::VisitBinary(BinaryExpr e) {
  switch (e.kind()) {
  case expr::ADD:
    Visit(e.lhs());
    Visit(e.rhs());
    break;
  case expr::SUB:
    Visit(e.lhs());
    VisitMultiplied(-1, e.rhs());
    break;
  case expr::MUL:
    if (auto n = Cast<NumericConstant>(e.lhs()))
      return VisitMultiplied(n.value(), e.rhs());
    if (auto n = Cast<NumericConstant>(e.rhs()))
      return VisitMultiplied(n.value(), e.lhs());
    throw MakeUnsupportedError("nonlinear expression");
    break;
  default:
    ExprVisitor<AffineExprExtractor, void>::VisitBinary(e);
  }
}
}  // namespace internal

void SPAdapter::GetRealizations(int con_index, CallExpr random,
                                int &arg_index) {
  // Add a random variable.
  auto &rv = rvs_[con_index];
  auto var = Cast<Reference>(random.arg(arg_index));
  int var_index = var.index();
  random_vars_.push_back(RandomVarInfo(
                           var_index, con_index, rv.num_elements()));
  int &core_var_index = var_orig2core_[var_index];
  if (core_var_index != 0) {
    throw Error("{}: variable {} used in multiple random vectors",
                lcon_name(con_index), var_index);
  }
  // Mark variable as random.
  core_var_index = -static_cast<int>(random_vars_.size());

  // Get realizatons.
  int start_arg_index = ++arg_index;
  for (; ; ++arg_index) {
    if (arg_index >= random.num_args()) break;
    auto constant = Cast<NumericConstant>(random.arg(arg_index));
    if (!constant) break;
    rv.Add(constant.value());
  }
  int num_realizations = arg_index - start_arg_index;
  if (rv.num_realizations() == 0) {
    rv.set_num_realizations(num_realizations);
  } else if (rv.num_realizations() != num_realizations) {
    throw Error("{}: inconsistent number of realizations",
                lcon_name(con_index));
  }
}

void SPAdapter::GetRandomVectors(const Problem &p) {
  var_orig2core_.resize(p.num_vars());
  int num_logical_cons = p.num_logical_cons();
  if (num_logical_cons == 0)
    return;
  rvs_.resize(num_logical_cons);
  for (int con_index = 0; con_index < num_logical_cons; ++con_index) {
    auto expr = p.logical_con(con_index).expr();
    if (expr.kind() != expr::NE)
      throw MakeUnsupportedError("logical constraint");
    auto relational = Cast<RelationalExpr>(expr);
    auto call = Cast<CallExpr>(relational.lhs());
    if (!call || call.function() != random_ || !IsZero(relational.rhs()))
      throw MakeUnsupportedError("logical constraint");
    int num_args = call.num_args();
    if (num_args == 0)
      continue;
    auto &rv = rvs_[con_index];
    int arg_index = 0;
    // Get probabilities.
    for (; arg_index < num_args; ++arg_index) {
      auto constant = Cast<NumericConstant>(call.arg(arg_index));
      if (!constant) break;
      double prob = constant.value();
      if (prob < 0 || prob > 1)
        throw Error("{}: invalid probability {}", lcon_name(con_index), prob);
      rv.AddProbability(prob);
    }
    // Get realizations.
    while (arg_index < num_args) {
      auto arg = call.arg(arg_index);
      if (arg.kind() != expr::VARIABLE)
        throw Error("{}: expected variable or constant", lcon_name(con_index));
      GetRealizations(con_index, call, arg_index);
    }
  }
}

template <typename Suffix>
void SPAdapter::GetVarStages(Suffix stage_suffix) {
  // Count the number of variables in each stage and temporarily store variable
  // stages in var_orig2core_.
  for (int i = 0, n = problem_.num_vars(); i < n; ++i) {
    if (var_orig2core_[i])
      continue;  // Skip random variables.
    int stage = std::max(stage_suffix.value(i) - 1, 0);
    var_orig2core_[i] = stage;
    if (stage >= static_cast<int>(num_stage_vars_.size()))
      num_stage_vars_.resize(stage + 1);
    ++num_stage_vars_[stage];
  }
}

// Update stages of all constraints containing this variable.
void SPAdapter::UpdateConStages(int var_index, int stage) {
  for (int k = problem_.col_start(var_index),
       end = problem_.col_start(var_index + 1); k != end; ++k) {
    int &con_stage = con_orig2core_[problem_.row_index(k)];
    con_stage = std::max(stage, con_stage);
  }
}

void SPAdapter::ProcessObjs() {
  int num_objs = problem_.num_objs();
  if (num_objs == 0)
    return;
  int num_stage1_vars = num_stage_vars_[0];
  nonlinear_objs_.reserve(num_objs);
  for (auto obj: problem_.objs())
    nonlinear_objs_.push_back(obj.nonlinear_expr());
  // Strip expectation from the first objective.
  auto obj = problem_.obj(0);
  int num_vars = this->num_vars();
  std::vector<double> core_obj(num_vars);
  for (auto term: obj.linear_expr()) {
    int core_var_index = var_orig2core_[term.var_index()];
    if (core_var_index >= num_stage1_vars && term.coef() != 0)
      throw Error("second-stage variable outside of expectation in objective");
    core_obj[core_var_index] = term.coef();
  }
  auto obj_expr = obj.nonlinear_expr();
  if (obj_expr) {
    auto call = Cast<CallExpr>(obj_expr);
    if (call && std::strcmp(call.function().name(), "expectation") == 0) {
      if (call.num_args() != 1 || !(obj_expr = Cast<NumericExpr>(call.arg(0))))
        throw Error("invalid arguments to expectation");
      nonlinear_objs_[0] = mp::NumericExpr();
    }
    internal::AffineExprExtractor extractor(*this, core_obj);
    extractor.Visit(obj_expr);
    if (auto constant = extractor.constant())
      nonlinear_objs_[0] = factory_.MakeNumericConstant(constant);
  }
  for (int i = 0; i < num_vars; ++i) {
    if (auto obj = core_obj[i])
      linear_obj_.AddTerm(i, obj);
  }
}

void SPAdapter::ProcessCons() {
  int num_cons = problem_.num_algebraic_cons();
  num_stage_cons_.resize(num_stages_);
  if (num_cons == 0)
    return;
  // Compute stage of each constraint as a maximum of stages of variables
  // in it and temporarily store stages in con_orig2core_.
  con_orig2core_.resize(num_cons);
  int core_var_index = num_stage_vars_[0];
  for (int stage = 1; stage < num_stages_; ++stage) {
    for (int n = core_var_index + num_stage_vars_[stage];
         core_var_index < n; ++core_var_index) {
      UpdateConStages(var_core2orig_[core_var_index], stage);
    }
  }
  // Constraints containing random variables should be at least in the
  // second stage.
  for (const auto &rv: random_vars_)
    UpdateConStages(rv.var_index, 1);
  for (int i = 0; i < num_cons; ++i) {
    int &stage = con_orig2core_[i];
    if (problem_.algebraic_con(i).nonlinear_expr())
      stage = std::max(con_orig2core_[i], 1);
    if (stage >= static_cast<int>(num_stage_cons_.size()))
      num_stage_cons_.resize(stage + 1);
    ++num_stage_cons_[stage];
  }
  if (num_stage_cons_.size() > num_stage_vars_.size()) {
    num_stages_ = num_stage_cons_.size();
    num_stage_vars_.resize(num_stages_);
  }
  // Reorder constraints by stages.
  con_core2orig_.resize(num_cons);
  std::vector<int> stage_offsets(num_stages_);
  for (int i = 1; i < num_stages_; ++i)
    stage_offsets[i] = stage_offsets[i - 1] + num_stage_cons_[i - 1];
  for (int i = 0; i < num_cons; ++i) {
    int stage = con_orig2core_[i];
    int core_con_index = stage_offsets[stage]++;
    con_core2orig_[core_con_index] = i;
    con_orig2core_[i] = core_con_index++;
  }
}

void SPAdapter::ExtractRandomTerms() {
  int num_stage1_cons = num_stage_cons_[0];
  int num_stage2_cons = num_stage_cons_[1];

  // A matrix containing linear terms involving random variables.
  // The major dimension is equal to the number of second-stage constraints.
  linear_random_.resize_major(num_stage2_cons);

  // Count random variables in second-stage constraints.
  for (const auto &rv: random_vars_) {
    for (int k = problem_.col_start(rv.var_index),
         end = problem_.col_start(rv.var_index + 1); k != end; ++k) {
       int core_con_index = con_orig2core_[problem_.row_index(k)];
       assert(core_con_index >= num_stage1_cons);
       ++linear_random_.start(core_con_index - num_stage1_cons + 1);
    }
  }
  // Acummulate counts to get vector starts.
  int start = 0;
  for (int i = 1; i <= num_stage2_cons; ++i) {
    int next = start + linear_random_.start(i);
    linear_random_.start(i) = start;
    start = next;
  }

  // Map second-stage constraints to random variables that appear linearly
  // in them.
  linear_random_.resize_elements(start);
  for (const auto &rv: random_vars_) {
    for (int k = problem_.col_start(rv.var_index),
         end = problem_.col_start(rv.var_index + 1); k != end; ++k) {
       int core_con_index = con_orig2core_[problem_.row_index(k)];
       assert(core_con_index >= num_stage1_cons);
       int index = core_con_index - num_stage1_cons + 1;
       int element_index = linear_random_.start(index)++;
       linear_random_.index(element_index) = rv.var_index;
       linear_random_.coef(element_index) = -problem_.value(k);
    }
  }

  // Get variables that appear in nonlinear parts of constraint expressions.
  VariableCollector collector(*this, vars_in_nonlinear_);
  collector.Collect();
}

SPAdapter::SPAdapter(const ColProblem &p): problem_(p), num_stages_(1) {
  // Find the random function.
  for (int i = 0, n = p.num_functions(); i < n; ++i) {
    auto function = p.function(i);
    if (std::strcmp(function.name(), "random") == 0) {
      random_ = function;
      break;
    }
  }

  GetRandomVectors(p);

  // Reorder variables by stages.
  if (IntSuffix stage_suffix = p.suffixes(suf::VAR).Find<int>("stage"))
    GetVarStages(stage_suffix);
  else
    GetVarStages(NullSuffix());
  if (num_stage_vars_.empty())
    num_stage_vars_.resize(1);
  num_stages_ = static_cast<int>(num_stage_vars_.size());
  int num_vars = p.num_vars();
  var_core2orig_.resize(num_vars - random_vars_.size());
  std::vector<int> stage_offsets(num_stages_);
  for (int i = 1; i < num_stages_; ++i)
    stage_offsets[i] = stage_offsets[i - 1] + num_stage_vars_[i - 1];
  for (int i = 0; i < num_vars; ++i) {
    int stage = var_orig2core_[i];
    if (stage < 0)
      continue;  // Skip random variables.
    int core_var_index = stage_offsets[stage]++;
    var_core2orig_[core_var_index] = i;
    var_orig2core_[i] = core_var_index;
  }

  ProcessObjs();
  ProcessCons();
  if (num_stages_ > 1)
    ExtractRandomTerms();
}

void SPAdapter::GetScenario(Scenario &s, int scenario_index) const {
  int num_stage1_cons = num_stage_cons_[0];
  int num_stage2_cons = num_stage_cons_[1];
  s.rhs_offsets_.assign(num_stage2_cons, 0);
  for (int stage2_con = 0; stage2_con < num_stage2_cons; ++stage2_con) {
    for (int i = linear_random_.start(stage2_con),
         end = linear_random_.start(stage2_con + 1); i < end; ++i) {
      auto random_var = this->random_var(linear_random_.index(i));
      assert(random_var);
      s.rhs_offsets_[stage2_con] +=
          linear_random_.coef(i) * random_var.realization(scenario_index);
    }
    int orig_con_index = con_core2orig_[num_stage1_cons + stage2_con];
    if (auto expr = problem_.algebraic_con(orig_con_index).nonlinear_expr()) {
      RandomAffineExprExtractor extractor(*this, scenario_index);
      s.rhs_offsets_[stage2_con] += extractor.Visit(expr);
      for (auto term: extractor.linear_expr()) {
        // TODO: get coefficient from the linear part
        fmt::print("Term: {} {}\n", term.coef(), term.var_index());
      }
    }
  }

  // TODO: remove if unused
  /*std::vector<int> nonzero_coef_indices;
  nonzero_coef_indices.reserve(p.num_algebraic_cons());
  for (int core_var_index = 0, n = var_core2orig_.size();
       core_var_index < n; ++core_var_index) {
    int var_index = var_core2orig_[core_var_index];
    // Write the core coefficients and store them in the core_coefs_ vector.
    for (int k = p.col_start(var_index), end = p.col_start(var_index + 1);
         k != end; ++k) {
      int con_index = p.row_index(k);
      int core_con_index = con_orig2core_[con_index];
      // TODO
      //handler.OnTerm(core_con_index, core_coefs_[k]);
    }
  }

  coefs.assign(problem_.values(),
               problem_.values() + problem_.col_start(problem_.num_vars()));
  // Handle random variables/parameters in nonlinear constraint expressions.
  for (int core_con_index = num_stage_cons_[0], num_cons = this->num_cons();
       core_con_index < num_cons; ++core_con_index) {
    int con_index = con_core2orig_[core_con_index];
    auto con = problem_.algebraic_con(con_index);
    auto expr = con.nonlinear_expr();
    if (!expr) continue;
    // Get constraint coefficients for scenario.
    RandomAffineExprExtractor extractor(*this, scenario);
    auto value = extractor.Visit(expr);
    for (auto term: extractor.linear_expr()) {
      int var_index = term.var_index();
      int index = FindTerm(problem_, con_index, var_index);
      if (index == -1)
        throw Error("cannot find term ({}, {})", con_index, var_index);
      // Add extracted term coefficient to the one in constraint matrix.
      coefs[index] += term.coef();
    }
  }*/
}
}  // namespace mp

// TODO: migrate
/*
// Finds term with the specified constraint and variable indices in the
// column-wise constraint matrix.
int FindTerm(const ColProblem &p, int con_index, int var_index) {
  for (int i = p.col_start(var_index),
       n = p.col_start(var_index + 1); i < n; ++i) {
    if (p.row_index(i) == con_index)
      return i;
  }
  return -1;
}*/
