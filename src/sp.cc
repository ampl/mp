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

#include <cstring>    // std::strcmp
#include <algorithm>  // std::max

namespace mp {
namespace {

// Returns the name of a logical constraint.
inline std::string lcon_name(int index) {
  return fmt::format("_slogcon[{}]", index + 1);
}

// Collects the list of variables that appear in a nonlinear expression.
class VariableCollector:
    public internal::BasicRandomAffineExprExtractor<VariableCollector> {
 private:
  SparseMatrix<double> &vars_in_nonlinear_;
  std::vector<bool> visited_vars_;

 public:
  VariableCollector(const SPAdapter &sp,
                    SparseMatrix<double> &vars_in_nonlinear):
    internal::BasicRandomAffineExprExtractor<VariableCollector>(sp, 0),
    vars_in_nonlinear_(vars_in_nonlinear), visited_vars_(sp.num_vars()) {}

  void Collect();

  double VisitNumericConstant(NumericConstant) { return 0; }

  void OnVar(int var_index) {
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
  for (int con_index = 0; con_index < num_stage2_cons; ++con_index) {
    auto expr = sp_.con(num_stage1_cons + con_index).nonlinear_expr();
    if (expr)
      Visit(expr);
    int num_elements = vars_in_nonlinear_.num_elements();
    vars_in_nonlinear_.start(con_index + 1) = num_elements;
    vars_in_nonlinear_.resize_elements(num_elements);
    for (int i = vars_in_nonlinear_.start(con_index),
         n = vars_in_nonlinear_.start(con_index + 1); i < n; ++i) {
      int core_var_index = vars_in_nonlinear_.index(i);
      double &coef = this->coef(core_var_index);
      vars_in_nonlinear_.value(i) = coef;
      coef = 0;
      visited_vars_[core_var_index] = false;
    }
  }
}

class NullSuffix {
 public:
  int value(int) const { return 0; }
};
}  // namespace

namespace internal {

void Transpose(SparseMatrix<double> &m, SparseMatrix<double*> &transposed,
               int major_size) {
  transposed.resize_major(major_size);
  transposed.resize_elements(m.num_elements());
  int orig_major_size = m.major_size();
  for (int i = 0; i < orig_major_size; ++i) {
    for (int j = m.start(i), n = m.start(i + 1); j < n; ++j)
      ++transposed.start(m.index(j) + 1);
  }
  int start = 0;
  for (int i = 1; i <= major_size; ++i) {
    double next = start + transposed.start(i);
    transposed.start(i) =  start;
    start = next;
  }
  for (int i = 0; i < orig_major_size; ++i) {
    for (int j = m.start(i), n = m.start(i + 1); j < n; ++j) {
      int index = transposed.start(m.index(j) + 1)++;
      transposed.index(index) = i;
      transposed.value(index) = &m.value(j);
    }
  }
}

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
    for (auto i = e.begin(), end = e.end(); i != end; ++i)
      Visit(*i);
  }
};

void AffineExprExtractor::VisitCommonExpr(Reference e) {
  auto common_expr = sp_.problem_.common_expr(e.index());
  const auto &linear = common_expr.linear_expr();
  for (auto i = linear.begin(), end = linear.end(); i != end; ++i)
    AddTerm(i->var_index(), i->coef());
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
  auto objs = problem_.objs();
  for (auto i = objs.begin(), end = objs.end(); i != end; ++i)
    nonlinear_objs_.push_back(i->nonlinear_expr());
  // Strip expectation from the first objective.
  auto obj = problem_.obj(0);
  int num_vars = this->num_vars();
  std::vector<double> core_obj(num_vars);
  const auto &linear = obj.linear_expr();
  for (auto i = linear.begin(), end = linear.end(); i != end; ++i) {
    int core_var_index = var_orig2core_[i->var_index()];
    if (core_var_index >= num_stage1_vars && i->coef() != 0)
      throw Error("second-stage variable outside of expectation in objective");
    core_obj[core_var_index] = i->coef();
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
       linear_random_.value(element_index) = -problem_.value(k);
    }
  }

  // Get variables that appear in nonlinear parts of constraint expressions.
  VariableCollector collector(*this, vars_in_nonlinear_);
  const double *coefs = problem_.values();
  collector.Collect();

  // Convert vars_in_nonlinear_ to column-major form.
  int num_vars = this->num_vars();
  SparseMatrix<double*> var2con;
  internal::Transpose(vars_in_nonlinear_, var2con, num_vars);

  // Combine the first scenario with core coefficients.
  core_coefs_.assign(coefs, coefs + problem_.col_start(problem_.num_vars()));
  for (int core_var_index = 0; core_var_index < num_vars; ++core_var_index) {
    int var_index = var_core2orig_[core_var_index];
    int elt_index = problem_.col_start(var_index);
    for (int j = var2con.start(core_var_index),
         n = var2con.start(core_var_index + 1); j < n; ++j) {
      int con_index = con_core2orig_[var2con.index(j) + num_stage1_cons];
      while (problem_.row_index(elt_index) != con_index)
        ++elt_index;
      assert(elt_index < problem_.col_start(var_index + 1));
      double coef = core_coefs_[elt_index] + *var2con.value(j);
      *var2con.value(j) = core_coefs_[elt_index];
      core_coefs_[elt_index] = coef;
    }
  }
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
}  // namespace mp
