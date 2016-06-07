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

// Extracts an affine expression from a nonlinear one.
class AffineExprExtractor : public ExprVisitor<AffineExprExtractor, void> {
 private:
  const Problem &problem_;
  LinearExpr linear_;
  double constant_;
  double coef_;

  void VisitMultiplied(double multiplier, NumericExpr expr) {
    double saved_coef = coef_;
    coef_ *= multiplier;
    Visit(expr);
    coef_ = saved_coef;
  }

 public:
  AffineExprExtractor(const Problem &p) : problem_(p), constant_(0), coef_(1) {}

  const LinearExpr &linear_expr() const { return linear_; }

  void VisitVariable(Variable v) {
    linear_.AddTerm(v.index(), coef_);
  }

  void VisitCommonExpr(CommonExpr e) {
    auto common_expr = problem_.common_expr(e.index());
    for (auto term: common_expr.linear_expr())
      linear_.AddTerm(term.var_index(), coef_ * term.coef());
    Visit(common_expr.nonlinear_expr());
  }

  void VisitBinary(BinaryExpr e);

  void VisitSum(SumExpr e) {
    for (auto arg: e)
      Visit(arg);
  }
};

void AffineExprExtractor::VisitBinary(BinaryExpr e) {
  switch (e.kind()) {
  case expr::ADD:
    Visit(e.lhs());
    Visit(e.rhs());
    break;
  case expr::SUB:
    Visit(e.lhs());
    coef_ = -coef_;
    Visit(e.rhs());
    coef_ = -coef_;
    break;
  case expr::MUL:
    if (auto n = Cast<NumericConstant>(e.lhs()))
      return VisitMultiplied(n.value(), e.rhs());
    if (auto n = Cast<NumericConstant>(e.rhs()))
      return VisitMultiplied(n.value(), e.lhs());
    throw UnsupportedError("nonlinear expression");
    break;
  default:
    ExprVisitor<AffineExprExtractor, void>::VisitBinary(e);
  }
}

template <typename Impl>
class RandomConstExprExtractor : public mp::ExprVisitor<Impl, double> {
 private:
  const SPAdapter &adapter_;
  int scenario_;

 public:
  explicit RandomConstExprExtractor(const SPAdapter &extractor, int scenario)
    : adapter_(extractor), scenario_(scenario) {}

  double VisitNumericConstant(NumericConstant n) {
    return n.value();
  }

  double VisitVariable(Reference v) {
    int rv_index = adapter_.GetRandVarIndex(v.index());
    return rv_index >= 0 ?
          adapter_.GetRealization(rv_index, scenario_) :
          mp::ExprVisitor<Impl, double>::VisitVariable(v);
  }

  /*double VisitCall(CallExpr e) {
    auto function = e.function();
    if (function != random_)
      throw Error("unsupported function: {}", function.name());
    // TODO: check the number of arguments
    return this->Visit(e.arg(scenario_));
  }*/

  // TODO
};

// Extracts an affine expression for a single scenario from an expression
// containing random variables.
class RandomAffineExprExtractor :
    public RandomConstExprExtractor<RandomAffineExprExtractor> {
 private:
  LinearExpr linear_;
  double coef_;

  typedef RandomConstExprExtractor<RandomAffineExprExtractor> Base;

  double ExtractTerm(Expr coef, Expr var) {
    RandomConstExprExtractor extractor(*this);
    linear_.AddTerm(Cast<Reference>(var).index(),
                    coef_ * extractor.Visit(coef));
    return 0;
  }

 public:
  explicit RandomAffineExprExtractor(const SPAdapter &extractor, int scenario)
    : Base(extractor, scenario), coef_(1) {}

  const LinearExpr &linear_expr() const { return linear_; }

  double VisitUnary(UnaryExpr e) {
    if (e.kind() != expr::MINUS)
      return Base::VisitUnary(e);
    double saved_coef = coef_;
    coef_ = -coef_;
    double result = -Visit(e.arg());
    coef_ = saved_coef;
    return result;
  }

  double VisitBinary(BinaryExpr e) {
    switch (e.kind()) {
    case expr::MUL:
      if (e.rhs().kind() == expr::VARIABLE)
        return ExtractTerm(e.lhs(), e.rhs());
      if (e.lhs().kind() == expr::VARIABLE)
        return ExtractTerm(e.rhs(), e.lhs());
      throw UnsupportedError("nonlinear expression");
      break;
    default:
      return Base::VisitBinary(e);
    }
  }

  // TODO
};

class NullSuffix {
 public:
  int value(int) const { return 0; }
};
}  // namespace

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
    int element_index = 0;
    // Get probabilities.
    for (; arg_index < num_args; ++arg_index) {
      auto arg = call.arg(arg_index);
      if (arg.kind() == expr::VARIABLE) {
        AddRVElement(arg, con_index, element_index);
        ++arg_index;
        break;
      }
      auto prob = Cast<NumericConstant>(arg);
      if (!prob)
        throw Error("{}: expected variable or constant", lcon_name(con_index));
      double p = prob.value();
      if (p < 0 || p > 1)
        throw Error("{}: invalid probability {}", lcon_name(con_index), p);
      rv.AddProbability(p);
    }
    // Get realizations.
    int num_realizations = rv.num_realizations();
    if (num_realizations == 0)
      num_realizations = -1;
    int realization_index = 0;
    for (; arg_index < num_args; ++arg_index) {
      auto arg = call.arg(arg_index);
      if (arg.kind() == expr::VARIABLE) {
        if (num_realizations == -1)
          num_realizations = realization_index;
        else if (realization_index != num_realizations)
          throw Error("RV {}: inconsistent number of realizations", con_index);
        ++element_index;
        realization_index = 0;
        AddRVElement(arg, con_index, element_index);
      } else if (auto constant = Cast<NumericConstant>(arg)) {
        rv.Add(constant.value());
        ++realization_index;
      } else {
        throw Error("RV {}: expected variable or constant", con_index);
      }
    }
    rv.set_num_realizations(num_realizations == -1 ?
                              realization_index : num_realizations);
  }
}

template <typename Suffix>
int SPAdapter::ProcessStage1Vars(Suffix stage) {
  int num_stage1_vars = 0;
  for (int i = 0, n = problem_.num_vars(); i < n; ++i) {
    if (var_orig2core_[i])
      continue;  // Skip random variables.
    int stage_plus_1 = stage.value(i);
    if (stage_plus_1 > 1) {
      num_stages_ = std::max(stage_plus_1, num_stages_);
    } else {
      var_core2orig_[num_stage1_vars] = i;
      var_orig2core_[i] = num_stage1_vars;
      ++num_stage1_vars;
    }
  }
  return num_stage1_vars;
}

void SPAdapter::ProcessObjs(int num_stage1_vars) {
  int num_objs = problem_.num_objs();
  if (num_objs == 0)
    return;
  nonlinear_objs_.reserve(num_objs);
  for (auto obj: problem_.objs())
    nonlinear_objs_.push_back(obj.nonlinear_expr());
  // Strip expectation from the first objective.
  auto obj = problem_.obj(0);
  for (auto term: obj.linear_expr()) {
    if (var_orig2core_[term.var_index()] >= num_stage1_vars && term.coef() != 0)
      throw Error("second-stage variable outside of expectation in objective");
  }
  auto obj_expr = obj.nonlinear_expr();
  if (!obj_expr)
    return;
  auto call = Cast<CallExpr>(obj_expr);
  if (!call || std::strcmp(call.function().name(), "expectation") != 0)
    return;
  NumericExpr arg;
  if (call.num_args() != 1 || !(arg = Cast<NumericExpr>(call.arg(0))))
    throw Error("invalid arguments to expectation");
  nonlinear_objs_[0] = arg;
}

int SPAdapter::ProcessCons(int num_stage1_vars) {
  int num_vars = problem_.num_vars();
  int num_cons = problem_.num_algebraic_cons();
  std::vector<int> con_stages(num_cons);
  if (num_stage1_vars != num_vars) {
    // Compute stage of each constraint as a maximum of stages of
    // variables in it.
    for (int j = 0; j < num_vars; ++j) {
      int core_var_index = var_orig2core_[j];
      if (core_var_index >= 0 && core_var_index < num_stage1_vars)
        continue;
      // Update stages of all constraints containing this variable.
      for (int k = problem_.col_start(j),
           end = problem_.col_start(j + 1); k != end; ++k) {
        con_stages[problem_.row_index(k)] = 1;
      }
    }
  }
  int num_stage1_cons = 0;
  for (int i = 0; i < num_cons; ++i) {
    if (problem_.algebraic_con(i).nonlinear_expr())
      con_stages[i] = 1;
    else if (con_stages[i] == 0)
      ++num_stage1_cons;
  }
  // Compute core indices for constraints.
  con_core2orig_.resize(num_cons);
  con_orig2core_.resize(num_cons);
  int stage1_index = 0, stage2_index = num_stage1_cons;
  for (int i = 0; i < num_cons; ++i) {
    int &index = con_stages[i] != 0 ? stage2_index : stage1_index;
    con_core2orig_[index] = i;
    con_orig2core_[i] = index++;
  }
  assert(stage1_index == num_stage1_cons && stage2_index == num_cons);
  if (num_stage1_cons != num_cons)
    num_stages_ = 2;
  return num_stage1_cons;
}

void SPAdapter::AddRVElement(Expr arg, int rv_index, int element_index) {
  auto var = Cast<Reference>(arg);
  int var_index = var.index();
  rv_info_.push_back(RVInfo(var_index, rv_index, element_index));
  int &core_var_index = var_orig2core_[var_index];
  if (core_var_index != 0)
    throw Error("RV {}: redefinition of variable {}", rv_index, var_index);
  // Mark variable as random.
  core_var_index = -static_cast<int>(rv_info_.size());
}

void SPAdapter::ExtractRandomTerms() {
  int num_stage1_cons = num_stage_cons_[0];
  int num_stage2_cons = num_stage_cons_[1];

  // A matrix containing linear terms involving random variables.
  // The major dimension is equal to the number of second-stage constraints.
  linear_random_.resize_major(num_stage2_cons);

  // Count random variables in second-stage constraints.
  for (auto info: rv_info_) {
    for (int k = problem_.col_start(info.var_index),
         end = problem_.col_start(info.var_index + 1); k != end; ++k) {
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
  for (auto info: rv_info_) {
    for (int k = problem_.col_start(info.var_index),
         end = problem_.col_start(info.var_index + 1); k != end; ++k) {
       int core_con_index = con_orig2core_[problem_.row_index(k)];
       assert(core_con_index >= num_stage1_cons);
       int index = core_con_index - num_stage1_cons + 1;
       int element_index = linear_random_.start(index)++;
       linear_random_.index(element_index) = info.var_index;
       linear_random_.coef(element_index) = -problem_.value(k);
    }
  }
}

SPAdapter::SPAdapter(const ColProblem &p)
  : problem_(p), num_stages_(1) {
  // Find the random function.
  for (int i = 0, n = p.num_functions(); i < n; ++i) {
    auto function = p.function(i);
    if (std::strcmp(function.name(), "random") == 0) {
      random_ = function;
      break;
    }
  }

  GetRandomVectors(p);

  var_core2orig_.resize(p.num_vars() - rv_info_.size());
  IntSuffix stage_suffix = p.suffixes(suf::VAR).Find<int>("stage");
  int num_stage1_vars = stage_suffix ?
        ProcessStage1Vars(stage_suffix) : ProcessStage1Vars(NullSuffix());
  if (num_stages_ > 2)
    throw Error("SP problems with more than 2 stages are not supported");

  ProcessObjs(num_stage1_vars);

  // Compute core indices for variables in later stages.
  int num_vars = p.num_vars();
  int stage2_index = num_stage1_vars;
  if (num_stages_ > 1) {
    // Temporary change mapping for the core variable 0, not to confuse it with
    // second-stage variables.
    if (num_stage1_vars > 0)
      var_orig2core_[var_core2orig_[0]] = 1;
    for (int i = 0; i < num_vars; ++i) {
      if (var_orig2core_[i] == 0) {
        var_core2orig_[stage2_index] = i;
        var_orig2core_[i] = stage2_index++;
      }
    }
    // Restore mapping for the core variable 0.
    if (num_stage1_vars > 0)
      var_orig2core_[var_core2orig_[0]] = 0;
    assert(static_cast<std::size_t>(num_vars) ==
           stage2_index + rv_info_.size());
  }

  int num_stage1_cons = ProcessCons(num_stage1_vars);

  num_stage_vars_.resize(num_stages_);
  num_stage_cons_.resize(num_stages_);
  num_stage_vars_[0] = num_stage1_vars;
  num_stage_cons_[0] = num_stage1_cons;
  if (num_stages_ > 1) {
    num_stage_vars_[1] = this->num_vars() - num_stage1_vars;
    num_stage_cons_[1] = this->num_cons() - num_stage1_cons;
    ExtractRandomTerms();
  }
}

std::vector<double> SPAdapter::core_obj() const {
  std::vector<double> obj(num_vars());
  if (problem_.num_objs() == 0)
    return obj;
  // Get objective coefficients in the core problem (first scenario).
  for (auto term: problem_.obj(0).linear_expr())
    obj[var_orig2core_[term.var_index()]] = term.coef();
  if (NumericExpr expr = nonlinear_objs_[0]) {
    AffineExprExtractor extractor(problem_);
    extractor.Visit(expr);
    for (auto term: extractor.linear_expr())
      obj[var_orig2core_[term.var_index()]] += term.coef();
  }
  return obj;
}

void SPAdapter::GetScenario(Scenario &s, int scenario_index) const {
  int num_stage1_cons = num_stage_cons_[0];
  int num_stage2_cons = num_stage_cons_[1];
  s.rhs_offsets_.assign(num_stage2_cons, 0);
  for (int stage2_con = 0; stage2_con < num_stage2_cons; ++stage2_con) {
    for (int i = linear_random_.start(stage2_con),
         end = linear_random_.start(stage2_con + 1); i < end; ++i) {
      int rv_index = GetRandVarIndex(linear_random_.index(i));
      assert(rv_index >= 0);
      s.rhs_offsets_[stage2_con] +=
          linear_random_.coef(i) * GetRealization(rv_index, scenario_index);
    }
    int orig_con_index = con_core2orig_[num_stage1_cons + stage2_con];
    if (auto expr = problem_.algebraic_con(orig_con_index).nonlinear_expr()) {
      RandomAffineExprExtractor extractor(*this, scenario_index);
      s.rhs_offsets_[stage2_con] += extractor.Visit(expr);
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
}

void SPAdapter::GetScenario(int scenario, std::vector<double> &coefs,
                            std::vector<Bounds> &rhs) {
  coefs.assign(problem_.values(),
               problem_.values() + problem_.col_start(problem_.num_vars()));
  // Handle random variables/parameters in the constraint matrix.
  for (const auto &info: rv_info_) {
    int var_index = info.var_index;
    for (int k = problem_.col_start(var_index),
         end = problem_.col_start(var_index + 1); k != end; ++k) {
      auto value = rvs_[info.rv_index].value(info.element_index, scenario);
      value *= coefs[k];
      auto &bounds = rhs[problem_.row_index(k)];
      bounds.lb -= value;
      bounds.ub -= value;
    }
  }
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
    auto &bounds = rhs[core_con_index];
    bounds.lb -= value;
    bounds.ub -= value;
    for (auto term: extractor.linear_expr()) {
      int var_index = term.var_index();
      int index = FindTerm(problem_, con_index, var_index);
      if (index == -1)
        throw Error("cannot find term ({}, {})", con_index, var_index);
      // Add extracted term coefficient to the one in constraint matrix.
      coefs[index] += term.coef();
    }
  }
}*/
