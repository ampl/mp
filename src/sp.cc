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
#include <vector>

namespace mp {
namespace {

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

  void VisitBinary(BinaryExpr e) {
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

  void VisitSum(SumExpr e) {
    for (auto arg: e)
      Visit(arg);
  }
};

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
    int rv_index = adapter_.GetRVIndex(v.index());
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

class NullSuffix {
 public:
  int value(int) const { return 0; }
};
}  // namespace

template <typename Suffix>
int SPAdapter::ProcessStage1Vars(const ColProblem &p, Suffix stage) {
  int num_stage1_vars = 0;
  for (int i = 0, n = p.num_vars(); i < n; ++i) {
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
        throw Error("RV {}: expected variable or constant", con_index);
      rv.AddProbability(prob.value());
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

void SPAdapter::GetScenario(const ColProblem &p, int scenario,
                              std::vector<double> &coefs,
                              std::vector<Bounds> &rhs) {
  coefs.assign(p.values(), p.values() + p.col_start(p.num_vars()));
  // Handle random variables/parameters in the constraint matrix.
  for (const auto &info: rv_info_) {
    int var_index = info.var_index;
    for (int k = p.col_start(var_index),
         end = p.col_start(var_index + 1); k != end; ++k) {
      auto value = rvs_[info.rv_index].value(info.element_index, scenario);
      value *= coefs[k];
      auto &bounds = rhs[p.row_index(k)];
      bounds.lb -= value;
      bounds.ub -= value;
    }
  }
  // Handle random variables/parameters in nonlinear constraint expressions.
  for (int core_con_index = num_stage_cons_[0], num_cons = this->num_cons();
       core_con_index < num_cons; ++core_con_index) {
    int con_index = con_core2orig_[core_con_index];
    auto con = p.algebraic_con(con_index);
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
      int index = FindTerm(p, con_index, var_index);
      if (index == -1)
        throw Error("cannot find term ({}, {})", con_index, var_index);
      // Add extracted term coefficient to the one in constraint matrix.
      coefs[index] += term.coef();
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
        ProcessStage1Vars(p, stage_suffix) : ProcessStage1Vars(p, NullSuffix());
  if (num_stages_ > 2)
    throw Error("SP problems with more than 2 stages are not supported");

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

  int num_cons = p.num_algebraic_cons();
  int num_stage1_cons = num_cons;
  con_core2orig_.resize(num_cons);
  con_orig2core_.resize(num_cons);
  if (num_stages_ > 1) {
    // Compute stage of each constraint as a maximum of stages of
    // variables in it.
    std::vector<int> con_stages(num_cons);
    for (int j = 0; j < num_vars; ++j) {
      if (var_orig2core_[j] < num_stage1_vars)
        continue;
      // Update stages of all constraints containing this variable.
      for (int k = p.col_start(j), end = p.col_start(j + 1); k != end; ++k) {
        int con_index = p.row_index(k);
        int &con_stage = con_stages[con_index];
        if (con_stage < 1) {
          if (con_stage == 0)
            --num_stage1_cons;
          con_stage = 1;
        }
      }
    }

    // Compute core indices for constraints.
    int stage1_index = 0, stage2_index = num_stage1_cons;
    for (int i = 0; i < num_cons; ++i) {
      int &index = con_stages[i] != 0 ? stage2_index : stage1_index;
      con_core2orig_[index] = i;
      con_orig2core_[i] = index++;
    }
    assert(stage1_index == num_stage1_cons && stage2_index == num_cons);

    // Get the number of scenarios from the expectation in the objective.
    CallExpr expr;
    if (p.num_objs() > 0) {
      auto obj = p.obj(0);
      for (auto term: obj.linear_expr()) {
        if (var_orig2core_[term.var_index()] >= num_stage1_vars &&
            term.coef() != 0) {
          throw mp::Error(
                "second-stage variable outside of expectation in objective");
        }
      }
      if (auto e = obj.nonlinear_expr()) {
        expr = Cast<CallExpr>(e);
        if (expr && std::strcmp(expr.function().name(), "expectation")) {
          // TODO: check that the number of arguments is 1.
          obj_expr_ = expr.arg(0);
        }
      }
    }
  } else {
    for (int i = 0; i < num_cons; ++i) {
      con_core2orig_[i] = i;
      con_orig2core_[i] = i;
    }
  }
  num_stage_vars_.resize(num_stages_);
  num_stage_cons_.resize(num_stages_);
  num_stage_vars_[0] = num_stage1_vars;
  num_stage_cons_[0] = num_stage1_cons;
  if (num_stages_ > 1) {
    num_stage_vars_[1] = this->num_vars() - num_stage1_vars;
    num_stage_cons_[1] = this->num_cons() - num_stage1_cons;
  }

  base_rhs_.resize(p.num_algebraic_cons());
  for (int i = 0; i < num_cons; ++i) {
    auto con = p.algebraic_con(con_core2orig_[i]);
    base_rhs_[i] = Bounds(con.lb(), con.ub());
  }

  core_rhs_ = base_rhs_;
  std::vector<double> core_coefs;
  GetScenario(p, 0, core_coefs, core_rhs_);

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

std::vector<double> SPAdapter::core_obj() const {
  std::vector<double> obj(num_vars());
  if (problem_.num_objs() == 0)
    return obj;
  // Get objective coefficients in the core problem (first scenario).
  for (auto term: problem_.obj(0).linear_expr())
    obj[var_orig2core_[term.var_index()]] = term.coef();
  if (obj_expr_) {
    AffineExprExtractor extractor(problem_);
    extractor.Visit(obj_expr_);
    for (auto term: extractor.linear_expr())
      obj[var_orig2core_[term.var_index()]] += term.coef();
  }
  return obj;
}
}  // namespace mp
