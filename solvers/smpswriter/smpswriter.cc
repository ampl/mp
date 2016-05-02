/*
 SMPS writer implemented as an AMPL solver.

 Copyright (C) 2013 - 2016 AMPL Optimization Inc

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

#include "smpswriter/smpswriter.h"
#include "mp/expr-visitor.h"

#include <cmath>
#include <cstdio>

namespace {

double GetConRHSAndType(const mp::Problem &p, int con_index, char &type) {
  mp::Problem::AlgebraicCon con = p.algebraic_con(con_index);
  double lb = con.lb(), ub = con.ub();
  double inf = std::numeric_limits<double>::infinity();
  if (lb <= -inf) {
    type = ub >= inf ? 'N' : 'L';
    return ub;
  }
  if (ub >= inf)
    type = 'G';
  else if (lb == ub)
    type = 'E';
  else
    throw mp::Error("SMPS writer doesn't support ranges");
  return lb;
}
}  // namespace

namespace mp {

class FileWriter {
 private:
  FILE *f_;
  FMT_DISALLOW_COPY_AND_ASSIGN(FileWriter);

 public:
  explicit FileWriter(fmt::CStringRef name)
    : f_(std::fopen(name.c_str(), "w")) {}
  ~FileWriter() { std::fclose(f_); }

  void Write(fmt::CStringRef format, const fmt::ArgList &args) {
    fmt::print(f_, format, args);
  }
  FMT_VARIADIC(void, Write, fmt::CStringRef)
};

SMPSWriter::SMPSWriter()
  : SolverImpl<ColProblem>("smpswriter", "SMPSWriter", 20160419),
    num_stage1_cons_(0) {
  AddSuffix("stage", 0, suf::VAR);
}

void SMPSWriter::AddRVElement(Expr arg, int rv_index, int element_index) {
  auto var = Cast<Reference>(arg);
  int var_index = var.index();
  rv_info_.push_back(RVInfo(var_index, rv_index, element_index));
  int &core_var_index = var_orig2core_[var_index];
  if (core_var_index != 0)
    throw Error("RV {}: redefinition of variable {}", rv_index, var_index);
  // Mark variable as random.
  core_var_index = -1;
}

void SMPSWriter::GetRandomVectors(const Problem &p) {
  var_orig2core_.resize(p.num_vars());
  int num_logical_cons = p.num_logical_cons();
  if (num_logical_cons == 0) {
    rvs_.resize(1);
    rvs_.back().set_num_realizations(1);
    return;
  }
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
    int num_realizations = -1;
    int realization_index = 0;
    int element_index = 0;
    auto arg = call.arg(0);
    if (arg.kind() != expr::VARIABLE)
      throw Error("RV {}: expected variable reference", con_index);
    AddRVElement(arg, con_index, element_index);
    for (int i = 1; i < num_args; ++i) {
      arg = call.arg(i);
      if (arg.kind() == expr::VARIABLE) {
        if (num_realizations == -1) {
          num_realizations = realization_index;
          rv.set_num_realizations(num_realizations);
          // TODO: get probabilties from the random function
          /*// Get probabilities.
          for (int i = 0; i < num_scenarios; ++i) {
            auto prob = Cast<mp::NumericConstant>(expr.arg(i));
            if (!prob)
              throw Error("probability is not constant");
            probabilities[i] = prob.value();
          }*/
        } else if (realization_index != num_realizations) {
          throw Error("RV {}: inconsistent number of realizations", con_index);
        }
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
  }
}

void SMPSWriter::WriteColumns(
    FileWriter &writer, const ColProblem &p, int num_core_cons,
    const std::vector<double> &core_obj_coefs,
    const std::vector<double> &coefs) {
  writer.Write("COLUMNS\n");
  std::vector<int> nonzero_coef_indices;
  nonzero_coef_indices.reserve(num_core_cons);
  int int_var_index = 0;
  bool integer_block = false;
  for (int core_var_index = 0, n = var_core2orig_.size();
       core_var_index < n; ++core_var_index) {
    int var_index = var_core2orig_[core_var_index];
    if (p.var(var_index).type() == var::CONTINUOUS) {
      if (integer_block) {
        writer.Write(
            "    INT{:<5}    'MARKER'      'INTEND'\n", int_var_index);
        integer_block = false;
      }
    } else if (!integer_block) {
      writer.Write(
          "    INT{:<5}    'MARKER'      'INTORG'\n", ++int_var_index);
      integer_block = true;
    }

    if (double obj_coef = core_obj_coefs[var_index])
      writer.Write("    C{:<7}  OBJ       {}\n", core_var_index + 1, obj_coef);

    // Write the core coefficients and store them in the core_coefs vector.
    for (int k = p.col_start(var_index), end = p.col_start(var_index + 1);
         k != end; ++k) {
      int con_index = p.row_index(k);
      int core_con_index = con_orig2core_[con_index];
      writer.Write("    C{:<7}  R{:<7}  {}\n",
          core_var_index + 1, core_con_index + 1, coefs[k]);
    }
  }
  if (integer_block)
    writer.Write("    INT{:<5}    'MARKER'      'INTEND'\n", int_var_index);
}

// Extracts an affine expression from a nonlinear one.
class AffineExprExtractor : public mp::ExprVisitor<AffineExprExtractor, void> {
 private:
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
  AffineExprExtractor() : constant_(0), coef_(1) {}

  const LinearExpr &linear_expr() const { return linear_; }

  void VisitVariable(Variable v) {
    linear_.AddTerm(v.index(), coef_);
  }

  void VisitBinary(BinaryExpr e) {
    switch (e.kind()) {
    case expr::ADD:
      Visit(e.lhs());
      Visit(e.rhs());
      break;
    case expr::MUL:
      if (auto n = Cast<mp::NumericConstant>(e.lhs()))
        return VisitMultiplied(n.value(), e.rhs());
      if (auto n = Cast<mp::NumericConstant>(e.rhs()))
        return VisitMultiplied(n.value(), e.lhs());
      throw UnsupportedError("nonlinear expression");
      break;
    default:
      mp::ExprVisitor<AffineExprExtractor, void>::VisitBinary(e);
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
  mp::Function random_;
  int scenario_;

 public:
  explicit RandomConstExprExtractor(mp::Function random, int scenario)
    : random_(random), scenario_(scenario) {}

  double VisitNumericConstant(NumericConstant n) {
    return n.value();
  }

  double VisitCall(CallExpr e) {
    auto function = e.function();
    if (function != random_)
      throw Error("unsupported function: {}", function.name());
    // TODO: check the number of arguments
    return this->Visit(e.arg(scenario_));
  }

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
  explicit RandomAffineExprExtractor(mp::Function random, int scenario)
    : Base(random, scenario), coef_(1) {}

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
int FindTerm(ColProblem &p, int con_index, int var_index) {
  for (int i = p.col_start(var_index),
       n = p.col_start(var_index + 1); i < n; ++i) {
    if (p.row_index(i) == con_index)
      return i;
  }
  return -1;
}

void SMPSWriter::GetScenario(ColProblem &p, int scenario,
                             std::vector<double> &coefs,
                             std::vector<double> &rhs) {
  coefs.assign(p.values(), p.values() + p.col_start(p.num_vars()));
  // Handle random variables/parameters in the constraint matrix.
  for (const auto &info: rv_info_) {
    int var_index = info.var_index;
    for (int k = p.col_start(var_index),
         end = p.col_start(var_index + 1); k != end; ++k) {
      auto value = rvs_[info.rv_index].value(info.element_index, scenario);
      fmt::print("RV: {} {} {}\n", info.element_index, scenario, value);
      rhs[p.row_index(k)] -= coefs[k] * value;
    }
  }
  // Handle random variables/parameters in nonlinear constraint expressions.
  for (int core_con_index = num_stage1_cons_, num_cons = p.num_algebraic_cons();
       core_con_index < num_cons; ++core_con_index) {
    int con_index = con_core2orig_[core_con_index];
    auto con = p.algebraic_con(con_index);
    auto expr = con.nonlinear_expr();
    if (!expr) continue;
    // Get constraint coefficients for scenario.
    RandomAffineExprExtractor extractor(random_, scenario);
    rhs[core_con_index] -= extractor.Visit(expr);
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

void SMPSWriter::Solve(ColProblem &p, SolutionHandler &) {
  // Find the random function.
  for (int i = 0, n = p.num_functions(); i < n; ++i) {
    auto function = p.function(i);
    if (std::strcmp(function.name(), "random") == 0) {
      random_ = function;
      break;
    }
  }

  GetRandomVectors(p);

  // Count the number of stages, the number of variables in the first stage and
  // compute core indices for the first-stage variables.
  int num_vars = p.num_vars();
  int num_stage1_vars = num_vars;
  int num_stages = 1;
  var_core2orig_.resize(p.num_vars() - rv_info_.size());
  if (IntSuffix stage_suffix = p.suffixes(suf::VAR).Find<int>("stage")) {
    num_stage1_vars = 0;
    for (int i = 0; i < num_vars; ++i) {
      if (var_orig2core_[i])
        continue;  // Skip random variables.
      int stage_plus_1 = stage_suffix.value(i);
      if (stage_plus_1 > 1) {
        num_stages = std::max(stage_plus_1, num_stages);
      } else {
        var_core2orig_[num_stage1_vars] = i;
        var_orig2core_[i] = num_stage1_vars;
        ++num_stage1_vars;
      }
    }
  }
  if (num_stages > 2)
    throw Error("SMPS writer doesn't support problems with more than 2 stages");

  // Compute core indices for variables in later stages.
  int stage2_index = num_stage1_vars;
  if (num_stages > 1) {
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

  std::string smps_basename = p.name();
  std::string::size_type ext_pos = smps_basename.rfind('.');
  if (ext_pos != std::string::npos)
    smps_basename.resize(ext_pos);

  int num_cons = p.num_algebraic_cons();
  num_stage1_cons_ = num_cons;
  Expr obj_expr;
  con_core2orig_.resize(num_cons);
  con_orig2core_.resize(num_cons);
  if (num_stages > 1) {
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
            --num_stage1_cons_;
          con_stage = 1;
        }
      }
    }

    // Compute core indices for constraints.
    int stage1_index = 0, stage2_index = num_stage1_cons_;
    for (int i = 0; i < num_cons; ++i) {
      int &index = con_stages[i] != 0 ? stage2_index : stage1_index;
      con_core2orig_[index] = i;
      con_orig2core_[i] = index++;
    }
    assert(stage1_index == num_stage1_cons_ && stage2_index == num_cons);

    // Get the number of scenarios from the expectation in the objective.
    CallExpr expr;
    if (p.num_objs() > 0) {
      if (auto e = p.obj(0).nonlinear_expr())
        expr = Cast<CallExpr>(e);
    }
    if (!expr || std::strcmp(expr.function().name(), "expectation") != 0)
      throw Error("objective doesn't contain expectation");
    // TODO: check that the number of arguments is 1.
    obj_expr = expr.arg(0);
    // TODO: check that second-stage variables only occur in expectation
  } else {
    for (int i = 0; i < num_vars; ++i) {
      var_core2orig_[i] = i;
      var_orig2core_[i] = i;
    }
    for (int i = 0; i < num_cons; ++i) {
      con_core2orig_[i] = i;
      con_orig2core_[i] = i;
    }
  }

  // Write the .tim file.
  {
    FileWriter writer(smps_basename + ".tim");
    writer.Write(
      "TIME          PROBLEM\n"
      "PERIODS\n"
      "    C1        OBJ                      T1\n");
    if (num_stages > 1) {
      writer.Write("    C{:<7}  R{:<7}                 T2\n",
          num_stage1_vars + 1, num_stage1_cons_ + 1);
    }
    writer.Write("ENDATA\n");
  }

  // Write the .cor file.
  std::vector<double> core_coefs;
  std::vector<double> base_rhs(p.num_algebraic_cons()), core_rhs;
  {
    FileWriter writer(smps_basename + ".cor");
    writer.Write(
      "NAME          PROBLEM\n"
      "ROWS\n"
      " N  OBJ\n");
    for (int i = 0; i < num_cons; ++i) {
      char type = 0;
      base_rhs[i] = GetConRHSAndType(p, con_core2orig_[i], type);
      writer.Write(" {}  R{}\n", type, i + 1);
    }

    std::vector<double> core_obj_coefs(num_vars);
    if (p.num_objs() != 0) {
      // Get objective coefficients in the core problem (first scenario).
      for (auto term: p.obj(0).linear_expr())
        core_obj_coefs[term.var_index()] = term.coef();
      if (obj_expr) {
        AffineExprExtractor extractor;
        extractor.Visit(obj_expr);
        for (auto term: extractor.linear_expr())
          core_obj_coefs[term.var_index()] += term.coef();
      }
    }

    core_rhs = base_rhs;
    GetScenario(p, 0, core_coefs, core_rhs);
    WriteColumns(writer, p, num_cons, core_obj_coefs, core_coefs);

    writer.Write("RHS\n");
    for (int i = 0; i < num_cons; ++i)
      writer.Write("    RHS1      R{:<7}  {}\n", i + 1, core_rhs[i]);

    writer.Write("BOUNDS\n");
    double inf = std::numeric_limits<double>::infinity();
    for (int core_var_index = 0, n = var_core2orig_.size();
         core_var_index < n; ++core_var_index) {
      auto var = p.var(var_core2orig_[core_var_index]);
      double lb = var.lb(), ub = var.ub();
      if (lb != 0)
        writer.Write(" LO BOUND1      C{:<7}  {}\n", core_var_index + 1, lb);
      if (ub < inf)
        writer.Write(" UP BOUND1      C{:<7}  {}\n", core_var_index + 1, ub);
    }

    writer.Write("ENDATA\n");
  }

  // Write the .sto file.
  {
    const auto &rv = rvs_[0];
    FileWriter writer(smps_basename + ".sto");
    writer.Write(
      "STOCH         PROBLEM\n"
      "SCENARIOS     DISCRETE\n");
    if (num_stages > 1) {
      std::vector<double> coefs, rhs;
      writer.Write(" SC SCEN1     'ROOT'    {:<12}   T1\n", rv.probability(0));
      for (size_t s = 1, num_scen = rv.num_realizations(); s < num_scen; ++s) {
        writer.Write(" SC SCEN{:<4}  SCEN1     {:<12}   T2\n",
            s + 1, rv.probability(s));
        rhs = base_rhs;
        GetScenario(p, s, coefs, rhs);
        // Compare to the core and write differences.
        for (int j = 0, num_vars = p.num_vars(); j < num_vars; ++j) {
          for (int k = p.col_start(j), n = p.col_start(j + 1); k != n; ++k) {
            double coef = coefs[k];
            if (coef == core_coefs[k]) continue;
            int core_con_index = con_orig2core_[p.row_index(k)];
            writer.Write("    C{:<7}  R{:<7}  {}\n",
                var_orig2core_[j] + 1, core_con_index + 1, coef);
          }
        }
        for (int i = num_stage1_cons_, n = p.num_algebraic_cons(); i < n; ++i) {
          double value = rhs[i];
          if (core_rhs[i] != value)
            writer.Write("    RHS1      R{:<7}  {}\n", i + 1, value);
        }
      }
    }
    writer.Write("ENDATA\n");
  }
}

SolverPtr create_smpswriter(const char *) {
  return SolverPtr(new SMPSWriter());
}
}  // namespace mp
