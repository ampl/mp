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

void SMPSWriter::WriteColumns(
    FileWriter &writer, const ColProblem &p, int num_core_cons,
    const std::vector<double> &core_obj_coefs,
    const std::vector<double> &coefs) {
  writer.Write("COLUMNS\n");
  std::vector<double> core_coefs(num_core_cons);
  std::vector<int> nonzero_coef_indices;
  nonzero_coef_indices.reserve(num_core_cons);
  int int_var_index = 0;
  bool integer_block = false;
  for (int core_var_index = 0, n = p.num_vars();
       core_var_index < n; ++core_var_index) {
    int var_index = var_indices_[core_var_index];
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
      int core_con_index = core_con_indices_[con_index];
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
class RandomConstantExprExtractor : public mp::ExprVisitor<Impl, double> {
 protected:
  int scenario_;

 public:
  explicit RandomConstantExprExtractor(int scenario) : scenario_(scenario) {}

  double VisitNumericConstant(NumericConstant n) {
    return n.value();
  }

  double VisitCall(CallExpr e) {
    // TODO: check the number of arguments
    if (std::strcmp(e.function().name(), "random") == 0)
      return this->Visit(e.arg(scenario_));
    return mp::ExprVisitor<Impl, double>::VisitCall(e);
  }

  // TODO
};

// Extracts an affine expression for a single scenario from an expression
// containing random variables.
class RandomAffineExprExtractor :
    public RandomConstantExprExtractor<RandomAffineExprExtractor> {
 private:
  LinearExpr linear_;
  double coef_;

  typedef RandomConstantExprExtractor<RandomAffineExprExtractor> Base;

  double ExtractTerm(Expr coef, Expr var) {
    RandomConstantExprExtractor extractor(scenario_);
    linear_.AddTerm(Cast<Reference>(var).index(),
                    coef_ * extractor.Visit(coef));
    return 0;
  }

 public:
  explicit RandomAffineExprExtractor(int scenario) : Base(scenario), coef_(1) {}

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
  for (int core_con_index = num_stage1_cons_, num_cons = p.num_algebraic_cons();
       core_con_index < num_cons; ++core_con_index) {
    int con_index = con_indices_[core_con_index];
    auto con = p.algebraic_con(con_index);
    auto expr = con.nonlinear_expr();
    if (!expr) continue;
    // Get constraint coefficients for scenario.
    RandomAffineExprExtractor extractor(scenario);
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
  // Count the number of stages, the number of variables in the first stage and
  // compute core indices for the first-stage variables.
  int num_vars = p.num_vars();
  int num_stage1_vars = num_vars;
  int num_stages = 1;
  var_indices_.resize(num_vars);
  core_var_indices_.resize(num_vars, -1);
  if (IntSuffix stage_suffix = p.suffixes(suf::VAR).Find<int>("stage")) {
    num_stage1_vars = 0;
    for (int i = 0; i < num_vars; ++i) {
      int stage_plus_1 = stage_suffix.value(i);
      if (stage_plus_1 > 1) {
        num_stages = std::max(stage_plus_1, num_stages);
      } else {
        var_indices_[num_stage1_vars] = i;
        core_var_indices_[i] = num_stage1_vars;
        ++num_stage1_vars;
      }
    }
  }
  if (num_stages > 2)
    throw Error("SMPS writer doesn't support problems with more than 2 stages");

  // Compute core indices for variables in later stages.
  int stage2_index = num_stage1_vars;
  if (num_stages > 1) {
    for (int i = 0; i < num_vars; ++i) {
      if (core_var_indices_[i] == -1) {
        var_indices_[stage2_index] = i;
        core_var_indices_[i] = stage2_index++;
      }
    }
    assert(stage2_index == num_vars);
  }

  std::string smps_basename = p.name();
  std::string::size_type ext_pos = smps_basename.rfind('.');
  if (ext_pos != std::string::npos)
    smps_basename.resize(ext_pos);

  int num_cons = p.num_algebraic_cons();
  num_stage1_cons_ = num_cons;
  std::vector<double> probabilities;
  mp::Expr obj_expr;
  con_indices_.resize(num_cons);
  core_con_indices_.resize(num_cons);
  if (num_stages > 1) {
    // Compute stage of each constraint as a maximum of stages of
    // variables in it.
    std::vector<int> con_stages(num_cons);
    for (int j = 0; j < num_vars; ++j) {
      if (core_var_indices_[j] < num_stage1_vars)
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
      con_indices_[index] = i;
      core_con_indices_[i] = index++;
    }
    assert(stage1_index == num_stage1_cons_ && stage2_index == num_cons);

    // Get the number of scenarios from the expectation in the objective.
    CallExpr expr;
    if (p.num_objs() > 0) {
      if (auto e = p.obj(0).nonlinear_expr())
        expr = Cast<CallExpr>(e);
    }
    if (!expr)
      throw Error("objective doesn't contain expectation");
    int num_scenarios = expr.num_args() - 1;
    probabilities.resize(num_scenarios);
    // Get probabilities.
    for (int i = 0; i < num_scenarios; ++i) {
      auto prob = Cast<mp::NumericConstant>(expr.arg(i));
      if (!prob)
        throw Error("probability is not constant");
      probabilities[i] = prob.value();
    }
    obj_expr = expr.arg(num_scenarios);
    // TODO: check that second-stage variables only occur in expectation
  } else {
    for (int i = 0; i < num_vars; ++i) {
      var_indices_[i] = i;
      core_var_indices_[i] = i;
    }
    for (int i = 0; i < num_cons; ++i) {
      con_indices_[i] = i;
      core_con_indices_[i] = i;
    }
    probabilities.push_back(1);
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
      base_rhs[i] = GetConRHSAndType(p, con_indices_[i], type);
      writer.Write(" {}  R{}\n", type, i + 1);
    }

    std::vector<double> core_obj_coefs(num_vars);
    if (p.num_objs() != 0) {
      // Get objective coefficients in the core problem (first scenario).
      AffineExprExtractor extractor;
      extractor.Visit(obj_expr);
      for (auto term: p.obj(0).linear_expr())
        core_obj_coefs[term.var_index()] = term.coef();
      for (auto term: extractor.linear_expr())
        core_obj_coefs[term.var_index()] += term.coef();
    }

    core_rhs = base_rhs;
    GetScenario(p, 0, core_coefs, core_rhs);
    WriteColumns(writer, p, num_cons, core_obj_coefs, core_coefs);

    writer.Write("RHS\n");
    for (int i = 0; i < num_cons; ++i)
      writer.Write("    RHS1      R{:<7}  {}\n", i + 1, core_rhs[i]);

    writer.Write("BOUNDS\n");
    double inf = std::numeric_limits<double>::infinity();
    for (int i = 0; i < num_vars; ++i) {
      auto var = p.var(i); // TODO: reorder by stages
      double lb = var.lb(), ub = var.ub();
      if (lb != 0)
        writer.Write(" LO BOUND1      C{:<7}  {}\n", i + 1, lb);
      if (ub < inf)
        writer.Write(" UP BOUND1      C{:<7}  {}\n", i + 1, ub);
    }

    writer.Write("ENDATA\n");
  }

  // Write the .sto file.
  {
    FileWriter writer(smps_basename + ".sto");
    writer.Write(
      "STOCH         PROBLEM\n"
      "SCENARIOS     DISCRETE\n");
    if (num_stages > 1) {
      std::vector<double> coefs, rhs;
      writer.Write(" SC SCEN1     'ROOT'    {:<12}   T1\n", probabilities[0]);
      for (size_t s = 1, num_scen = probabilities.size(); s < num_scen; ++s) {
        writer.Write(" SC SCEN{:<4}  SCEN1     {:<12}   T2\n",
            s + 1, probabilities[s]);
        rhs = base_rhs;
        GetScenario(p, s, coefs, rhs);
        // Compare to the core and write differences.
        for (int j = 0, num_vars = p.num_vars(); j < num_vars; ++j) {
          for (int k = p.col_start(j), n = p.col_start(j + 1); k != n; ++k) {
            double coef = coefs[k];
            if (coef == core_coefs[k]) continue;
            int core_con_index = core_con_indices_[p.row_index(k)];
            writer.Write("    C{:<7}  R{:<7}  {}\n",
                core_var_indices_[j] + 1, core_con_index + 1, coef);
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
