/*
 A C++ interface to an AMPL problem.

 Copyright (C) 2012 AMPL Optimization LLC

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization LLC disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Victor Zverovich
 */

#include "solvers/util/problem.h"

using ampl::NumericConstant;
using ampl::RelationalExpr;
using ampl::UnaryExpr;

namespace {

// Writes a linear term.
void WriteTerm(fmt::Writer &w, double coef, unsigned var) {
  if (coef != 1)
    w << coef << " * ";
  w << "x" << (var + 1);
}

class ExprPrinter : public ampl::ExprVisitor<ExprPrinter, void, void> {
 private:
  fmt::Writer &writer_;

 private:
  void Print(RelationalExpr e, const char *op) {
    Visit(e.lhs());
    writer_ << ' ' << op << ' ';
    Visit(e.rhs());
  }

 public:
  ExprPrinter(fmt::Writer &w) : writer_(w) {}

  /*void VisitPlus(BinaryExpr e) {
     // TODO
  }

  void VisitMinus(BinaryExpr e) {
     // TODO
  }

  void VisitMult(BinaryExpr e) {
     // TODO
  }

  void VisitDiv(BinaryExpr e) {
     // TODO
  }

  void VisitRem(BinaryExpr e) {
     // TODO
  }

  void VisitPow(BinaryExpr e) {
     // TODO
  }

  void VisitNumericLess(BinaryExpr e) {
     // TODO
  }

  void VisitMin(VarArgExpr e) {
     // TODO
  }

  void VisitMax(VarArgExpr e) {
     // TODO
  }

  void VisitFloor(UnaryExpr e) {
     // TODO
  }

  void VisitCeil(UnaryExpr e) {
     // TODO
  }

  void VisitAbs(UnaryExpr e) {
     // TODO
  }*/

  void VisitUnaryMinus(UnaryExpr e) {
    writer_ << '-';
    Visit(e.arg());
  }

  void VisitIf(ampl::IfExpr e) {
    writer_ << "if ";
    Visit(e.condition());
    writer_ << " then ";
    Visit(e.true_expr());
    ampl::NumericExpr false_expr = e.false_expr();
    NumericConstant c = ampl::Cast<NumericConstant>(false_expr);
    if (!c || c.value() != 0) {
      writer_ << " else ";
      Visit(false_expr);
    }
  }

  /*void VisitTanh(UnaryExpr e) {
     // TODO
  }

  void VisitTan(UnaryExpr e) {
     // TODO
  }

  void VisitSqrt(UnaryExpr e) {
     // TODO
  }

  void VisitSinh(UnaryExpr e) {
     // TODO
  }

  void VisitSin(UnaryExpr e) {
     // TODO
  }

  void VisitLog10(UnaryExpr e) {
     // TODO
  }

  void VisitLog(UnaryExpr e) {
     // TODO
  }

  void VisitExp(UnaryExpr e) {
     // TODO
  }

  void VisitCosh(UnaryExpr e) {
     // TODO
  }

  void VisitCos(UnaryExpr e) {
     // TODO
  }

  void VisitAtanh(UnaryExpr e) {
     // TODO
  }

  void VisitAtan2(BinaryExpr e) {
     // TODO
  }

  void VisitAtan(UnaryExpr e) {
     // TODO
  }

  void VisitAsinh(UnaryExpr e) {
     // TODO
  }

  void VisitAsin(UnaryExpr e) {
     // TODO
  }

  void VisitAcosh(UnaryExpr e) {
     // TODO
  }

  void VisitAcos(UnaryExpr e) {
     // TODO
  }*/

  void VisitSum(ampl::SumExpr e) {
    writer_ << "sum(";
    ampl::SumExpr::iterator i = e.begin(), end = e.end();
    if (i != end) {
      Visit(*i);
      for (++i; i != end; ++i) {
        writer_ << ", ";
        Visit(*i);
      }
    }
    writer_ << ")";
  }

  /*void VisitIntDiv(BinaryExpr e) {
     // TODO
  }

  void VisitPrecision(BinaryExpr e) {
     // TODO
  }

  void VisitRound(BinaryExpr e) {
     // TODO
  }

  void VisitTrunc(BinaryExpr e) {
     // TODO
  }

  void VisitCount(CountExpr e) {
     // TODO
  }

  void VisitNumberOf(NumberOfExpr e) {
     // TODO
  }

  void VisitPLTerm(PiecewiseLinearTerm t) {
     // TODO
  }

  void VisitPowConstExp(BinaryExpr e) {
     // TODO
  }

  void VisitPow2(UnaryExpr e) {
     // TODO
  }

  void VisitPowConstBase(BinaryExpr e) {
     // TODO
  }

  void VisitCall(CallExpr e) {
     // TODO
  }*/

  void VisitNumericConstant(NumericConstant c) { writer_ << c.value(); }
  void VisitVariable(ampl::Variable v) { writer_ << 'x' << (v.index() + 1); }

  /*void VisitOr(BinaryLogicalExpr e) {
     // TODO
  }

  void VisitAnd(BinaryLogicalExpr e) {
     // TODO
  }*/

  void VisitLess(RelationalExpr e) { Print(e, "<"); }
  void VisitLessEqual(RelationalExpr e) { Print(e, "<="); }
  void VisitEqual(RelationalExpr e) { Print(e, "="); }
  void VisitGreaterEqual(RelationalExpr e) { Print(e, ">="); }
  void VisitGreater(RelationalExpr e) { Print(e, ">"); }
  void VisitNotEqual(RelationalExpr e) { Print(e, "!="); }

  /*void VisitNot(NotExpr e) {
     // TODO
  }

  void VisitAtLeast(LogicalCountExpr e) {
     // TODO
  }

  void VisitAtMost(LogicalCountExpr e) {
     // TODO
  }

  void VisitExactly(LogicalCountExpr e) {
     // TODO
  }

  void VisitNotAtLeast(LogicalCountExpr e) {
     // TODO
  }

  void VisitNotAtMost(LogicalCountExpr e) {
     // TODO
  }

  void VisitNotExactly(LogicalCountExpr e) {
     // TODO
  }

  void VisitForAll(IteratedLogicalExpr e) {
     // TODO
  }

  void VisitExists(IteratedLogicalExpr e) {
     // TODO
  }

  void VisitImplication(ImplicationExpr e) {
     // TODO
  }

  void VisitIff(BinaryLogicalExpr e) {
     // TODO
  }

  void VisitAllDiff(AllDiffExpr e) {
     // TODO
  }

  void VisitLogicalConstant(LogicalConstant c) {
     // TODO
  }*/
};

template <typename LinearExpr>
void WriteExpr(fmt::Writer &w, LinearExpr linear, ampl::NumericExpr nonlinear) {
  bool has_terms = false;
  typedef typename LinearExpr::iterator Iterator;
  for (Iterator i = linear.begin(), e = linear.end(); i != e; ++i) {
    double coef = i->coef();
    if (coef != 0) {
      if (has_terms)
        w << " + ";
      else
        has_terms = true;
      WriteTerm(w, coef, i->var_index());
    }
  }
  if (!has_terms)
    w << "0";
  if (!nonlinear)
    return;
  NumericConstant c = ampl::Cast<NumericConstant>(nonlinear);
  if (c && c.value() == 0)
    return;
  w << " + ";
  ExprPrinter(w).Visit(nonlinear);
}
}

namespace ampl {

Solution::Solution()
: solve_code_(-1), num_vars_(0), num_cons_(0), values_(0), dual_values_(0) {}

Solution::~Solution() {
  free(values_);
  free(dual_values_);
}

void Solution::Swap(Solution &other) {
  std::swap(solve_code_, other.solve_code_);
  std::swap(num_vars_, other.num_vars_);
  std::swap(num_cons_, other.num_cons_);
  std::swap(values_, other.values_);
  std::swap(dual_values_, other.dual_values_);
}

void Solution::Read(fmt::StringRef stub, int num_vars, int num_cons) {
  // Allocate filename large enough to hold stub, ".sol" and terminating zero.
  std::vector<char> filename(stub.size() + 5);
  std::strcpy(&filename[0], stub.c_str());
  ASL asl = {};
  asl.i.n_var_ = num_vars;
  asl.i.n_con_ = num_cons;
  asl.i.ASLtype = 1;
  asl.i.filename_ = &filename[0];
  asl.i.stub_end_ = asl.i.filename_ + stub.size();
  Solution sol;
  sol.num_vars_ = num_vars;
  sol.num_cons_ = num_cons;
  char *message = read_sol_ASL(&asl, &sol.values_, &sol.dual_values_);
  if (!message)
    throw Error("Error reading solution file");
  free(message);
  Swap(sol);
  solve_code_ = asl.p.solve_code_;
}

void Problem::Free() {
  if (var_capacity_) {
    delete [] asl_->i.LUv_;
    delete [] asl_->i.Uvx_;
    delete [] var_types_;
    asl_->i.LUv_ = asl_->i.Uvx_ = 0;
    var_capacity_ = 0;
  }
  if (obj_capacity_) {
    delete [] asl_->I.obj_de_;
    delete [] asl_->i.objtype_;
    delete [] asl_->i.Ograd_;
    asl_->I.obj_de_ = 0;
    asl_->i.objtype_ = 0;
    asl_->i.Ograd_ = 0;
    obj_capacity_ = 0;
  }
  if (logical_con_capacity_) {
    delete [] asl_->I.lcon_de_;
    asl_->I.lcon_de_ = 0;
    logical_con_capacity_ = 0;
  }
}

Problem::Problem()
: asl_(reinterpret_cast<ASL_fg*>(ASL_alloc(ASL_read_fg))),
  var_capacity_(0), obj_capacity_(0), logical_con_capacity_(0), var_types_(0) {
}

Problem::~Problem() {
  Free();
  ASL_free(reinterpret_cast<ASL**>(&asl_));
}

fmt::Writer &operator<<(fmt::Writer &w, const Problem &p) {
  // Write variables.
  int num_vars = p.num_vars();
  for (int i = 0; i < num_vars; ++i) {
    w << "var x" << (i + 1);
    double lb = p.var_lb(i), ub = p.var_ub(i);
    if (lb == ub) {
      w << " = " << lb;
    } else {
      if (lb != -Infinity)
        w << " >= " << lb;
      if (ub != Infinity)
        w << " <= " << ub;
    }
    w << ";\n";
  }

  // Write objectives.
  for (int i = 0, n = p.num_objs(); i < n; ++i) {
    w << (p.obj_type(i) == MIN ? "minimize" : "maximize") << " o: ";
    WriteExpr(w, p.linear_obj_expr(i), p.nonlinear_obj_expr(i));
    w << ";\n";
  }

  // Write constraints.
  for (int i = 0, n = p.num_cons(); i < n; ++i) {
    w << "s.t. c" << (i + 1) << ": ";
    double lb = p.con_lb(i), ub = p.con_ub(i);
    if (lb != ub && lb != -Infinity && ub != Infinity)
      w << lb << " <= ";
    WriteExpr(w, p.linear_con_expr(i), p.nonlinear_con_expr(i));
    if (lb == ub)
      w << " = " << lb;
    else if (ub != Infinity)
      w << " <= " << ub;
    else if (lb != -Infinity)
      w << " >= " << lb;
    w << ";\n";
  }
  return w;
}

// A manager of temporary files.
class TempFiles : Noncopyable {
 private:
  char *name_;

 public:
  TempFiles() : name_(tempnam(0, 0)) {}
  ~TempFiles() {
    std::remove(c_str(fmt::Format("{}.nl") << name_));
    std::remove(c_str(fmt::Format("{}.sol") << name_));
    free(name_);
  }

  const char *stub() const { return name_; }
};

void Problem::AddVar(double lb, double ub, VarType type) {
  int &num_vars = asl_->i.n_var_;
  if (num_vars >= var_capacity_) {
    IncreaseCapacity(num_vars, var_capacity_);
    Grow(asl_->i.LUv_, num_vars, var_capacity_);
    Grow(asl_->i.Uvx_, num_vars, var_capacity_);
    if (var_types_)
      Grow(var_types_, num_vars, var_capacity_);
  }
  if (type != CONTINUOUS) {
    // Allocate var_types_ if this is the first integer variable added
    // after continuous.
    int num_integer_vars = Problem::num_integer_vars();
    if (!var_types_ && num_vars != num_integer_vars) {
      var_types_ = new VarType[var_capacity_];
      std::fill_n(var_types_, num_integer_vars, INTEGER);
      std::fill(var_types_ + num_integer_vars,
          var_types_ + num_vars, CONTINUOUS);
    }
    ++asl_->i.niv_;
  }
  asl_->i.LUv_[num_vars] = lb;
  asl_->i.Uvx_[num_vars] = ub;
  if (var_types_)
    var_types_[num_vars] = type;
  ++num_vars;
}

void Problem::AddObj(ObjType type, NumericExpr expr) {
  int &num_objs = asl_->i.n_obj_;
  if (num_objs >= obj_capacity_) {
    IncreaseCapacity(num_objs, obj_capacity_);
    Grow(asl_->I.obj_de_, num_objs, obj_capacity_);
    Grow(asl_->i.objtype_, num_objs, obj_capacity_);
    Grow(asl_->i.Ograd_, num_objs, obj_capacity_);
  }
  cde e = {expr.expr_};
  asl_->I.obj_de_[num_objs] = e;
  asl_->i.objtype_[num_objs] = type;
  asl_->i.Ograd_[num_objs] = 0;
  ++num_objs;
  ++asl_->i.nlo_;
}

void Problem::AddCon(LogicalExpr expr) {
  int &num_logical_cons = asl_->i.n_lcon_;
  if (num_logical_cons >= logical_con_capacity_) {
    IncreaseCapacity(num_logical_cons, logical_con_capacity_);
    Grow(asl_->I.lcon_de_, num_logical_cons, logical_con_capacity_);
  }
  cde e = {expr.expr_};
  asl_->I.lcon_de_[num_logical_cons] = e;
  ++num_logical_cons;
}

void Problem::Read(fmt::StringRef stub) {
  Free();
  ASL *asl = reinterpret_cast<ASL*>(asl_);
  FILE *nl = jac0dim_ASL(asl, const_cast<char*>(stub.c_str()),
      static_cast<ftnlen>(stub.size()));
  efunc *r_ops_int[N_OPS];
  for (int i = 0; i < N_OPS; ++i)
    r_ops_int[i] = reinterpret_cast<efunc*>(i);
  asl_->I.r_ops_ = r_ops_int;
  asl_->p.want_derivs_ = 0;
  fg_read_ASL(asl, nl, ASL_allow_CLP | ASL_sep_U_arrays);
  asl_->I.r_ops_ = 0;
}

void Problem::WriteNL(fmt::StringRef stub, ProblemChanges *pc, unsigned flags) {
  int nfunc = asl_->i.nfunc_;
  if ((flags & IGNORE_FUNCTIONS) != 0)
    asl_->i.nfunc_ = 0;
  int result = fg_write_ASL(reinterpret_cast<ASL*>(asl_),
      stub.c_str(), pc ? pc->vco() : 0, ASL_write_ASCII);
  asl_->i.nfunc_ = nfunc;
  if (result)
    throw Error("Error writing .nl file");
}

void Problem::Solve(fmt::StringRef solver_name,
    Solution &sol, ProblemChanges *pc, unsigned flags) {
  TempFiles temp;
  WriteNL(temp.stub(), pc, flags);
  // Run the solver and read the solution file.
  int exit_code = std::system(
      c_str(fmt::Format("{} {} -AMPL") << solver_name.c_str() << temp.stub()));
  if (exit_code != 0) {
    throw Error(fmt::Format("Error running solver {}, exit code = {}")
        << solver_name.c_str() << exit_code);
  }
  sol.Read(temp.stub(), num_vars() + (pc ? pc->num_vars() : 0),
      num_cons() + (pc ? pc->num_cons() : 0));
}

NewVCO *ProblemChanges::vco() {
  static double dummy;
  vco_.nnv = static_cast<int>(var_lb_.size());
  vco_.nnc = static_cast<int>(cons_.size());
  vco_.nno = static_cast<int>(objs_.size());
  if (!var_lb_.empty()) {
    vco_.LUnv = &var_lb_[0];
    vco_.Unv = &var_ub_[0];
  } else {
    vco_.LUnv = &dummy;
  }
  if (!cons_.empty()) {
    vco_.LUnc = &con_lb_[0];
    vco_.Unc = &con_ub_[0];
    vco_.newc = &cons_[0];
  } else {
    vco_.LUnc = &dummy;
  }
  if (!objs_.empty()) {
    vco_.newo = &objs_[0];
    vco_.ot = &obj_types_[0];
  }
  return &vco_;
}

void ProblemChanges::AddObj(
    ObjType type, unsigned size, const double *coefs, const int *vars) {
  std::size_t start = obj_terms_.size();
  obj_terms_.resize(start + size);
  ograd dummy;
  ograd *prev = &dummy;
  for (unsigned i = 0; i < size; ++i) {
    ograd &term = obj_terms_[start + i];
    term.coef = coefs[i];
    term.varno = vars[i];
    prev->next = &term;
    prev = &term;
  }
  objs_.push_back(&obj_terms_[start]);
  obj_types_.push_back(type);
}

void ProblemChanges::AddCon(const double *coefs, double lb, double ub) {
  con_lb_.push_back(lb);
  con_ub_.push_back(ub);
  std::size_t start = con_terms_.size();
  std::size_t num_vars = problem_->num_vars() + var_lb_.size();
  con_terms_.resize(start + num_vars);
  ograd dummy;
  ograd *prev = &dummy;
  for (std::size_t i = 0; i < num_vars; ++i) {
    ograd &term = con_terms_[start + i];
    term.coef = coefs[i];
    term.varno = i;
    prev->next = &term;
    prev = &term;
  }
  cons_.push_back(&con_terms_[start]);
}
}
