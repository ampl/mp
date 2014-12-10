/*
 CPLEX solver with SOCP transformations

 Copyright (C) 2014 AMPL Optimization Inc

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

#include "mp/nl.h"
#include "mp/problem.h"
#include "asl/aslbuilder.h"
#include "asl.h"

// Adapts Problem interface for use with .nl reader.
class ProblemBuilder : public mp::Problem {
 public:
  typedef mp::Function Function;
  typedef mp::Expr Expr;
  typedef mp::NumericExpr NumericExpr;
  typedef mp::LogicalExpr LogicalExpr;
  typedef mp::CountExpr CountExpr;
  typedef mp::Variable Variable;

  typedef IteratedExprBuilder NumericExprBuilder;
  typedef IteratedExprBuilder VarArgExprBuilder;
  typedef IteratedExprBuilder NumberOfExprBuilder;

  VarArgExprBuilder BeginVarArg(mp::expr::Kind kind, int num_args) {
    return BeginIterated(kind, num_args);
  }
  NumericExpr EndVarArg(VarArgExprBuilder builder) {
    return EndIterated(builder);
  }

  NumericExprBuilder BeginSum(int num_args) {
    return BeginIterated(mp::expr::SUM, num_args);
  }
  NumericExpr EndSum(NumericExprBuilder builder) {
    return EndIterated(builder);
  }

  struct ColumnSizeHandler {
    void Add(int) {
      // Ignore column sizes as the constraints are stored row-wise.
    }
  };

  // Returns a handler that receives column sizes in Jacobian.
  ColumnSizeHandler GetColumnSizeHandler() {
    return ColumnSizeHandler();
  }
};

// Detects if a problem is convertible to an SOCP.
class SOCPDetector {
 // TODO
};

class SOCPConverter {
 private:
  ProblemBuilder problem_;

  typedef mp::asl::internal::ASLBuilder ASLBuilder;
  ASLBuilder builder_;

 public:
  explicit SOCPConverter(ASL *asl) : builder_(asl) {}

  void Run(const char *stub);

  // Converts the problem into ASL format.
  void ConvertToASL();
};

void SOCPConverter::Run(const char *stub) {
  mp::ProblemBuilderToNLAdapter<ProblemBuilder> adapter(problem_);
  ReadNLFile(fmt::format("{}.nl", stub), adapter);
  if (!problem_.HasComplementarity()) {
    // TODO
    // 1. check if the problem can be converted to SOCP
    // 2. convert to SOCP
  }
  mp::ProblemInfo info = mp::ProblemInfo();
  info.num_vars = problem_.num_vars();
  info.num_objs = problem_.num_objs();
  info.num_algebraic_cons = problem_.num_algebraic_cons();
  for (int i = 0, n = problem_.num_objs(); i < n; ++i)
    info.num_obj_nonzeros += problem_.obj(i).linear_expr().num_terms();
  // TODO: convert all problem info
  builder_.SetInfo(info);
  builder_.set_stub(stub);
}

void SOCPConverter::ConvertToASL() {
  // Convert variables.
  int num_vars = problem_.num_vars();
  for (int i = 0; i < num_vars; ++i) {
    mp::Problem::Variable var = problem_.var(i);
    builder_.AddVar(var.lb(), var.ub(), var.type());
  }

  // Convert objectives.
  for (int i = 0, n = problem_.num_objs(); i < n; ++i) {
    mp::Problem::Objective obj = problem_.obj(i);
    mp::LinearExpr expr = obj.linear_expr();
    ASLBuilder::LinearObjBuilder obj_builder =
        builder_.AddObj(mp::obj::MIN, mp::asl::NumericExpr(), expr.num_terms());
    for (mp::LinearExpr::iterator i = expr.begin(), e = expr.end(); i != e; ++i)
      obj_builder.AddTerm(i->var_index(), i->coef());
    // TODO: handle nonlinear part of objective expression
  }

  // Convert constraints.
  for (int i = 0, n = problem_.num_algebraic_cons(); i < n; ++i) {
    mp::Problem::AlgebraicCon con = problem_.algebraic_con(i);
    builder_.AddCon(con.lb(), con.ub(), mp::asl::NumericExpr(), 0);
    // TODO: handle nonlinear part of constraint expression
  }

  ASLBuilder::ColumnSizeHandler cols = builder_.GetColumnSizeHandler();
  for (int i = 1; i < num_vars; ++i)
    cols.Add(0);
  builder_.EndBuild();
  // TODO
}

extern "C" void *socp_jac0dim(ASL *asl, const char *stub, ftnlen) {
  SOCPConverter *converter = new SOCPConverter(asl); // TODO: smart pointer
  converter->Run(stub);
  return converter;
}

extern "C" int socp_qp_read(ASL *, void *converter_ptr, int) {
  SOCPConverter *converter = static_cast<SOCPConverter*>(converter_ptr);
  converter->ConvertToASL();
  return 0;
}

extern "C" void socp_write_sol(
    ASL *asl, const char *msg, double *x, double *y, Option_Info *oi) {
  // TODO: convert and write solution
  write_sol_ASL(asl, msg, x, y, oi);
}
