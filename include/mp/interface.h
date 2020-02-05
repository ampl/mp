/*
 An abstract interface between modeling system and solver backend.

 Copyright (C) 2020 AMPL Optimization Inc

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

 Author: Gleb Belov <gleb.belov@monash.edu>
 */

#ifndef INTERFACE_H_
#define INTERFACE_H_


namespace mp {

template <class Problem, class Interface>
class ProblemToInterfaceFeeder {
  const Problem& problem_;
  Interface& interface_;
public:
  ProblemToInterfaceFeeder(const Problem& p, Interface& i)
    : problem_(p), interface_(i)  { }

  void PushWholeProblem() {
    InitProblemModificationPhase();
    PushVariables();
    PushCommonSubExpr();
    PushObjectives();
    PushAlgebraicConstraints();
    PushLogicalConstraints();
    FinishProblemModificationPhase();
  }

  void InitProblemModificationPhase() {
    interface_.InitProblemModificationPhase(problem_);       // TODO remove problem_ here
  }

  void PushVariables() {
    int num_vars = problem_.num_vars();
    for (int j = 0; j < num_vars; ++j) {
      typename Problem::Variable var = problem_.var(j);
      double lb = var.lb();
      double ub = var.ub();
      var::Type ty = var.type();
      interface_.AddVariables(1, &lb, &ub, &ty);
    }
  }

  void PushCommonSubExpr() {
    int num_common_exprs = problem_.num_common_exprs();
    for (int i = 0; i < num_common_exprs; ++i) {
      typename Problem::CommonExpr expr = problem_.common_expr(i);
      interface_.AddCommonExpressions(1, &expr);
    }
  }

  void PushObjectives() {
    if (int num_objs = problem_.num_objs()) {
      for (int i = 0; i < num_objs; ++i) {
        typename Problem::Objective obj = problem_.obj(i);
        interface_.AddObjectives(1, &obj);
      }
    }
  }

  void PushAlgebraicConstraints() {
    if (int n_cons = problem_.num_algebraic_cons()) {
      for (int i = 0; i < n_cons; ++i) {
        typename Problem::AlgebraicCon con = problem_.algebraic_con(i);
        interface_.AddAlgebraicConstraints(1, &con);
      }
    }
  }

  void PushLogicalConstraints() {
    if (int n_lcons = problem_.num_logical_cons()) {
      for (int i = 0; i < n_lcons; ++i) {
        typename Problem::LogicalCon con = problem_.logical_con(i);
        interface_.AddLogicalConstraints(1, &con);
      }
    }
  }

  void FinishProblemModificationPhase() {
    interface_.FinishProblemModificationPhase();
  }
};

}  // namespace mp

#endif  // INTERFACE_H_
