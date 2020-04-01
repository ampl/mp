/*
 Abstract solver backend wrapper and
 interfaces between converters and solver backends.

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

#ifndef BACKEND_H_
#define BACKEND_H_

#include "mp/problem.h"
#include "mp/convert/constraint.h"

namespace mp {

/// Converting linear expr to 2 vectors.
/// TODO Change mp::LinearExpr to store them right away.
struct LinearExprUnzipper {
  std::vector<double> c_;
  std::vector<int> v_;
  LinearExprUnzipper(const LinearExpr& e) {
    for (LinearExpr::const_iterator it=e.begin(); it!=e.end(); ++it) {
      c_.push_back(it->coef());
      v_.push_back(it->var_index());
    }
  }
};

/// Basic backend wrapper.
/// Used by converter to directly access a solver.
/// The basic wrapper provides common functionality: option handling
/// and placeholders for solver API
template <class Impl>
class BasicBackend : public BasicConstraintAdder {
public:
  /// [[ TODO Incrementality ]]

  /// TODO InitProblemModificationPhase: demand redefinition in concrete backend?
  /// TODO Remove the function stubs, just don't compile when not defined?
  void InitProblemModificationPhase(const Problem& p) { }  // TODO Get rid of Problem here
  void FinishProblemModificationPhase() { }
  void AddVariables(int n, double* lbs, double* ubs, mp::var::Type* types) {
    throw MakeUnsupportedError("BasicBackend::AddVariables");
  }
  void AddCommonExpressions(int n, Problem::CommonExpr* cexprs) {
    throw MakeUnsupportedError("BasicBackend::AddCommonExpressions");
  }
  void AddLogicalConstraints(int n, Problem::LogicalCon* lcons) {
    throw MakeUnsupportedError("BasicBackend::AddLogicalConstraints");
  }

  /// We distinguish linear-only vs general objectives
  void AddObjectives(int n, Problem::Objective* objs) {
    for (int i=0; i!=n; ++i) {
      if (!objs[i].nonlinear_expr()) {
        LinearExprUnzipper leu(objs[i].linear_expr());
        /// Add to the concrete interface
        MP_DISPATCH( AddLinearObjective( objs[i].type(), leu.c_.size(),
                                         leu.c_.data(), leu.v_.data()) );
      } else {
        MP_DISPATCH( AddGeneralObjective( objs[i] ) );
      }
    }
  }
  void AddLinearObjective( obj::Type sense, int nnz,
                           const double* c, const int* v) {
    throw MakeUnsupportedError("BasicBackend::AddLinearObjective");
  }
  void AddGeneralObjective(Problem::Objective obj) {
    throw MakeUnsupportedError("BasicBackend::AddGeneralObjective");
  }

  /// We distinguish linear-only vs general algebraic constraint
  void AddAlgebraicConstraints(int n, Problem::AlgebraicCon* cons) {
    for (int i=0;i!=n;++i) {
      if (!cons[i].nonlinear_expr()) {
        LinearExprUnzipper leu(cons[i].linear_expr());
        /// Add to the concrete interface
        MP_DISPATCH( AddLinearConstraint(leu.c_.size(), leu.c_.data(), leu.v_.data(),
                                         cons[i].lb(), cons[i].ub()) );
      } else {
        MP_DISPATCH( AddGeneralConstraint( cons[i] ) );
      }
    }
  }
  /// TODO Do we need ability to add several at once?
  /// TODO Attributes (lazy/user cut, etc)
  void AddLinearConstraint(int nnz, const double* c, const int* v,
                           double lb, double ub) {
    throw MakeUnsupportedError("BasicBackend::AddLinearConstraint");
  }
  void AddGeneralConstraint(Problem::AlgebraicCon con) {
    throw MakeUnsupportedError("BasicBackend::AddGeneralConstraint");
  }
};


template <class Model, class Backend>
class ModelToBackendFeeder {
  const Model& model_;
  Backend& backend_;
public:
  ModelToBackendFeeder(const Model& p, Backend& i)
    : model_(p), backend_(i)  { }

  void PushWholeProblem_noCustomConstraints() {
    InitProblemModificationPhase();
    PushStandardItems();
    FinishProblemModificationPhase();
  }
  void PushWholeProblem() {
    InitProblemModificationPhase();
    PushStandardItems();
    PushCustomConstraints();
    FinishProblemModificationPhase();
  }

protected:
  void PushStandardItems() {
    PushVariables();
    PushCommonSubExpr();
    PushObjectives();
    PushAlgebraicConstraints();
    PushLogicalConstraints();
  }

  void InitProblemModificationPhase() {
    backend_.InitProblemModificationPhase(model_);       // TODO remove problem_ here
  }

  void PushVariables() {
    int num_vars = model_.num_vars();
    for (int j = 0; j < num_vars; ++j) {
      typename Model::Variable var = model_.var(j);
      double lb = var.lb();
      double ub = var.ub();
      var::Type ty = var.type();
      backend_.AddVariables(1, &lb, &ub, &ty);
    }
  }

  void PushCommonSubExpr() {
    int num_common_exprs = model_.num_common_exprs();
    for (int i = 0; i < num_common_exprs; ++i) {
      typename Model::CommonExpr expr = model_.common_expr(i);
      backend_.AddCommonExpressions(1, &expr);
    }
  }

  void PushObjectives() {
    if (int num_objs = model_.num_objs()) {
      for (int i = 0; i < num_objs; ++i) {
        typename Model::Objective obj = model_.obj(i);
        backend_.AddObjectives(1, &obj);
      }
    }
  }

  void PushAlgebraicConstraints() {
    if (int n_cons = model_.num_algebraic_cons()) {
      for (int i = 0; i < n_cons; ++i) {
        typename Model::AlgebraicCon con = model_.algebraic_con(i);
        backend_.AddAlgebraicConstraints(1, &con);
      }
    }
  }

  void PushLogicalConstraints() {
    if (int n_lcons = model_.num_logical_cons()) {
      for (int i = 0; i < n_lcons; ++i) {
        typename Model::LogicalCon con = model_.logical_con(i);
        backend_.AddLogicalConstraints(1, &con);
      }
    }
  }

  void PushCustomConstraints() {
    if (int n_ccons = model_.num_custom_cons()) {
      for (int i = 0; i < n_ccons; ++i) {
        model_.custom_con(i)->AddToBackend(backend_);
      }
    }
  }

  void FinishProblemModificationPhase() {
    backend_.FinishProblemModificationPhase();
  }
};

}  // namespace mp

#endif  // BACKEND_H_
