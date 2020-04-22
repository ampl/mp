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
#include "mp/convert/constraint_keeper.h"
#include "mp/convert/std_constr.h"

namespace mp {

/// Basic backend wrapper.
/// Used by converter to directly access a solver.
/// The basic wrapper provides common functionality: option handling
/// and placeholders for solver API
template <class Impl>
class BasicBackend : public BasicConstraintAdder {
public:
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

  ////////////////// Some basic custom constraints /////////////////
  USE_BASE_CONSTRAINT_HANDLERS(BasicConstraintAdder)

  /// Optionally exclude LFDs from being posted,
  /// then all those are converted to LinearConstraint's first
  ACCEPT_CONSTRAINT(LinearDefiningConstraint, NotAccepted)
  void AddConstraint(const LinearDefiningConstraint& ldc) {
    MP_DISPATCH( AddConstraint(ldc.to_linear_constraint()) );
  }

  ACCEPT_CONSTRAINT(LinearConstraint, Recommended)
  void AddConstraint(const LinearConstraint& ldc) {     // TODO make this primary
    MP_DISPATCH( AddLinearConstraint(ldc.nnz(), ldc.coefs(), ldc.vars(),
                                     ldc.lb(), ldc.ub()) );
  }


};


}  // namespace mp

#endif  // BACKEND_H_
