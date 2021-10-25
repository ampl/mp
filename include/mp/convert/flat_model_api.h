/*
 Flat model API.
 This is how Backend receives a flat model.
 TODO Separate NL entities from Flat objects.

 Copyright (C) 2021 AMPL Optimization Inc

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
*/
#ifndef FLAT_MODEL_API_H_
#define FLAT_MODEL_API_H_

#include "mp/convert/flat_model_api_basic.h"
#include "mp/convert/std_constr.h"
#include "mp/convert/std_obj.h"
#include "mp/convert/model.h"

namespace mp {

template <class Impl>
class FlatModelAPI : public BasicFlatModelAPI {
  ////////////////////////////////////////////////////////////////////////
  ///////////////////////////// MODEL MANIP //////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  /// TODO To be moved into Converter
public:
  using BaseBackend = BasicFlatModelAPI;
  using Model = BasicModel<>;
  using Variable = typename Model::Variable;

  void AddVariable(Variable ) {
    throw MakeUnsupportedError("BasicBackend::AddVariable");
  }
  void AddCommonExpression(Problem::CommonExpr ) {
    throw MakeUnsupportedError("BasicBackend::AddCommonExpressions");
  }
  void AddLogicalConstraint(Problem::LogicalCon ) {
    throw MakeUnsupportedError("BasicBackend::AddLogicalConstraints");
  }

  void AddObjective(typename Model::Objective obj) {
    if (obj.nonlinear_expr()) {
      MP_DISPATCH( AddGeneralObjective( obj ) );
    } else {
      LinearExprUnzipper leu(obj.linear_expr());
      LinearObjective lo { obj.type(),
            std::move(leu.c_), std::move(leu.v_) };
      if (nullptr==obj.p_extra_info()) {
        MP_DISPATCH( SetLinearObjective( obj.index(), lo ) );
      } else {
        auto qt = obj.p_extra_info()->qt_;
        assert(!qt.empty());
        MP_DISPATCH( SetQuadraticObjective( obj.index(),
                       QuadraticObjective{std::move(lo), std::move(qt)} ) );
      }
    }
  }
  void AddGeneralObjective(typename Model::Objective ) {
    throw MakeUnsupportedError("BasicBackend::AddGeneralObjective");
  }
  void SetLinearObjective( int, const LinearObjective& ) {
    throw MakeUnsupportedError("BasicBackend::AddLinearObjective");
  }
  void SetQuadraticObjective( int, const QuadraticObjective& ) {
    throw MakeUnsupportedError("BasicBackend::AddQuadraticObjective");
  }

  void AddAlgebraicConstraint(typename Model::AlgebraicCon con) {
    if (con.nonlinear_expr()) {
      MP_DISPATCH( AddGeneralConstraint( con ) );
    } else {
      LinearExprUnzipper leu(con.linear_expr());
      auto lc = LinearConstraint{
          std::move(leu.c_), std::move(leu.v_),
          con.lb(), con.ub() };
      if (nullptr==con.p_extra_info()) {
        MP_DISPATCH( AddConstraint( lc ) );
        orig_lin_constr_.push_back(n_alg_constr_);
      } else {
        auto qt = con.p_extra_info()->qt_;
        assert(!qt.empty());
        MP_DISPATCH( AddConstraint( QuadraticConstraint{std::move(lc), std::move(qt)} ) );
      }
      ++n_alg_constr_;
    }
  }

  void AddGeneralConstraint(typename Model::AlgebraicCon ) {
    throw MakeUnsupportedError("BasicBackend::AddGeneralConstraint");
  }

  ////////////////// Some basic custom constraints /////////////////
  USE_BASE_CONSTRAINT_HANDLERS(BasicFlatModelAPI)

  /// Optionally exclude LDCs from being posted,
  /// then all those are converted to LinearConstraint's first
  ACCEPT_CONSTRAINT(LinearDefiningConstraint, NotAccepted)
  void AddConstraint(const LinearDefiningConstraint& ldc) {
    MP_DISPATCH( AddConstraint(ldc.to_linear_constraint()) );
  }

  ACCEPT_CONSTRAINT(LinearConstraint, Recommended)
  /// TODO Attributes (lazy/user cut, etc)
  void AddConstraint(const LinearConstraint& ) {
    throw MakeUnsupportedError("BasicBackend::AddLinearConstraint");
  }


  /// Convenience method
  /// Gurobi reports duals separately for linear and QCP constraints
  /// We rely on QCP ones coming first in NL
  static
      std::vector<double> MakeDualsFromLPAndQCPDuals(
        std::vector<double> pi, std::vector<double> qcpi) {
    qcpi.insert(qcpi.end(), pi.begin(), pi.end());
    return qcpi;
  }

  /// Gurobi handles linear constraints as a separate class.
  /// AMPL provides suffixes for all constraints together.
  /// The method returns the indexes of linear constraints
  /// which have suffixes, in the overall constraints list.
  const std::vector<size_t>& GetIndexesOfLinearConstraintsWithSuffixes() const
  { return orig_lin_constr_; }

private:
  /// Indices of NL original linear constr in the total constr ordering
  std::vector<size_t> orig_lin_constr_;
  size_t n_alg_constr_=0;

};

} // namespace mp

#endif // FLAT_MODEL_API_H_
