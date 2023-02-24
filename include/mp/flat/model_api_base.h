/*
 Basic flat model API definitions.

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
#ifndef FLAT_MODEL_API_BASE_H_
#define FLAT_MODEL_API_BASE_H_

/**
  * Basic definitions for a FlatModelAPI
  */

#include <string>

#include "mp/arrayref.h"
#include "mp/common.h"
#include "mp/flat/constr_base.h"
#include "mp/flat/obj_std.h"
#include "mp/flat/model_info.h"

namespace mp {

/// Define an array of variables
class VarArrayDef {
  ArrayRef<double> lbs_;
  ArrayRef<double> ubs_;
  ArrayRef<var::Type> types_;
  ArrayRef<const char*> names_;
public:
  VarArrayDef() = default;
  template <class BndVec, class TypeVec,
            class NameVec=ArrayRef<const char*> >
  VarArrayDef(BndVec&& lbs, BndVec&& ubs,
              TypeVec&& tys, NameVec&& nms={}) :
    lbs_(std::forward<BndVec>(lbs)),
    ubs_(std::forward<BndVec>(ubs)),
    types_(std::forward<TypeVec>(tys)),
    names_(std::forward<NameVec>(nms)) { }
  VarArrayDef(std::initializer_list<double> lbs,
              std::initializer_list<double> ubs,
              std::initializer_list<var::Type> tys,
              std::initializer_list<const char*> nms = {}) :
    lbs_((lbs)), ubs_((ubs)), types_((tys)), names_(nms) { }
  int size() const { assert(check()); return lbs_.size(); }
  const double* plb() const { return lbs_.data(); }
  const double* pub() const { return ubs_.data(); }
  const var::Type* ptype() const { return types_.data(); }
  void set_lb_ub_types(ArrayRef<double> lbs, ArrayRef<double> ubs,
                       ArrayRef<var::Type> types) {
    lbs_=lbs; ubs_=ubs; types_=types;
    assert(check());
  }
  const char* const* pnames() const { return names_.data(); }
  bool check() const
  { return ubs_.size()==lbs_.size() &&
        types_.size()==lbs_.size(); }
};


/// Level of acceptance of a constraint by a backend
enum ConstraintAcceptanceLevel {
  NotAccepted=0,
  AcceptedButNotRecommended=1,
  Recommended=2
};


/// Constraint groups
///
/// This is used to access constraint attributes (basis status, duals, ...)
/// Convenient when the solver accesses constraint attributes in groups
/// For example, Gurobi 9.5 has linear, quadratic, SOS, and general
enum ConstraintGroup {
  CG_Default,
  CG_Linear,
  CG_Quadratic,
  CG_General,
  CG_Piecewiselinear,
  CG_SOS,
  CG_SOS1,
  CG_SOS2,
  CG_Logical
};


/// ModelAPIs handling custom flat constraints should derive from
class BasicFlatModelAPI {
public:
  /// Placeholder for GetTypeName()
  static const char* GetTypeName()    { return "BasicBackendFlatModelAPI"; }
  /// Placeholder for GetLongName()
  static const char* GetLongName() { return nullptr; }

  /// Placeholder for InitCustomOptions()
  void InitCustomOptions() { }

  /// Pass on a FlatModelInfo object
  void PassFlatModelInfo(std::unique_ptr<FlatModelInfo>&& pfmi) {
    pfmi_ = std::move(pfmi);
  }

  /// Retrieve FlatModelInfo*
  const FlatModelInfo* GetFlatModelInfo() const { return pfmi_.get(); }

  /// Chance to prepare problem update,
  /// e.g., allocate storage
  void InitProblemModificationPhase(const FlatModelInfo*) {  }
  /// Chance to end problem update
  void FinishProblemModificationPhase() {  }

  ////////////////// Some standard items /////////////////

  /// Placeholder for SetLinearObjective()
  void SetLinearObjective(int , const LinearObjective& ) {
    MP_UNSUPPORTED("FlatModelAPI::SetLinearObjective()");
  }

  /// Whether accepting quadratic objectives:
  /// 0 - no, 1 - convex, 2 - nonconvex
  static int AcceptsQuadObj() { return 0; }

  /// Placeholder for SetQuadraticObjective()
  void SetQuadraticObjective(int , const QuadraticObjective& ) {
    MP_UNSUPPORTED("FlatModelAPI::SetQuadraticObjective()");
  }

  /// Placeholder for AddConstraint<>()
  template <class Constraint>
  void AddConstraint(const Constraint& ) {
    MP_RAISE(
          std::string("Not handling constraint type '") +
          Constraint::GetTypeName() +
          "'. Provide a handler or a converter method");
  }
  /// Derived backends have to tell C++ to use default handlers if they are needed
  /// when they overload AddConstraint(), due to C++ name hiding
#define USE_BASE_CONSTRAINT_HANDLERS(BaseBackend) \
  using BaseBackend::AddConstraint; \
  using BaseBackend::AcceptanceLevel; \
  using BaseBackend::GroupNumber;

  /// Default constraint group
  static constexpr ConstraintGroup GroupNumber(const BasicConstraint*) {
    return CG_Default;
  }

  /// By default, we say constraint XYZ is not accepted but...
  static constexpr ConstraintAcceptanceLevel AcceptanceLevel(const BasicConstraint*) {
    return NotAccepted;
  }

  /// Specifically, ask if the solver accepts non-convex quadratic constraints
  static constexpr bool AcceptsNonconvexQC() { return false; }


private:
  std::unique_ptr< FlatModelInfo > pfmi_ { nullptr };
};


/// ... then for a certain constraint it can be specified
#define ACCEPT_CONSTRAINT(ConstrType, level, con_grp) \
  static mp::ConstraintAcceptanceLevel \
    AcceptanceLevel(const ConstrType*) \
  { return (mp::ConstraintAcceptanceLevel)level; } \
  static constexpr int \
    GroupNumber(const ConstrType*) { return con_grp; }

} // namespace mp

#endif // FLAT_MODEL_API_BASE_H_
