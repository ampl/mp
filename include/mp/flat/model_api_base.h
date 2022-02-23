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

#include <string>
#include <stdexcept>

#include "mp/arrayref.h"
#include "mp/common.h"
#include "mp/flat/constraint_base.h"
#include "mp/flat/obj_std.h"

namespace mp {

/// Define an array of variables
class VarArrayDef {
  ArrayRef<double> lbs_;
  ArrayRef<double> ubs_;
  ArrayRef<var::Type> types_;
public:
  VarArrayDef() = default;
  template <class BndVec, class TypeVec>
  VarArrayDef(BndVec&& lbs, BndVec&& ubs, TypeVec&& tys) :
    lbs_(std::forward<BndVec>(lbs)),
    ubs_(std::forward<BndVec>(ubs)),
    types_(std::forward<TypeVec>(tys)) { }
  VarArrayDef(std::initializer_list<double> lbs,
              std::initializer_list<double> ubs,
              std::initializer_list<var::Type> tys) :
    lbs_((lbs)), ubs_((ubs)), types_((tys)) { }
  int size() const { assert(check()); return lbs_.size(); }
  const double* plb() const { return lbs_.data(); }
  const double* pub() const { return ubs_.data(); }
  const var::Type* ptype() const { return types_.data(); }
  void set_lb_ub_types(ArrayRef<double> lbs, ArrayRef<double> ubs,
                       ArrayRef<var::Type> types) {
    lbs_=lbs; ubs_=ubs; types_=types;
    assert(check());
  }
  bool check() const
  { return ubs_.size()==lbs_.size() && types_.size()==lbs_.size(); }
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

/// Backends handling custom flat constraints should derive from
class BasicFlatModelAPI {
public:
  /// Placeholder
  static const char* GetName()    { return "BasicFlatModelAPI"; }
  /// Placeholder
  static const char* GetLongName() { return nullptr; }

  /// Placeholder
  void InitOptions() { }

  /// Chance to prepare problem update
  void InitProblemModificationPhase() {  }
  /// Chance to end problem update
  void FinishProblemModificationPhase() {  }

  ////////////////// Some standard items /////////////////
  void SetLinearObjective(int , const LinearObjective& ) {
    throw MakeUnsupportedError("FlatBackend::SetLinearObjective()");
  }

  void SetQuadraticObjective(int , const QuadraticObjective& ) {
    throw MakeUnsupportedError("FlatBackend::SetQuadraticObjective()");
  }

  template <class Constraint>
  void AddConstraint(const Constraint& ) {
    throw std::logic_error(
          std::string("Not handling constraint '") +
          Constraint::GetName() +
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
};

/// ... then for a certain constraint it can be specified
#define ACCEPT_CONSTRAINT(ConstrType, level, con_grp) \
  mp::ConstraintAcceptanceLevel \
    AcceptanceLevel(const ConstrType*) const \
  { return (mp::ConstraintAcceptanceLevel)level; } \
  static constexpr int \
    GroupNumber(const ConstrType*) { return con_grp; }


} // namespace mp

#endif // FLAT_MODEL_API_BASE_H_
