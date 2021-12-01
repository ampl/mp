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
#ifndef FLAT_MODEL_API_BASIC_H_
#define FLAT_MODEL_API_BASIC_H_

#include <string>
#include <stdexcept>

#include "mp/arrayref.h"
#include "mp/common.h"
#include "mp/flat/basic_constr.h"

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
  NotAccepted,
  AcceptedButNotRecommended,
  Recommended
};

/// Backends handling custom flat constraints should derive from
class BasicFlatBackend {
public:
  template <class Constraint>
  void AddConstraint(const Constraint& ) {
    throw std::logic_error(
          std::string("Not handling constraint ") +
          Constraint::GetConstraintName() +
          ". Provide a handler or a converter method");
  }
  /// Derived backends have to tell C++ to use default handlers if they are needed
  /// when they overload AddConstraint(), due to C++ name hiding
#define USE_BASE_CONSTRAINT_HANDLERS(BaseBackend) \
  using BaseBackend::AddConstraint; \
  using BaseBackend::AcceptanceLevel;
  /// By default, we say constraint XYZ is not accepted but...
  static constexpr ConstraintAcceptanceLevel AcceptanceLevel(const BasicConstraint*) {
    return NotAccepted;
  }
};

/// ... then for a certain constraint it can be specified
#define ACCEPT_CONSTRAINT(ConstrType, level) \
  static constexpr mp::ConstraintAcceptanceLevel \
    AcceptanceLevel(const ConstrType*) { return level; }


} // namespace mp

#endif // FLAT_MODEL_API_BASIC_H_
