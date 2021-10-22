#ifndef CONSTRAINT_ADDER_H
#define CONSTRAINT_ADDER_H

#include "mp/convert/basic_constr.h"

namespace mp {

/// Level of acceptance of a constraint by a backend
enum ConstraintAcceptanceLevel {
  NotAccepted,
  AcceptedButNotRecommended,
  Recommended
};

/// Backends handling custom constraints should derive from
class BasicConstraintAdder {
public:
  template <class Constraint>
  void AddConstraint(const Constraint& ) {
    throw std::logic_error(
          std::string("Not handling constraint ") + Constraint::GetConstraintName());
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

#endif // CONSTRAINT_ADDER_H
