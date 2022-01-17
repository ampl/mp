#ifndef BASE_H
#define BASE_H

#include "mp/presolve.h"

namespace mp {

/// A base class for specific Converters, such as
/// individual constraint converters.
/// Template parameter ModelConverter should have GetPresolver() method.
template <class ModelConverter>
class BasicConverter {
public:
  /// Constructor
  BasicConverter(ModelConverter& mc) : mdl_cvt_(mc) { }

protected:
  /// Access ModelConverter
  ModelConverter& GetMC() { return mdl_cvt_; }
  /// Access Presolver
  pre::Presolver& GetPresolver() { return GetMC().GetPresolver(); }
  /// Access specific constraint's ValueNode
  template <class Constraint>
  pre::ValueNode& GetValueNode(Constraint* pc)
  { return GetMC().GetConstraintKeeper(pc).GetValueNode(); }

private:
  ModelConverter& mdl_cvt_;
};

/// To be used by descendants of BasiccOnverter
#define GET_CONSTRAINT_VALUE_NODE(con_type) \
  this->GetMC().GetValueNode((con_type*)nullptr)

} // namespace mp

#endif // BASE_H
