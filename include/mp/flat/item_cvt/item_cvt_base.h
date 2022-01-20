#ifndef BASE_H
#define BASE_H

#include "mp/presolve.h"

namespace mp {

/// A base class for specific Converters, such as
/// individual constraint converters.
/// Template parameter ModelConverter should have a GetPresolver() method.
template <class ModelConverter>
class BasicItemConverter {
public:
  /// Constructor
  BasicItemConverter(ModelConverter& mc) : mdl_cvt_(mc) { }

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

/// In the ModelConverter: to use a specific item_cvt_type<>
/// Assumes Impl is the final ModelConverter type
#define INSTALL_ITEM_CONVERTER(item_cvt_type) \
  item_cvt_type<Impl> item_cvt__ ## item_cvt_type ## _ \
   { *static_cast<Impl*>(this) }; \
  void Convert(const typename \
      item_cvt_type<Impl>::ItemType& con, int i) { \
    item_cvt__ ## item_cvt_type ## _ . Convert(con, i); \
  }


} // namespace mp

#endif // BASE_H
