#ifndef NOT_H
#define NOT_H

#include "mp/flat/redef/redef_base.h"
#include "mp/flat/constraints_std.h"

namespace mp {

/// Converts NOT for MIP
template <class ModelConverter>
class NotConverter_MIP :
    public BasicFuncConstrCvt<
      NotConverter_MIP<ModelConverter>, ModelConverter> {
public:
  /// Base class
  using Base = BasicFuncConstrCvt<
    NotConverter_MIP<ModelConverter>, ModelConverter>;
  /// Constructor
  NotConverter_MIP(ModelConverter& mc) : Base(mc) { }
  /// Converted item type
  using ItemType = NotConstraint;

  /// Convert in both contexts (full reification)
  void Convert(const ItemType& nc, int ) {
    GetMC().AddConstraint(LinearFunctionalConstraint(
      nc.GetResultVar(), {{{-1.0}, {nc.GetArguments()[0]}}, 1.0}));
  }

  /// Reuse the stored ModelConverter
  using Base::GetMC;
};

} // namespace mp

#endif // NOT_H
