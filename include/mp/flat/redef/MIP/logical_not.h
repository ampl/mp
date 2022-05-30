#ifndef LOGICAL_NOT_H
#define LOGICAL_NOT_H

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
    /// Obtain negation variable via map
    int var_res_lin = GetMC().AssignResultVar2Args(
          LinearFunctionalConstraint(
            {{{-1.0}, {nc.GetArguments()[0]}}, 1.0}));
    GetMC().AddConstraint(LinConEQ{
                            { {-1.0, 1.0},
                              {nc.GetResultVar(), var_res_lin} },
                            {0.0}});
  }

  /// Reuse the stored ModelConverter
  using Base::GetMC;
};

} // namespace mp

#endif // LOGICAL_NOT_H
