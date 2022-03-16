#ifndef NUMBEROF_VAR_H
#define NUMBEROF_VAR_H

#include "mp/flat/redef/redef_base.h"
#include "mp/flat/constraints_std.h"

namespace mp {

/// Converts NOV for MIP
template <class ModelConverter>
class NumberofVarConverter_MIP :
    public BasicFuncConstrCvt<
      NumberofVarConverter_MIP<ModelConverter>, ModelConverter> {
public:
  /// Base class
  using Base = BasicFuncConstrCvt<
    NumberofVarConverter_MIP<ModelConverter>, ModelConverter>;
  /// Constructor
  NumberofVarConverter_MIP(ModelConverter& mc) : Base(mc) { }
  /// Converted item type
  using ItemType = NumberofVarConstraint;

  /// Convert in any context.
  /// Very basic, could be improved
  void Convert(const ItemType& novc, int ) {
    const auto& args = novc.GetArguments();
    std::vector<double> coefs(args.size(), 1.0);
    coefs.front() = -1.0;
    std::vector<int> flags(args.size(), novc.GetResultVar());
    for (size_t ivar = 1; ivar < args.size(); ++ivar) {
      flags[ivar] = GetMC().AssignResultVar2Args(   // flag = (args[i]==args[0])
            EQ0Constraint(
                 { { {1.0, -1.0}, {args[ivar], args[0]} }, 0.0 } ) );
    }
    GetMC().AddConstraint( LinConEQ( {coefs, flags}, {0.0} ) );
  }

  /// Reuse the stored ModelConverter
  using Base::GetMC;
};

} // namespace mp

#endif // NUMBEROF_VAR_H
