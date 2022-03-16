#ifndef NUMBEROF_CONST_H
#define NUMBEROF_CONST_H

#include "mp/flat/redef/redef_base.h"
#include "mp/flat/constraints_std.h"

namespace mp {

/// Converts NOC for MIP
template <class ModelConverter>
class NumberofConstConverter_MIP :
    public BasicFuncConstrCvt<
      NumberofConstConverter_MIP<ModelConverter>, ModelConverter> {
public:
  /// Base class
  using Base = BasicFuncConstrCvt<
    NumberofConstConverter_MIP<ModelConverter>, ModelConverter>;
  /// Constructor
  NumberofConstConverter_MIP(ModelConverter& mc) : Base(mc) { }
  /// Converted item type
  using ItemType = NumberofConstConstraint;

  /// Convert in any context
  void Convert(const ItemType& nocc, int ) {
    const auto& args = nocc.GetArguments();
    const double k = nocc.GetParameters()[0];
    std::vector<double> coefs(args.size()+1, 1.0);
    std::vector<int> flags(args.size()+1, nocc.GetResultVar());
    for (size_t ivar = 0; ivar < args.size(); ++ivar) {
      flags[ivar] = GetMC().AssignResultVar2Args(  // flag = (args[i]==k)
            EQ0Constraint( { {{1.0}, {args[ivar]}}, -k } ) );
    }
    coefs.back() = -1.0;
    GetMC().AddConstraint( LinConEQ( {coefs, flags}, {0.0} ) );
  }

  /// Reuse the stored ModelConverter
  using Base::GetMC;
};

} // namespace mp

#endif // NUMBEROF_CONST_H
