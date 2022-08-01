#ifndef COUNT_H
#define COUNT_H

#include "mp/flat/redef/redef_base.h"
#include "mp/flat/constr_std.h"

namespace mp {

/// Converts count() for MIP
template <class ModelConverter>
class CountConverter_MIP :
    public BasicFuncConstrCvt<
      CountConverter_MIP<ModelConverter>, ModelConverter> {
public:
  /// Base class
  using Base = BasicFuncConstrCvt<
    CountConverter_MIP<ModelConverter>, ModelConverter>;
  /// Constructor
  CountConverter_MIP(ModelConverter& mc) : Base(mc) { }
  /// Converted item type
  using ItemType = CountConstraint;

  /// Convert in any context
  void Convert(const ItemType& cc, int ) {
    const auto& args = cc.GetArguments();
    std::vector<double> coefs(args.size()+1, 1.0);
    coefs.back() = -1.0;
    std::vector<int> flags(args.size()+1, cc.GetResultVar());
    for (size_t ivar = 0; ivar < args.size(); ++ivar) {
      flags[ivar] = args[ivar];
      /// Force booleanize: reify if we have a "!=0" expression
      if ( !GetMC().is_binary_var(args[ivar]) ) {
        auto feq0 = GetMC().AssignResultVar2Args(   // feq0 = (args[i]==0)
            CondLinConEQ( { {{1.0}, {args[ivar]}}, 0.0 } ) );
        flags[ivar] = GetMC().AssignResultVar2Args(   // flag = (args[i]!=0)
            NotConstraint( {feq0} ));
      }
    }
    GetMC().AddConstraint( LinConEQ( {coefs, flags}, 0.0 ) );
  }

  /// Reuse the stored ModelConverter
  using Base::GetMC;
};

} // namespace mp

#endif // COUNT_H
