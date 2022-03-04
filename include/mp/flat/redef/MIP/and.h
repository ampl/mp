#ifndef AND_H
#define AND_H

#include "mp/flat/redef/redef_base.h"
#include "mp/flat/constraints_std.h"

namespace mp {

/// Converts And/Forall for MIP
template <class ModelConverter>
class AndConverter_MIP :
    public BasicFuncConstrCvt<
      AndConverter_MIP<ModelConverter>, ModelConverter> {
public:
  /// Base class
  using Base = BasicFuncConstrCvt<
    AndConverter_MIP<ModelConverter>, ModelConverter>;
  /// Constructor
  AndConverter_MIP(ModelConverter& mc) : Base(mc) { }
  /// Converted item type
  using ItemType = AndConstraint;

  /// Convert in positive context
  void ConvertCtxPos(const ItemType& conj, int ) {
    std::array<double, 2> coefs{-1.0, 1.0};
    std::array<int, 2> vars{-1, conj.GetResultVar()};
    for (auto arg: conj.GetArguments()) {   // res <= arg[i] in CTX+
      vars[0] = arg;
      GetMC().AddConstraint(
                     LinConLE({coefs, vars}, {0.0} ));
    }
  }

  /// Convert in negative context
  void ConvertCtxNeg(const ItemType& conj, int ) {
    const auto& args = conj.GetArguments();
    auto flags = args;
    flags.push_back(conj.GetResultVar());
    std::vector<double> ones(args.size(), 1.0); // res+n-1 >= sum(args) in CTX-
    ones.push_back(-1.0);
    GetMC().AddConstraint(
                   LinConLE({ones, flags},
                               {(double)args.size()-1} ));
  }

  /// Reuse the stored ModelConverter
  using Base::GetMC;
};

} // namespace mp

#endif // AND_H
