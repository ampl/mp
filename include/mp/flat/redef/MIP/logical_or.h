#ifndef LOGICAL_OR_H
#define LOGICAL_OR_H

#include "mp/flat/redef/redef_base.h"
#include "mp/flat/constr_std.h"

namespace mp {

/// Converts Or/Exists for MIP
template <class ModelConverter>
class OrConverter_MIP :
    public BasicFuncConstrCvt<
      OrConverter_MIP<ModelConverter>, ModelConverter> {
public:
  /// Base class
  using Base = BasicFuncConstrCvt<
    OrConverter_MIP<ModelConverter>, ModelConverter>;
  /// Constructor
  OrConverter_MIP(ModelConverter& mc) : Base(mc) { }
  /// Converted item type
  using ItemType = OrConstraint;

  /// Convert in positive context
  void ConvertCtxPos(const ItemType& disj, int ) {
    const auto& args = disj.GetArguments();
    auto flags = args;
    flags.push_back(disj.GetResultVar());
    std::vector<double> ones(args.size(), 1.0);  // res <= sum(args) in CTX+
    ones.push_back(-1.0);
    GetMC().AddConstraint(
                   LinConGE({ones, flags}, {0.0} ));
  }

  /// Convert in negative context
  void ConvertCtxNeg(const ItemType& disj, int ) {
    std::array<double, 2> coefs{1.0, -1.0};
    std::array<int, 2> vars{-1, disj.GetResultVar()};
    for (auto arg: disj.GetArguments()) {        // res >= arg[i] in CTX-
      vars[0] = arg;
      GetMC().AddConstraint(
                     LinConLE({coefs, vars}, {0.0} ));
    }
  }

  /// Reuse the stored ModelConverter
  using Base::GetMC;
};

} // namespace mp

#endif // LOGICAL_OR_H
