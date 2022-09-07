#ifndef SOS2_H
#define SOS2_H

#include <cmath>

#include "mp/flat/redef/redef_base.h"
#include "mp/flat/constr_std.h"
#include "mp/common.h"

namespace mp {

/// Converts SOS2 for MIP
template <class ModelConverter>
class SOS2Converter_MIP :
    public BasicFuncConstrCvt<
      SOS2Converter_MIP<ModelConverter>, ModelConverter> {
public:
  /// Base class
  using Base = BasicFuncConstrCvt<
    SOS2Converter_MIP<ModelConverter>, ModelConverter>;
  /// Constructor
  SOS2Converter_MIP(ModelConverter& mc) : Base(mc) { }
  /// Converted item type
  using ItemType = SOS2Constraint;

  /// Convert in any context
  void Convert(const ItemType& cc, int ) {
    if (SOSExtraInfo::Bounds{1.0, 1.0} == cc.get_sum_of_vars_range())
      ConvertSOS2FromPL(cc);
    else
      ConvertGeneralSOS2(cc);
  }


protected:
  void ConvertGeneralSOS2(const ItemType& ) {
    MP_RAISE("MIP conversion of general SOS2 (not used to "
             "encode a pl): not implemented");
  }
  void ConvertSOS2FromPL(const ItemType& cc) {
    const int d = int(cc.get_vars().size()-1);
    const int r = int(std::ceil(std::log2(d)));
    auto lambda = cc.get_vars();
    lambda.push_back(-1);             // reserve for 1 extra var
    auto y = GetMC().AddVars_returnIds(r, 0.0, 1e100, var::INTEGER);
    for (int k=1; k<=r; ++k) {
      lambda.back() = y[k-1];
      auto Ckr_0_d = GetMC().GetZZIExtendedColumn(r, k, 0, d);
      Ckr_0_d.push_back(-1.0);
      GetMC().AddConstraint(
            LinConLE({std::move(Ckr_0_d), lambda}, {0.0}) );
      auto Ckr_1_dp1 = GetMC().GetZZIExtendedColumn(r, k, 1, d+1);
      Ckr_1_dp1.push_back(-1.0);
      GetMC().AddConstraint(
            LinConGE({std::move(Ckr_1_dp1), lambda}, {0.0}) );
    }
  }

  /// Reuse the stored ModelConverter
  using Base::GetMC;
};

} // namespace mp

#endif // SOS2_H
