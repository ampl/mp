#ifndef ABS_H
#define ABS_H

#include "mp/flat/redef/redef_base.h"
#include "mp/flat/constraints_std.h"

namespace mp {

/// Converts abs for MIP
template <class ModelConverter>
class AbsConverter_MIP :
    public BasicFuncConstrCvt<
      AbsConverter_MIP<ModelConverter>, ModelConverter> {
public:
  /// Base class
  using Base = BasicFuncConstrCvt<
    AbsConverter_MIP<ModelConverter>, ModelConverter>;
  /// Constructor
  AbsConverter_MIP(ModelConverter& mc) : Base(mc) { }
  /// Converted item type
  using ItemType = AbsConstraint;

  /// Convert in positive context
  void ConvertCtxPos(const ItemType& ac, int ) {
    const int arg = ac.GetArguments()[0];
    const int res = ac.GetResultVar();
    const int flag = GetMC().AddVar(0.0, 1.0, var::INTEGER);
    GetMC().AddConstraint(
          IndicatorConstraintLinLE(flag, 1, {{{1.0, 1.0}, {res, arg}}, 0.0}));
    GetMC().AddConstraint(
          IndicatorConstraintLinLE(flag, 0, {{{1.0, -1.0}, {res, arg}}, 0.0}));
  }

  /// Convert in negative context
  void ConvertCtxNeg(const ItemType& ac, int ) {
    const int arg = ac.GetArguments()[0];
    const int res = ac.GetResultVar();
    GetMC().AddConstraint(LinConGE({{1.0, 1.0}, {res, arg}}, {0.0}));
    GetMC().AddConstraint(LinConGE({{1.0, -1.0}, {res, arg}}, {0.0}));
  }

  /// Reuse the stored ModelConverter
  using Base::GetMC;
};

} // namespace mp

#endif // ABS_H
