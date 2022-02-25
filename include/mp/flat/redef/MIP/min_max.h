#ifndef MIN_MAX_H
#define MIN_MAX_H

/*
 * Convert Min and Max for MIP
 */

#include "mp/flat/redef/redef_base.h"
#include "mp/flat/constraints_std.h"

namespace mp {

/// Type parameterization
template <int sense>
struct MaxOrMinConstraint { };
/// MinConstraint
template <>
struct MaxOrMinConstraint<-1> { using ItemType = MinimumConstraint; };
/// MinConstraint
template <>
struct MaxOrMinConstraint< 1> { using ItemType = MaximumConstraint; };


/// Convert Min and Max for MIP
template <class ModelConverter, int sense>
class MinOrMaxConverter_MIP :
    public BasicFuncConstrCvt<
      MinOrMaxConverter_MIP<ModelConverter, sense>, ModelConverter> {
public:
  /// Base class
  using Base = BasicFuncConstrCvt<
    MinOrMaxConverter_MIP<ModelConverter, sense>, ModelConverter>;
  /// Constructor
  MinOrMaxConverter_MIP(ModelConverter& mc) : Base(mc) { }
  /// ItemType
  using ItemType = typename MaxOrMinConstraint<sense>::ItemType;

  /// Convert in positive context
  template <class MinOrMaxConstraint>
  void ConvertCtxPos(const MinOrMaxConstraint& mc, int ) {
    if (sense>0)
      ConvertNonConvexPart(mc);
    else
      ConvertConvexPart(mc);
  }

  /// Convert in negative context
  template <class MinOrMaxConstraint>
  void ConvertCtxNeg(const MinOrMaxConstraint& mc, int ) {
    if (sense>0)
      ConvertConvexPart(mc);
    else
      ConvertNonConvexPart(mc);
  }

  /// Convert the convex part
  template <class MinOrMaxConstraint>
  void ConvertConvexPart(const MinOrMaxConstraint& mc) {
    const auto& args = mc.GetArguments();
    const std::size_t nargs = args.size();
    const auto resvar = mc.GetResultVar();
    for (size_t i=0; i<nargs; ++i) {
      GetMC().AddConstraint(LinConLE(
                         { {1.0*sense, -1.0*sense},
                           {args[i], resvar} }, 0.0));
    }
  }

  /// Convert the non-convex part
  template <class MinOrMaxConstraint>
  void ConvertNonConvexPart(const MinOrMaxConstraint& mc) {
    const auto& args = mc.GetArguments();
    const std::size_t nargs = args.size();
    const auto flags =                            // binary flags
        GetMC().AddVars_returnIds(nargs, 0.0, 1.0, var::Type::INTEGER);
    GetMC().AddConstraint(LinConGE(               // sum of the flags >= 1
               { std::vector<double>(nargs, 1.0),
                                   flags }, 1.0));
    const auto resvar = mc.GetResultVar();
    for (size_t i=0; i<nargs; ++i) {
      GetMC().AddConstraint(
                     IndicatorConstraintLinLE{flags[i], 1,
                         { { {1.0*sense, -1.0*sense},
                             {resvar, args[i]} }, 0.0 }});
    }
  }

  /// Reuse the stored ModelConverter
  using Base::GetMC;
};

/// Typedef of MinConverter
template <class ModelConverter>
using MinConverter_MIP = MinOrMaxConverter_MIP<ModelConverter, -1>;
/// Typedef of MaxConverter
template <class ModelConverter>
using MaxConverter_MIP = MinOrMaxConverter_MIP<ModelConverter,  1>;

} // namespace mp

#endif // MIN_MAX_H
