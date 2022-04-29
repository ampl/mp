#ifndef COND_LT_H
#define COND_LT_H

#include "mp/flat/redef/MIP/cond_le.h"

namespace mp {

/// Converts conditional strict < for MIP
template <class ModelConverter, class AlgConBody>
class CondLTConverter_MIP :
    public BasicFuncConstrCvt<
      CondLTConverter_MIP<ModelConverter, AlgConBody>,
      ModelConverter> {
public:
  /// Base class
  using Base = BasicFuncConstrCvt<
    CondLTConverter_MIP<ModelConverter, AlgConBody>,
    ModelConverter>;

  /// Constructor
  CondLTConverter_MIP(ModelConverter& mc) : Base(mc) { }

  /// Converted item type
  using ItemType = ConditionalConstraint<
      AlgebraicConstraint< AlgConBody, AlgConRhs<-2> > >;

  /// Convert in positive context
  void ConvertCtxPos(const ItemType& le0c, int ) {
    auto bNt = GetMC().ComputeBoundsAndType(le0c.GetArguments());
    double cmpEps = GetMC().ComparisonEps( bNt.get_result_type() );
    ConvertLELT_MIP_CtxPos(GetMC(), le0c, cmpEps);
  }

  /// Convert in negative context.
  /// resvar==0 --> (body) >= d
  void ConvertCtxNeg(const ItemType& le0c, int ) {
    ConvertLELT_MIP_CtxNeg(GetMC(), le0c, 0.0);
  }

  /// Reuse the stored ModelConverter
  using Base::GetMC;
};


/// Specialize for linear constraint
template <class MC>
using CondLinLTConverter_MIP = CondLTConverter_MIP<MC, LinTerms>;

/// Specialize for quadratic constraint
template <class MC>
using CondQuadLTConverter_MIP = CondLTConverter_MIP<MC, QuadAndLinTerms>;

} // namespace mp

#endif // COND_LT_H
