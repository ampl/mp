#ifndef COND_CMP_H
#define COND_CMP_H

#include "mp/flat/redef/redef_base.h"
#include "mp/flat/constraints_std.h"

namespace mp {

/// Converts conditional <=, <, >, >= for MIP.
/// The comparison kind is given by con.GetConstraint().kind().
template <class ModelConverter, class AlgCon>
class Cond_LE_LT_GT_GE_Converter_MIP :
    public BasicFuncConstrCvt<
      Cond_LE_LT_GT_GE_Converter_MIP<ModelConverter, AlgCon>,
      ModelConverter> {
public:
  /// Base class
  using Base = BasicFuncConstrCvt<
    Cond_LE_LT_GT_GE_Converter_MIP<ModelConverter, AlgCon>,
    ModelConverter>;

  /// Constructor
  Cond_LE_LT_GT_GE_Converter_MIP(ModelConverter& mc) : Base(mc) { }

  /// Converted item type
  using ItemType = ConditionalConstraint< AlgCon >;

  /// Comparison kind
  static constexpr int kind_input = AlgCon::kind();

  /// Convert in positive context:
  /// resvar==1 --> <body> (cmp) <rhs>
  void ConvertCtxPos(const ItemType& cc, int ) {
    constexpr auto kind_output =               // output comparison is <= or >=
        ( kind_input>0 ) ? 1 : -1;             // even for < or >
    auto eps_output =
        ( kind_input==1 || kind_input==-1 ) ?  // non-strict comparison?
          0.0 :
          kind_output *
            GetMC().ComparisonEps(   // 1 for integer expressions
              GetMC().ComputeBoundsAndType(cc.GetArguments()).
              get_result_type() );
    ConvertCondIneq<kind_output>(cc, 1, eps_output);
  }

  /// Convert in negative context:
  /// resvar==0 --> (body) (!cmp) <rhs>
  void ConvertCtxNeg(const ItemType& cc, int ) {
    constexpr auto kind_output =               // output comparison is <= or >=
        ( kind_input>0 ) ? -1 : 1;             // even for < or >
    auto eps_output =
        ( kind_input==1 || kind_input==-1 ) ?  // non-strict comparison?
          kind_output *
            GetMC().ComparisonEps(   // 1 for integer expressions
              GetMC().ComputeBoundsAndType(cc.GetArguments()).
              get_result_type() ) :
          0.0;
    ConvertCondIneq<kind_output>(cc, 0, eps_output);
  }


protected:
  /// Convert conditional inequality to indicator
  /// res==value ==> Body (Cmp<kind>) Rhs+eps
  template <int kind>
  void ConvertCondIneq(const ItemType& cc, int value, double eps) {
    const auto& con = cc.GetConstraint();
    const auto res = cc.GetResultVar();
    using AlgConOutput = AlgebraicConstraint<
      typename AlgCon::BodyType, AlgConRhs<kind> >;
    if (con.empty()) {                    // empty body
      if (kind*(con.rhs() + eps) > 0.0)   // TODO option to switch off
        GetMC().NarrowVarBounds(res, !value, !value);    // fix result
    } else {
      if (GetMC().is_fixed(res)) {
        if (value==GetMC().fixed_value(res)) { // rsult already fixed,
          GetMC().AddConstraint(               // add static constraint
                AlgConOutput{ con.GetBody(), con.rhs()+eps } );
        }                                      // otherwise, forget
      } else {
        GetMC().AddConstraint( IndicatorConstraint< AlgConOutput >(
                                res, value,
                                { con.GetBody(), con.rhs()+eps }));
      }
    }
  }

  /// Reuse the stored ModelConverter
  using Base::GetMC;
};


/// Specialize for linear constraint
template <class MC>
using CondLinLEConverter_MIP = Cond_LE_LT_GT_GE_Converter_MIP<MC, LinConLE>;

/// Specialize for quadratic constraint
template <class MC>
using CondQuadLEConverter_MIP = Cond_LE_LT_GT_GE_Converter_MIP<MC, QuadConLE>;

/// Specialize for linear constraint
template <class MC>
using CondLinGEConverter_MIP = Cond_LE_LT_GT_GE_Converter_MIP<MC, LinConGE>;

/// Specialize for quadratic constraint
template <class MC>
using CondQuadGEConverter_MIP = Cond_LE_LT_GT_GE_Converter_MIP<MC, QuadConGE>;



/// Specialize for linear constraint
template <class MC>
using CondLinLTConverter_MIP = Cond_LE_LT_GT_GE_Converter_MIP<MC, LinConLT>;

/// Specialize for quadratic constraint
template <class MC>
using CondQuadLTConverter_MIP = Cond_LE_LT_GT_GE_Converter_MIP<MC, QuadConLT>;


/// Specialize for linear constraint
template <class MC>
using CondLinGTConverter_MIP = Cond_LE_LT_GT_GE_Converter_MIP<MC, LinConGT>;

/// Specialize for quadratic constraint
template <class MC>
using CondQuadGTConverter_MIP = Cond_LE_LT_GT_GE_Converter_MIP<MC, QuadConGT>;

} // namespace mp

#endif // COND_CMP_H
