#ifndef COND_INEQ_H
#define COND_INEQ_H

#include "mp/flat/redef/redef_base.h"
#include "mp/flat/constr_std.h"
#include "mp/flat/redef/MIP/cond_literal.h"

namespace mp {

/// Strict comparison offset added to the rhs of the
/// output inequality (<= or >=) of an indicator implying
/// conditional inequality \a cc, or its negation if \a f_neg.
/// Strict outputs exclude equality by ComparisonEps():
/// 1 for integer bodies, cvt:mip:eps otherwise.
template <class MC, class AlgCon>
double CondIneqRhsOffset(
    MC& mc, const ConditionalConstraint<AlgCon>& cc, bool f_neg) {
  constexpr int kind_input = AlgCon::kind();
  const bool f_strict_output =                 // < or >, or negated <=, >=
      (kind_input==1 || kind_input==-1) == f_neg;
  if (!f_strict_output)
    return 0.0;
  const int kind_output = (kind_input>0) != f_neg ? 1 : -1;
  return kind_output * mc.ComparisonEps(
        mc.ComputeBoundsAndType(cc.GetArguments()).get_result_type() );
}

/// Add (lit) ==> body (cmp) rhs+offset, where (cmp) is <= or >=.
/// Handles an empty body and a fixed literal.
template <int kind, class MC, class AlgCon>
void AddIneqIndicator(MC& mc, const ConditionalConstraint<AlgCon>& cc,
                      CondLiteral lit, double offset) {
  using AlgConOutput = AlgebraicConstraint<
    typename AlgCon::BodyType, AlgConRhs<kind> >;
  const auto& con = cc.GetConstraint();
  AlgConOutput con_out { con.GetBody(), con.rhs()+offset };
  if (con.empty()) {                         // empty body: 0 (cmp) rhs
    if (kind*(con.rhs() + offset) > 0.0)     // violated
      ForbidTrigger(mc, lit);
  } else if (SurelyTriggers(mc, lit))
    mc.AddConstraint( std::move(con_out) );  // static constraint
  else if (CanTrigger(mc, lit))
    mc.AddConstraint( IndicatorConstraint< AlgConOutput >(
                        lit.var_, lit.val_, std::move(con_out) ));
}

/// Add the indicator (lit) ==> conditional inequality \a cc,
/// or (lit) ==> its negation if \a f_neg.
/// The output comparison is <= or >=, even for < or >.
/// Shared by Cond_LE_LT_GT_GE_Converter_MIP
/// (literal = result variable) and by disjunction absorption.
template <class MC, class AlgCon>
void AddCondIneqAsIndicator(
    MC& mc, const ConditionalConstraint<AlgCon>& cc,
    CondLiteral lit, bool f_neg) {
  constexpr int kind_pos = AlgCon::kind()>0 ? 1 : -1;
  const auto offset = CondIneqRhsOffset(mc, cc, f_neg);
  if (f_neg)
    AddIneqIndicator<-kind_pos>(mc, cc, lit, offset);
  else
    AddIneqIndicator<kind_pos>(mc, cc, lit, offset);
}


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

  /// Convert in positive context:
  /// resvar==1 --> (body) (cmp) (rhs)
  void ConvertCtxPos(const ItemType& cc, int ) {
    AddCondIneqAsIndicator(GetMC(), cc, {cc.GetResultVar(), 1}, false);
  }

  /// Convert in negative context:
  /// resvar==0 --> (body) (!cmp) (rhs)
  void ConvertCtxNeg(const ItemType& cc, int ) {
    AddCondIneqAsIndicator(GetMC(), cc, {cc.GetResultVar(), 0}, true);
  }


protected:
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

#endif // COND_INEQ_H
