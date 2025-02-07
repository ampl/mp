#ifndef COND_EQ_H
#define COND_EQ_H

#include "mp/common.h"
#include "mp/flat/redef/redef_base.h"
#include "mp/flat/constr_std.h"

namespace mp {

/// Converts conditional equality for MIP
template <class ModelConverter, class AlgConBody>
class CondEQConverter_MIP :
    public BasicFuncConstrCvt<
      CondEQConverter_MIP<ModelConverter, AlgConBody>,
      ModelConverter> {
public:
  /// Base class
  using Base = BasicFuncConstrCvt<
    CondEQConverter_MIP<ModelConverter, AlgConBody>,
    ModelConverter>;

  /// Constructor
  CondEQConverter_MIP(ModelConverter& mc) : Base(mc) { }

  /// Underlying algebraic constraint,
  /// parameterized by comparison type
  template <int kind>
  using AlgCon =
    AlgebraicConstraint< AlgConBody, AlgConRhs<kind> >;

  /// Converted item type
  using ItemType = ConditionalConstraint< AlgCon<0> >;

  /// Reuse
  using Base::IfDelayConversion;

  /// Generic check whether the constraint
  /// needs to be skipped from conversion,
  /// despite being not accepted by ModelAPI.
  /// Reimplementing for the linear case.
  bool IfDelayConversion(const CondLinConEQ& eq0c, int ) {
    const auto& args = eq0c.GetArguments();
    return 1==args.size() &&               // 1 variable
        GetMC().IfMightUseEqualityEncodingForVar(
            args.var(0));
  } // If yes, might use unary encoding whose flags are,
    // in the fixed case, fixed by PropagateResult()

  /// Convert in positive context
  /// resvar==1 => body==d
  void ConvertCtxPos(const ItemType& eq0c, int ) {
    const auto& con = eq0c.GetConstraint();
    const auto res = eq0c.GetResultVar();
    if (con.empty()) {
      bool valid = !con.rhs();
      GetMC().NarrowVarBounds(res, valid, valid);
    } else {
      if (GetMC().is_fixed(res)) {
        if (GetMC().fixed_value(res)) {  // fixed to 1
          GetMC().AddConstraint( con );
        } // else, skip
      } else
        GetMC().AddConstraint(IndicatorConstraint< AlgCon<0> >(
                              res, 1, con));
    }
  }

  /// Convert in negative context
  /// resvar==0 ==> body!=d
  void ConvertCtxNeg(const ItemType& eq0c, int ) {
    const auto& con = eq0c.GetConstraint();
    const auto res = eq0c.GetResultVar();
    if (con.empty()) {
      bool valid = !con.rhs();
      GetMC().NarrowVarBounds(res, valid, valid);
    } else if ( !GetMC().is_fixed(res) ||   // not fixed, or
                !GetMC().fixed_value(res) ) // fixed to 0
    {
      auto bNt = GetMC().ComputeBoundsAndType(con.GetBody());
      double cmpEps = GetMC().ComparisonEps( bNt.get_result_type() );
#define COND_LIN_EQ__USE_FLAT_ALGEBRA
#ifndef COND_LIN_EQ__USE_FLAT_ALGEBRA
      // Not doing this.
      // Reason: the below OR propagates CTX+ into res,
      // which is the expression we are redefining,
      // but this is not necessary. See #248.
      // Happens on eqVarConst01.mod.
      /// res3 ==> (res || body <= rhs-eps || body >= rhs+eps)
      auto res3 = GetMC().AssignResultVar2Args(
          OrConstraint{ {
            res,
            GetMC().AssignResultVar2Args(   // <=
                ConditionalConstraint< AlgCon<-1> >
                { { con.GetBody(), con.rhs() - cmpEps } }),
            GetMC().AssignResultVar2Args(   // >=
                ConditionalConstraint< AlgCon<1> >
                { { con.GetBody(), con.rhs() + cmpEps } })
          } });
      GetMC().FixAsTrue(res3);
#else  // COND_LIN_EQ__USE_FLAT_ALGEBRA
      // Old way: straight to algebra and indicators.
      // This could duplicate the below indicators,
      // but simpler as of now:
      // the indicators redefine conditional inequalities
      // which, in the above variant, could have been redefined
      // in CTX- before, see int_ne_05_redef_ctx.mod,
      // thus losing the CTX+ redefinition, see #248.
      if (GetMC().is_fixed(res)  // This prepro would be unnec with logicals
          && !GetMC().fixed_value(res)) {  // preprocess: fixed to 0
        if (bNt.lb() > con.rhs() - cmpEps) { // body > value-eps
          GetMC().AddConstraint(
              AlgCon<1>{ con.GetBody(), con.rhs() + cmpEps } );
          return;
        }
        if (bNt.ub() < con.rhs() + cmpEps) { // body < value+eps
          GetMC().AddConstraint(
              AlgCon<-1>{ con.GetBody(), con.rhs() - cmpEps } );
          return;
        }
      }
      auto newvars = GetMC().AddVars_returnIds(2, 0.0, 1.0, var::INTEGER);
      newvars.push_back( res );
      GetMC().AddConstraint( LinConGE(   // b1+b2+resvar >= 1
          {{1.0, 1.0, 1.0}, newvars},
          1.0 ) );
      {
        GetMC().AddConstraint(IndicatorConstraint< AlgCon<-1> >(
            newvars[0], 1,
            { con.GetBody(),
             con.rhs() - cmpEps }));
      }
      GetMC().AddConstraint(IndicatorConstraint< AlgCon<1> >(
          newvars[1], 1,
          { con.GetBody(),
           con.rhs() + cmpEps }));
#endif  // USE_FLAT_ALGEBRA
    } // else, skip
  }

  /// Reuse the stored ModelConverter
  using Base::GetMC;
};


/// Specialize for linear constraint
template <class MC>
using CondLinEQConverter_MIP = CondEQConverter_MIP<MC, LinTerms>;

/// Specialize for quadratic constraint
template <class MC>
using CondQuadEQConverter_MIP = CondEQConverter_MIP<MC, QuadAndLinTerms>;

} // namespace mp

#endif // COND_EQ_H
