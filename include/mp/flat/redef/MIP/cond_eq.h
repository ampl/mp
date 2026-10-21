#ifndef COND_EQ_H
#define COND_EQ_H

#include "mp/common.h"
#include "mp/flat/redef/redef_base.h"
#include "mp/flat/constr_std.h"
#include "mp/flat/redef/MIP/cond_literal.h"

namespace mp {

/// Conditional equality body == rhs with body type \a Body
template <class Body>
using CondEQ = ConditionalConstraint<
  AlgebraicConstraint< Body, AlgConRhs<0> > >;

/// Add (lit) ==> body == rhs for a non-empty body.
/// Shared by CondEQConverter_MIP (literal = result variable)
/// and by disjunction absorption.
template <class MC, class Body>
void AddCondEqAsIndicator(
    MC& mc, const CondEQ<Body>& eq0c, CondLiteral lit) {
  const auto& con = eq0c.GetConstraint();
  if (SurelyTriggers(mc, lit))
    mc.AddConstraint( con );                  // static constraint
  else if (CanTrigger(mc, lit))
    mc.AddConstraint( IndicatorConstraint<
                     AlgebraicConstraint< Body, AlgConRhs<0> > >(
                        lit.var_, lit.val_, con) );
}

/// If the body's bounds leave just one side of body != rhs
/// (body >= rhs+eps or body <= rhs-eps), add it statically.
/// @return whether added
template <class MC, class Body, class BndNType>
bool AddNeAsStaticIfOneSided(MC& mc, const CondEQ<Body>& eq0c,
                             const BndNType& bnt, double eps) {
  const auto& con = eq0c.GetConstraint();
  if (bnt.lb() > con.rhs() - eps) {          // body > rhs-eps
    mc.AddConstraint( AlgebraicConstraint< Body, AlgConRhs<1> >{
                        con.GetBody(), con.rhs() + eps } );
    return true;
  }
  if (bnt.ub() < con.rhs() + eps) {          // body < rhs+eps
    mc.AddConstraint( AlgebraicConstraint< Body, AlgConRhs<-1> >{
                        con.GetBody(), con.rhs() - eps } );
    return true;
  }
  return false;
}

/// Add one side of body != rhs as an indicator on a new binary
/// q: q==1 ==> body (cmp) rhs_side, where (cmp) is <= or >=
/// for \a kind -1 or 1. Add q to the disjunction \a disj.
template <int kind, class MC, class Body>
void AddNeSide(MC& mc, const CondEQ<Body>& eq0c,
               double rhs_side, LinTerms& disj) {
  const int q = (int)mc.AddVar(0.0, 1.0, var::INTEGER);
  disj.add_term(1.0, q);
  mc.AddConstraint( IndicatorConstraint<
                   AlgebraicConstraint< Body, AlgConRhs<kind> > >(
                      q, 1, { eq0c.GetConstraint().GetBody(), rhs_side }) );
}

/// Add (lit) ==> body <= rhs-eps || body >= rhs+eps
/// as  [lit does not trigger] + q1 + q2 >= 1
/// with indicators on q1, q2, for the sides
/// allowed by the body's bounds.
template <class MC, class Body, class BndNType>
void AddNeDisjunction(MC& mc, const CondEQ<Body>& eq0c,
                      CondLiteral lit, const BndNType& bnt, double eps) {
  const auto rhs = eq0c.GetConstraint().rhs();
  auto non_trigger = NonTriggerExpr(lit);
  LinTerms disj = non_trigger.GetBody();
  if (bnt.lb() <= rhs - eps)
    AddNeSide<-1>(mc, eq0c, rhs - eps, disj);
  if (bnt.ub() >= rhs + eps)
    AddNeSide<1>(mc, eq0c, rhs + eps, disj);
  if (disj.size() > non_trigger.GetBody().size())
    mc.AddConstraint( LinConGE(
                        std::move(disj), 1.0 - non_trigger.constant_term() ) );
  else
    ForbidTrigger(mc, lit);                  // body == rhs for sure
}

/// Add (lit) ==> body != rhs for a non-empty body,
/// directly in flat algebra and indicators.
/// Shared by CondEQConverter_MIP (literal: result var == 0)
/// and by disjunction absorption.
///
/// @note An alternative would be the logical constraint
///   res || body <= rhs-eps || body >= rhs+eps.
///   It is not used: it would propagate CTX+ into the result
///   variable being redefined, and the new conditional
///   inequalities could have been redefined in CTX- before,
///   losing the CTX+ redefinition (#248, eqVarConst01.mod,
///   int_ne_05_redef_ctx.mod).
template <class MC, class Body>
void AddCondNeAsIndicators(
    MC& mc, const CondEQ<Body>& eq0c, CondLiteral lit) {
  if (!CanTrigger(mc, lit))
    return;
  const auto bnt = mc.ComputeBoundsAndType(eq0c.GetConstraint().GetBody());
  const double eps = mc.ComparisonEps( bnt.get_result_type() );
  if (SurelyTriggers(mc, lit)
      && AddNeAsStaticIfOneSided(mc, eq0c, bnt, eps))
    return;
  AddNeDisjunction(mc, eq0c, lit, bnt, eps);
}


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

  /// Converted item type
  using ItemType = CondEQ<AlgConBody>;

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
    if (!FixResultIfEmpty(eq0c))
      AddCondEqAsIndicator(GetMC(), eq0c, {eq0c.GetResultVar(), 1});
  }

  /// Convert in negative context
  /// resvar==0 ==> body!=d
  void ConvertCtxNeg(const ItemType& eq0c, int ) {
    if (!FixResultIfEmpty(eq0c))
      AddCondNeAsIndicators(GetMC(), eq0c, {eq0c.GetResultVar(), 0});
  }


protected:
  /// For an empty body (0 == rhs), fix the result.
  /// @return whether the body is empty
  bool FixResultIfEmpty(const ItemType& eq0c) {
    const auto& con = eq0c.GetConstraint();
    if (!con.empty())
      return false;
    const bool valid = !con.rhs();
    GetMC().NarrowVarBounds(eq0c.GetResultVar(), valid, valid);
    return true;
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
