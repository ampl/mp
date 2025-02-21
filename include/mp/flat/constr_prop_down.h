#ifndef CONSTR_PROP_DOWN_H
#define CONSTR_PROP_DOWN_H

/*
 * Propagate flat constraints from result (result bounds & context)
 * "down", i.e., to the arguments.
 *
 * In the below methods, for functional constraints,
 * lb and ub are bounds on the result.
 *
 * For static constraints, they may be the constraint range.
 */

#include "mp/common.h"
#include "mp/flat/constr_std.h"
#include "mp/flat/nl_expr/constr_nl.h"

namespace mp {

/// A mix-in base class
/// providing "down propagators" of flat constraints, i.e.,
/// from result bounds & context to arguments.
template <class Impl>
class ConstraintPropagatorsDown {
public:

  /// By default, add mixed context for argument variables
  template <class Constraint>
  void PropagateResult(Constraint& con, double lb, double ub, Context ctx) {
    internal::Unused(&con, lb, ub, ctx);
    con.AddContext(ctx);           // merge context
    PropagateResult2Args(con.GetArguments(),     // we don't know the constraint
                         MPD( MinusInfty() ), MPD( Infty() ), Context::CTX_MIX);
  }

  void PropagateResult(LinearFunctionalConstraint& con, double lb, double ub,
                       Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    PropagateResult2LinTerms(con.GetAffineExpr(),   // @todo better in special cases
                             MPD( MinusInfty() ), MPD( Infty() ), +ctx);
  }

  void PropagateResult(QuadraticFunctionalConstraint& con, double lb, double ub,
                       Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    const auto& args = con.GetArguments();
    PropagateResult2LinTerms(args.GetLinTerms(),
                             MPD( MinusInfty() ), MPD( Infty() ), +ctx);
    PropagateResult2QuadTerms(args.GetQPTerms(),
                              MPD( MinusInfty() ), MPD( Infty() ), +ctx);
  }

  /// Propagate a root algebraic range constraint
  template <class Body, class RangeOrRhs>
  void PropagateResult(const AlgebraicConstraint<Body, RangeOrRhs>& con) {
    PropagateResult(con, Context::CTX_ROOT);
  }

  /// Conditional algebraic con
  template <class Body, class RngOrRhs>
  void PropagateResult(const AlgebraicConstraint< Body, RngOrRhs >& con,
                       Context ctx) {
    /// Distinguish bounds' finiteness for context
    auto ctx_new = con.lb()<=MPCD( PracticallyMinusInf() )
                   ? -ctx : (con.ub()>=MPCD( PracticallyInf() )
                          ? +ctx : Context::CTX_MIX);
    auto LB = MPCD( MinusInfty() );
    auto UB = MPCD( Infty() );
    if (ctx.IsRoot()) {             // not for conditional comparisons
      LB = con.lb();
      UB = con.ub();
    }
    PropagateResult2Args(con.GetBody(), LB, UB, ctx_new);
  }

  template <class Body, int sens>
  void PropagateResult(IndicatorConstraint<
                         AlgebraicConstraint< Body, AlgConRhs<sens> > >& con,
                       double lb, double ub, Context ctx) {
    internal::Unused(lb, ub);
    MPD( PropagateResultOfInitExpr(con.get_binary_var(),
                              MPD( MinusInfty() ), MPD( Infty() ),
                              1==con.get_binary_value() ?  // b==1 means b in CTX_NEG
                                Context::CTX_NEG : Context::CTX_POS) );
    PropagateResult(con.get_constraint(), Context::CTX_POS);      // implication
  }

  template <int type>
  void PropagateResult(SOS_1or2_Constraint<type>& con, double lb, double ub, Context ) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    MPD( PropagateResult2Vars(con.get_vars(),
                        MPD( MinusInfty() ), MPD( Infty() ), Context::CTX_MIX) );
  }

  /// Propagate root complementarity constraint
  template <class ExprBody>
  void PropagateResult(ComplementarityConstraint<ExprBody>& con) {
    MPD( PropagateResult(con,
                         MPD(MinusInfty()), MPD(Infty()), Context::CTX_MIX) );
  }

  /// Not used?
  template <class ExprBody>
  void PropagateResult(ComplementarityConstraint<ExprBody>& con,
                       double lb, double ub, Context ctx) {
    internal::Unused(ctx);
    MPD( PropagateResult2Args(con.GetExpression().GetBody(),
                         lb, ub, Context::CTX_MIX) );
    MPD( PropagateResultOfInitExpr(con.GetVariable(), lb, ub, Context::CTX_MIX) );
  }

  void PropagateResult(NotConstraint& con, double lb, double ub, Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    if (lb==ub) {                         // result fixed
      if (!lb && ctx.HasNegative()) {       // result==0 && ctx-
        MPD( PropagateResultOfInitExpr(con.GetArguments()[0], 1.0, 1.0, -ctx) );
        MPD( DecrementVarUsage(con.GetResultVar()) );
        return;
      } else if (lb && ctx.HasPositive()) { // result==1 && ctx+
        MPD( PropagateResultOfInitExpr(con.GetArguments()[0], 0.0, 0.0, -ctx) );
        MPD( DecrementVarUsage(con.GetResultVar()) );
        return;
      }
    }
    MPD( PropagateResultOfInitExpr(con.GetArguments()[0], 0.0, 1.0, -ctx) );
  }

  void PropagateResult(AndConstraint& con, double lb, double ub, Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    if (lb>0.5 && ctx.HasPositive()) {                  // Remove, arguments are fixed
      MPD( PropagateResult2Vars(con.GetArguments(), lb, 1.0, +ctx) );
      MPD( DecrementVarUsage(con.GetResultVar()) );    // Or, remove completely?
    } else
      MPD( PropagateResult2Vars(con.GetArguments(), 0.0, 1.0, +ctx) );  // in any ctx??
  }

  void PropagateResult(OrConstraint& con, double lb, double ub, Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    if (ub<0.5 && ctx.HasNegative()) {                 // Remove, arguments are fixed
      MPD( PropagateResult2Vars(con.GetArguments(), 0.0, ub, +ctx) );
      MPD( DecrementVarUsage(con.GetResultVar()) );
    } else
      MPD( PropagateResult2Vars(con.GetArguments(), 0.0, 1.0, +ctx) );
  }

  void PropagateResult(IfThenConstraint& con, double lb, double ub, Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    auto& args = con.GetArguments();
    MPD( PropagateIfThenResultIntoCondition(args, ctx) );
    MPD( PropagateResultOfInitExpr(args[1], MPD( MinusInfty() ), MPD( Infty() ), +ctx) );
    MPD( PropagateResultOfInitExpr(args[2], MPD( MinusInfty() ), MPD( Infty() ), +ctx) );
  }

  /// Context of the condition in IfThen.
  /// @param args: [condition, then, else] result variables
  /// @param ctx: context of the overall expression
  template <class Array3>
  void PropagateIfThenResultIntoCondition(Array3 args, Context ctx) {
    Context ctx_cond = Context::CTX_MIX;
    if (ctx.IsPositive() || ctx.IsNegative()) {
      if (MPCD( lb(args[1]) ) >= MPCD( ub(args[2]) ))
        ctx_cond = +ctx;
      else if (MPCD( lb(args[2]) ) >= MPCD( ub(args[1]) ))
        ctx_cond = -ctx;
    }
    MPD( PropagateResultOfInitExpr(args[0], 0.0, 1.0, ctx_cond) );
  }

  void PropagateResult(ImplicationConstraint& con, double lb, double ub, Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    auto& args = con.GetArguments();
    MPD( PropagateResultOfInitExpr(args[0], 0.0, 1.0, Context::CTX_MIX) );
    MPD( PropagateResultOfInitExpr(args[1], 0.0, 1.0, +ctx) );
    MPD( PropagateResultOfInitExpr(args[2], 0.0, 1.0, +ctx) );
  }

  void PropagateResult(AllDiffConstraint& con, double lb, double ub, Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    MPD( PropagateResult2Vars(con.GetArguments(), MPD( MinusInfty() ), MPD( Infty() ),
                         Context::CTX_MIX) );
  }

  /// @todo Propagate CTX+/- into the conditional equalities
  void PropagateResult(NumberofConstConstraint& con, double lb, double ub, Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    MPD( PropagateResult2Vars(con.GetArguments(), MPD( MinusInfty() ), MPD( Infty() ),
                         Context::CTX_MIX) );
  }

  void PropagateResult(NumberofVarConstraint& con, double lb, double ub, Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    MPD( PropagateResult2Vars(con.GetArguments(), MPD( MinusInfty() ), MPD( Infty() ),
                         Context::CTX_MIX) );
  }

  void PropagateResult(CountConstraint& con, double lb, double ub, Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    MPD( PropagateResult2Vars(con.GetArguments(),       // forward same context
                              MPD( MinusInfty() ), MPD( Infty() ), ctx) );
  }


  //////////////////////////// NONLINEAR /////////////////////////////

  /// \brief PropagateResult(PowConstraint)
  void PropagateResult(
      PowConstraint& con, double , double , Context ctx) {
    con.AddContext(ctx);           // merge context
    auto x = con.GetArguments()[0];
    auto y = con.GetArguments()[1];
    Context ctx_x, ctx_y;
    if (MPCD( lb(x) )>=0.0 ) {     // x >= 0
      ctx_x = (MPCD( lb(y) ) >= 0.0)
                  ? +ctx
                  : MPCD( ub(y) ) <= 0.0
                        ? -ctx: Context::CTX_MIX;
      ctx_y = (MPCD( ub(x) ) <= 1.0)  // 0<=x<=1
                  ? -ctx
                  : MPCD( lb(x) ) >= 1.0  // x>=1
                        ? +ctx
                        : Context::CTX_MIX;
    } else {
      ctx_x = ctx_y = Context::CTX_MIX;
    }
    VarArray1 xx {x};
    PropagateResult2Args(xx,
                         MPD( MinusInfty() ), MPD( Infty() ),
                         ctx_x);
    VarArray1 yy {y};
    PropagateResult2Args(yy,
                         MPD( MinusInfty() ), MPD( Infty() ),
                         ctx_y);
  }


  void PropagateResult(
      PowConstExpConstraint& con, double , double , Context ctx) {
    con.AddContext(ctx);           // merge context
    auto pwr = con.GetParameters()[0];
    auto arg = con.GetArguments()[0];
    Context ctx_new;
    bool is_pow_int = (MPD( is_integer_value(pwr) ));
    bool is_pow_odd = is_pow_int
        && !MPD( is_integer_value(pwr / 2.0) );
    bool is_pow_odd_pos = (pwr >= 0.0 && is_pow_odd);
    if (is_pow_odd_pos) {            // some monotone cases
      ctx_new = +ctx;
    } else
      if (MPD( lb(arg) )>=0.0 ) {     // arg >= 0
        ctx_new = (pwr >= 0.0) ? +ctx : -ctx;
      } else
        if (MPD( ub(arg) )<=0.0 && is_pow_int) {
          if (pwr >= 0.0)
            ctx_new = -ctx;          // arg<=0, pwr even >=0
          else                       // pwr < 0
            ctx_new = is_pow_odd ?
                          -ctx : +ctx;
        } else
          ctx_new.Add(Context::CTX_MIX);  // undecidable
    PropagateResult2Args(con.GetArguments(),
                         MPD( MinusInfty() ), MPD( Infty() ),
                         ctx_new);
  }

  void PropagateResult(LogConstraint& con, double , double , Context ctx) {
    con.AddContext(ctx);           // merge context
    PropagateResult2Args(con.GetArguments(),     // monotone
                         MPD( MinusInfty() ), MPD( Infty() ), +ctx);
  }

  void PropagateResult(LogAConstraint& con, double , double , Context ctx) {
    con.AddContext(ctx);           // merge context
    auto ctx_new = (con.GetParameters()[0]>=0.0) ? +ctx : -ctx;
    PropagateResult2Args(con.GetArguments(),     // monotone
                         MPD( MinusInfty() ), MPD( Infty() ), ctx_new);
  }

  void PropagateResult(ExpConstraint& con, double , double , Context ctx) {
    con.AddContext(ctx);           // merge context
    PropagateResult2Args(con.GetArguments(),     // monotone
                         MPD( MinusInfty() ), MPD( Infty() ), +ctx);
  }

  void PropagateResult(ExpAConstraint& con, double , double , Context ctx) {
    con.AddContext(ctx);           // merge context
    auto ctx_new = (con.GetParameters()[0]>=0.0) ? +ctx : -ctx;
    PropagateResult2Args(con.GetArguments(),     // monotone
                         MPD( MinusInfty() ), MPD( Infty() ), ctx_new);
  }


  //////////////////////////// ALGEBRAIC /////////////////////////////

  template <class Body, int kind>
  void PropagateResult(
      ConditionalConstraint<
        AlgebraicConstraint< Body, AlgConRhs<kind> > >& con,
      double lb, double ub, Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    if (lb>0 && ctx.HasPositive()) {              // Is true
      if constexpr (kind*kind<=1) {               // == or <= or >=
        MPD(AddConstraint_AS_ROOT(con.GetConstraint()));
        MPD( DecrementVarUsage(con.GetResultVar()) );
      }
    } else {
      MPD( PropagateResult(con.GetConstraint(), ctx) );
    }
  }


  /// Propagate given bounds & context into arguments of a constraint.
  /// The default template assumes it just a vector of variables.
  /// @param lb, ub: bounds for each variable
  template <class Args>
  void PropagateResult2Args(
      const Args& vars, double lb, double ub, Context ctx) {
    PropagateResult2Vars(vars, lb, ub, ctx);
  }

  /// Specialize: propagate result into LinTerms
  void PropagateResult2Args(
      const LinTerms& lint, double lb, double ub, Context ctx) {
    PropagateResult2LinTerms(lint, lb, ub, ctx);
  }

  /// Specialize: propagate result into QuadAndLinTerms
  void PropagateResult2Args(
      const QuadAndLinTerms& qlt, double lb, double ub, Context ctx) {
    PropagateResult2QuadAndLinTerms(qlt, lb, ub, ctx);
  }

  /// Propagate result into QuadAndLinTerms
  void PropagateResult2QuadAndLinTerms(
      const QuadAndLinTerms& qlt, double lb, double ub, Context ctx) {
    auto LB_lt = MPD( MinusInfty() );
    auto UB_lt = MPD( Infty() );
    if (qlt.GetQPTerms().empty()) {
      LB_lt = lb;
      UB_lt = ub;
    }
    PropagateResult2LinTerms(qlt.GetLinTerms(),
                             LB_lt, UB_lt, ctx);
    PropagateResult2QuadTerms(qlt.GetQPTerms(),
                              MPD( MinusInfty() ), MPD( Infty() ), ctx);
  }

  /// Propagate result into LinTerms
  void PropagateResult2LinTerms(const LinTerms& lint,
                                double lb, double ub, Context ctx) {
    auto LB = MPD( MinusInfty() );
    auto UB = MPD( Infty() );
    if (1==lint.size()) {
      auto coef = lint.coef(0);
      if (coef > 0.0) {
        LB = lb / coef;
        UB = ub / coef;
      } else if (coef < 0.0) {
        LB = ub / coef;
        UB = lb / coef;
      }
    }
    for (auto i=lint.size(); i--; ) {
      if (0.0!=std::fabs(lint.coef(i))) {
        MPD( PropagateResultOfInitExpr(lint.var(i), LB, UB,
                                (lint.coef(i)>=0.0) ? +ctx : -ctx) );
      }
    }
  }

  /// Propagate given bounds & context into a vector of variables
  /// @param lb, ub: bounds for each variable
  template <class Vec>
  void PropagateResult2Vars(const Vec& vars, double lb, double ub, Context ctx) {
    for (auto v: vars) {
      MPD( PropagateResultOfInitExpr(v, lb, ub, ctx) );
    }
  }

  /// Propagate result into QuadTerms
  void PropagateResult2QuadTerms(const QuadTerms& quadt, double , double , Context ctx) {
    for (auto i=quadt.size(); i--; ) {
      double coef_i = quadt.coef(i);
      if (0.0!=std::fabs(coef_i)) {
        // Propagate context in some cases.
        auto var1 = quadt.var1(i), var2 = quadt.var2(i);
        auto ctx12 = (coef_i>0.0) ? +ctx : -ctx;
        if (MPD( lb(var1) ) >= 0.0 && MPD( lb(var2) ) >= 0.0) {
          // leave as is
        } else if (MPD( ub(var1) ) <= 0.0 && MPD( ub(var2) ) <= 0.0) {
          ctx12 = -ctx12;
        } else // Propagate mixed if not decidable, otherwise we miss some cases
          ctx12 = Context::CTX_MIX;
        MPD( PropagateResultOfInitExpr(var1, ctx12) );
        if (var1!=var2)
          MPD( PropagateResultOfInitExpr(var2, ctx12) );
      }
    }
  }

  /// Do nothing
  inline void PropagateResult(NLConstraint& , double , double , Context ) { }
  /// Do nothing
  inline void PropagateResult(NLLogical& , double , double , Context ) { }
  /// Do nothing
  inline void PropagateResult(NLReifEquiv& , double , double , Context ) { }
  /// Do nothing
  inline void PropagateResult(NLReifImpl& , double , double , Context ) { }
  /// Do nothing
  inline void PropagateResult(NLReifRimpl& , double , double , Context ) { }

};

} // namespace mp

#endif // CONSTR_PROP_DOWN_H
