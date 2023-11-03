#ifndef CONSTR_PROP_DOWN_H
#define CONSTR_PROP_DOWN_H

/**
 * Propagate flat constraints from result (result bounds & context)
 * "down", i.e., to the arguments
 */

#include "mp/common.h"
#include "mp/flat/constr_std.h"

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
    PropagateResult2LinTerms(con.GetAffineExpr(),
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
  template <class Body>
  void PropagateResult(AlgebraicConstraint<Body, AlgConRange>& con) {
    /// Distinguish bounds' finiteness for context
		auto ctx = con.lb()<=MPD( PracticallyMinusInf() ) ?
					Context::CTX_NEG : con.ub()>=MPD( PracticallyInf() ) ?
            Context::CTX_POS : Context::CTX_MIX;
    PropagateResult2Args(con.GetBody(), con.lb(), con.ub(), ctx);
  }

  template <class Body, int sens>
  void PropagateResult(IndicatorConstraint<
                         AlgebraicConstraint< Body, AlgConRhs<sens> > >& con,
                       double lb, double ub, Context ctx) {
    internal::Unused(lb, ub, ctx);
    MPD( PropagateResultOfInitExpr(con.get_binary_var(),
                              MPD( MinusInfty() ), MPD( Infty() ),
                              1==con.get_binary_value() ?  // b==1 means b in CTX_NEG
                                Context::CTX_NEG : Context::CTX_POS) );
    PropagateResult2Args(con.get_constraint().GetBody(),   // Assume Con::BodyType is handled
                             MPD( MinusInfty() ), MPD( Infty() ),
                             0==sens ? Context::CTX_MIX :
                                       0<sens ? +ctx : -ctx);
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
    internal::Unused(lb, ub, ctx);
    MPD( PropagateResult2Args(con.GetExpression().GetBody(),
                         lb, ub, Context::CTX_MIX) );
    MPD( PropagateResultOfInitExpr(con.GetVariable(), lb, ub, Context::CTX_MIX) );
  }

  void PropagateResult(NotConstraint& con, double lb, double ub, Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    MPD( PropagateResultOfInitExpr(con.GetArguments()[0], 1.0-ub, 1.0-lb, -ctx) );
  }

  void PropagateResult(AndConstraint& con, double lb, double ub, Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    MPD( PropagateResult2Vars(con.GetArguments(), lb, 1.0, +ctx) );
    if (lb>0.5)                                 // Remove, arguments are fixed
      MPD( DecrementVarUsage(con.GetResultVar()) );
  }

  void PropagateResult(OrConstraint& con, double lb, double ub, Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    MPD( PropagateResult2Vars(con.GetArguments(), 0.0, ub, +ctx) );
    if (ub<=0.5)                                 // Remove, arguments are fixed
      MPD( DecrementVarUsage(con.GetResultVar()) );
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
  /// @args: [condition, then, else] result variables
  /// @ctx: context of the overall expression
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


  //////////////////////////// NONLINEAR /////////////////////////////

  void PropagateResult(PowConstraint& con, double , double , Context ctx) {
    con.AddContext(ctx);           // merge context
    auto pwr = con.GetParameters()[0];
    auto ctx_new = (pwr>=0.0 &&    // some monotone cases
                    MPD( is_integer_value(pwr) ) &&
                    !MPD( is_integer_value(pwr / 2.0) )) ||
        MPD( lb(con.GetArguments()[0])>=0.0 ) ? ctx : Context::CTX_MIX;
    PropagateResult2Args(con.GetArguments(),
                         MPD( MinusInfty() ), MPD( Infty() ), ctx_new);
  }

  void PropagateResult(LogConstraint& con, double , double , Context ctx) {
    con.AddContext(ctx);           // merge context
    PropagateResult2Args(con.GetArguments(),     // monotone
                         MPD( MinusInfty() ), MPD( Infty() ), ctx);
  }

  void PropagateResult(LogAConstraint& con, double , double , Context ctx) {
    con.AddContext(ctx);           // merge context
    auto ctx_new = (con.GetParameters()[0]>=0.0) ? ctx : -ctx;
    PropagateResult2Args(con.GetArguments(),     // monotone
                         MPD( MinusInfty() ), MPD( Infty() ), ctx_new);
  }

  void PropagateResult(ExpConstraint& con, double , double , Context ctx) {
    con.AddContext(ctx);           // merge context
    PropagateResult2Args(con.GetArguments(),     // monotone
                         MPD( MinusInfty() ), MPD( Infty() ), ctx);
  }

  void PropagateResult(ExpAConstraint& con, double , double , Context ctx) {
    con.AddContext(ctx);           // merge context
    PropagateResult2Args(con.GetArguments(),     // monotone
                         MPD( MinusInfty() ), MPD( Infty() ), ctx);
  }


  //////////////////////////// ALGEBRAIC /////////////////////////////

  void PropagateResult(CondLinConEQ& con, double lb, double ub, Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    MPD( PropagateResult2LinTerms(con.GetConstraint().GetBody(),
                         lb, ub, Context::CTX_MIX) );
  }

  void PropagateResult(CondQuadConEQ& con, double lb, double ub, Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    MPD( PropagateResult2QuadAndLinTerms(con.GetConstraint().GetBody(),
                         lb, ub, Context::CTX_MIX) );
  }

  template <class Body, int kind>
  void PropagateResult(
      ConditionalConstraint<
        AlgebraicConstraint< Body, AlgConRhs<kind> > >& con,
      double lb, double ub, Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    auto ctx_new = kind>0 ? ctx : kind<0 ? -ctx : Context::CTX_MIX;
    MPD( PropagateResult2Args(con.GetConstraint().GetBody(), lb, ub,
                             ctx_new) );
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
    PropagateResult2LinTerms(qlt.GetLinTerms(), lb, ub, ctx);
    PropagateResult2QuadTerms(qlt.GetQPTerms(), lb, ub, ctx);
  }

  /// Propagate result into LinTerms
  void PropagateResult2LinTerms(const LinTerms& lint, double , double , Context ctx) {
    for (auto i=lint.size(); i--; ) {
      if (0.0!=std::fabs(lint.coef(i))) {
        MPD( PropagateResultOfInitExpr(lint.var(i),
                                MPD( MinusInfty() ), MPD( Infty() ),
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
      if (0.0!=std::fabs(quadt.coef(i))) {
        // Propagate context in some cases.
        auto var1 = quadt.var1(i), var2 = quadt.var2(i);
        auto ctx12 = ctx;
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

};

} // namespace mp

#endif // CONSTR_PROP_DOWN_H
