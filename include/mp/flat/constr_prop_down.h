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

  /// By default, set mixed context for argument variables
  template <class Constraint>
  void PropagateResult(Constraint& con, double lb, double ub, Context ctx) {
    internal::Unused(&con, lb, ub, ctx);
    con.SetContext(ctx);
    PropagateResult2Args(con.GetArguments(),
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

  void PropagateResult(QuadConRange& con, double lb, double ub,
                       Context ctx) {
    internal::Unused(lb, ub, ctx);
    PropagateResult2LinTerms(con.GetLinTerms(), // TODO sense dep. on bounds
                             MPD( MinusInfty() ), MPD( Infty() ), ctx);
    PropagateResult2QuadTerms(con.GetQPTerms(), // TODO bounds?
                              MPD( MinusInfty() ), MPD( Infty() ), ctx);
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

  void PropagateResult(ComplementarityLinear& con, double lb, double ub,
                       Context ctx) {
    internal::Unused(lb, ub, ctx);
    MPD( PropagateResult2LinTerms(con.GetExpression().GetLinTerms(),
                         lb, ub, Context::CTX_MIX) );
    MPD( PropagateResultOfInitExpr(con.GetVariable(), lb, ub, Context::CTX_MIX) );
  }

  void PropagateResult(ComplementarityQuadratic& con, double lb, double ub,
                       Context ctx) {
    internal::Unused(lb, ub, ctx);
    MPD( PropagateResult2LinTerms(con.GetExpression().GetLinTerms(),
                    lb, ub, Context::CTX_MIX) );
    MPD( PropagateResult2QuadTerms(con.GetExpression().GetQPTerms(),
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
  }

  void PropagateResult(OrConstraint& con, double lb, double ub, Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    MPD( PropagateResult2Vars(con.GetArguments(), 0.0, ub, +ctx) );
  }

  void PropagateResult(IfThenConstraint& con, double lb, double ub, Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    auto& args = con.GetArguments();
    /// TODO consider bounds for then/else for the context:
    MPD( PropagateResultOfInitExpr(args[0], 0.0, 1.0, Context::CTX_MIX) );
    MPD( PropagateResultOfInitExpr(args[1], MPD( MinusInfty() ), MPD( Infty() ), +ctx) );
    MPD( PropagateResultOfInitExpr(args[2], MPD( MinusInfty() ), MPD( Infty() ), -ctx) );
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
    MPD( PropagateResult2Args(con.GetConstraint().GetBody(), lb, ub,
                             kind>0 ? ctx : -ctx) );
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
      MPD( PropagateResultOfInitExpr(lint.var(i),      /// TODO bounds as well
                                MPD( MinusInfty() ), MPD( Infty() ),
                                (lint.coef(i)>=0.0) ? +ctx : -ctx) );
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
  void PropagateResult2QuadTerms(const QuadTerms& quadt, double , double , Context ) {
    for (auto i=quadt.size(); i--; ) {             /// TODO context for special cases
      MPD( PropagateResultOfInitExpr(quadt.var1(i),     /// TODO bounds as well
                                MPD( MinusInfty() ), MPD( Infty() ), Context::CTX_MIX) );
      MPD( PropagateResultOfInitExpr(quadt.var2(i),
                                MPD( MinusInfty() ), MPD( Infty() ), Context::CTX_MIX) );
    }
  }

};

} // namespace mp

#endif // CONSTR_PROP_DOWN_H
