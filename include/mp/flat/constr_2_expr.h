#ifndef CONSTR_2_EXPR_H
#define CONSTR_2_EXPR_H

#include <functional>
#include <type_traits>

#include "mp/flat/constr_keeper.h"
#include "mp/flat/nl_expr/constr_nl.h"
#include "mp/valcvt-link.h"

namespace mp {

/// A mix-in base class
/// facilitating "inlining" of some functional constraints
/// into expression trees.
template <class Impl>
class Constraints2Expr {
public:
  /// Mark some functional constraints as expressions.
  /// Then convert some objectives and constraints
  /// to "NL...." items using the expressions.
  /// The auxiliary result variables of the functional
  /// constraints which were marked as expressions,
  /// will be handled in a special way.
  ///
  /// @note Beyond the initial marking stage,
  /// the result variables of any new constraints/expressions,
  /// and any other new variables should be marked too.
  void Convert2NL() {
    MPD( MarkExpressions() );
    /// Objectives before constraints,
    /// because if we produce new expressions and their result variables
    /// need to be explicit, this is done in constraints conversion.
    MPD( ConvertObjectivesWithExpressions() );
    stage_cvt2expr_ = 1;               // static cons
                               // -> new static cons + func cons
    MPD( GetModel() ).ConvertAllWithExpressions(*(Impl*)this);
    stage_cvt2expr_ = 2;               // func cons -> explicifiers
    MPD( GetModel() ).ConvertAllWithExpressions(*(Impl*)this);
    MPD( EliminateExprResultVars() );  // In the very end
    MPD( PassHooksToModelAPI() );
  }

  /// Mark which functional constraints to be used as expressions,
  /// vs assigning their result to a variable
  void MarkExpressions() {
    MPD( MarkAllResultVarsAsVars() );      // tentatively
    MPD( GetModel() ).MarkExprResultVars(*(Impl*)this); // As expressions
    MPD( GetModel() ).MarkArguments(*(Impl*)this);   // As vars if required
  }

  /// Consider marking the result variables as
  /// possible expressions
  template <class Con>
  void ConsiderMarkingResultVar(
      const Con& con, int i, ExpressionAcceptanceLevel eal) {
    assert(ExpressionAcceptanceLevel::NotAccepted!=eal);
    if (con.HasResultVar()) {         // A functional constraint
      assert(                     // Check: the result var has \a con as the init expr
          MPD( template GetInitExpressionOfType<Con>(con.GetResultVar()) )
          == &con);
      if (con.IsLogical()         // Fixed logical results handled differently
          || !MPD( IfVarBoundsStrongerThanInitExpr(con.GetResultVar()) ) )
        MPD( MarkAsExpression(con.GetResultVar()) );   // can be changed later
    }
  }

  /// Consider marking the argument variables as
  /// "explicit variables" (not expressions.)
  /// Generic request.
  /// For func cons, once not accepted as expressions, consider their args
  /// as variables.
  /// For static cons, always.
  template <class Con>
  void ConsiderMarkingArguments(
      const Con& con, int i, ExpressionAcceptanceLevel eal) {
    bool fMarkArgs = false;
    if (con.HasResultVar())    // func cons: those not accepted as expr
      fMarkArgs = (ExpressionAcceptanceLevel::NotAccepted==eal);
    else
      fMarkArgs = true;        // static cons: all non-algebraic by default
    if (fMarkArgs)
      MPD( DoMarkArgsAsVars(con, i) );
  }


  /// Convert algebraic constraint for use with expressions
  /// @param i: constraint index
  /// @param cal: (user chosen) acceptance level for \a con.
  ///   While it is recommended that flat algebraic constraints
  ///   are still accepted along with NLConstraint,
  ///   user might set acc:linle=0 etc.
  /// @return true if this constraint has been eliminated/replaced.
  template <class Body, class RhsOrRange>
  bool ConvertWithExpressions(
      const AlgebraicConstraint<Body, RhsOrRange>& con,
      int i,
      ConstraintAcceptanceLevel cal, ExpressionAcceptanceLevel ) {
    assert(stage_cvt2expr_>0);
    /// Replace \a con by a NLConstraint,
    /// if either the ModelAPI does not accept it,
    /// or the linear/quadratic terms have expressions
    /// (and then they are non-flat.)
    if (1==stage_cvt2expr_
        && (ConstraintAcceptanceLevel::Recommended != cal
            || HasExpressionArgs(con.GetBody()))) {
      auto alscope = MPD( MakeAutoLinker( con, i ) );   // link from \a con
      if (HasLogicalExpressionArgs(con.GetBody()))
        if (HandleLogicalArgs(con, i))
          return true;                     // to remove the original \a con
      return ConvertToNLCon(con, i);       // convert if expr's, or not acc
    }
    return false;
  }

  /// Convert conditional algebraic constraint for use with expressions.
  /// Basically we need the LHS to be an expression.
  /// This can be simplified if we switch to do it earlier in flat conversions.
  /// @return true if this constraint has been eliminated/replaced.
  template <class Body, class RhsOrRange>
  bool ConvertWithExpressions(
      const ConditionalConstraint< AlgebraicConstraint<Body, RhsOrRange> >& con,
      int i,
      ConstraintAcceptanceLevel , ExpressionAcceptanceLevel eal) {
    assert(stage_cvt2expr_>0 && stage_cvt2expr_<=2);
    assert(!con.GetContext().IsNone());
    // See if the item is going into an expr.
    // Otherwise it's a flat con.
    if (ExpressionAcceptanceLevel::NotAccepted != eal) {
      if (1==stage_cvt2expr_) {
        HandleLogicalArgs(con, i);         // expicify logical args
        if (!con.GetConstraint().GetBody().is_variable()) {  // already a variable
          ConvertConditionalConLHS(con, i);
          return true;
        }
      }
      else   // if (2==stage_cvt2expr_)
        ConsiderExplicifyingExpression(con, i);        // this is a func con too
    }
    return false;
  }

  /// Convert other functional constraints for use with expressions.
  /// Basically, see if the result variable is explicit.
  /// @return true if this constraint has been eliminated/replaced.
  template <class FuncCon,
           std::enable_if_t<
               std::is_base_of_v<FunctionalConstraint, FuncCon>, bool > = true >
  bool ConvertWithExpressions(
      const FuncCon& con, int i,
      ConstraintAcceptanceLevel , ExpressionAcceptanceLevel eal) {
    if (ExpressionAcceptanceLevel::NotAccepted != eal) {     // going into an expr
      if constexpr (std::is_base_of_v<NumericFunctionalConstraintTraits, FuncCon>)
        if (1==stage_cvt2expr_)
          HandleLogicalArgs(con, i);
      if (2==stage_cvt2expr_)
        ConsiderExplicifyingExpression(con, i);
    }
    return false;                          // leave it active
  }

  /// Special handling for LinearFunctionalConstraint
  bool ConvertWithExpressions(
      const LinearFunctionalConstraint& con, int i,
      ConstraintAcceptanceLevel , ExpressionAcceptanceLevel eal) {
    if (1==stage_cvt2expr_)
      HandleLogicalArgs(con, i);         // explicify logical args
    if (2==stage_cvt2expr_) {
      return ConsiderExplicifyingAlgebraic(con, i);
    }
    return false;
  }

  /// Special handling for LinearFunctionalConstraint
  bool ConvertWithExpressions(
      const QuadraticFunctionalConstraint& con, int i,
      ConstraintAcceptanceLevel , ExpressionAcceptanceLevel eal) {
    if (1==stage_cvt2expr_)
      HandleLogicalArgs(con, i);         // explicify logical args
    if (2==stage_cvt2expr_) {
      return ConsiderExplicifyingAlgebraic(con, i);
    }
    return false;
  }

  /// Convert complementarity constraint for use with expressions.
  /// Similarly to Conditional, we need the expression part to be an NL expression.
  /// @return true if this constraint has been eliminated/replaced.
  template <class Expr>
  bool ConvertWithExpressions(
      const ComplementarityConstraint<Expr>& con,
      int i,
      ConstraintAcceptanceLevel , ExpressionAcceptanceLevel ) {
    if (1==stage_cvt2expr_)
      HandleLogicalArgs(con, i);         // explicify logical args
    if (false                      // TODO check acc for NLCompl
        && 1==stage_cvt2expr_
        && !con.GetExpression().is_variable()) {             // already a variable
      ConvertComplementarityExpr(con, i);
      return true;
    }
    return false;
  }

  /// Template for Indicators: do nothing.
  /// But check that they are flat?
  template <class SubCon>
  bool ConvertWithExpressions(
      const IndicatorConstraint<SubCon>& , int ,
      ConstraintAcceptanceLevel , ExpressionAcceptanceLevel ) {
    return false;
  }

  /// SOS1: do nothing.
  /// But check that they are flat?
  bool ConvertWithExpressions(
      const SOS1Constraint& con, int ,
      ConstraintAcceptanceLevel , ExpressionAcceptanceLevel ) {
    if (2==stage_cvt2expr_)
    for (int v: con.GetArguments()) {
      assert(MPCD( IsProperVar(v) ));
    }
    return false;
  }

  /// SOS1: do nothing.
  /// But check that they are flat?
  bool ConvertWithExpressions(
      const SOS2Constraint& con, int ,
      ConstraintAcceptanceLevel , ExpressionAcceptanceLevel ) {
    if (2==stage_cvt2expr_)
    for (int v: con.GetArguments()) {
      assert(MPCD( IsProperVar(v) ));
    }
    return false;
  }

  /// NLConstraint: just produced.
  bool ConvertWithExpressions(
      const NLConstraint& , int ,
      ConstraintAcceptanceLevel , ExpressionAcceptanceLevel ) {
    return false;
  }
  /// NLAssignEQ: just produced.
  bool ConvertWithExpressions(
      const NLAssignEQ& , int ,
      ConstraintAcceptanceLevel , ExpressionAcceptanceLevel ) {
    return false;
  }
  /// NLAssignLE: just produced.
  bool ConvertWithExpressions(
      const NLAssignLE& , int ,
      ConstraintAcceptanceLevel , ExpressionAcceptanceLevel ) {
    return false;
  }
  /// NLAssignGE: just produced.
  bool ConvertWithExpressions(
      const NLAssignGE& , int ,
      ConstraintAcceptanceLevel , ExpressionAcceptanceLevel ) {
    return false;
  }

  /// NLLogical: just produced.
  bool ConvertWithExpressions(
      const NLLogical& , int ,
      ConstraintAcceptanceLevel , ExpressionAcceptanceLevel ) {
    return false;
  }
  /// NLBaseReif: just produced.
  template <int sense>
  bool ConvertWithExpressions(
      const NLBaseReif<sense>& , int ,
      ConstraintAcceptanceLevel , ExpressionAcceptanceLevel ) {
    return false;
  }

  /// Any other static con.
  template <class A, class P, class I>
  bool ConvertWithExpressions(
      const CustomStaticConstraint<A, P, I>& , int ,
      ConstraintAcceptanceLevel , ExpressionAcceptanceLevel ) {
    return false;
  }

  /// Convert objectives
  void ConvertObjectivesWithExpressions() {
    auto& objs_emulated = MPD( get_emulated_objectives() );
    auto& objs_original = MPD( get_objectives() );
    auto& objs
        = objs_emulated.size() ? objs_emulated : objs_original;
    for (size_t iobj=0; iobj<objs.size(); ++iobj) {
      HandleLogicalArgs(objs[iobj], iobj);
      Convert1ObjWithExpressions(iobj, objs[iobj]);
    }
  }

  /// Mark expr result vars for elimination
  void EliminateExprResultVars() {
    for (auto i = MPCD(num_vars()); i--; )
      if (!MPCD( IsProperVar(i) ))
        MPD( MarkVarAsEliminated(i) );
  }

  /// Pass some infos and callbacks to the ModelAPI
  void PassHooksToModelAPI() {
    MPD( GetModelAPI() ).PassVarProperFlags(
        MPCD( GetVarProperFlags() ));
    MPD( GetModelAPI() ).PassInitExprGetter(
        [this](int i_res_var, void* pexpr) {
      auto cloc = MPD( GetInitExpression(i_res_var) );
      cloc.StoreSolverExpression(
          MPD(GetModelAPI()), pexpr);
    });
  }


protected:
  /// Algebraic cons: no marking (even when NLConstraint not accepted.)
  /// @todo Do we need to consider NLConstraint / NLObjective at this step?
  /// Are they added during marking?
  template <class Body, class RhsOrRange>
  void DoMarkArgsAsVars(   // needs to appear before the most generic template
      const AlgebraicConstraint<Body, RhsOrRange>& , int ) { }

  /// Complementarity: only mark the var
  /// (NL accepts expressions for the expression part)
  template <class Expr>
  void DoMarkArgsAsVars(   // needs to appear before the most generic template
      const ComplementarityConstraint<Expr>& cc, int ) {
    MPD( MarkAsResultVar(cc.GetVariable()) );
  }

  /// Generic arguments marking call
  template <class Con>
  void DoMarkArgsAsVars(const Con& con, int ) {
    VisitArguments(con, MarkVar_);
  }


  /// Check if the algebraic body has expressions
  bool HasExpressionArgs(const QuadAndLinTerms& qlt) const {
    return HasExpressionArgs(qlt.GetLinTerms())
    || HasExpressionArgs(qlt.GetQPTerms());
  }

  bool HasExpressionArgs(const LinTerms& lt) const {
    for (auto v: lt.vars())
      if (!MPCD( IsProperVar(v) )) {
        return true;
      }
    return false;
  }

  bool HasExpressionArgs(const QuadTerms& qt) const {
    for (auto v: qt.get_folded())
      if (!MPCD( IsProperVar(v.first.first) )
          || !MPCD( IsProperVar(v.first.second) ))
        return true;
    return false;
  }

  /// Check if the algebraic body has logical expressions.
  /// This is bad for MP2NL.
  bool HasLogicalExpressionArgs(const QuadAndLinTerms& qlt) const {
    return HasLogicalExpressionArgs(qlt.GetLinTerms())
    || HasLogicalExpressionArgs(qlt.GetQPTerms());
  }

  bool HasLogicalExpressionArgs(const LinTerms& lt) const {
    for (auto v: lt.vars())
      if (!MPCD( IsProperVar(v) )             // still an expression
          && MPCD( IsInitExprLogical(v) )) {
        return true;
      }
    return false;
  }

  bool HasLogicalExpressionArgs(const QuadTerms& qt) const {
    for (auto v: qt.get_folded())
      if ((!MPCD( IsProperVar(v.first.first) )
           && MPCD( IsInitExprLogical(v.first.first) ))
          || (!MPCD( IsProperVar(v.first.second) )
              && MPCD( IsInitExprLogical(v.first.second) )))
        return true;
    return false;
  }

  /// Handle logical expressions in an algebraic con
  /// @return whether to remove the original \a con.
  template <class ConObj>
  bool HandleLogicalArgs(const ConObj& con, int ) {
    VisitArguments(con, MarkVarIfLogical_);          // Mark as proper vars
    return false;                                    // don't remove immediately
  }

  /// Handle logical expressions in an IfThen
  /// @return whether to remove the original \a con.
  template <typename fake> // to avoid "specialization" not declared at namespace
  // with certain versions of gcc
  bool HandleLogicalArgs(const IfThenConstraint& con, int ) {
    MarkVarIfLogical_(con.GetArguments()[1]);    // then part
    MarkVarIfLogical_(con.GetArguments()[2]);    // else part
    return false;
  }

  /// Special linear cases.
  /// Not doing any more because these simplifications
  /// interfere with result variable marking.
  /// These simplifications are general presolve
  /// and should better have been done in normal conversion stage.
  /// @todo atleast, atmost, exactly
  template <class RhsOrRange>
  bool HandleLogicalArgs_SpecialCases(
      const AlgebraicConstraint<LinTerms, RhsOrRange>& con, int ) {
    /*
    const auto& body = con.GetBody();
    if (!con.lb() && !con.ub()          // == 0.0
        && 2==body.size()) {            // 2 terms
      if (MPCD( IsInitExprLogical(body.var(0)) )
          && !MPCD( IsProperVar(body.var(0)) )
          && MPCD( IsInitExprLogical(body.var(1)) )
          && !MPCD( IsProperVar(body.var(1)) )
          && 1.0==std::fabs(body.coef(0))
          && 1.0==std::fabs(body.coef(1))
          && body.coef(0) == -body.coef(0)
          && MPCD( template ModelAPIOk<EquivalenceConstraint>() )) {
        int resvar = MPD( AssignResultVar2Args(
            EquivalenceConstraint({body.var(0), body.var(1)})) );
        MPD( FixAsTrue(resvar) );
        return true;
      }
    }
    else if (1==body.size()
             && MPCD( IsInitExprLogical(body.var(0)) )
             && !MPCD( IsProperVar(body.var(0)) )
             && body.coef(0))                    // != 0
    {
      if (!con.lb() && !con.ub()                 // == 0.0
          && MPCD( template ModelAPIOk<NotConstraint>() ) ) {
        int resvar = MPD( AssignResultVar2Args(
            NotConstraint({body.var(0)})) );
        MPD( FixAsTrue(resvar) );
        return true;
      }
      else if (con.lb() && con.ub() && con.lb()==con.ub()
               && con.lb()==body.coef(0) ) {
        MPD( FixAsTrue(body.var(0)) );
        return true;
      }
    } */
    return false;
  }

  /// Convert algebraic con to a \a NLConstraint,
  /// if \a NLConstraint's are accepted. Otherwise,
  /// explicify the expression and convert to
  /// \a LinConLE/EQ/GE/Range.
  /// @return true iff the original constraint should be deleted.
  /// @todo leave QCP terms here if accepted, even if other expr terms?
  template <class Body, class RhsOrRange>
  bool ConvertToNLCon(
      const AlgebraicConstraint<Body, RhsOrRange>& con, int ) {
    LinTerms lt;
    /// exprTerm will be a LinearFunctionalConstraint or a Quadratic...
    auto exprTerm = ExtractLinAndExprArgs(con.GetBody(), lt);
    assert(0.0 == exprTerm.GetArguments().constant_term());
    AlgConRange rng { con.GetRhsOrRange().lb(), con.GetRhsOrRange().ub() };
    /// Store full LFC only if it is not 1.0*var
    int exprResVar = -1;
    if (exprTerm.GetArguments().is_variable()) {
      exprResVar = exprTerm.GetArguments().get_representing_variable();
    } else if ( !exprTerm.GetArguments().empty() ) {  // has more terms, or coef != 1.0
      exprTerm.AddContext(                            // Context is compulsory
          rng.lb() > MPCD( PracticallyMinusInf() )    // Should not need to propagate
          ? (rng.ub() < MPCD( PracticallyInf() )
                ? Context::CTX_MIX : Context::CTX_POS)
              : Context::CTX_NEG);
      exprResVar = MPD( AssignResultVar2Args(std::move(exprTerm)) );
    }
    bool need_nlc {false};
    if (exprResVar >= 0) {                            // Some expressions are there
      if (!MPCD(VarHasMarking(exprResVar)))             // mark as expr if new
        MPD( MarkAsExpression(exprResVar) );
      if ( !MPCD( UserAcceptsAndRecommends((const NLConstraint*)nullptr) ) )
        MPD( MarkAsResultVar(exprResVar) );
      /// Exists and marked a variable
      if (MPCD( IsProperVar(exprResVar) )) {            // Not an expression after all
        lt.add_term(1.0, exprResVar);        // @todo When exprTerm was originally a var,
        lt.sort_terms();                     // this would reproduce the original con.
        if (lt.size()>1) {                              // ... has other variables
          if (MPCD( UserAcceptsAndRecommends(       // Accepts LinCon..
                  (const AlgebraicConstraint<LinTerms, RhsOrRange>*)nullptr) )) {
            AlgebraicConstraint<LinTerms, RhsOrRange> lc {lt, con.GetRhsOrRange(), false};
            MPD( AddConstraint( std::move(lc) ) );
            return true;
          }
          need_nlc = true;
        } else {         // single variable, its expression will be explicified
          MPD( IncrementVarUsage(exprResVar) );  // Because removed from top-level con #201
          MPD( NarrowVarBounds(exprResVar, rng.lb(), rng.ub()) );
          return true;
        }
        exprResVar = -1;                                // no expression
      }
    }
    if (0<=exprResVar                                   // either: have expression
        || !MPCD( UserAcceptsAndRecommends(         // or, not accepts source \a con
            (const AlgebraicConstraint<Body, RhsOrRange>*)nullptr) )
        || need_nlc) {                                  // or, other reason
      assert( MPCD( UserAcceptsAndRecommends((const NLConstraint*)nullptr) ) );
      NLConstraint nlc{lt, exprResVar, rng, false};     // false: no sort any more
      MPD( AddConstraint( std::move(nlc) ) );
      return true;
    }
    return false;
  }

  /// Convert (if needed) an objective to NLObjective.
  /// Logic similar to NLConstraint.
  void Convert1ObjWithExpressions(int iobj, QuadraticObjective& qobj) {
    LinTerms lt_varsonly;
    LinTerms lt_in_expr = SplitLinTerms(qobj.GetLinTerms(), lt_varsonly);
    // Have expression(s) or QP terms?
    // Might need to hide them into the expression part.
    // @todo leave QP terms here if accepted, even if other expr terms?
    if (lt_in_expr.size()
        || (qobj.GetQPTerms().size()
            && (!MPCD(IfPassQuadObj())         // cannot or want not
                || HasExpressionArgs(qobj.GetQPTerms())))) {
      MPD( UncountArgRefs(qobj) );
      int exprResVar = -1;
      if (lt_in_expr.is_variable() && qobj.GetQPTerms().empty()) {
        exprResVar = lt_in_expr.get_representing_variable();
      } else {                      // We need a new expression
        // Set up AutoLink
        auto obj_src =              // source value node for this obj
            MPD( GetValuePresolver() ).GetSourceNodes().GetObjValues()().Select(iobj);
        pre::AutoLinkScope<Impl> auto_link_scope{ *(Impl*)this, obj_src };
        if (qobj.GetQPTerms().empty())
          exprResVar = MPD( AssignResultVar2Args(
              LinearFunctionalConstraint{ {lt_in_expr, 0.0} } ) );
        else {                       // Move QP terms into the expr
          exprResVar = MPD( AssignResultVar2Args(
              QuadraticFunctionalConstraint
              { {{lt_in_expr, std::move(qobj.GetQPTerms())}, 0.0} } ) );
          qobj.GetQPTerms().clear();           // std::move() does not clear
        }
        MPD( AddInitExprContext(exprResVar,             // Context is compulsory
                               obj::MAX==qobj.obj_sense_true()    // no need to propagate
                                   ? Context::CTX_POS : Context::CTX_NEG) );
      }
      if ( !MPCD(VarHasMarking(exprResVar) ))         // mark as expr if new
        MPD( MarkAsExpression(exprResVar) );
      if ( !MPCD( GetModelAPI() ).AcceptsNLObj() )    // But as var if NLObj not accepted
        MPD( MarkAsResultVar(exprResVar) );
      if ( MPCD( IsProperVar(exprResVar) ) ) {        // Not an expression after all
        lt_varsonly.add_term(1.0, exprResVar);
        exprResVar = -1;                              // no expression
        lt_varsonly.sort_terms();
      }
      qobj.GetLinTerms() = lt_varsonly;
      if (exprResVar>=0)
        qobj.SetExprIndex(exprResVar);
      MPD( CountArgRefs(qobj) );
    }
  }

  /// Extract linear and expression args
  /// from Quad+Linear terms
  /// @param ltout: pure linear part (coefs*vars)
  ///   for NLConstraint
  /// @return expression part (QFC)
  QuadraticFunctionalConstraint ExtractLinAndExprArgs(
      const QuadAndLinTerms& qltin,
      LinTerms& ltout) {
    /// Create full QFC because we have QP terms
    LinTerms lt_in_expr = SplitLinTerms(qltin.GetLinTerms(), ltout);
    assert(qltin.GetQPTerms().size());
    return
        {{{std::move(lt_in_expr), qltin.GetQPTerms()}, 0.0}};
  }

  /// Extract linear and expression args
  /// from Linear terms.
  /// @param ltout: pure linear part (coefs*vars)
  ///   for NLConstraint
  /// @return expression part (LFC)
  LinearFunctionalConstraint ExtractLinAndExprArgs(
      const LinTerms& ltin,
      LinTerms& ltout) {
    LinTerms lt_in_expr = SplitLinTerms(ltin.GetLinTerms(), ltout);
    return {{std::move(lt_in_expr), 0.0}};
  }

  /// Split linterms into pure-var and expressions
  /// @return the expressions part
  LinTerms SplitLinTerms(const LinTerms& ltin, LinTerms& lt_out_vars) {
    LinTerms result;
    int nvars=0;
    for (int v: ltin.vars())
      if (MPCD( IsProperVar(v) ))
        ++nvars;
    lt_out_vars.reserve(nvars);
    result.reserve(ltin.size() - nvars);
    int v=0;
    double c=0.0;
    for (size_t i=0; i<ltin.size(); ++i) {
      if ((c = ltin.coef(i))) {                       // non-0
        if (MPCD( IsProperVar(v = ltin.var(i)) )) {
          lt_out_vars.add_term(c, v);
        } else {
          result.add_term(c, v);
        }
      }
    }
    return result;
  }

  /// Convert the LHS of a conditional con
  template <class Body, class RhsOrRange>
  void ConvertConditionalConLHS(
      const ConditionalConstraint< AlgebraicConstraint<Body, RhsOrRange> >& con,
      int i) {
    auto alscope = MPD( MakeAutoLinker( con, i ) );       // link from \a con
    /// Create a functional constraint from the LHS
    auto fc = MakeFunctionalConstraint(
        AlgebraicExpression{con.GetConstraint().GetArguments(), 0.0});
    auto resvar = MPD( AssignResultVar2Args(std::move(fc)) );
    if ( !MPCD(VarHasMarking(resvar) ))         // mark as expr if new
      MPD( MarkAsExpression(resvar) );
    /// resvar can be a proper variable - ModelAPI should flexibly handle this
    LinTerms lt { {1.0}, {resvar} };
    ConditionalConstraint< AlgebraicConstraint<LinTerms, RhsOrRange> >
        ccnew { { std::move(lt), con.GetConstraint().GetRhsOrRange() } };
    MPD( RedefineVariable(con.GetResultVar(), std::move(ccnew)) );  // Use new CondCon
    MPD( PropagateResultOfInitExpr(con.GetResultVar(), con.GetContext()) ); // context
  }

  /// Convert the expression part of complementarity.
  /// Similar to the argument of a conditional con.
  template <class Expr>
  void ConvertComplementarityExpr(
      const ComplementarityConstraint<Expr>& con,
      int i) {
    auto alscope = MPD( MakeAutoLinker( con, i ) );       // link from \a con
    /// Create a functional constraint from the LHS
    auto fc = MakeFunctionalConstraint(con.GetExpression());
    fc.SetContext(Context::CTX_MIX);                      // need context
    auto resvar = MPD( AssignResultVar2Args(std::move(fc)) );
    if ( !MPCD(VarHasMarking(resvar) ))         // mark as expr if new
      MPD( MarkAsExpression(resvar) );
    /// resvar can be a proper variable - ModelAPI should flexibly handle this
    LinTerms lt { {1.0}, {resvar} };
    ComplementarityConstraint< AlgebraicExpression<LinTerms> >
        ccnew { AffineExpr{ std::move(lt), 0.0 }, con.GetVariable() };
    MPD( AddConstraint(std::move(ccnew)) );  // Use new CondCon
  }

  /// Consider explicifying an expression
  template <class FuncCon,
           std::enable_if_t<
               std::is_base_of_v<          // functional cons only
                   FunctionalConstraint, FuncCon>, bool > = true >
  void ConsiderExplicifyingExpression(const FuncCon& con, int i) {
    if (MPCD( IsProperVar(con.GetResultVar()) ))
      DoExplicify(con, i);
    else {
      if (con.IsLogical()) {
        assert(con.HasResultVar());        // is a func con
        auto resvar = con.GetResultVar();
        if (MPCD( is_fixed(resvar) )) {    // we don't explicify just for this
          auto val = MPCD( fixed_value(resvar) );
          if constexpr (std::is_same_v<FuncCon, NotConstraint>) {
            assert(MPCD( is_fixed(con.GetArguments()[0]) ));
            assert(val != MPCD( fixed_value(con.GetArguments()[0]) ));
            // that's it, leave it because
            // we should handle the arg expr separately as below
          } else {
            if ((val && con.GetContext().HasPositive())
                || (!val && con.GetContext().HasNegative())) {
              MPD( AddConstraint(NLLogical(resvar, val)) );  // static con
            }  // else, skip
          }
        }
      }
    }
  }

  /// Add expr = var assignment for algebraic expression
  template <class FuncCon,
           std::enable_if_t<
               std::is_base_of_v<
                   NumericFunctionalConstraintTraits, FuncCon>, bool > = true >
  void DoExplicify(const FuncCon& con, int i) {
    auto alscope = MPD( MakeAutoLinker( con, i ) );       // link from \a con
    assert(!con.GetContext().IsNone());
    auto resvar = con.GetResultVar();
    if (con.GetContext().IsMixed())
      MPD( AddConstraint(NLAssignEQ(resvar)) );
    else if (con.GetContext().HasPositive())
      MPD( AddConstraint(NLAssignLE(resvar)) );
    else
      MPD( AddConstraint(NLAssignGE(resvar)) );
  }

  /// Add expr = var assignment (NLReifEquiv, NLReifImpl, NLRImpl)
  /// for logical expression.
  /// Even if this expression is root-level (is true).
  template <class FuncCon,
           std::enable_if_t<
               std::is_base_of_v<
                   LogicalFunctionalConstraintTraits, FuncCon>, bool > = true >
  void DoExplicify(const FuncCon& con, int i) {
    auto alscope = MPD( MakeAutoLinker( con, i ) );       // link from \a con
    auto resvar = con.GetResultVar();
    assert(!con.GetContext().IsNone());
    if (con.GetContext().IsMixed())
      MPD( AddConstraint(NLReifEquiv(resvar)) );
    else if (con.GetContext().HasPositive())
      MPD( AddConstraint(NLReifImpl(resvar)) );
    else
      MPD( AddConstraint(NLReifRimpl(resvar)) );
  }

  /// Special handling for algebraic functional constraints (LFC, QFC)
  /// @return whether the \a con should be deleted
  template <class AlgFuncCon>
  bool ConsiderExplicifyingAlgebraic(const AlgFuncCon& con, int i) {
    if (MPCD( IsProperVar(con.GetResultVar()) )) {
      using TargetCon = AlgebraicConstraint<
          std::decay_t<decltype(con.GetArguments().GetBody())>,
          AlgConRhs<0> >;  // @todo can be ,=, >=
      if (!MPCD( template ModelAPIOk< TargetCon >() )
          || HasExpressionArgs(con.GetArguments())) {
        DoExplicify(con, i);          // as other explicified expressions
        return false;
      }
      auto& ck = GET_CONSTRAINT_KEEPER(AlgFuncCon);
      const auto& ie = MPD( GetInitExpression(con.GetResultVar()) );
      ck.ConvertConstraint(ie.GetIndex());
      return true;
    }
    return false;
  }

  /// Recompute implicit aux vars
  /// (those corresponding to expressions.)
  /// needed for MO emulator and sol checker.
  void RecomputeNLAuxVars(pre::ModelValuesDbl& sol) {
    if (MPCD( IfWantNLOutput() )) {
      auto& xx = sol.GetVarValues()();
      if (xx.size()) {                    // solution available
        const auto& var_is_proper = MPCD( GetVarProperFlags() );
        xx = MPD( RecomputeAuxVars(xx, var_is_proper) );
      }
    }
  }

private:
  /// (Argument) variable visitor: mark var as proper
  std::function<void( int )> MarkVar_ = [this](int v){
    MPD( MarkAsResultVar(v) );       // no recursion
  };

  /// (Argument) variable visitor: mark var as proper,
  /// if it represents a logical expression
  std::function<void( int )> MarkVarIfLogical_
      = [this](int v){
          if (MPCD( IsInitExprLogical(v) ))
            MPD( MarkAsResultVar(v) ); // no recursion
  };

  /// NL conversion stage
  int stage_cvt2expr_ {-1};
};

}  // namespace mp

#endif // CONSTR_2_EXPR_H
