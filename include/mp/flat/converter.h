#ifndef CONVERTER_FLAT_H
#define CONVERTER_FLAT_H

#include <unordered_map>
#include <map>
#include <cmath>
#include <utility>
#include <cstdio>
#include <cassert>

#include "mp/env.h"
#include "mp/format.h"
#include "mp/solver-base.h"
#include "mp/suffix.h"
#include "mp/flat/converter_info.h"
#include "mp/flat/converter_model.h"
#include "mp/flat/convert_functional.h"
#include "mp/flat/constr_keeper.h"
#include "mp/flat/constr_std.h"
#include "mp/flat/expr_bounds.h"
#include "mp/flat/constr_prepro.h"
#include "mp/flat/constr_prop_down.h"
#include "mp/flat/expr_alg_inline.h"
#include "mp/flat/converter_multiobj.h"
#include "mp/flat/constr_2_expr.h"
#include "mp/flat/sol_check.h"
#include "mp/valcvt.h"
#include "mp/flat/redef/std/range_con.h"
#include "mp/flat/redef/conic/cones.h"
#include "mp/flat/redef/conic/qcones2qc.h"
#include "mp/flat/redef/SDP/sdp.h"
#include "mp/ampls-ccallbacks.h"
#include "mp/utils-misc.h"

namespace mp {

/// FlatConverter: preprocesses and manages flat constraints.
/// Such constraints might need to be converted to others, which is
/// handled by overloaded methods in derived classes.
/// @param Impl: the final CRTP class
/// @param ModelAPI: the solver's model API wrapper
/// @param FlatModel: internal representation of a flat model
template <class Impl, class ModelAPI,
          class FlatModel = FlatModel< > >
class FlatConverter :
                      public ConverterInfo,
                      public BasicFlatConverter,
                      public FlatModel,
                      public BoundComputations<Impl>,
                      public ConstraintPreprocessors<Impl>,
                      public ConstraintPropagatorsDown<Impl>,
                      public AlgebraicExpressionInliner<Impl>,
                      public MOManager<Impl>,
                      public Constraints2Expr<Impl>,
                      public SolutionChecker<Impl>,
                      public EnvKeeper
{
public:
  /// Class name
  static const char* GetTypeName() { return "FlatConverter"; }

  /// Construct with Env&
  FlatConverter(Env& e) : EnvKeeper(e), modelapi_(e) {
    this->AddConversionAction(
        [this](BasicFlatConverter& cvt) {
          MP_ASSERT_ALWAYS(&cvt == this, "Bad ptr");
          return this->InlineAlgSubexpr();
        }, 3099);     // 3099: before LinFuncCon/QuadFuncCon
  }

  /// Converter info
  const ConverterInfo* GetConverterInfo() const override
  { return this; }

  /// Trying to use 'Var' instead of bare 'int'
  using Var = typename FlatModel::Var;

  /// 'Invalid' var id
  static constexpr Var VoidVar() { return FlatModel::VoidVar(); }

  /// Array of variable Id's
  using VarArray = std::vector<int>;

protected:
  using ClassType = FlatConverter<Impl, ModelAPI, FlatModel>;
  using BaseConverter = BasicFlatConverter;
  using BaseFlatModel = FlatModel;


  //////////////////////////// CONVERTERS OF STANDARD MP ITEMS //////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////
public:
  /// Fix the resulting variable of a logical expression as true
  /// and propagate positive ctx.
  /// Currently this happens for all root-context logical constraints,
  /// i.e., we create an auxiliary variable which is later fixed to 1.
  void FixAsTrue(int resvar) {
    IncrementVarUsage(resvar);
    PropagateResultOfInitExpr(resvar, 1.0, 1.0, +Context());  // afterwards #201 #266
  }


public:
  /// Reverse propagate result variable of an expression
  void PropagateResultOfInitExpr(int var, Context ctx) {
    PropagateResultOfInitExpr(var, lb(var), ub(var), ctx);
  }

  /// Reverse propagate result variable of an expression
  void PropagateResultOfInitExpr(int var, double lb, double ub, Context ctx) {
    assert(!ctx.IsNone());
    bool tighterBounds = (lb > MPCD(lb(var)) || ub < MPCD(ub(var)));
    lb = std::max(lb, MPCD(lb(var)));
    ub = std::min(ub, MPCD(ub(var)));
    if (tighterBounds)
      NarrowVarBestBounds(var, lb, ub);       // cvt:pre:boundsbest
    if (HasInitExpression(var)) {
      const auto& ckid = GetInitExpression(var);
      const auto ctx_old = ckid.GetCK()->GetContext(ckid.GetIndex());
      if (tighterBounds
          || !ctx.IsSubsetOf(ctx_old)) {      // new context
        ckid.GetCK()->PropagateResult(*this, ckid.GetIndex(), lb, ub, ctx);
      }
    }
  }

  /// Propagate objective contexts
  void PropagateObjContexts() {
    auto& objs = MPD( get_objectives() );
    const auto objwgt = MPD( GetMOWeights() );
    if (GetEnv().multiobj())               // only in obj:multi mode
      assert(objs.size() == objwgt.size());
    for (size_t i=0; i<objs.size(); ++i) {
      auto isMax = obj::MAX==objs[i].obj_sense();
      if (GetEnv().multiobj() && objwgt[i]<0.0) {  // only in obj:multi
        isMax = !isMax;
        objs[i].set_sense_true(isMax ? obj::MAX : obj::MIN);
      }
      auto ctx = isMax ? Context::CTX_POS : Context::CTX_NEG;
      MPD( PropagateResult2LinTerms(objs[i].GetLinTerms(),
                          MPD( MinusInfty() ), MPD( Infty() ), ctx) );
      MPD( PropagateResult2QuadTerms(objs[i].GetQPTerms(),
                          MPD( MinusInfty() ), MPD( Infty() ), ctx) );
    }
  }

public:
  //////////////////////////////////// VISITOR ADAPTERS /////////////////////////////////////////
  /// These are called to transform expressions, either by FlatCvt itself,
  /// or when flattening NL model

  /// From an affine expression:
  /// Adds a result variable r and constraint r == expr
  int Convert2Var(AffineExpr ee) {
    if (ee.is_variable())
      return ee.get_representing_variable();
    if (ee.is_constant())
      return int( MakeFixedVar(ee.constant_term()) );
    ee.sort_terms();
    ee.shrink_to_fit();
    return AssignResultVar2Args(
            LinearFunctionalConstraint(std::move(ee)));
  }

  /// From a quadratic expression:
  /// Adds a result variable r and constraint r == expr
  int Convert2Var(QuadraticExpr ee) {
    if (ee.is_variable())
      return ee.get_representing_variable();
    if (ee.is_constant())
      return int(MakeFixedVar(ee.constant_term()));
    ee.sort_terms();
    ee.GetLinTerms().shrink_to_fit();
    ee.GetQPTerms().shrink_to_fit();
    if (ee.is_affine())
      return AssignResultVar2Args(
            LinearFunctionalConstraint(
              MoveOutAffineExpr(std::move(ee))));
    return AssignResultVar2Args(
        QuadraticFunctionalConstraint(std::move(ee)));
  }


  /// Take FuncConstraint with arguments
  ///
  /// Prefer this over AddConstraint() for mapped functional
  /// constraints.
  /// If the result of the function can be presolved or
  /// is known via map, return it.
  /// Otherwise, create a result variable and add the constraint.
  /// @return VarOrConst
  template <class FuncConstraint>
  typename FCC<Impl, FuncConstraint>::VarOrConst
  AssignResult2Args(FuncConstraint&& fc) {
    auto fcc = MakeFuncConstrConverter<Impl, FuncConstraint>(
          *this, std::forward<FuncConstraint>(fc));
		return fcc.Convert();
  }

  /// Same, but always return a variable
  /// (a fixed varible if the result is a constant).
  template <class FuncConstraint>
  typename FCC<Impl, FuncConstraint>::Var
  AssignResultVar2Args(FuncConstraint&& fc) {
    auto vc = AssignResult2Args(std::move(fc));
    if (vc.is_const())
      return int( MPD( MakeFixedVar(vc.get_const()) ) );
    return vc.get_var();
  }

  /// Same, but for constant result, still add the full
  /// expression.
  /// This is necessary for complementarity constraints
  /// in NL expression output where the constant part
  /// needs to be an actual expression.
  template <class FuncConstraint>
  typename FCC<Impl, FuncConstraint>::Var
  AssignResult2Args__FullExpression(FuncConstraint&& fc) {
    auto fcc = MakeFuncConstrConverter<Impl, FuncConstraint>(
        *this, std::forward<FuncConstraint>(fc));
    return fcc.Convert(true).get_var();
  }

	/// Typedef ConInfo; constraint location
	using ConInfo = AbstractConstraintLocation;

	/// Replace functional expression defining a given variable.
  /// After this, if the arguments are any new expressions,
  /// the original context should be propagated into \a res_var:
  /// PropagateResultOfInitExpr(res_var, ctx);
  template <class FuncConstraint>
  void RedefineVariable(int res_var, FuncConstraint&& fc) {
    assert( MPD( HasInitExpression(res_var) ) );
		auto ci_old = MPD( GetInitExpression(res_var) );
    fc.SetResultVar(res_var);
		// If this expression exists, use it.
		// TODO make sure any new context is re-converted
		// if necessary.
		auto i = MPD( MapFind(fc) );
    // TODO preprocess, try map again, and use the result.
    // assert(i<0);
    auto& ck = GET_CONSTRAINT_KEEPER( FuncConstraint );
    if (i<0)
      i = int( MPD( AddFunctionalConstraint(std::move(fc)) ) );
    else {       // #270 see redefvar_02.mod. Just add a new copy
      i = ck.AddConstraint(constr_depth_, std::move(fc));
      AutoLink(ck.SelectValueNodeRange(i));
    }
    ConInfo ci{&ck, i};
    ReplaceInitExpression(res_var, ci);
    MarkAsUsed(ci);          // Now manually #201 #266
    MarkAsBridged(ci_old);
  }


	/// Variables' reference counting ///////////////////////////////////
	/// Currently only for defined variables ////////////////////////////

  /// Use "+1" a variable.
  /// When changing from 0 to 1,
  /// mark as "used" if not redefined.
  void IncrementVarUsage(int v) {
    if (1==++VarUsageRef(v)) {
      if (HasInitExpression(v)) {
        auto& ci = GetInitExpression(v);
        if (ci.GetCK()->IsUnused(ci.GetIndex())
            && !ci.GetCK()->IsBridged(ci.GetIndex())) {
          MarkAsUsed(ci);
        }
      }
    }
    // Not catching reuse after redef:
    // @todo check new context in context propagation.
#ifdef CATCH_REUSE_AFTER_REDEF
    // If unused, no reformulation tried,
    // currently no repetition of reformulation cycle.
    // @todo could allow if redefined in CTX_MIX.
    MP_ASSERT_ALWAYS(!IsUnused(GetInitExpression(v))
                     || IsBridgingToBeConsidered(GetInitExpression(v)),
                     "An expression's redefinition\n"
                     "could be lost. Please contact\n"
                     "AMPL customer support.");
#endif
  }

	/// Unuse result variable.
  /// Actually this is to 'unuse' the init expression
  /// - might change naming.
  /// Throw if already not used.
  /// When changing from 1 to 0,
  /// mark "unused" if not already and not redefined
  void DecrementVarUsage(int v) {
    assert(VarUsageRef(v)>0);
    if (VarUsageRef(v)>0)          // in Release build
      if (! (--VarUsageRef(v))) {
        if (HasInitExpression(v)) {
          auto& ci = GetInitExpression(v);
          if (IsConActive(ci)) { // used && !bridged
            MarkAsUnused(ci);
          }
        }
      }
  }

  /// Count argument references.
  /// @todo currently called manually for objectives
  template <class ConObj>
  void CountArgRefs(const ConObj& con) {
    VisitArguments(con,
                   [this](int v) {
      IncrementVarUsage(v);
    });
  }

  /// Uncount argument references
  /// @todo currently called manually for objectives
  template <class ConObj>
  void UncountArgRefs(const ConObj& con) {
    VisitArguments(con,
                   [this](int v) {
      DecrementVarUsage(v);
    });
  }

  /// Mark unused defined vars for elimination.
	/// Normally should delete them.
  void EliminateUnusedDefinedVars() {
		for (auto i=num_vars(); i--; ) {
			if (HasInitExpression(i) &&
					! VarUsageRef(i)) {
        MPD( MarkVarAsEliminated(i) );
			}
		}
	}

  /// Var usage
  int VarUsage(int i) const {
    return (i<(int)refcnt_vars_.size())
               ? refcnt_vars_[i] : 0;
  }

protected:
	int& VarUsageRef(int i) {
		assert(i>=0 && i<num_vars());
		if ((size_t)i>=refcnt_vars_.size())
			refcnt_vars_.resize(
						std::max((size_t)num_vars(),
										 (size_t)(refcnt_vars_.size()*1.4)));
		return refcnt_vars_[i];
	}


public:
  /// Signal if the model needs a solve.
  /// This should be used in particular with MO emulation.
  /// @return true if the next solve should be executed
  /// and its results passed via ProcessSolveIterationSolution().
  bool PrepareNextSolveIteration(
      std::function<sol::Status(void)> get_stt, std::function<Solution(void)> get_sol) {
    if (MPCD( IsMOActive() ))
      return MPD( PrepareMOIteration(get_stt, get_sol) );
    return !(n_solve_iter_++);
  }

  /// Objective weights, adapted according to obj:multi:weight
  ArrayRef<double> GetObjWeightsAdapted() { return MPD( GetMOWeightsLegacy() ); }


protected:
  //////////////////////////// CUSTOM CONSTRAINTS CONVERSION ////////////////////////////
  ///
  //////////////////////////// THE CONVERSION LOOP: BREADTH-FIRST ///////////////////////
  void ConvertItems() {
    try {
      MPD( OutputModelInfo("AMPL MP initial flat model", 0, "flat0_"); );
      MPD( Convert2SDP(); );                   // Scan SDP before cones?
      MPD( Convert2Cones(); );                 // sweep before other conversions
      MP_DISPATCH( ConvertAllConstraints() );
      // MP_DISPATCH( PreprocessIntermediate() );     // preprocess after each level
      constr_depth_ = 1;  // Workaround. TODO have maps as special constraints
			MP_DISPATCH( ConvertMaps() );
      MP_DISPATCH( PreprocessFlatFinal() );           // final flat model prepro
      MP_DISPATCH( ConsiderEmulatingMultiobj() );     // Before NL conversion
      if constexpr (IfAcceptingNLOutput()) {
        if (IfWantNLOutput()) {
          MPD( Convert2NL() );                        // Possibly for emulated objectives
          MPD( PreprocessNLFinal() );
        }
      }
    } catch (const ConstraintConversionFailure& cff) {
      MP_RAISE(cff.message());
    }
  }

  /// Print model info and/or output such suffixes
  void OutputModelInfo(const char* header, bool aux_vars,
      const char* suf_prefix) {
    if (GetEnv().verbose_mode() || GetEnv().debug_mode()) {
      MPD( CreateFlatModelInfo(GetModelAPI()) );
      if (GetEnv().debug_mode())
        ReportModelInfoSuffixes(
            *MPCD( GetModelInfo() ), suf_prefix, suf_get_set_);
      if (GetEnv().verbose_mode()) {
        const auto* fmi = MPCD( GetModelInfo() );
        if (!aux_vars) {                      // 1st output
          fmt::print("\n");
            modelinfo_flat0_vars_ = fmi->GetVarInfo();
          modelinfo_flat0_objs_ = fmi->GetObjInfo();
          modelinfo_flat0_cons_ = fmi->GetConstraintTypes();
          PrintModelInfo(
              *fmi, header, aux_vars);
        } else {                              // 2nd output
          if (fmi->GetVarInfo() != modelinfo_flat0_vars_
              || fmi->GetObjInfo() != modelinfo_flat0_objs_
              || fmi->GetConstraintTypes() != modelinfo_flat0_cons_)
            PrintModelInfo(
                *fmi, header, aux_vars);
          else
            fmt::print("AMPL MP did not modify the model.\n\n");
          fmt::print("\n");
          std::fflush(stdout);
        }
      }
    }
  }

  void OpenGraphExporter() {
    if (graph_export_file().size()) {
      if (!GetModel().GetFileAppender().Open(
            graph_export_file().c_str(), true))
        MP_RAISE("Failed to open the graph export file.");
    }
  }

  /// Offload the SDP logic to a functor
  void Convert2SDP() {
    sdp_cvt_.Run();
  }

  /// Offload the conic logic to a functor
  void Convert2Cones() {
    conic_cvt_.Run();
  }

  /// Can be called from ConvertMaps()
  void ConvertAllConstraints() {
    GetModel().ConvertAllConstraints(*this);
  }

  /// Inline algebraic subexpr.
  /// @return true iff anything changed.
  bool InlineAlgSubexpr() {
    auto preu = MPCD( IfPreproUnnest() );
    return MPD(
        ConsiderInliningAlgExpr(preu & 2, preu & 4) );
  }

  /// Default map conversions. Currently empty
  void ConvertMaps() { }

  /// Acceptance level for all constr/expr,
  /// if provided (>=0)
  int AcceptanceLevelCommon() const { return options_.accAll_; }

public:
  /// Option to actually use expressions if available
  bool IfWantNLOutput() const { return options_.accExpr_==1; }

  /// Whether solver CAN accept expressions
  static constexpr bool IfAcceptingNLOutput()
  { return
        ExpressionAcceptanceLevel::AcceptedButNotRecommended
        == ModelAPI::ExpressionInterfaceAcceptanceLevel()
        ||
           ExpressionAcceptanceLevel::Recommended
               == ModelAPI::ExpressionInterfaceAcceptanceLevel(); }

protected:
  /// Finish exporting the reformulation graph
  void CloseGraphExporter() {
    value_presolver_.FinishExportingLinkEntries();
    GetModel().GetFileAppender().Close();
  }

  //////////////////////// WHOLE-MODEL PREPROCESSING /////////////////////////
  void PreprocessIntermediate() { }
  void PreprocessFlatFinal() { }
  void PreprocessNLFinal() { }


  //////////////////////////// CONSTRAINT PROPAGATORS ///////////////////////////////////

  /// Allow FCC to access Preprocess methods
  template <class Impl1, class Converter, class Constraint>
  friend class BasicFCC;

  //////////////////////////// SPECIFIC CONSTRAINT RESULT-TO-ARGUMENTS PROPAGATORS //////
  /// Currently we should propagate to all arguments, be it always the CTX_MIX.

  /// Allow ConstraintKeeper to PropagateResult(), use GetModelAPI() etc
  template <class , class , class >
  friend class ConstraintKeeper;


  //////////////////////////// CUSTOM CONSTRAINTS CONVERSION ////////////////////////////
  ///
public: // for ConstraintKeeper
  /// RunConversion() of a constraint in the flat phase:
  /// Assume mixed context if not set.
  /// @return Context used for redefinition
  /// @note Do not use directly. Call via
  ///   ConstraintKeeper.ConvertConstraint().
  template <class Constraint>
  Context RunConversion(const Constraint& con, int i, int depth) {
    assert(
        !GET_CONSTRAINT_KEEPER(Constraint).IsRedundant(i));
    constr_depth_ = depth+1;
    if (con.UsesContext())              // If context relevant,
      if (con.GetContext().IsNone())    // ensure we have context, mixed if none
        con.SetContext(Context::CTX_MIX);
    pre::AutoLinkScope<Impl> auto_link_scope{
      *static_cast<Impl*>(this),
      GET_CONSTRAINT_KEEPER(Constraint).SelectValueNodeRange(i)
    };
    return MP_DISPATCH(Convert(con, i));
  }

  /// Query the number of addable constraints of type.
  template <class Con>
  int GetNumberOfAddable(Con* ) const {
    return GET_CONST_CONSTRAINT_KEEPER(Con).GetNumberOfAddable();
  }

  /// Query if the constraint type
  /// can be converted.
  /// This method should not be redefined;
  /// specialize IfHasCvt_impl instead.
  template <class Constraint>
  bool IfHasConversion(const Constraint* c) {
    return MPD( IfHasCvt_impl(c) );
  }

  /// Generic query if a constraint type can be converted.
  /// Specialize this method, but normally it's specialized
  /// by INSTALL_CONSTRAINT_CONVERTER.
  template <class Constraint>
  bool IfHasCvt_impl(const Constraint* ) {
    return false;
  }

  /// Query if the specific item of the constraint
  /// needs to be converted,
  /// despite being accepted by the ModelAPI.
  /// For example, Gurobi only accepts Pow with non-negative
  /// argument.
  /// This method should not be redefined;
  /// specialize IfNeedsCvt_impl instead.
  template <class Constraint>
  bool IfNeedsConversion(const Constraint& con, int i) {
    return MPD( IfNeedsCvt_impl(con, i) );
  }

  /// Generic query if a constraint needs to be converted,
  /// despite being accepted by the ModelAPI.
  /// Specialize this method, or even
  /// ConstraintConverter::IfNeedsConversion
  /// (see class PowConstExponentConverter_MIP).
  template <class Constraint>
  bool IfNeedsCvt_impl(const Constraint& , int ) {
    return false;
  }

  /// Query if the specific item of the constraint
  /// needs to be skipped in regular conversion,
  /// despite being not accepted by the ModelAPI.
  /// For example, some conditional comparisons.
  /// This method should not be redefined;
  /// specialize IfNeedsCvt_impl instead.
  template <class Constraint>
  bool IfDelayConversion(const Constraint& con, int i) {
    return MPD( IfDelayCvt_impl(con, i) );
  }

  /// Generic query if a constraint needs to be
  /// skipped from regular conversion,
  /// despite being not accepted by the ModelAPI.
  /// Specialize this method, or even
  /// ConstraintConverter::IfDelayConversion
  /// (see class CondEQConverter_MIP).
  template <class Constraint>
  bool IfDelayCvt_impl(const Constraint& , int ) {
    return false;
  }


  ////////////////// Constraint/expression accpetance ////////////////

  /// Check whether ModelAPI (not user) accept and recommend the constraint
  template <class Constraint>
  bool ModelAPIOk() const {
    return ModelAPIAcceptsAndRecommends((const Constraint*)0);
  }

  /// Whether user (ONLY) recommends either the constraint or corr. expr,
  /// for the expression also check if expression output is desired
  template <class Constraint>
  bool UserAcceptsConOrExpr() const {
    return
        MPCD( UserAcceptsAndRecommends((const Constraint*)0) )
           || (MPCD( IfWantNLOutput() )
               && MPCD( template UserAcceptsExprForCon<Constraint>()));
  }


  /// Check whether ModelAPI and user accept and recommend the constraint
  template <class Constraint>
  bool UserAcceptsAndRecommends(const Constraint* pcon) const {
    return ConstraintAcceptanceLevel::Recommended ==
           GetConstraintAcceptance_USER(pcon);
  }

  /// Check whether ModelAPI accept and recommend the constraint
  template <class Constraint>
  bool ModelAPIAcceptsAndRecommends(const Constraint* pcon) const {
    return ConstraintAcceptanceLevel::Recommended ==
           GetConstraintAcceptance_DEFAULT(pcon);
  }

  /// Check whether user accepts and recommends
  /// the expression corresponding to constraint
  template <class Constraint>
  bool UserAcceptsExprForCon() const {
    return UserAcceptsAndRecommendsEXPR((
        const ExprWrapper<Constraint>*)0);
  }

  /// Check whether ModelAPI and user accept and recommend the expression
  template <class Expression>
  bool UserAcceptsAndRecommendsEXPR(const Expression* pcon) const {
    return ExpressionAcceptanceLevel::Recommended ==
           GetConstraintAcceptanceEXPR_USER(pcon);
  }

  /// Check whether ModelAPI accept and recommend the expression
  template <class Expression>
  bool ModelAPIAcceptsAndRecommendsEXPR(const Expression* pcon) const {
    return ExpressionAcceptanceLevel::Recommended ==
           GetConstraintAcceptanceEXPR_DEFAULT(pcon);
  }

  /// Query if a constraint type
  /// is natively accepted by the solver (and user setting).
  /// The parameter is only needed for type.
  template <class Con>
  ConstraintAcceptanceLevel GetConstraintAcceptance_USER(Con* ) const {
    return GET_CONST_CONSTRAINT_KEEPER(Con).GetChosenAcceptanceLevel();
  }

  /// Query if an expression type
  /// is natively accepted by the solver (and user setting).
  /// The parameter is only needed for type.
  template <class Con>
  ExpressionAcceptanceLevel GetConstraintAcceptanceEXPR_USER(
      const ExprWrapper< Con >* ) const {
    return GET_CONST_CONSTRAINT_KEEPER(Con).GetChosenAcceptanceLevelEXPR();
  }

  /// Query if a constraint type
  /// is natively accepted by the solver.
  /// The parameter is only needed for type.
  template <class Con>
  ConstraintAcceptanceLevel GetConstraintAcceptance_DEFAULT(Con* ) const {
    return GET_CONST_CONSTRAINT_KEEPER(Con).GetModelAPIAcceptance();
  }

  /// Query if an expression type
  /// is natively accepted by the solver.
  /// The parameter is only needed for type.
  template <class Con>
  ExpressionAcceptanceLevel GetConstraintAcceptanceEXPR_DEFAULT(
      const ExprWrapper< Con >* ) const {
    return GET_CONST_CONSTRAINT_KEEPER(Con).GetModelAPIAcceptanceEXPR();
  }


  /// Generic adapter for old non-bridged Convert() methods
  ///
  /// New way is to use the \a i parameter for bridging
  template <class Constraint>
  Context Convert(const Constraint& con, int ) {
    return MPD( Convert(con) );
  }

  /// By default, we complain about someone trying to
  /// convert an unknown constraint
  template <class Constraint>
  Context Convert(const Constraint& ) {
    MP_RAISE(
          std::string("Constraint type '") +
            Constraint::GetTypeName() +
            "' is neither accepted by '" +
            ModelAPI::GetTypeName() +
            "', nor is conversion implemented");
    return Context::CTX_NONE;
  }

  //////////////////////////// SOME SPECIFIC CONSTRAINT CONVERTERS
  /// ///////////////////////////////////// ///////////////////////////

  /// If backend does not like LFC, we redefine it here
  Context Convert(const LinearFunctionalConstraint& ldc) {
    MPD( AddConstraint(ldc.to_linear_constraint()) );
    return Context::CTX_MIX;
  }
  /// Say we can (for acc:_all=0)
  bool IfHasCvt_impl(const LinearFunctionalConstraint* ) {
    return true;
  }

  /// If backend does not like QFC, we redefine it
  Context Convert(const QuadraticFunctionalConstraint& qdc) {
    return qdc.AddQuadraticConstraint(*(Impl*)this);
  }
  /// Say we can
  bool IfHasCvt_impl(const QuadraticFunctionalConstraint* ) {
    return true;
  }


public:
  /// Add objective.
  ///
  /// Currently handling quadratic objectives.
  /// Its quadratic terms will be empty for linear objectives.
  ///
  /// Linking NL objectives straight into solver's objectives.
  /// If any conversions are performed, need to have intermediate nodes,
  /// as for constraints
  pre::NodeRange AddObjective(QuadraticObjective&& qo) {
    CountArgRefs(qo);     // no "used" flag for objectives
    GetModel().AddObjective( std::move(qo) );
    /// Temporarily removing AutoLinking for objectives
    // return AutoLink( GetObjValueNode().Add() );
    return GetObjValueNode().Select(-1);
  }

  /// ADD STATIC CONSTRAINT.
  ///
  /// Does not propagate result
  /// (use AddConstraint_AS_ROOT() otherwise).
  ///
  /// Counts argument references.
  ///
  /// This method is enough
  /// (instead of the _AS_ROOT() version)
  /// if the arguments already have contexts.
  ///
  /// Use only for non-mapped constraints. For functional constraints
  /// stored __WITH_MAP, use AssignResult(Var)2Args().
  /// Non-functional constraints cannot be unified currently.
  /// Takes ownership.
  ///
  /// @return Node reference for the stored constraint
  template <class Constraint>
  pre::NodeRange AddConstraint(Constraint con) {
    assert(!con.HasResultVar());
    if (MPD( PreprocessStaticConstraint(con) ))
      return {};  // we should not need the presolver nodes
    auto node_range =
        AddConstraintAndTryNoteResultVariable( std::move(con) );
    auto& ck = GET_CONSTRAINT_KEEPER( Constraint );
    ConInfo ci{&ck, int(node_range)};
    MarkAsUsed(ci);      // this also counts arg refs #266
    return AutoLink( node_range );
  }

  /// ADD STATIC CONSTRAINT and propagate root-ness
  /// (use AddConstraint() otherwise).
  ///
  /// Use only for non-mapped constraints. For functional constraints
  /// stored __WITH_MAP, use AssignResult(Var)2Args().
  /// Takes ownership.
  /// @return Node reference for the stored constraint
  template <class Constraint>
  pre::NodeRange AddConstraint_AS_ROOT(Constraint con) {
    auto nr = AddConstraint( std::move(con) );
    MPD( PropagateResult(             // after AddConstraint() #201 #266
           GetConstraint<Constraint>(int(nr))) );
    return nr;
  }

  /// ADD FUNCTIONAL CONSTRAINT.
  ///
  /// Do not use directly. For functional constraints
  /// stored __WITH_MAP, use AssignResult(Var)2Args().
  /// Takes ownership.
  ///
  /// @note Does not propagate result
  ///   (use PropagateResult()).
  ///
  /// @return Node reference for the stored constraint
  template <class Constraint>
  pre::NodeRange AddFunctionalConstraint(Constraint con) {
    assert(con.HasResultVar());
    auto node_range =
        AddConstraintAndTryNoteResultVariable( std::move(con) );
    return AutoLink( node_range );
  }

  /// Retrieve constraint of specified type at location \a ci.
  template <class Constraint>
	const Constraint& GetConstraint(const ConInfo& ci) const {
		assert(MPCD(template IsConInfoType<Constraint>(ci) ));
		return GET_CONST_CONSTRAINT_KEEPER(Constraint).
				GetConstraint(ci.GetIndex());
  }

	/// Retrieve constraint of specified type at index \a i.
	template <class Constraint>
	const Constraint& GetConstraint(int i) const {
    return
        GET_CONST_CONSTRAINT_KEEPER(Constraint).GetConstraint(i);
	}

  /// Mark constraint as reformulated
  void MarkAsBridged(const ConInfo& ci) {
		ci.GetCK()->MarkAsBridged(ci.GetIndex());
	}

  /// Mark constraint as unused
  void MarkAsUnused(const ConInfo& ci) {
    ci.GetCK()->MarkAsUnused(ci.GetIndex());
  }

  /// Mark constraint as unused.
  /// Do not propagate to arguments.
  void MarkAsUnused_ThisOnly(const ConInfo& ci) {
    ci.GetCK()->MarkAsUnused_ThisOnly(ci.GetIndex());
  }

  /// Mark constraint as used
  void MarkAsUsed(const ConInfo& ci) {
    ci.GetCK()->MarkAsUsed(ci.GetIndex());
  }

  /// Is bridging of constraint \a i
  /// to be considered yet?
  bool IsBridgingToBeConsidered(const ConInfo& ci) const {
    return ci.GetCK()->IsBridgingToBeConsidered(ci.GetIndex());
  }


  /// Is constraint reformulated?
  bool IsBridged(const ConInfo& ci) const {
    return ci.GetCK()->IsBridged(ci.GetIndex());
  }

  /// Is constraint unused?
  bool IsUnused(const ConInfo& ci) const {
    return ci.GetCK()->IsUnused(ci.GetIndex());
  }


protected:
  USE_BASE_MAP_FINDERS( BaseConverter )

  template <class Constraint>
  pre::NodeRange AddConstraintAndTryNoteResultVariable(Constraint&& con) {
    const auto resvar = con.GetResultVar();
    auto& ck = GET_CONSTRAINT_KEEPER( Constraint );
    auto i = ck.AddConstraint(constr_depth_, std::move(con));
    ConInfo ci{&ck, i};
    if (resvar>=0)
      AddInitExpression(resvar, ci);
    /// Can also cache non-functional constraints,
    /// but then implement checking before
    if (! MP_DISPATCH( MapInsert(
                         MPD(template GetConstraint<Constraint>(i)), i ) ))
      MP_RAISE("Trying to MapInsert() duplicated constraint: " +
                             ck.GetDescription());
    pre::NameChunk old_nc {};
    if (name_chunk_)
      old_nc = ck.GetValueNode().SetNameChunk(name_chunk_);
    auto result = ck.SelectValueNodeRange(i);
    if (name_chunk_)
      ck.GetValueNode().SetNameChunk(old_nc);
    return result;
  }


public:
	/// Select value node \a i for constraint type \a Con.
	template <class Constraint>
	pre::NodeRange SelectValueNode(int i) {
		auto& ck = GET_CONSTRAINT_KEEPER( Constraint );
		return ck.SelectValueNodeRange(i);
	}

  /// Handle start of model input
  void StartModelInput() {
    MPD( OpenGraphExporter() );
  }

  /// Handle end of model input
  void FinishModelInput() {
    MPD( ConvertModel() );
    if (relax())
      GetModel().RelaxIntegrality();
    EliminateUnusedDefinedVars();       // Until we have proper var deletion
    CheckLinearCons();
    PresolveNames();
    GetModel().PrepareVariables();      // Should this be Converter's task?
    MPD( OutputModelInfo("AMPL MP final model", 1, "flat1_"); );
    GetModel().PushModelTo(GetModelAPI());
    MPD( CloseGraphExporter() );
    if (value_presolver_.GetExport())
      assert( value_presolver_.AllEntriesExported() );
//  Printing always.  if (GetEnv().verbose_mode())
    GetEnv().PrintWarnings();
  }

  /// Check linear constraints.
  /// Solvers complain about close-to-0 coefficients,
  /// so fail on this in the debug build.
  void CheckLinearCons() {
#ifndef NDEBUG
    CheckLinearConType< AlgConRange >();
    CheckLinearConType< AlgConRhs<-1> >();
    CheckLinearConType< AlgConRhs<0> >();
    CheckLinearConType< AlgConRhs<1> >();
#endif
  }

  template <class Rhs>
  void CheckLinearConType() {
    using LinConType = AlgebraicConstraint< LinTerms, Rhs >;
    auto& ck = GetConstraintKeeper(
          (LinConType*)nullptr );
    ck.ForEachActive( [](const LinConType& lc, int ){
      const auto& lt = lc.GetBody();
      for (auto i=lt.size(); i--; ) {
        assert(0.0 != std::fabs(lt.coef(i)) &&
            "Most solvers don't like near-zero coefficients");
      }
      return false;     // don't delete
    } );
  }

  /// Reset local item counters, for redefinitions
  void ResetLocalCounters() {
    ConstraintManager::ResetLocalCounters();
    GetVarValueNode().ResetLocalCounter();
  }

  /// Presolve item names
  void PresolveNames() {
    if (var_names_.size()) {
      /// Check that constr / obj names are present too?
      GetValuePresolver().CleanUpNameNodes();
      // They are at top level of the reformulation tree
      TransferNames2Node((SOS1Constraint*)nullptr);
      TransferNames2Node((SOS2Constraint*)nullptr);
      auto vm = GetValuePresolver().
          PresolveNames({
                          {var_names_},
                          {con_names_},
                          {obj_names_}
                        });
      const auto& vcs = vm.GetVarValues()();    // vars
      std::vector<std::string> vs(vcs.begin(), vcs.end());
      BaseFlatModel::AddVarNames(vs);
      const auto& ocs = vm.GetObjValues()();    // objs
      auto& obj = BaseFlatModel::get_objectives();
      assert(obj.size() == ocs.size());
      for (auto i=obj.size(); i--;)
        obj[i].set_name(ocs[i]);
      ConstraintManager::CopyNamesFromValueNodes();  // cons
    }
  }

  template <class Con>
  void TransferNames2Node(Con* pcon) {
    auto& ck = GetConstraintKeeper(pcon);
    ck.CopyNames2ValueNodes();
  }

  /// Fill model traits for license check.
  /// To be called after ConvertModel().
  /// KEEP THIS UP2DATE.
  void FillModelTraits(AMPLS_ModelTraits& mt) {
    const auto& fmi = GetModelAPI().GetFlatModelInfo();
    mt.n_vars = num_vars();
    mt.n_quad_con =
        fmi->GetNumberOfConstraints(typeid(QuadConRange))
        + fmi->GetNumberOfConstraints(typeid(QuadConGE))
        + fmi->GetNumberOfConstraints(typeid(QuadConEQ))
        + fmi->GetNumberOfConstraints(typeid(QuadConLE));
    mt.n_conic_con =
        fmi->GetNumberOfConstraints(typeid(QuadraticConeConstraint))
        + fmi->GetNumberOfConstraints(typeid(RotatedQuadraticConeConstraint))
        + fmi->GetNumberOfConstraints(typeid(ExponentialConeConstraint))
        + fmi->GetNumberOfConstraints(typeid(PowerConeConstraint))
        + fmi->GetNumberOfConstraints(typeid(GeometricConeConstraint));
    mt.n_alg_con =
        fmi->GetNumberOfConstraints(typeid(LinConRange))
        + fmi->GetNumberOfConstraints(typeid(LinConGE))
        + fmi->GetNumberOfConstraints(typeid(LinConEQ))
        + fmi->GetNumberOfConstraints(typeid(LinConLE))
        + mt.n_quad_con
        + fmi->GetNumberOfConstraints(typeid(ComplementarityLinear))
                   + fmi->GetNumberOfConstraints(typeid(ComplementarityQuadratic))
                   + fmi->GetNumberOfConstraints(typeid(NLComplementarity))
        ;
    mt.n_log_con =
        fmi->GetNumberOfConstraints(typeid(AndConstraint))
        + fmi->GetNumberOfConstraints(typeid(OrConstraint))
        + fmi->GetNumberOfConstraints(typeid(MaxConstraint))
        + fmi->GetNumberOfConstraints(typeid(MinConstraint))
        + fmi->GetNumberOfConstraints(typeid(IndicatorConstraintLinGE))
        + fmi->GetNumberOfConstraints(typeid(IndicatorConstraintLinEQ))
        + fmi->GetNumberOfConstraints(typeid(IndicatorConstraintLinLE))
        + fmi->GetNumberOfConstraints(typeid(IndicatorConstraintQuadGE))
        + fmi->GetNumberOfConstraints(typeid(IndicatorConstraintQuadEQ))
        + fmi->GetNumberOfConstraints(typeid(IndicatorConstraintQuadLE));
  }


protected:
  void ConvertModel() {
    MPD( PrepareConversion() );
    MPD( ConvertItems() );
    MPD( WindupConversion() );
  }

  void PrepareConversion() {
  }

  void WindupConversion() {
  }


  //////////////////////////// UTILITIES /////////////////////////////////
  ///
public:
  /// Expose abstract Backend
  const ModelAPI& GetModelAPI() const { return modelapi_; }
  ModelAPI& GetModelAPI() { return modelapi_; }

  /// Expose ValuePresolver
  const pre::ValuePresolver& GetValuePresolver() const { return value_presolver_; }
  pre::ValuePresolver& GetValuePresolver() { return value_presolver_; }


private:
  std::unordered_map<double, int> map_fixed_vars_;


public:
  //////////////////////////// CREATE OR FIND A FIXED VARIABLE //////////////////////////////
  pre::NodeRange MakeFixedVar(double value) {
    auto it = map_fixed_vars_.find(value);
    if (map_fixed_vars_.end()!=it)
      return AutoLink( GetVarValueNode().Select( it->second ) );
    auto v = MPD( DoAddVar(value, value) );
    map_fixed_vars_[value] = (int)v;
    return GetVarValueNode().Select( (int)v );  // no autolink, done in DoAddVar()
  }

  /// Create or find a fixed variable
  pre::NodeRange AddVar(double lb, double ub, var::Type type = var::CONTINUOUS) {
    if (lb!=ub)
      return DoAddVar(lb, ub, type);
    return MakeFixedVar(lb);
  }

  /// Add several variables once.
  /// @note this is only to be called once for the original vars.
  /// @return value node range for them
  pre::NodeRange AddVars(const typename BaseFlatModel::VarBndVec& lbs,
               const typename BaseFlatModel::VarBndVec& ubs,
               const typename BaseFlatModel::VarTypeVec& types) {
    assert(0==BaseFlatModel::num_vars());                     // allow this only once
    BaseFlatModel::AddVars__basic(lbs, ubs, types);
    return AutoLink( GetVarValueNode().Add( lbs.size() ) );
  }

  void AddVarNames(std::vector<std::string> names) {
    var_names_ = std::move(names);
  }
  void AddConNames(std::vector<std::string> names) {
    con_names_ = std::move(names);
  }
  void AddObjNames(std::vector<std::string> names) {
    obj_names_ = std::move(names);
  }
  /// Reuse ValuePresolver's target nodes for all variables
  pre::ValueNode& GetVarValueNode()
  { return GetValuePresolver().GetTargetNodes().GetVarValues().MakeSingleKey(); }

  /// Constraint type's Value Node
  template <class Constraint>
  pre::ValueNode& GetValueNode(Constraint*)
  { return GET_CONSTRAINT_KEEPER(Constraint).GetValueNode(); }

  /// Reuse ValuePresolver's source nodes for all objectives
  pre::ValueNode& GetObjValueSourceNode()
  { return GetValuePresolver().GetSourceNodes().GetObjValues().MakeSingleKey(); }

  /// Reuse ValuePresolver's target nodes for all objectives
  pre::ValueNode& GetObjValueNode()
  { return GetValuePresolver().GetTargetNodes().GetObjValues().MakeSingleKey(); }


public:
  /// Shortcut num_vars()
  int num_vars() const { return MPCD(GetModel()).num_vars(); }
  /// Shortcut is_var_original()
  int is_var_original(int i) const
  { return MPCD(GetModel()).is_var_original(i); }
  /// Shortcut lb(var)
  double lb(int var) const { return this->GetModel().lb(var); }
  /// Shortcut ub(var)
  double ub(int var) const { return this->GetModel().ub(var); }
  /// Shortcut lb_hard(var)
  double lb_hard(int var) const { return this->GetModel().lb_hard(var); }
  /// Shortcut ub_hard(var)
  double ub_hard(int var) const { return this->GetModel().ub_hard(var); }
  /// lb_array().
  /// @todo best-known bounds currently. ?
  template <class VarArray>
  double lb_array(const VarArray& va) const
  { return this->GetModel().lb_array(va); }
  /// ub_array()
  template <class VarArray>
  double ub_array(const VarArray& va) const
  { return this->GetModel().ub_array(va); }
  /// Does the variable have stronger solver-submitted bounds
  /// than its init expression?
  /// We might extend this to automatically recompute bounds from bottom up.
  /// @note this was done when creating the functional constraint,
  ///   but we might obtain stronger new bounds,
  ///   such as in redefinition of complementarity.
  /// @note Considers option *cvt:pre:boundsbest*.
  /// @note Should not be used directly,
  ///   use CanBeEliminated().
  bool IfSubmittedVarBoundsStrongerThanInitExpr(int res_var) const {
    if (MPCD( HasInitExpression(res_var) )) {
      if (lb(res_var)>MPCD( MinusInfty() )
          || ub(res_var)<MPCD( Infty() )) {
        const auto& cloc = MPCD( GetInitExpression(res_var) );
        PreprocessInfoStd preinfo;
        cloc.GetCK()->PreprocessConstraint(cloc.GetIndex(), preinfo);
        if (lb_hard(res_var) > preinfo.lb()     // If some hard bound better:
            || ub_hard(res_var) < preinfo.ub())
          return true;
        return (lb(res_var) > preinfo.lb() // If some best-known bound better:
                || ub(res_var) < preinfo.ub()) ?
                   (GetModel().if_submit_best_known_bounds() || is_fixed(res_var)) :
                   is_fixed(res_var);         // Only if fixed by default
      }
    }
    return false;
  }
  /// Set lb(var)
  void set_var_lb(int var, double lb) { this->GetModel().set_lb(var, lb); }
  /// Set ub(var)
  void set_var_ub(int var, double ub) { this->GetModel().set_ub(var, ub); }
  /// Set lb(var), propagate context if functional result
  void set_var_lb_context(int var, double lb, Context ctx) {
    set_var_lb(var, lb);   // Because PropResult() only hint bounds
    PropagateResultOfInitExpr(var, lb, ub(var), ctx);
  }
  /// Set ub(var), propagate context
  void set_var_ub_context(int var, double ub, Context ctx) {
    set_var_ub(var, ub);   // Because PropResult() only hint bounds
    PropagateResultOfInitExpr(var, lb(var), ub, ctx);
  }
  /// Set bounds(var), propagate context
  void set_var_bounds_context(int var, double lb, double ub, Context ctx) {
    NarrowVarBounds(var, lb, ub);       // Because PropResult() only hint bounds
    PropagateResultOfInitExpr(var, lb, ub, ctx);
  }

  /// Narrow "best-known" variable domain range
  /// @todo we could automatically decrement var usage
  ///   in appropriate context when the variable is fixed
  ///   at one of the original/upwards-implied bounds.
  /// Then remove manual DecrementVarUsage()'s
  void NarrowVarBestBounds(int var, double lb, double ub) {
    auto& m = GetModel();
    m.set_best_lb(var, lb);
    m.set_best_ub(var, ub);
    if (m.lb(var)>m.ub(var))
      CheckVarConDomain(m.lb(var), m.ub(var), "_svar", var, true);
  }

  /// Narrow variable domain range
  /// @todo we could automatically decrement var usage
  ///   in appropriate context when the variable is fixed
  ///   at one of the original/upwards-implied bounds.
  /// Then remove manual DecrementVarUsage()'s
  void NarrowVarBounds(int var, double lb, double ub) {
    auto& m = GetModel();
    m.set_lb(var, lb);
    m.set_ub(var, ub);
    if (m.lb(var)>m.ub(var))
      CheckVarConDomain(m.lb(var), m.ub(var), "_svar", var);
  }

  /// Check var/con domain
  bool CheckVarConDomain(
      double lb, double ub, const char* kind, int i, bool bestb=false) {
    if (lb>ub
        && lb-ub > MPCD( model_feas_tol() )
        && lb-ub
               > std::max(std::abs(lb), std::abs(ub))
                     * MPCD( model_feas_tol_rel() )) {
      GetEnv().AddWarning(
          std::string(kind) + (bestb ? " best-known" : "") + " bounds",
          fmt::format("Bounds [{:.17}, {:.17}]\nof {}[{}] "
                      "contradict pre:eps and pre:epsrel.\n"
                      "Model can be infeasible",
                      lb, ub, kind, i+1).c_str());
      return true;
    }
    return true;
  }

  /// var_type()
  var::Type var_type(int var) const { return this->GetModel().var_type(var); }
  /// is_fixed(), uses best-known bounds
  bool is_fixed(int var) const { return this->GetModel().is_fixed(var); }
  /// fixed_value()
  double fixed_value(int var) const
  { assert(is_fixed(var)); return this->GetModel().fixed_value(var); }

  /// MakeComplementVar()
  int MakeComplementVar(int bvar) {
    if ( !(lb_hard(bvar)==0.0 && ub_hard(bvar)==1.0) ) {
      // Should be hard-fixed at 0 or 1
      MP_ASSERT_ALWAYS( ((!lb_hard(bvar) && !ub_hard(bvar))
                        || (1.0==lb_hard(bvar) && 1.0==ub_hard(bvar))),
                "Asked to complement variable with bounds "
                    + std::to_string(lb_hard(bvar))
                    + ".." + std::to_string(ub_hard(bvar)));
    }
    /// Algebraic way: AffineExpr ae({{-1.0}, {bvar}}, 1.0);
    /// return MP_DISPATCH( Convert2Var(std::move(ae)) );
    return
        AssignResultVar2Args( NotConstraint{{bvar}} );
  }

  /// Add vector of variables. Type: var::CONTINUOUS by default
  /// @return vector of the Ids of the new vars
  std::vector<int> AddVars_returnIds(std::size_t nvars,
                           double lb=MinusInfty(), double ub=Infty(),
                           var::Type type = var::CONTINUOUS) {
    std::vector<int> newVars(nvars);
    for (std::size_t  i=0; i<nvars; ++i)
      newVars[i] = int( AddVar(lb, ub, type) );
    return newVars;
  }

  bool is_var_integer(int var) const
  { return MPCD( GetModel() ).is_integer_var(var); }


private:
  std::vector<ConInfo> var_info_;

  /// If non-0, use this chunk when adding cons
  const char* name_chunk_ {nullptr};


public:
  /// Pass 0 to use default for constraints
  void SetNameChunk(const char* nc) { name_chunk_ = nc; }


protected:
  /// Add variable. Type: var::CONTINUOUS by default
  pre::NodeRange DoAddVar(double lb=MinusInfty(), double ub=Infty(),
             var::Type type = var::CONTINUOUS) {
    int v = GetModel().AddVar__basic(lb, ub, type);
    return AutoLink( GetVarValueNode().Select( v ) );
  }

  /// Add init expr for \a var
  void AddInitExpression(int var, const ConInfo& vi) {
    if (var_info_.size() <= (size_t)var)
      var_info_.resize(((size_t)(var+1)*2));
    var_info_[var] = vi;
  }

  /// Replace init expression for \a var
  void ReplaceInitExpression(int var, const ConInfo& vi) {
    var_info_.at(var) = vi;
  }


public:
  /// Variable has an init expr?
  bool HasInitExpression(int var) const {
    return int(var_info_.size())>var && var_info_[var].HasId();
  }

  /// Get the init expr
  const ConInfo& GetInitExpression(int var) const {
    assert(HasInitExpression(var));
		return var_info_.at(var);
  }

  /// The variable has an init expr,
  /// and the init expr is a logical constraint?
  bool IsInitExprLogical(int var) const {
    if (!MPCD(HasInitExpression(var)))
      return false;
    const auto& ie = GetInitExpression(var);
    return ie.GetCK()->IsLogical();
  }

  /// Do the LinTerms have an arg with init expr?
  bool HasInitExpression(const LinTerms& lt) {
    for (auto v: lt.vars())
      if (HasInitExpression(v))
        return true;
    return false;
  }

  /// Do the QuadTerms have an arg with init expr?
  bool HasInitExpression(const QuadTerms& qt) {
    for (auto e: qt.get_folded())
      if (HasInitExpression(e.first.first)
          || HasInitExpression(e.first.second))
        return true;
    return false;
  }

  /// Get func con context.
  /// Does not check if that's a func con,
  /// user check for != CTX_NONE
  Context GetInitExprContext(int var) const {
    const auto& ie = GetInitExpression(var);
    return ie.GetCK()->GetContext(ie.GetIndex());
  }

  /// Set func expr context.
  /// @warning 1 level only, no propagation.
  ///   Use PropagateResultOfInitExpr() otherwise.
  void SetInitExprContext_NoProp(int var, Context ctx) {
    const auto& ie = GetInitExpression(var);
    ie.GetCK()->SetContext(ie.GetIndex(), ctx);
  }

  /// Add func expr context.
  /// @warning 1 level only, no propagation.
  ///   Use PropagateResultOfInitExpr() otherwise.
  void AddInitExprContext_NoProp(int var, Context ctx) {
    const auto& ie = GetInitExpression(var);
    ie.GetCK()->AddContext(ie.GetIndex(), ctx);
  }

  /// Get the init expression pointer.
	/// @return nullptr if no init expr or not this type
	template <class ConType>
  const ConType* GetInitExpressionOfType(int var) const {
		if (MPCD( HasInitExpression(var) )) {
      const auto& ci0 = MPCD( GetInitExpression(var) );
			if (IsConInfoType<ConType>(ci0)) {
				const auto& con =
						GetConstraint<ConType>(ci0);
				assert(&con);
				return &con;
			}
		}
		return nullptr;
	}

  /// Get the init expression pointer.
  /// @return nullptr if no init expr,
  ///   or not this type, or redefined/eliminated.
  template <class ConType>
  const ConType* GetActiveInitExpressionOfType(int var) const {
    if (MPCD( HasInitExpression(var) )) {
      const auto& ci0 = MPCD( GetInitExpression(var) );
      if (IsConInfoType<ConType>(ci0)
              && IsConActive(ci0)) {
        const auto& con =
            GetConstraint<ConType>(ci0);
        assert(&con);
        return &con;
      }
    }
    return nullptr;
  }

  /// Check if the constraint location points to the
	/// constraint keeper used for this ConType.
	template <class ConType>
	bool IsConInfoType(const ConInfo& ci) const {
		return &(BasicConstraintKeeper&)
				(GET_CONST_CONSTRAINT_KEEPER(ConType))
				== ci.GetCK();
	}

  /// Check if \a ci points to an active constraint
  bool IsConActive(const ConInfo& ci) const {
    return !ci.GetCK()->IsRedundant(ci.GetIndex());
  }


  /////////////////////// AUTO LINKING ////////////////////////////

  /// Make AutoLinker.
  /// @param con: source constraint
  /// @param i: \a con's index
  /// @return An AutoLinkScope object for \a con
  ///   which is about to be converted.
  template <class Con>
  pre::AutoLinkScope<Impl> MakeAutoLinker(const Con& , int i) {
    return {
        *(Impl*)this,
        MPD( template SelectValueNode<Con>(i) )
    };
  }

  /// Make an empty autolinker.
  /// We'd link manually but at least we care.
  pre::AutoLinkScope<Impl> MakeEmptyLinker(
      pre::NodeRange src) {
    pre::AutoLinkScope<Impl> result
        {
                  *(Impl*)this,   // 1-index source allowed only
                  {src.GetValueNode(), src.GetIndexRange().beg_}
        };
    TurnOffAutoLinking(false);
    return result;
  }

  /// Auto link node range \a nr.
  /// The nodes of \a nr will be autolinked with \a auto_link_src_item_.
  /// Means, a link is created automatically, without the
  /// conversion/flattening code doing anything.
  /// This is used to propagate values via flattened expression trees
  /// and conversions, as well as to export the conversion tree.
  pre::NodeRange AutoLink(pre::NodeRange nr) {
    if (DoingAutoLinking()) {
      if (auto_link_targ_items_.empty() ||
          !auto_link_targ_items_.back().TryExtendBy(nr))
        auto_link_targ_items_.push_back(nr);
    } else {
      // If not proper autolinking,
      // check that we have called MakeEmptyLinker
      // (hoping we care to link manually)
      // or called TurnAutoLinkingOff(false)
      assert(is_autolinking_requested_);
    }
    return nr;
  }

  /// Whether we should auto-link new items
  bool DoingAutoLinking() const
  { return auto_link_src_item_.IsValid(); }

  /// Turn off auto-linking for current conversion.
  /// @param f_full: if it's full stop,
  ///   otherwise we continue conversions but
  ///   add manual links.
  ///   This mechanism is necessary to check
  ///   that we always care about linking.
  ///
  ///   So when switching to manual linking,
  ///   call with \a f_full=false,
  ///   e.g., from MakeEmptyLinker().
  void TurnOffAutoLinking(bool f_full=true) {
    auto_link_src_item_.Invalidate();
    auto_link_targ_items_.clear();
    is_autolinking_requested_ = !f_full;
  }

  /// Get autolink source node range
  pre::NodeRange GetAutoLinkSource() const
  { return auto_link_src_item_; }

  /// Set autolink source node range
  void SetAutoLinkSource(pre::NodeRange nr)
  { assert(nr.IsSingleIndex()); auto_link_src_item_=nr; }

  /// Get autolink target node ranges
  const std::vector<pre::NodeRange>& GetAutoLinkTargets() const
  { return auto_link_targ_items_; }


public:
  /// The internal flat model type
  using ModelType = FlatModel;
  /// The internal flat model object, const ref
  const ModelType& GetModel() const { return *this; }
  /// The internal flat model object, ref
  ModelType& GetModel() { return *this; }


  ///////////////////////////////////////////////////////////////////////
  /////////////////////// OPTIONS /////////////////////////
public:
  /// Whether the ModelAPI accepts quadratic objectives
  static bool ModelAPIAcceptsQuadObj() {
    return 0 < ModelAPI::AcceptsQuadObj();
  }

  /// Whether the ModelAPI accepts quadratic constraints
  bool ModelAPIAcceptsQC() const {
    return ModelAPIAcceptsAndRecommends(
          (const QuadConLE*)nullptr)     // if accepts QuadConLE
        && ModelAPIAcceptsAndRecommends(
          (const QuadConEQ*)nullptr)
        && ModelAPIAcceptsAndRecommends(
          (const QuadConGE*)nullptr);
  }

  /// Whether the ModelAPI accepts nonconvex QC
  static bool ModelAPIAcceptsNonconvexQC() {
    return ModelAPI::AcceptsNonconvexQC();
  }

  /// Whether the ModelAPI recommends
  /// logicalizing products of 2 binaries
  static bool ModelAPIWantsLogicalProd2Bins() {
    return ModelAPI::WantLogicalizedProd2Bin();
  }

  /// Whether the ModelAPI recommends
  /// recognizing signpow()
  static bool ModelAPIWantsSignpow() {
    return ModelAPI::WantSignPow();
  }

  /// Whether the ModelAPI recommends
  /// recognizing logistic()
  static bool ModelAPIWantsLogistic() {
    return ModelAPI::WantLogistic();
  }

  /// Ask if the solver can recognize SOCP corner cases
  /// (non-std representations such as xy>=1, see tests)
  /// from quadratic representations
  static bool ModelAPICanSOCPCornerCasesFromQC() {
    return ModelAPI::CanSOCPCornerCasesFromQC();
  }

  /// Whether the solver can mix conic quadratic
  /// (entered via dedicated API)
  /// and direct quadratic constraints
  static bool ModelAPICanMixConicQCAndQC() {
    return ModelAPI::CanMixConicQCAndQC();
  }

  /// Whether the ModelAPI accepts quadratic cones
  int ModelAPIAcceptsQuadraticCones() const {
		return
        0 != std::max(
          (int)GetConstraintAcceptance_DEFAULT(
            (QuadraticConeConstraint*)nullptr),
          (int)GetConstraintAcceptance_DEFAULT(
            (RotatedQuadraticConeConstraint*)nullptr));
	}

  /// Number of QC -> SOCP conversions
  void IncQC2SOCPAttempted() { ++nQC2SOCPAttempted_; }
  void IncQC2SOCPSucceeded() { ++nQC2SOCPSucceeded_; }
  int NumQC2SOCPAttempted() const { return nQC2SOCPAttempted_; }
  int NumQC2SOCPSucceeded() const { return nQC2SOCPSucceeded_; }

  /// Number of exp cones recognized
  void IncExpConeCounter() { ++nExpConesRecognized_; }
  int NumExpConesRecognized() const { return nExpConesRecognized_; }

	/// Whether the ModelAPI accepts exp cones
	int ModelAPIAcceptsExponentialCones() {
		return
        (int)ModelAPIAcceptsAndRecommends(
          (ExponentialConeConstraint*)nullptr);
	}


private:
  struct Options {
    std::string file_graph_export_;
    int preprocessAnything_ = 1;
    int preprocessEqualityResultBounds_ = 1;
    int preprocessEqualityBvar_ = 1;
    int preprocessInequalityRhs_ = 1;
    int preprocessInequalityResultBounds_ = 1;
    int preprocessIneq2BndEq_ = 1;
    int preprocessIneq2Related_ = 1;

    int preproUnnest_ = 15;
    int preproSortUnify_  = 1;
    int boundLogArg_ = 0;

    int propCtxIneq_ = 1;
    int propCtxBndEq_ = 1;
    int propCtxCountNumberof_ = 7;

    int passQuadObj_ = ModelAPIAcceptsQuadObj();
    int passQuadCon_ = 1;
    int useQP2Pass_ = 1;
    double QPMultOutCard_ = 1e9;
    int passSOCPCones_ = 0;
    int passSOCP2QC_ = 0;
    int passExpCones_ = 0;

    int accAll_ = -1;
    int accExpr_ = static_cast<
                       std::underlying_type_t<ExpressionAcceptanceLevel> >
                   (ModelAPI::ExpressionInterfaceAcceptanceLevel())
                   -1;               // If available, 0 or 1

    int relax_ = 0;

    double modelfeastol_ = 1e-6;
    double modelfeastolrel_ = 1e-6;

    int mo_options_ {1};

    int solcheckmode_ = 1+2+512;
    bool solcheckinfeas_ = false;
    bool solcheckfail_ = false;
    double solfeastol_ = 1e-6;
    double solfeastolrel_ = 1e-6;
    double solinttol_ = 1e-5;
    int sol_round_ = 100;
    int sol_prec_ = 100;
    int solchkoutlev_ = 0;

    int nlassign_lev_ = ModelAPI::NLAssignLevelDefault();
    int nlreif_lev_ = ModelAPI::NLReifLevelDefault();
  };
  Options options_;


public:             // public for CRTP
  /// Graph export file
  const std::string& graph_export_file() const
  { return options_.file_graph_export_; }

  /// Whether to parse QP expressions in 2 passes
  int IfParseQPIn2Passes() const { return options_.useQP2Pass_; }

  /// Up to which QP matrix cardinality should we multiply out
  double QPMultOutCard() const { return options_.QPMultOutCard_; }

  /// Whether we should relax integrality
  int relax() const { return options_.relax_; }

  /// Propagate context into conditional inequalities?
  int IfPropCtxCondIneq() const { return options_.propCtxIneq_; }
  /// Propagate context into conditional equalities-to-bound?
  int IfPropCtxCondEqBnd() const { return options_.propCtxBndEq_; }
  /// Propagate context into count/numberof?
  int IfPropCtxCountNumberof() const { return options_.propCtxCountNumberof_; }
  /// Bound argument of logarithm?
  bool IfBoundLogArg() const { return options_.boundLogArg_; }

  /// Model checking options
  double model_feas_tol() const { return options_.modelfeastol_; }
  double model_feas_tol_rel() const { return options_.modelfeastolrel_; }

  /// Use multiobjective options?
  int multiobj_options() const { return options_.mo_options_; }

  /// Solution checking options
  int sol_check_mode() const { return options_.solcheckmode_; }
  bool sol_check_infeas() const { return options_.solcheckinfeas_; }
  bool sol_check_fail() const { return options_.solcheckfail_; }
  double sol_feas_tol() const { return options_.solfeastol_; }
  double sol_feas_tol_rel() const { return options_.solfeastolrel_; }
  double sol_int_tol() const { return options_.solinttol_; }
  int sol_round() const { return options_.sol_round_; }
  int sol_prec() const { return options_.sol_prec_; }
  int sol_check_outlev() const { return options_.solchkoutlev_; }

  /// Option *cvt:expr:nlassign*
  int NLAssignLevel() const { return options_.nlassign_lev_; }
  /// Option *cvt:expr:nlreif*
  int NLReifLevel() const { return options_.nlreif_lev_; }
  /// Implement for ConverterInfo
  int RefCountMaxAlgebraic() const override { return NLAssignLevel(); }

public:
  /// Init FlatConverter options
  void InitOptions() {
    InitOwnOptions();
    GetModelAPI().InitStandardOptions();
    GetModelAPI().InitCustomOptions();
  }


private:
  const std::string solchkfailtext_ {
    "Fail on MP solution check violations, with solve result "
    + std::to_string(sol::MP_SOLUTION_CHECK) + '.'
  };

  int DefaultSOCPMode() const {
    return
        !ModelAPIAcceptsQC() && !ModelAPIAcceptsQuadraticCones()
        ? 0
        : ModelAPICanSOCPCornerCasesFromQC() ? 1
                                             : 2;
  }
  int DefaultSOCP2QCMode() const {
    return
        ((!ModelAPIAcceptsQC() || ModelAPICanMixConicQCAndQC())
         && ModelAPIAcceptsQuadraticCones())
        ? 0
        : (!ModelAPICanMixConicQCAndQC()
           && ModelAPIAcceptsQuadraticCones()) ? 1
                                             : 2;
  }
  std::string socp_mode_text_;
  std::string socp2qc_mode_text_;
  const mp::OptionValueInfo socp_values_[3] = {
    { "0", "Do not recognize SOCP forms", 0},
    { "1", "Recognize from non-quadratic expressions only (sqrt, abs)", 1},
    { "2",
      "Recognize from quadratic and non-quadratic SOCP forms. "
      "Helpful if the solver does not recognize non-standard forms", 2}
  };
  const mp::OptionValueInfo socp2qc_values_[3] = {
    { "0", "Do not convert", 0},
    { "1", "Convert if no other cone types found, and "
      "not all original quadratics could be recognized as SOC, "
      "in particular if the objective is quadratic", 1},
    { "2", "Always convert", 2}
  };
  const mp::OptionValueInfo values_allexpr_acceptance_[2] = {
      { "0", "Not accepted, all expressions will be treated as flat constraints, "
            "or redefined", 0},
      { "1", "Accepted. See the individual acc:... options", 1}
  };


  void InitOwnOptions() {
    /// Should be called after adding all constraint keepers
    FlatModel::ConsiderItemTypeOptions(*this, GetModelAPI(), GetEnv());

    GetEnv().AddStoredOption("tech:writegraph cvt:writegraph writegraph exportgraph",
        "File to export conversion graph. Format: JSON Lines.",
        options_.file_graph_export_);
    GetEnv().AddOption("cvt:pre:all",
        "0/1*: Set to 0 to disable most presolve in the flat converter.",
        options_.preprocessAnything_, 0, 1);
    GetEnv().AddOption("cvt:pre:eqresult",
        "0/1*: Preprocess reified equality comparison's decidable cases.",
        options_.preprocessEqualityResultBounds_, 0, 1);
    GetEnv().AddOption("cvt:pre:eqbinary",
        "0/1*: Preprocess reified equality comparison with a binary variable.",
        options_.preprocessEqualityBvar_, 0, 1);
    GetEnv().AddOption("cvt:pre:ineqresult",
                       "0/1*: Preprocess reified inequality comparison's decidable cases.",
                       options_.preprocessInequalityResultBounds_, 0, 1);
    GetEnv().AddOption("cvt:pre:ineqrhs",
                       "0/1*: Preprocess reified inequality comparison's right-hand sides "
                       "(round for integer expression body).",
                       options_.preprocessInequalityRhs_, 0, 1);
    GetEnv().AddOption("cvt:pre:ineq2bndeq ineq2bndeq",
                       "0/1*: Preprocess reified inequality expr <(=) c, where "
                       "c <)=( lb(expr)+cvt:mip:eps, into expr == lb(expr), "
                       "which works better on some benchmarks/solvers.",
                       options_.preprocessIneq2BndEq_, 0, 1);
    GetEnv().AddOption("cvt:pre:ineq2related ineq2related ineq2rel",
                       "0/1*: Unify related reified inequalities: "
                       "<=c, <c+cvt:mip:eps, >c, >=c+cvt:mip:eps.",
                       options_.preprocessIneq2Related_, 0, 1);

    GetEnv().AddOption("cvt:pre:unnest cvt:unnest cvt:pre:inline cvt:inline",
        "Inline nested expressions. Bitwise OR of the following values:\n"
                       "\n"
                       "|  1 - AND/FORALL and OR/EXISTS expressions\n"
                       "|  2 - Linear subexpressions\n"
                       "|  4 - Quadratic subexpressions\n"
                       "|  8 - MIN/MAX.\n"
                       "\n"
                       "See also option cvt:dvelim concerning only the input model. "
                       "Default 15.",
        options_.preproUnnest_, 0, 15);
    GetEnv().AddOption("cvt:pre:sort cvt:sort",
                       "0/1*: Sort and eliminate duplicates in arguments "
                       "of AND, OR, MIN, MAX. Sort arguments of "
                       "COUNT, ATLEAST, EXACTLY, ATMOST, NUMBEROF, ALLDIFF. "
                       "Can be necessary for some solvers.",
                       options_.preproSortUnify_, 0, 1);

    GetEnv().AddOption("cvt:pre:ctx2ineq ctx2ineq",
                       "0/1*: Propagate exact context into conditional inequalities, "
                       "vs always mixed. See #267.\n"
                       "\n"
                       "Finer control provided by cvt:pre:ctx:cond...(le/ge) options.",
                       options_.propCtxIneq_, 0, 1);
    GetEnv().AddOption("cvt:pre:ctx2bndeq ctx2bndeq",
                       "0/1*: Propagate exact context into conditional "
                       "(dis)equalities-to-bound, vs always mixed. "
                       "Can be affected by cvt:pre:ineq2bndeq. See #267.",
                       options_.propCtxBndEq_, 0, 1);
    GetEnv().AddOption("cvt:pre:ctx2count ctx2count",
                       "DEPRECATED. Use ctx2bndeq. NEW DEFAULT.\n\n"
                       "Propagate exact context into atleast/atmost/exactly, "
                       "count and numberof expressions, "
                       "vs always mixed. Bitwise OR of the following values:\n"
                       "\n"
                       "|  1 - atleast/atmost/exactly, count\n"
                       "|  2 - numberof with constant reference value\n"
                       "|  4 - numberof with variable reference value.\n"
                       "\n"
                       "Default 7, see #267.",
                       options_.propCtxCountNumberof_, 0, 7);

    GetEnv().AddOption("cvt:pre:boundsbest boundsbest",
                       "0*/1: Submit best-known variable bounds to the solver. "
                       "Can inhibit its presolve.\n"
                       "\n"
                       "Note: when a variable can be fixed, the stronger bounds "
                       "are always submitted.",
                       GetModel().if_submit_best_known_bounds(), 0, 1);

    GetEnv().AddOption("cvt:pre:continuous_fixed_vars continuous_fixed_vars ctg_fixed",
                       "0/1*: Make fixed variables continuous, "
                       "to avoid fake MIPs.",
                       GetModel().if_make_fixed_vars_continuous(), 0, 1);

    GetEnv().AddOption("cvt:pre:boundlogarg boundlogarg",
                       "0*/1: Bound logarithm arguments to nonnegative.",
                       options_.boundLogArg_, 0, 1);

    GetEnv().AddOption("cvt:quadobj passquadobj",
                       ModelAPIAcceptsQuadObj() ?
        "0/1*: Pass quadratic objective terms to the solver. "
        "When 0, if the solver accepts quadratic constraints, "
        "such a constraint will be created with those, "
        "otherwise linearly approximated."
                       :
        "0*/1: Pass quadratic objective terms to the solver. "
        "When 0, if the solver accepts quadratic constraints, "
        "such a constraint will be created with those, "
        "otherwise linearly approximated.",
        options_.passQuadObj_, 0, 1);
    GetEnv().AddOption("cvt:quadcon passquadcon",
                       "0/1*: set to 0 to disable quadratic constraints. "
                       "Synonym for acc:quad..=0. "
                       "Setting to 0 disables out-multiplication "
                       "of quadratic terms, then they are linearized.",
                       options_.passQuadCon_, 0, 1);
    GetEnv().AddOption("cvt:qp2passes cvt:qp2pass qp2passes qp2pass",
                       "0/1*: Parse sums of QP expressions in 2 passes. "
                       "Usually faster.",
                       options_.useQP2Pass_, 0, 1);
    GetEnv().AddOption("cvt:multoutcard multoutcard",
                       "Up to which (estimated) QP matrix cardinality "
                       "should a product of 2 linear expressions "
                       "be multiplied out. Default 1e9.\n"
                       "\n"
                       "Low value can speed up model input, but prone to "
                       "numerical issues.",
                       options_.QPMultOutCard_, 0.0, 1e20);


    GetEnv().AddOption("cvt:expcones expcones",
                       ModelAPIAcceptsExponentialCones() ?
                         "0/1*: Recognize exponential cones." :
                         "0*/1: Recognize exponential cones.",
                       options_.passExpCones_, 0, 1);
    options_.passExpCones_ = ModelAPIAcceptsExponentialCones();
    // Should be after construction
    socp_mode_text_ =
      "Second-Order Cone recognition mode:\n"
      "\n.. value-table::\n"
      "Recognized SOCP forms can be further converted to "
      "(SOCP-standardized) quadratic constraints, see cvt:socp2qc. "
      "Default: " + std::to_string(DefaultSOCPMode()) + ".";
    GetEnv().AddStoredOption("cvt:socp socpmode socp",
                       socp_mode_text_.c_str(),
                       options_.passSOCPCones_, socp_values_);
    options_.passSOCPCones_ = DefaultSOCPMode();
    socp2qc_mode_text_ =
      "Mode to convert recognized SOCP forms to "
      "SOCP-standardized quadratic constraints:\n"
      "\n.. value-table::\n"
      "Such conversion can be necessary "
      "if the solver does not accept "
      "a mix of conic and quadratic constraints/objectives. "
      "Default: " + std::to_string(DefaultSOCP2QCMode()) + ".";
    GetEnv().AddStoredOption("cvt:socp2qc socp2qcmode socp2qc",
                       socp2qc_mode_text_.c_str(),
                       options_.passSOCP2QC_, socp2qc_values_);
    options_.passSOCP2QC_ = DefaultSOCP2QCMode();

    GetEnv().AddStoredOption("acc:_all",
                             "Solver acceptance level for all "
                             "constraints and expressions. "
                             "Value meaning: as described in the specific "
                             "acc:... options.\n"
                             "\n"
                             "Can be useful to disable all reformulations (acc:_all=2), "
                             "or force linearization (acc:_all=0.)",
                             options_.accAll_, 0, 4);

    if constexpr (IfAcceptingNLOutput())
      GetEnv().AddStoredOption("acc:_expr",
                        fmt::format(
                            "Solver acceptance level for all expressions, "
                            "default {}:\n\n.. value-table::",
                                   options_.accExpr_).c_str(),
                        options_.accExpr_, values_allexpr_acceptance_);
    else
      GetEnv().AddStoredOption("acc:_expr", "HIDDEN", options_.accExpr_, 0, 1);

    GetEnv().AddOption("alg:relax relax",
        "0*/1: Whether to relax integrality of variables.",
        options_.relax_, 0, 1);

    GetEnv().AddOption("obj:multi:options multiobjoptions",
                       "0/1*: Regard multiobjective option suffixes "
                       "which are objective suffixes beginning with option_. "
                       "Example: suffix option_timelim; let _obj[2].option_timelim:=15;",
                       options_.mo_options_, 0, 1);

    GetEnv().AddOption("cvt:pre:feastol pre:feastol pre:eps pre:feastolabs pre:epsabs",
                       "Absolute tolerance to check variable "
                       "and constraint bound contraditions. "
                       "Only warns if also pre:feastolrel is violated. "
                       "See also sol:chk:feastol. "
                       "Default 1e-6.",
                       options_.modelfeastol_, 0.0, 1e100);
    GetEnv().AddOption("cvt:pre:feastolrel pre:feastolrel pre:epsrel",
                       "Relative tolerance to check variable "
                       "and constraint bound contradictions. "
                       "Only warns if also pre:feastol is violated. "
                       "See also sol:chk:feastol. "
                       "Default 1e-6.",
                       options_.modelfeastolrel_, 0.0, 1e100);

    GetEnv().AddStoredOption(
          "sol:chk:mode solcheck checkmode chk:mode",
        "Solution checking mode. "
        "Sum of a subset of the following bits:\n"
        "\n"
        "| 1 - Check variable bounds and integrality.\n"
        "| 2 - Check original model constraints, as well as "
        "      any non-linear expression values "
        "      reported by the solver.\n"
        "| 4 - Check intermediate auxiliary constraints "
        "      (i.e., those which were reformulated further).\n"
        "| 8 - Check final auxiliary constraints sent to solver.\n"
        "| 16 - Check objective values.\n"
        "| 32, 64, 128, 256, 512 - similar, but "
        "      non-linear expressions are recomputed "
        "      (vs using their values reported by the solver.) "
        "      *Experimental.* This is an idealistic check, because "
        "      it does not consider possible tolerances "
        "      applied by the solver when computing "
        "      expression values.\n"
                             "\n"
                             "Default: 1+2+512.",
        options_.solcheckmode_, 0, 1023);
    GetEnv().AddOption("sol:chk:feastol sol:chk:eps chk:eps chk:feastol",
        "Absolute tolerance to check objective values', variable "
        "and constraint bounds' violations. "
                       "Only triggers if also sol:chk:feastolrel is violated. "
                       "See also pre:feastol. "
                       "Default 1e-6.",
        options_.solfeastol_, 0.0, 1e100);
    GetEnv().AddOption("sol:chk:feastolrel sol:chk:epsrel chk:epsrel chk:feastolrel",
        "Relative tolerance to check objective values', variable "
        "and constraint bounds' violations. "
                       "Only triggers if also sol:chk:feastol is violated. "
                       "See also pre:feastol. "
                       "Default 1e-6.",
        options_.solfeastolrel_, 0.0, 1e100);
    GetEnv().AddOption("sol:chk:inttol sol:chk:inteps sol:inteps chk:inttol",
        "Solution checking tolerance for variables' integrality. "
        "Default 1e-5.",
        options_.solinttol_, 0.0, 1e100);
    GetEnv().AddOption("sol:chk:infeas chk:infeas checkinfeas",
                       "Check even infeasible solution condidates, "
                       "whenever solver reports them.",
                       options_.solcheckinfeas_, false, true);
    GetEnv().AddOption("sol:chk:fail chk:fail checkfail",
                       solchkfailtext_.c_str(),
                       options_.solcheckfail_, false, true);
    GetEnv().AddOption("sol:chk:round chk:round chk:rnd",
        "AMPL solution_round option when checking: "
                       "round to this number of decimals after comma "
                       "(before comma if negative.)",
                       options_.sol_round_, -1000, 1000);
    GetEnv().AddOption("sol:chk:prec chk:prec chk:precision",
        "AMPL solution_precision option when checking: "
                       "number of significant digits.",
                       options_.sol_prec_, -1000, 1000);

    GetEnv().AddOption("cvt:expr:nlassign expr:nlassign",
                       fmt::format("Above which reference count, "
                                   "an algebraic formula node should be assigned to a variable "
                                   "(see acc: options). 0 means all nodes assigned. "
                                   "Default {}.", options_.nlassign_lev_).c_str(),
                       options_.nlassign_lev_, 0, INT_MAX);

    GetEnv().AddOption("cvt:expr:nlreif expr:nlreif expr:nlreify",
                       fmt::format("Above which reference count, "
                                   "a logical formula node should be assigned (reified) "
                                   " to a variable "
                                   "(see acc: options). 0 means all nodes reified. "
                                   "Default {}.", options_.nlreif_lev_).c_str(),
                       options_.nlreif_lev_, 0, INT_MAX);

    ////////////////////// Solve result codes ////////////////////////
    GetEnv().AddSolveResults({
                               {sol::MP_SOLUTION_CHECK,
                                "solved? MP solution check failed "
                                "(option sol:chk:fail) "}
                             });
  }


public:
  /// Wrapper about a specific preprocess option:
  /// checks whether \a preprocessAnything_ is on.
  int CanPreprocess(int f) const {
    return options_.preprocessAnything_ ? f : 0;
  }

  /// Whether preprocess equality result bounds
  bool IfPreproEqResBounds() const
  { return MPCD( CanPreprocess(options_.preprocessEqualityResultBounds_) ); }

  /// Whether preprocess conditional equality of a binary variable
  bool IfPreproEqBinVar() const
  { return MPCD( CanPreprocess(options_.preprocessEqualityBvar_) ); }

  /// Whether preprocess inequality result bounds
  bool IfPreproIneqResBounds() const
  { return MPCD( CanPreprocess(options_.preprocessInequalityResultBounds_) ); }

  /// Whether preprocess inequality rhs
  bool IfPreproIneqRHS() const
  { return MPCD( CanPreprocess(options_.preprocessInequalityRhs_) ); }

  /// Whether preprocess inequality into ==LB/UB
  bool IfPreproIneq2BndEq() const
  { return MPCD( CanPreprocess(options_.preprocessIneq2BndEq_) ); }

  /// Whether to unify related inequalities
  bool IfPreproIneq2Related() const
  { return MPCD( CanPreprocess(options_.preprocessIneq2Related_) ); }

  /// Whether inline nested forall, exists, lin/quad expr
  int IfPreproUnnest() const
  { return MPCD( CanPreprocess(options_.preproUnnest_) ); }

  /// Whether sort and elim duplicates in argument lists
  int IfPreproSortUnify() const
  { return MPCD( CanPreprocess(options_.preproSortUnify_) ); }


  /// Whether we pass quad obj terms to the solver without linearization
  bool IfPassQuadObj() const { return options_.passQuadObj_; }

  /// Whether we pass quad con terms to the solver without linearization
  bool IfPassQuadCon() const
  { return options_.passQuadCon_ && ModelAPIAcceptsQC(); }

  /// Whether to quadratize pow(..., const_pos_int).
  /// The fact that we use IfPassQuadCon()
  /// is much Gurobi-biased: v9.5 does not PL-linearize Pow
  /// for negative arguments
  bool IfQuadratizePowConstPosIntExp() const
  { return IfPassQuadCon(); }

  /// Recognition mode for SOCP cones
  int IfPassSOCPCones() const { return options_.passSOCPCones_; }

  /// Mode for SOCP -> QC conversion
  int SOCP2QCMode() const { return options_.passSOCP2QC_; }

  /// Decide to convert SOCP -> QC
  void Setup2ConvertSOCP2QC() { ifCvtSOCP2QC_=true; }
  /// If decided to convert SOCP -> QC
  bool IfConvertSOCP2QC() const { return ifCvtSOCP2QC_; }

  /// Recognition mode for exp cones
  int IfPassExpCones() const { return options_.passExpCones_; }


public:
  /// Typedef ModelAPIType. For tests
  using ModelAPIType = ModelAPI;

  /// AddWarning.
  /// @param key: warning category
  /// @param msg: detailed message
  void AddWarning(
      std::string key, std::string msg, bool replace=false) {
    GetEnv().AddWarning(
          std::move(key), std::move(msg), replace);
  }

  /// Provide suffix getters and setters
  void SetSuffixManip(SuffixGetterSetter sgs)
  { suf_get_set_ = sgs; }


public:
  /// Read int suffix
  ArrayRef<int> ReadIntSuffix(const SuffixDef<int>& sd)
  { assert(suf_get_set_.sgi_); return suf_get_set_.sgi_(sd); }

  /// Read double suffix
  ArrayRef<double> ReadDblSuffix(const SuffixDef<double>& sd)
  { assert(suf_get_set_.sgd_); return suf_get_set_.sgd_(sd); }

  const SuffixSet& Suffixes(suf::Kind kind) {
    { assert(suf_get_set_.ss_); return suf_get_set_.ss_(kind); }
  }

protected:
  void CheckNumVars(pre::ModelValuesDbl& sol) {
    auto& xx = sol.GetVarValues()();
    if (xx.size()) {                    // solution available
      if (xx.size() < MPCD( num_vars() )) {
        MPCD( GetEnv() ).AddWarning(
            "FewerVariables",
            fmt::format("Solver reported {} variables,\n"
                        "less than {} received.\n"
                        "Solution might be incorrect.",
                        xx.size(), MPCD( num_vars() )));
        xx.resize( MPCD( num_vars() ) );
      }
    }
  }

  /// Recompute implicit aux vars
  /// (those corresponding to expressions
  ///   and/or eliminated functional constraints).
  /// Needed for MO emulator and sol checker.
  void RecomputeUnusedAuxVars(pre::ModelValuesDbl& sol) {
    auto& xx = sol.GetVarValues()();
    if (xx.size()) {                    // solution available
      auto var_is_used = MPCD( GetVarElimFlags() );
      var_is_used.flip();
      xx = MPD( RecomputeAuxVars(xx, var_is_used) );
    }
  }


private:
  /// We store ModelApi in the converter for speed.
  /// Should be before constraints
  ModelAPIType modelapi_;
  /// solve iteration
  int n_solve_iter_ {0};
  /// Suffix getters and setters
  SuffixGetterSetter suf_get_set_;
  /// ValuePresolver: should be init before constraint keepers
  /// and links
  pre::ValuePresolver value_presolver_
  {
      GetModel(), GetEnv(), (BasicLogger&)GetModel().GetFileAppender(),
      [this](                           // Solution checker
          ArrayRef<double> x,
          const pre::ValueMapDbl& y,
          ArrayRef<double> obj,
          void* p_extra) -> bool {
        return !this->options_.solcheckmode_   // not desired
               || MPD( CheckSolution(x, y, obj, p_extra) );
      },
      [this](pre::ModelValuesDbl& sol)  // Solution pre-postsolver
      {
        MPD( CheckNumVars(sol) );       // XPRESS 44.01.04
        MPD( RecomputeUnusedAuxVars(sol) );
        MPD( ProcessMOIterationUnpostsolvedSolution(sol) );
      }
  };
  pre::CopyLink copy_link_ { GetValuePresolver() }; // the copy links
  pre::One2ManyLink one2many_link_ { GetValuePresolver() }; // the 1-to-many links
  pre::NodeRange auto_link_src_item_;   // the source item for autolinking
  std::vector<pre::NodeRange> auto_link_targ_items_;
  /// This is to check that we always care
  /// to set up autolinking,
  /// or declare it empty
  bool is_autolinking_requested_ {false};

  pre::Many2OneLink many2one_link_ { GetValuePresolver() }; // the many-to-one links

	ConicConverter<Impl> conic_cvt_ { *static_cast<Impl*>(this) };
  int nQC2SOCPAttempted_= 0;
  int nQC2SOCPSucceeded_= 0;
  int nExpConesRecognized_ = 0;
  bool ifCvtSOCP2QC_ = 0;

  SDPConverter<Impl> sdp_cvt_ { *static_cast<Impl*>(this) };

	std::vector<int> refcnt_vars_;
  int constr_depth_ = 0;    // tree depth of new constraints

  /// AMPL item names
  std::vector<std::string>
  var_names_,
  con_names_,   // no SOS here, they go directly into the SOS
  obj_names_;

  /// Model stats
  FlatModelInfo::VarInfo modelinfo_flat0_vars_;
  FlatModelInfo::ObjInfo modelinfo_flat0_objs_;
  FlatModelInfo::ConstrTypeMapByName modelinfo_flat0_cons_;

protected:
  /////////////////////// CONSTRAINT KEEPERS /////////////////////////
  /// Constraint keepers and converters should be initialized after
  ///  \a value_presolver_

  /// Define constraint keepers for all constraint types.
  /// No maps for static constraints.
  /// 2nd parameter: solver options for this constraint,
  /// in case it is accepted by the solver natively and
  /// is convertible by us.
  /// 3rd parameter: reformulation priority (double).
  /// Can be changed in a derived class by ConstraintCvtPriority().
  /// NOTE: The reformulation meta-graph should be acyclic #248.

  /// Flattened NL expressions
  STORE_CONSTRAINT_TYPE__WITH_MAP(AbsConstraint, "acc:abs", 100)
  STORE_CONSTRAINT_TYPE__WITH_MAP(AllDiffConstraint, "acc:alldiff", 200)

  STORE_CONSTRAINT_TYPE__NO_MAP(
      ComplementarityQuadratic, "acc:complquad", 300)
  STORE_CONSTRAINT_TYPE__NO_MAP(
      ComplementarityLinear, "acc:compl acc:compllin", 350)
  STORE_CONSTRAINT_TYPE__NO_MAP(
      NLComplementarity, "acc:nlcompl", 360)

  STORE_CONSTRAINT_TYPE__WITH_MAP(CountConstraint, "acc:count", 400)

  STORE_CONSTRAINT_TYPE__WITH_MAP(DivConstraint, "acc:div", 600)
  STORE_CONSTRAINT_TYPE__WITH_MAP(IfThenConstraint, "acc:ifthen", 700)
  STORE_CONSTRAINT_TYPE__WITH_MAP(ImplicationConstraint, "acc:impl", 800)

  STORE_CONSTRAINT_TYPE__WITH_MAP(LogisticConstraint,
                                  "acc:logi acc:logistic", 880)
  STORE_CONSTRAINT_TYPE__WITH_MAP(SignpowConstExpConstraint,
                                  "acc:signpowc acc:signpowconstexp", 890)

  STORE_CONSTRAINT_TYPE__WITH_MAP(PowConstExpConstraint, "acc:powc acc:powconstexp", 900)

  STORE_CONSTRAINT_TYPE__WITH_MAP(PowConstraint, "acc:pow", 950) // -> exp, log
  STORE_CONSTRAINT_TYPE__WITH_MAP(ExpConstraint, "acc:exp", 1000)
  STORE_CONSTRAINT_TYPE__WITH_MAP(ExpAConstraint, "acc:expa acc:expA", 1002)
  STORE_CONSTRAINT_TYPE__WITH_MAP(LogConstraint, "acc:log", 1004)
  STORE_CONSTRAINT_TYPE__WITH_MAP(LogAConstraint, "acc:loga acc:logA", 1006)
  STORE_CONSTRAINT_TYPE__WITH_MAP(SinConstraint, "acc:sin", 1008)
  STORE_CONSTRAINT_TYPE__WITH_MAP(CosConstraint, "acc:cos", 1010)
  STORE_CONSTRAINT_TYPE__WITH_MAP(TanConstraint, "acc:tan", 1012)
  STORE_CONSTRAINT_TYPE__WITH_MAP(AsinConstraint, "acc:asin", 1014)
  STORE_CONSTRAINT_TYPE__WITH_MAP(AcosConstraint, "acc:acos", 1016)
  STORE_CONSTRAINT_TYPE__WITH_MAP(AtanConstraint, "acc:atan", 1018)
  STORE_CONSTRAINT_TYPE__WITH_MAP(SinhConstraint, "acc:sinh", 1020)
  STORE_CONSTRAINT_TYPE__WITH_MAP(CoshConstraint, "acc:cosh", 1022)
  STORE_CONSTRAINT_TYPE__WITH_MAP(TanhConstraint, "acc:tanh", 1024)
  STORE_CONSTRAINT_TYPE__WITH_MAP(AsinhConstraint, "acc:asinh", 1026)
  STORE_CONSTRAINT_TYPE__WITH_MAP(AcoshConstraint, "acc:acosh", 1028)
  STORE_CONSTRAINT_TYPE__WITH_MAP(AtanhConstraint, "acc:atanh", 1030)

  STORE_CONSTRAINT_TYPE__WITH_MAP(CallConstraint, "acc:call", 1040)

  STORE_CONSTRAINT_TYPE__WITH_MAP(SDPDotProdConstraint, "acc:sdpdotprod acc:sdpdot", 1045)

  STORE_CONSTRAINT_TYPE__WITH_MAP(MaxConstraint, "acc:max", 1100)
  STORE_CONSTRAINT_TYPE__WITH_MAP(MinConstraint, "acc:min", 1200)
  STORE_CONSTRAINT_TYPE__WITH_MAP(NumberofConstConstraint,
                                  "acc:numberofconst", 1300)
  STORE_CONSTRAINT_TYPE__WITH_MAP(NumberofVarConstraint,
                                  "acc:numberofvar", 1350)
  STORE_CONSTRAINT_TYPE__WITH_MAP(PLConstraint,
      "acc:pl acc:pwl acc:piecewise", 1500)
  STORE_CONSTRAINT_TYPE__NO_MAP(SOS1Constraint, "acc:sos1", 1600)
  STORE_CONSTRAINT_TYPE__NO_MAP(SOS2Constraint, "acc:sos2", 1700)

  /// Store UEncConstr
  STORE_CONSTRAINT_TYPE__NO_MAP(
      UnaryEncodingConstraint, "acc:uenc", 1800)
  /// Dummy conversion for UEncConstr
  Context Convert(const UnaryEncodingConstraint& )
  { return Context::CTX_ROOT; }
  /// Say we can (for acc:_all=0)
  bool IfHasCvt_impl(const UnaryEncodingConstraint* ) {
    return true;
  }

  STORE_CONSTRAINT_TYPE__WITH_MAP(CondQuadConEQ, "acc:condquadeq", 1900)
  STORE_CONSTRAINT_TYPE__WITH_MAP(CondLinConEQ, "acc:condlineq", 1950)

  STORE_CONSTRAINT_TYPE__WITH_MAP(CondQuadConLE, "acc:condquadle", 2000)
  STORE_CONSTRAINT_TYPE__WITH_MAP(CondQuadConLT, "acc:condquadlt", 2010)
  STORE_CONSTRAINT_TYPE__WITH_MAP(CondQuadConGE, "acc:condquadge", 2020)
  STORE_CONSTRAINT_TYPE__WITH_MAP(CondQuadConGT, "acc:condquadgt", 2030)

  STORE_CONSTRAINT_TYPE__WITH_MAP(CondLinConLE, "acc:condlinle", 2050)
  STORE_CONSTRAINT_TYPE__WITH_MAP(CondLinConLT, "acc:condlinlt", 2060)
  STORE_CONSTRAINT_TYPE__WITH_MAP(CondLinConGE, "acc:condlinge", 2070)
  STORE_CONSTRAINT_TYPE__WITH_MAP(CondLinConGT, "acc:condlingt", 2080)

  STORE_CONSTRAINT_TYPE__WITH_MAP(NotConstraint, "acc:not", 2100)

  /// Used only for expression output,
  /// flat model keeps this in algebraic form
  STORE_CONSTRAINT_TYPE__WITH_MAP(
      EquivalenceConstraint, "acc:equiv acc:equivalence", 2150)

  STORE_CONSTRAINT_TYPE__WITH_MAP(AndConstraint,
                                  "acc:and acc:forall", 2200)
  STORE_CONSTRAINT_TYPE__WITH_MAP(OrConstraint,
                                  "acc:or acc:exists", 2300)

  /// No maps for static constra ints
  STORE_CONSTRAINT_TYPE__NO_MAP(
      IndicatorConstraintQuadEQ, "acc:indquadeq", 2400)
  STORE_CONSTRAINT_TYPE__NO_MAP(
      IndicatorConstraintQuadLE, "acc:indquadle", 2410)
  STORE_CONSTRAINT_TYPE__NO_MAP(
      IndicatorConstraintQuadGE, "acc:indquadge", 2420)
  STORE_CONSTRAINT_TYPE__NO_MAP(
      IndicatorConstraintLinEQ, "acc:indeq acc:indlineq", 2450)
  STORE_CONSTRAINT_TYPE__NO_MAP(
      IndicatorConstraintLinLE, "acc:indle acc:indlinle", 2460)
  STORE_CONSTRAINT_TYPE__NO_MAP(
      IndicatorConstraintLinGE, "acc:indge acc:indlinge", 2470)

  STORE_CONSTRAINT_TYPE__NO_MAP(
      PowerConeConstraint, "acc:powercone", 3000)
  STORE_CONSTRAINT_TYPE__NO_MAP(
      RotatedQuadraticConeConstraint, "acc:rotatedquadcone", 3001)
  STORE_CONSTRAINT_TYPE__NO_MAP(
      QuadraticConeConstraint, "acc:quadcone", 3002)
  STORE_CONSTRAINT_TYPE__NO_MAP(
      ExponentialConeConstraint, "acc:expcone", 3010)
  STORE_CONSTRAINT_TYPE__NO_MAP(
      GeometricConeConstraint, "acc:geomcone", 3020)

  /// Static algebraic cons
  STORE_CONSTRAINT_TYPE__NO_MAP(QuadConRange,  // Before LinConRange
                                "acc:quadrange acc:quadrng", 3090)
  STORE_CONSTRAINT_TYPE__NO_MAP(LinConRange,   // before QuadFuncCon
                                "acc:linrange acc:linrng", 3091)

  /// Our own functional constraints: LFC, QFC.
  /// We'll also add inlining with priority 3099,
  /// see AddConversionAction() in the constructor #266.
  /// Inline algebraic subexpressions in algebraic
  /// constraints and objectives.
  /// This includes algebralized indicators.
  STORE_CONSTRAINT_TYPE__WITH_MAP(
      QuadraticFunctionalConstraint, "acc:quadfn acc:quadfunccon", 3100)
  STORE_CONSTRAINT_TYPE__WITH_MAP(
      LinearFunctionalConstraint, "acc:linfn acc:linfunccon", 3200)

  STORE_CONSTRAINT_TYPE__NO_MAP(QuadConLE, "acc:quadle", 4100)
  STORE_CONSTRAINT_TYPE__NO_MAP(QuadConEQ, "acc:quadeq", 4200)
  STORE_CONSTRAINT_TYPE__NO_MAP(QuadConGE, "acc:quadge", 4300)

  STORE_CONSTRAINT_TYPE__NO_MAP(LinConLE, "acc:linle", 5100)
  STORE_CONSTRAINT_TYPE__NO_MAP(LinConEQ, "acc:lineq", 5200)
  STORE_CONSTRAINT_TYPE__NO_MAP(LinConGE, "acc:linge", 5300)

  ////////////////////// NL constraints & expressions ///////////////////////
  STORE_CONSTRAINT_TYPE__NO_MAP(
      NLConstraint, "acc:nlcon acc:nlalgcon", 10000)
  STORE_CONSTRAINT_TYPE__NO_MAP(
      NLAssignEQ, "acc:nlassigneq", 10100)
  STORE_CONSTRAINT_TYPE__NO_MAP(
      NLAssignLE, "acc:nlassignle", 10200)
  STORE_CONSTRAINT_TYPE__NO_MAP(
      NLAssignGE, "acc:nlassignge", 10300)

  STORE_CONSTRAINT_TYPE__NO_MAP(
      NLLogical, "acc:nllogcon acc:nllogical", 12000)
  STORE_CONSTRAINT_TYPE__NO_MAP(
      NLReifEquiv, "acc:nlreifequiv", 12100)
  STORE_CONSTRAINT_TYPE__NO_MAP(
      NLReifImpl, "acc:nlreifimpl", 12200)
  STORE_CONSTRAINT_TYPE__NO_MAP(
      NLReifRimpl, "acc:nlreifrimpl", 12300)


  protected:
  ////////////////////// Default map accessors /////////////////////////
  /// Constraints without map should overload these by empty methods ///

  /// MapFind.
  /// Can be overloaded for more complex behavior.
  /// @param con: constraint reference
  /// @return constraint index, or -1
  template <class Constraint>
  int MapFind(const Constraint& con) {
    return MPD( MapFind__Impl(con) );
  }

  /// MapInsert.
  /// Can be overloaded for more complex behavior.
  /// @param con: the constraint
  /// @param i: ConstraintKeeper index
  /// @return false when inserted a duplicate (should not happen)
  template <class Constraint>
  bool MapInsert(const Constraint& con, int i) {
    return MPD( MapInsert__Impl(con, i) );
  }


  /// MapFind__Impl.
  /// Default version for functional constraints with a map.
  /// @param con: constraint reference
  /// @return constraint index, or -1
  template <class Constraint>
  int MapFind__Impl(const Constraint& con) {
    const auto& map = GET_CONST_CONSTRAINT_MAP(Constraint);
    auto it = map.find( con );
    return (map.end() != it) ? it->second : -1;
  }

  /// MapInsert__Impl.
  /// Default version for functional constraints with a map.
  /// @param con: the constraint
  /// @param i: ConstraintKeeper index
  /// @return false when inserted a duplicate (should not happen)
  template <class Constraint>
  bool MapInsert__Impl(const Constraint& con, int i) {
    auto& map = GET_CONSTRAINT_MAP(Constraint);
    auto result = map.insert( { con, i } );
    return result.second;
  }


  /////////////////////// CONSTRAINT CONVERTERS /////////////////////////
  /// Constraint keepers and converters should be initialized after \a presolver_

  /// Convert linear range constraints, if not accepted by ModelAPI
  INSTALL_ITEM_CONVERTER(RangeLinearConstraintConverter)
  /// Convert quadratic range constraints, if necessary
  INSTALL_ITEM_CONVERTER(RangeQuadraticConstraintConverter)

  /// Convert quadratic cones, if necessary
  INSTALL_ITEM_CONVERTER(QConeConverter)
  /// Convert rotated quadratic cones, if necessary
  INSTALL_ITEM_CONVERTER(RQConeConverter)


public:
  /// ValuePresolve link copying values 1:1 between model items
  pre::CopyLink& GetCopyLink() { return copy_link_; }

  /// ValuePresolve link copying values 1:many
  pre::One2ManyLink& GetOne2ManyLink() { return one2many_link_; }

  /// ValuePresolve link copying values many:1
  pre::Many2OneLink& GetMany2OneLink() { return many2one_link_; }
};


/// A 'final' flat converter in a CRTP hierarchy
template <template <typename, typename, typename> class FlatCvt,
          class Backend, class Model = FlatModel< > >
class FlatCvtImpl :
    public FlatCvt<FlatCvtImpl<FlatCvt, Backend, Model>, Backend, Model> {
public:
  /// Base type
  using Base = FlatCvt<FlatCvtImpl<FlatCvt, Backend, Model>, Backend, Model>;

  /// Construct
  FlatCvtImpl(Env& e) : Base(e) { }
};

} // namespace mp

#endif // CONVERTER_FLAT_H
