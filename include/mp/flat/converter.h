#ifndef CONVERTER_FLAT_H
#define CONVERTER_FLAT_H

#include <unordered_map>
#include <map>
#include <cstdio>
#include <cmath>
#include <utility>
#include <cassert>

#include "mp/env.h"
#include "mp/format.h"
#include "mp/solver-base.h"
#include "mp/flat/converter_model.h"
#include "mp/flat/convert_functional.h"
#include "mp/flat/constr_keeper.h"
#include "mp/flat/constr_std.h"
#include "mp/flat/expr_bounds.h"
#include "mp/flat/constr_prepro.h"
#include "mp/flat/constr_prop_down.h"
#include "mp/valcvt.h"
#include "mp/flat/redef/std/range_con.h"
#include "mp/flat/redef/conic/cones.h"
#include "mp/flat/redef/conic/qcones2qc.h"
#include "mp/utils-file.h"
#include "mp/ampls-ccallbacks.h"

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
    public BasicFlatConverter,
    public FlatModel,
    public BoundComputations<Impl>,
    public ConstraintPreprocessors<Impl>,
    public ConstraintPropagatorsDown<Impl>,
    public EnvKeeper
{
public:
  /// Class name
  static const char* GetTypeName() { return "FlatConverter"; }

  /// Construct with Env&
  FlatConverter(Env& e) : EnvKeeper(e), modelapi_(e) { }

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
    PropagateResultOfInitExpr(resvar, 1.0, 1.0, +Context());
  }


public:
  /// Reverse propagate result variable of an expression
  void PropagateResultOfInitExpr(int var, Context ctx) {
    PropagateResultOfInitExpr(var, lb(var), ub(var), ctx);
  }

  /// Reverse propagate result variable of an expression
  void PropagateResultOfInitExpr(int var, double lb, double ub, Context ctx) {
    NarrowVarBounds(var, lb, ub);
    if (HasInitExpression(var)) {
      const auto& ckid = GetInitExpression(var);
      ckid.GetCK()->PropagateResult(*this, ckid.GetIndex(), lb, ub, ctx);
    }
  }


public:
  //////////////////////////////////// VISITOR ADAPTERS /////////////////////////////////////////
  /// These are called to transform expressions, either by FlatCvt itself,
  /// or when flattening NL model

  /// From an affine expression:
  /// Adds a result variable r and constraint r == expr
  int Convert2Var(AffineExpr&& ee) {
    if (ee.is_variable())
      return ee.get_representing_variable();
    if (ee.is_constant())
      return MakeFixedVar(ee.constant_term());
    return AssignResultVar2Args(
            LinearFunctionalConstraint(std::move(ee)));
  }

  /// From a quadratic expression:
  /// Adds a result variable r and constraint r == expr
  int Convert2Var(QuadraticExpr&& ee) {
    if (ee.is_variable())
      return ee.get_representing_variable();
    if (ee.is_constant())
      return MakeFixedVar(ee.constant_term());
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
  template <class FuncConstraint>
  typename FCC<Impl, FuncConstraint>::Var
  AssignResultVar2Args(FuncConstraint&& fc) {
    auto vc = AssignResult2Args(std::move(fc));
    if (vc.is_const())
      return MPD( MakeFixedVar(vc.get_const()) );
    return vc.get_var();
  }

	/// Typedef ConInfo; constraint location
	using ConInfo = AbstractConstraintLocation;

	/// Replace functional expression defining a given variable.
  template <class FuncConstraint>
  void RedefineVariable(int res_var, FuncConstraint&& fc) {
    assert( MPD( HasInitExpression(res_var) ) );
		auto ci_old = MPD( GetInitExpression(res_var) );
    fc.SetResultVar(res_var);
		// If this expression exists, use it.
		// TODO make sure any new context is re-converted
		// if necessary.
		auto i = MPD( MapFind(fc) );
		if (i<0)
			i = MPD( AddConstraint(std::move(fc)) );
		auto& ck = GET_CONSTRAINT_KEEPER( FuncConstraint );
    ConInfo ci{&ck, i};
    ReplaceInitExpression(res_var, ci);
    MarkAsBridged(ci_old);
  }


	/// Variables' reference counting ///////////////////////////////////
	/// Currently only for defined variables ////////////////////////////

	/// Use "+1" a variable
	void IncrementVarUsage(int v) {
		++VarUsageRef(v);
	}

	/// Unuse result variable.
  /// Actually this is to 'unuse' the init expression
  /// - might change naming.
	/// Throw if already not used.
	void DecrementVarUsage(int v) {
		assert(VarUsageRef(v)>0);
		if (! (--VarUsageRef(v))) {
			if (HasInitExpression(v))
        MarkAsUnused(GetInitExpression(v));
		}
	}

	/// Fix unused defined vars.
	/// Normally should delete them.
	void FixUnusedDefinedVars() {
		for (auto i=num_vars(); i--; ) {
			if (HasInitExpression(i) &&
					! VarUsageRef(i)) {
				set_var_lb(i, 0.0);      // fix to 0
				set_var_ub(i, 0.0);
			}
		}
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

  //////////////////////////// CUSTOM CONSTRAINTS CONVERSION ////////////////////////////
  ///
  //////////////////////////// THE CONVERSION LOOP: BREADTH-FIRST ///////////////////////
  void ConvertItems() {
    try {
			MPD( Convert2Cones(); );                 // sweep before other conversions
      MP_DISPATCH( ConvertAllConstraints() );
      // MP_DISPATCH( PreprocessIntermediate() );     // preprocess after each level
      constr_depth_ = 1;  // Workaround. TODO have maps as special constraints
			MP_DISPATCH( ConvertMaps() );
      MP_DISPATCH( PreprocessFinal() );               // final prepro
    } catch (const ConstraintConversionFailure& cff) {
      MP_RAISE(cff.message());
    }
  }

  void OpenGraphExporter() {
    if (graph_export_file().size()) {
      if (!graph_exporter_app_->Open(graph_export_file().c_str(), true))
        MP_RAISE("Failed to open the graph export file.");
      value_presolver_.SetExport(true);
    }
  }

	/// Offload the conic logic to a functor
	void Convert2Cones() {
		conic_cvt_.Run();
	}

  void ConvertAllConstraints() {
    GetModel().ConvertAllConstraints(*this);
  }

  /// Default map conversions. Currently empty
  void ConvertMaps() { }

  void CloseGraphExporter() {
    value_presolver_.FinishExportingLinkEntries();
    graph_exporter_app_->Close();
  }

  //////////////////////// WHOLE-MODEL PREPROCESSING /////////////////////////
  void PreprocessIntermediate() { }
  void PreprocessFinal() { }


  //////////////////////////// CONSTRAINT PROPAGATORS ///////////////////////////////////

  /// Allow FCC to access Preprocess methods
  template <class Impl1, class Converter, class Constraint>
  friend class BasicFCC;

  //////////////////////////// SPECIFIC CONSTRAINT RESULT-TO-ARGUMENTS PROPAGATORS //////
  /// Currently we should propagate to all arguments, be it always the CTX_MIX.

  /// Allow ConstraintKeeper to PropagateResult(), use GetBackend() etc
  template <class , class , class >
  friend class ConstraintKeeper;


  //////////////////////////// CUSTOM CONSTRAINTS CONVERSION ////////////////////////////
  ///
public: // for ConstraintKeeper
  /// RunConversion() of a constraint:
  /// Assume mixed context if not set.
  template <class Constraint>
  void RunConversion(const Constraint& con, int i, int depth) {
    constr_depth_ = depth+1;
    if (con.UsesContext())              // If context relevant,
      if (con.GetContext().IsNone())    // ensure we have context, mixed if none
        con.SetContext(Context::CTX_MIX);
    pre::AutoLinkScope<Impl> auto_link_scope{
      *static_cast<Impl*>(this),
      GET_CONSTRAINT_KEEPER(Constraint).SelectValueNodeRange(i)
    };
    MP_DISPATCH(Convert(con, i));
  }

  /// Query if a constraint type
  /// is natively accepted by the solver.
  /// The parameter is only needed for type.
  template <class Con>
  ConstraintAcceptanceLevel GetConstraintAcceptance(Con* ) const {
    return GET_CONST_CONSTRAINT_KEEPER(Con).GetChosenAcceptanceLevel();
  }

  /// Query the nuber of addable constraints of type.
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

  /// Check whether ModelAPI accepts and recommends the constraint
  template <class Constraint>
  static bool ModelAPIAcceptsAndRecommends(const Constraint* pcon) {
    return ConstraintAcceptanceLevel::Recommended ==
        ModelAPI::AcceptanceLevel(pcon);
  }

  /// Generic adapter for old non-bridged Convert() methods
  ///
  /// New way is to use the \a i parameter for bridging
  template <class Constraint>
  void Convert(const Constraint& con, int ) {
    MPD( Convert(con) );
  }

  /// By default, we complain about someone trying to convert an unknown constraint
  template <class Constraint>
  void Convert(const Constraint& ) {
    MP_RAISE(
          std::string("Constraint type '") +
            Constraint::GetTypeName() +
            "' is neither accepted by '" +
            ModelAPI::GetTypeName() +
            "', nor is conversion implemented");
  }

  //////////////////////////// SOME SPECIFIC CONSTRAINT CONVERTERS
  /// ///////////////////////////////////// ///////////////////////////

  /// If backend does not like LFC, we redefine it here
  void Convert(const LinearFunctionalConstraint& ldc) {
    MPD( AddConstraint(ldc.to_linear_constraint()) );
  }

  /// If backend does not like QFC, we redefine it
  void Convert(const QuadraticFunctionalConstraint& qdc) {
    qdc.AddQuadraticConstraint(*(Impl*)this);
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
    GetModel().AddObjective( std::move(qo) );
    /// Temporarily removing AutoLinking for objectives
    // return AutoLink( GetObjValueNode().Add() );
    return GetObjValueNode().Select(-1);
  }

  /// ADD CUSTOM CONSTRAINT, does not propagate result
  /// (use AddConstraint_AS_ROOT() otherwise).
  ///
  /// Use only for non-mapped constraints. For functional constraints
  /// stored __WITH_MAP, use AssignResult(Var)2Args().
  /// Non-functional constraints cannot be unified currently.
  /// Takes ownership.
  /// @return Node reference for the stored constraint
  template <class Constraint>
  pre::NodeRange AddConstraint(Constraint con) {
    auto node_range =
        AddConstraintAndTryNoteResultVariable( std::move(con) );
    return AutoLink( node_range );
  }

  /// ADD CUSTOM CONSTRAINT and propagate root-ness
  /// (use AddConstraint() otherwise).
  ///
  /// Use only for non-mapped constraints. For functional constraints
  /// stored __WITH_MAP, use AssignResult(Var)2Args().
  /// Takes ownership.
  /// @return Node reference for the stored constraint
  template <class Constraint>
  pre::NodeRange AddConstraint_AS_ROOT(Constraint con) {
    MPD( PropagateResult(con) );
    return AddConstraint( std::move(con) );
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
		return GET_CONST_CONSTRAINT_KEEPER(Constraint).GetConstraint(i);
	}

  /// Mark constraint as reformulated
  void MarkAsBridged(const ConInfo& ci) {
		ci.GetCK()->MarkAsBridged(ci.GetIndex());
	}

  /// Mark constraint as unused
  void MarkAsUnused(const ConInfo& ci) {
    ci.GetCK()->MarkAsUnused(ci.GetIndex());
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
    return ck.SelectValueNodeRange(i);
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
		FixUnusedDefinedVars();       // Until we have proper var deletion
    CheckLinearCons();
    PresolveNames();
    GetModel().PushModelTo(GetModelAPI());
    MPD( CloseGraphExporter() );
    if (value_presolver_.GetExport())
      assert( value_presolver_.AllEntriesExported() );
//    if (GetEnv().verbose_mode())
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

  /// Presolve item names
  void PresolveNames() {
    if (var_names_.size()) {
      /// Check that constr / obj names are present too?
      GetValuePresolver().CleanUpNameNodes();
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
        + fmi->GetNumberOfConstraints(typeid(ComplementarityQuadratic));
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


  /// Check unpostsolved solution
  /// in various ways
  bool CheckSolution(
      ArrayRef<double> x,
      const pre::ValueMapDbl& duals,
      ArrayRef<double> obj) {
    bool result = true;
    std::string err_msg;
    try {            // protect
      std::vector<double> x_back; // to extract x used
      if (options_.solcheckmode_ & (1+2+4+8+16)) {
        if (!DoCheckSol(x, duals, obj, {}, x_back, false))
          result = false;
      }
      if (options_.solcheckmode_ & (32+64+128+256+512)) {
        auto x1 = RecomputeAuxVars(x);
        if (!DoCheckSol(x1, duals, obj, x_back, x_back, true))
          result = false;
      }
    } catch (const mp::Error& err) {
      if (options_.solcheckfail_)
        throw;
      err_msg = err.what();
    } catch (const std::exception& exc) {
      err_msg = exc.what();
    } catch (...) {
      err_msg = "unknown error";
    }
    if (err_msg.size()) {
      if (options_.solcheckfail_)
        MP_RAISE("Solution check aborted: " + err_msg);
      AddWarning("Solution check aborted", err_msg);
      result = false;
    }
    return result;
  }

  /// Functor to recompute auxiliary var \a i
  VarsRecomputeFn recomp_fn
  = [this](int i, const VarInfoRecomp& x) {
    if (MPCD( HasInitExpression(i) )) {
      const auto& iexpr = MPCD( GetInitExpression(i) );
      auto pCK = iexpr.GetCK();
      auto resvar = pCK->GetResultVar(iexpr.GetIndex());
      assert(resvar == i);
      if (!pCK->IsUnused(iexpr.GetIndex()))
        return pCK->ComputeValue(iexpr.GetIndex(), x);
    }
    return x.get_x().get_x()[i];  // no recomputation
  };

  /// Recompute auxiliary variables
  ArrayRef<double> RecomputeAuxVars(ArrayRef<double> x) {
    VarInfoRecomp vir {
      options_.solfeastol_,
          true,        // currently not relevant for recomputation
      {x, recomp_fn},
      {},              // now raw values
      GetModel().var_type_vec(),
          GetModel().var_lb_vec(),
          GetModel().var_ub_vec(),
          nullptr, nullptr
    };
    vir.get_x().set_p_var_info(&vir);
    for (auto i=vir.size(); i--; )
      vir[i];         // touch the variable to be recomputed
    return std::move(vir.get_x().get_x());
  }

  /// Check single unpostsolved solution
  bool DoCheckSol(
      ArrayRef<double> x,
      const pre::ValueMapDbl& duals,
      ArrayRef<double> obj,
      ArrayRef<double> x_raw,
      std::vector<double>& x_back,
      bool if_recomp_vals) {
    SolCheck chk(x, duals, obj, x_raw,
                 GetModel().var_type_vec(),
                 GetModel().var_lb_vec(),
                 GetModel().var_ub_vec(),
                 options_.solfeastol_,
                 options_.solinttol_,
                 options_.dont_use_sol_round_
                   ? "" : std::getenv("solution_round"),
                 options_.dont_use_sol_prec_
                   ? "" : std::getenv("solution_precision"),
                 if_recomp_vals,
                 if_recomp_vals
                 ? (options_.solcheckmode_ >> 5)
                 : (options_.solcheckmode_));  // for raw values
    if (chk.check_mode() & 1)
      CheckVars(chk);
    if (chk.check_mode() & (2+4+8))
      CheckCons(chk);
    if (chk.check_mode() & 16)
      CheckObjs(chk);
    GenerateViolationsReport(chk);
    // What if this is an intermediate solution?
    // Should be fine - warning by default,
    // fail if requested explicitly.
    // If warning, we should add the report
    // to that solutions' solve message, and
    // a summary in the final solve message.
    // For now, do this via warnings?
    if (chk.HasAnyViols()) {
      if (options_.solcheckfail_)
        MP_RAISE_WITH_CODE(520,    // numeric error
                           chk.GetReport());
      else
        AddWarning(
              GetEnv().GetSolCheckWarningKey(if_recomp_vals),
              chk.GetReport(),
              true);  // replace for multiple solutions
    }
    x_back = chk.x_ext().get_x();
    return !chk.HasAnyViols();
  }

  void CheckVars(SolCheck& chk) {
    for (auto i=num_vars(); i--; ) {
      auto x = chk.x(i);
      bool aux = !MPCD( is_var_original(i) );
      // no aux vars for idealistic mode
      if (!aux || !chk.if_recomputed()) {
        chk.VarViolBnds().at(aux).CheckViol(
              MPCD( lb(i) ) - x,
              options_.solfeastol_,
              GetModel().var_name(i));
        chk.VarViolBnds().at(aux).CheckViol(
              x - MPCD( ub(i) ),
              options_.solfeastol_,
              GetModel().var_name(i));
        if (is_var_integer(i))
          chk.VarViolIntty().at(aux).CheckViol(
                std::fabs(x - std::round(x)),
                options_.solinttol_,
                GetModel().var_name(i));
      }
    }
  }

  /// Includes logical constraints.
  void CheckCons(SolCheck& chk) {
    GetModel().ComputeViolations(chk);
  }

  void CheckObjs(SolCheck& chk) {
    const auto& objs = GetModel().get_objectives();
    // Solvers might have dummy obj.
    // Unbounded problems might have no obj value.
    for (auto i
         =std::min(objs.size(), chk.obj_vals().size());
         i--; ) {
      chk.ObjViols().CheckViol(
            std::fabs(chk.obj_vals()[i]
                      - ComputeValue(objs[i], chk.x_ext())),
            options_.solfeastol_,
            objs[i].name());
    }
  }

  void GenerateViolationsReport(SolCheck& chk) {
    fmt::MemoryWriter wrt;
    if (chk.HasAnyConViols()) {
      wrt.write(
            "Constraint violations\n"
            "    (sol:chk:feastol={}, sol:chk:inttol={},\n"
            "     solution_round='{}', solution_precision='{}'):\n",
            options_.solfeastol_, options_.solinttol_,
            chk.x_ext().solution_round(),
            chk.x_ext().solution_precision());
      Gen1Viol(chk.VarViolBnds().at(0), wrt,
               "  - {} original variable(s) violate bounds,\n"
               "      max by {}");
      Gen1Viol(chk.VarViolBnds().at(1), wrt,
               "  - {} auxiliary variable(s) violate bounds,\n"
               "      max by {}");
      Gen1Viol(chk.VarViolIntty().at(0), wrt,
               "  - {} original variable(s) violate integrality,\n"
               "      max by {}");
      Gen1Viol(chk.VarViolIntty().at(1), wrt,
               "  - {} auxiliary variable(s) violate integrality,\n"
               "      max by {}");
    }
    GenConViol(chk.ConViolAlg(), wrt, "Algebraic");
    GenConViol(chk.ConViolLog(), wrt, "Logical");
    if (chk.HasAnyObjViols()) {
      wrt.write("Objective value violations"
                "    (sol:chk:feastol={},\n"
                "     solution_round='{}', solution_precision='{}'):\n",
                options_.solfeastol_,
                chk.x_ext().solution_round(),
                chk.x_ext().solution_precision());
      Gen1Viol(chk.ObjViols(), wrt,
               "  - {} objective value(s) violated,\n      max by {}");
    }
    chk.SetReport( wrt.str() );
  }

  void Gen1Viol(
      const ViolSummary& vs, fmt::MemoryWriter& wrt,
      const std::string& format) {
    if (vs.N_) {
      wrt.write(format, vs.N_, vs.epsMax_);
      if (vs.name_ && *vs.name_ != '\0')
        wrt.write(" (item '{}')", vs.name_);
      wrt.write("\n");
    }
  }

  void GenConViol(
      const std::map< std::string, ViolSummArray<3> >& cvmap,
      fmt::MemoryWriter& wrt, const std::string& classnm) {
    if (cvmap.size()) {
      wrt.write(classnm + " constraint violations:\n");
      for (const auto& cva: cvmap) {
        Gen1Viol(cva.second.at(0), wrt,
                 "  - {} original constraint(s) of type '"
                 + std::string(cva.first)
                 + "' violate bounds,\n      max by {}");
        Gen1Viol(cva.second.at(1), wrt,
                 "  - {} intermediate auxiliary constraint(s) of type '"
                 + std::string(cva.first)
                 + "' violate bounds,\n      max by {}");
        Gen1Viol(cva.second.at(2), wrt,
                 "  - {} final auxiliary constraint(s) of type '"
                 + std::string(cva.first)
                 + "' violate bounds,\n      max by {}");
      }
    }
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
    map_fixed_vars_[value] = v;
    return GetVarValueNode().Select( v );  // no autolink, done in DoAddVar()
  }

  /// Create or find a fixed variable
  pre::NodeRange AddVar(double lb, double ub, var::Type type = var::CONTINUOUS) {
    if (lb!=ub)
      return DoAddVar(lb, ub, type);
    return MakeFixedVar(lb);
  }

  /// Add several variables once.
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
  /// lb_array()
  template <class VarArray>
  double lb_array(const VarArray& va) const
  { return this->GetModel().lb_array(va); }
  /// ub_array()
  template <class VarArray>
  double ub_array(const VarArray& va) const
  { return this->GetModel().ub_array(va); }
  /// Set lb(var)
  void set_var_lb(int var, double lb) { this->GetModel().set_lb(var, lb); }
  /// Set ub(var)
  void set_var_ub(int var, double ub) { this->GetModel().set_ub(var, ub); }
  /// Set lb(var), propagate context if functional result
  void set_var_lb_context(int var, double lb, Context ctx) {
    set_var_lb(var, lb);
    PropagateResultOfInitExpr(var, ctx);
  }
  /// Set ub(var), propagate context
  void set_var_ub_context(int var, double ub, Context ctx) {
    set_var_ub(var, ub);
    PropagateResultOfInitExpr(var, ctx);
  }
  /// Set bounds(var), propagate context
  void set_var_bounds_context(int var, double lb, double ub, Context ctx) {
    set_var_lb(var, lb);
    set_var_ub(var, ub);
    PropagateResultOfInitExpr(var, ctx);
  }

  /// Narrow variable domain range
  void NarrowVarBounds(int var, double lb, double ub) {
    auto& m = GetModel();
    m.set_lb(var, std::max(m.lb(var), lb));
    m.set_ub(var, std::min(m.ub(var), ub));
    if (m.lb(var)>m.ub(var))
      MP_INFEAS("empty variable domain");
  }

  /// var_type()
  var::Type var_type(int var) const { return this->GetModel().var_type(var); }
  /// is_fixed()
  bool is_fixed(int var) const { return this->GetModel().is_fixed(var); }
  /// fixed_value()
  double fixed_value(int var) const
  { assert(is_fixed(var)); return this->GetModel().fixed_value(var); }

  /// MakeComplementVar()
  int MakeComplementVar(int bvar) {
    if (! (lb(bvar)==0.0 && ub(bvar)==1.0) )
      MP_RAISE("Asked to complement variable with bounds "
                             + std::to_string(lb(bvar)) + ".." + std::to_string(ub(bvar)));
    AffineExpr ae({{-1.0}, {bvar}}, 1.0);
    return MP_DISPATCH( Convert2Var(std::move(ae)) );
  }

  /// Add vector of variables. Type: var::CONTINUOUS by default
  /// @return vector of the Ids of the new vars
  std::vector<int> AddVars_returnIds(std::size_t nvars,
                           double lb=MinusInfty(), double ub=Infty(),
                           var::Type type = var::CONTINUOUS) {
    std::vector<int> newVars(nvars);
    for (std::size_t  i=0; i<nvars; ++i)
      newVars[i] = AddVar(lb, ub, type);
    return newVars;
  }

  bool is_var_integer(int var) const
  { return MPCD( GetModel() ).is_integer_var(var); }


private:
  std::vector<ConInfo> var_info_;


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
		return var_info_.at(var);
  }

	/// Get the init expression pointer.
	/// @return nullptr if no init expr or not this type
	template <class ConType>
	const ConType* GetInitExpressionOfType(int var) {
		if (MPCD( HasInitExpression(var) )) {
			auto ci0 = MPCD( GetInitExpression(var) );
			if (IsConInfoType<ConType>(ci0)) {
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

  /////////////////////// AUTO LINKING ////////////////////////////

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
    }
    return nr;
  }

  /// Whether we should auto-link new items
  bool DoingAutoLinking() const
  { return auto_link_src_item_.IsValid(); }

  /// Turn off auto-linking for current conversion
  void TurnOffAutoLinking() {
    auto_link_src_item_.Invalidate();
    auto_link_targ_items_.clear();
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
  static bool ModelAPIAcceptsQC() {
    return ModelAPIAcceptsAndRecommends(
          (const QuadConLE*)nullptr);   // if accepts QuadConLE
  }

  /// Whether the ModelAPI accepts nonconvex QC
  static bool ModelAPIAcceptsNonconvexQC() {
    return ModelAPI::AcceptsNonconvexQC();
  }

  /// Whether the ModelAPI accepts quadratic cones
	int ModelAPIAcceptsQuadraticCones() {
		return
				std::max(
					(int)GetConstraintAcceptance((QuadraticConeConstraint*)nullptr),
					(int)GetConstraintAcceptance((RotatedQuadraticConeConstraint*)nullptr));
	}

	/// Whether the ModelAPI accepts exp cones
	int ModelAPIAcceptsExponentialCones() {
		return
				(int)GetConstraintAcceptance((ExponentialConeConstraint*)nullptr);
	}


private:
  struct Options {
    std::string file_graph_export_;
    int preprocessAnything_ = 1;
    int preprocessEqualityResultBounds_ = 1;
    int preprocessEqualityBvar_ = 1;

    int passQuadObj_ = ModelAPIAcceptsQuadObj();
		int passQuadCon_ = ModelAPIAcceptsQC();
		int passSOCPCones_ = 0;
		int passExpCones_ = 0;

    int relax_ = 0;

    int solcheckmode_ = 1+2+16+512;
    bool solcheckfail_ = false;
    double solfeastol_ = 1e-6;
    double solinttol_ = 1e-5;
    bool dont_use_sol_round_ = false;
    bool dont_use_sol_prec_ = false;
    int solchkoutlev_ = 0;
  };
  Options options_;


protected:
  /// Graph export file
  const std::string& graph_export_file() const
  { return options_.file_graph_export_; }

  /// Whether we should relax integrality
  int relax() const { return options_.relax_; }


public:
  /// Init FlatConverter options
  void InitOptions() {
    InitOwnOptions();
    GetModelAPI().InitCustomOptions();
  }


private:
  void InitOwnOptions() {
    /// Should be called after adding all constraint keepers
    FlatModel::ConsiderAcceptanceOptions(*this, GetModelAPI(), GetEnv());

    GetEnv().AddStoredOption("tech:writegraph writegraph exportgraph",
        "File to export conversion graph. Format: JSON Lines.",
        options_.file_graph_export_);
    GetEnv().AddOption("cvt:pre:all",
        "0/1*: Set to 0 to disable most presolve in the flat converter.",
        options_.preprocessAnything_, 0, 1);
    GetEnv().AddOption("cvt:pre:eqresult",
        "0/1*: Preprocess reified equality comparison's boolean result bounds.",
        options_.preprocessEqualityResultBounds_, 0, 1);
    GetEnv().AddOption("cvt:pre:eqbinary",
        "0/1*: Preprocess reified equality comparison with a binary variable.",
        options_.preprocessEqualityBvar_, 0, 1);
    GetEnv().AddOption("cvt:quadobj passquadobj",
                       ModelAPIAcceptsQuadObj() ?
        "0/1*: Multiply out and pass quadratic objective terms to the solver, "
        "vs. linear approximation."
                       :
        "0*/1: Multiply out and pass quadratic objective terms to the solver, "
                         "vs. linear approximation.",
        options_.passQuadObj_, 0, 1);
    GetEnv().AddOption("cvt:quadcon passquadcon",
                       ModelAPIAcceptsQC() ?
        "0/1*: Multiply out and pass quadratic constraint terms to the solver, "
        "vs. linear approximation."
                       :
        "0*/1: Multiply out and pass quadratic constraint terms to the solver, "
                         "vs. linear approximation.",
        options_.passQuadCon_, 0, 1);
		if (ModelAPIAcceptsExponentialCones())
			GetEnv().AddOption("cvt:expcones expcones",
												 ModelAPIAcceptsExponentialCones()>1 ?
														 "0/1*: Recognize exponential cones." :
														 "0*/1: Recognize exponential cones.",
												 options_.passExpCones_, 0, 1);
		options_.passExpCones_ = ModelAPIAcceptsExponentialCones()>1;
		if (ModelAPIAcceptsQuadraticCones())
			GetEnv().AddOption("cvt:socp passsocp socp",
												 ModelAPIAcceptsQuadraticCones()>1 ?
                           "0/1*: Recognize quadratic cones vs passing them "
                           "as pure quadratic constraints." :
                           "0*/1: Recognize quadratic cones vs passing them "
                           "as pure quadratic constraints.",
					options_.passSOCPCones_, 0, 1);
		options_.passSOCPCones_ = ModelAPIAcceptsQuadraticCones()>1;
		GetEnv().AddOption("alg:relax relax",
        "0*/1: Whether to relax integrality of variables.",
        options_.relax_, 0, 1);
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
                             "Default: 1+2+16+512.",
        options_.solcheckmode_, 0, 1024);
    GetEnv().AddOption("sol:chk:feastol sol:chk:eps sol:eps chk:eps",
        "Solution checking tolerance for objective values, variable "
        "and constraint bounds. Default 1e-6.",
        options_.solfeastol_, 0.0, 1e100);
    GetEnv().AddOption("sol:chk:inttol sol:chk:inteps sol:inteps chk:inteps",
        "Solution checking tolerance for variables' integrality. "
        "Default 1e-5.",
        options_.solinttol_, 0.0, 1e100);
    GetEnv().AddOption("sol:chk:fail chk:fail checkfail",
        "Fail on constraint violations.",
        options_.solcheckfail_, false, true);
    GetEnv().AddOption("sol:chk:noround chk:noround chk:no_solution_round",
        "Don't use AMPL solution_round option when checking.",
        options_.dont_use_sol_round_, false, true);
    GetEnv().AddOption("sol:chk:noprec chk:noprec chk:no_solution_precision",
        "Don't use AMPL solution_precision option when checking.",
        options_.dont_use_sol_prec_, false, true);
  }


public:
  /// Wrapper about a specific preprocess option:
  /// checks whether \a preprocessAnything_ is on.
  bool CanPreprocess(int f) const {
    return 0!=options_.preprocessAnything_ && 0!=f;
  }

  /// Whether preprocess equality result bounds
  bool IfPreproEqResBounds() const
  { return MPCD( CanPreprocess(options_.preprocessEqualityResultBounds_) ); }

  /// Whether preprocess conditional equality of a binary variable
  bool IfPreproEqBinVar() const
  { return MPCD( CanPreprocess(options_.preprocessEqualityBvar_) ); }

  /// Whether we pass quad obj terms to the solver without linearization
  bool IfPassQuadObj() const { return options_.passQuadObj_; }

  /// Whether we pass quad con terms to the solver without linearization
  bool IfPassQuadCon() const { return options_.passQuadCon_; }

  /// Whether to quadratize pow(..., const_pos_int).
  /// The fact that we use the passQuadCon_ flag
  /// is much Gurobi-biased: v9.5 does not PL-linearize Pow
  /// for negative arguments
  bool IfQuadratizePowConstPosIntExp() const
  { return options_.passQuadCon_; }

	/// Whether we pass SOCP cones
	bool IfPassSOCPCones() const { return options_.passSOCPCones_; }

	/// Whether we pass exp cones
	bool IfPassExpCones() const { return options_.passExpCones_; }


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


private:
  /// We store ModelApi in the converter for speed.
  /// Should be before constraints
  ModelAPIType modelapi_;
  /// Conversion graph exporter file appender
  std::unique_ptr<BasicFileAppender> graph_exporter_app_{MakeFileAppender()};
  /// Conversion graph exporter functor
  pre::ValuePresolver::ExporterFn graph_exporter_fn_{
    [this](const char* s){
      graph_exporter_app_->Append(s);
    }
  };
  /// ValuePresolver: should be init before constraint keepers
  /// and links
  pre::ValuePresolver value_presolver_
  {
    GetModel(), GetEnv(), graph_exporter_fn_,
        [this](
        ArrayRef<double> x,
        const pre::ValueMapDbl& y,
        ArrayRef<double> obj) -> bool {
      return !this->options_.solcheckmode_   // not desired
          || this->CheckSolution(x, y, obj);
    }
  };
  pre::CopyLink copy_link_ { GetValuePresolver() }; // the copy links
  pre::One2ManyLink one2many_link_ { GetValuePresolver() }; // the 1-to-many links
  pre::NodeRange auto_link_src_item_;   // the source item for autolinking
  std::vector<pre::NodeRange> auto_link_targ_items_;

	ConicConverter<Impl> conic_cvt_ { *static_cast<Impl*>(this) };

	std::vector<int> refcnt_vars_;
  int constr_depth_ = 0;    // tree depth of new constraints

  /// AMPL item names
  std::vector<std::string>
  var_names_,
  con_names_,   // no SOS here, they go directly into the SOS
  obj_names_;


protected:
  /////////////////////// CONSTRAINT KEEPERS /////////////////////////
  /// Constraint keepers and converters should be initialized after
  ///  \a value_presolver_

  /// Define constraint keepers for all constraint types.
  /// No maps for static constraints.
  /// 2nd parameter: solver options for this constraint,
  /// in case it is accepted by the solver natively and
  /// is convertible by us.
  STORE_CONSTRAINT_TYPE__NO_MAP(LinConRange,
                                "acc:linrange acc:linrng")
  STORE_CONSTRAINT_TYPE__NO_MAP(LinConLE, "acc:linle")
  STORE_CONSTRAINT_TYPE__NO_MAP(LinConEQ, "acc:lineq")
  STORE_CONSTRAINT_TYPE__NO_MAP(LinConGE, "acc:linge")

  STORE_CONSTRAINT_TYPE__NO_MAP(QuadConRange,
                                "acc:quadrange acc:quadrng")
  STORE_CONSTRAINT_TYPE__NO_MAP(QuadConLE, "acc:quadle")
  STORE_CONSTRAINT_TYPE__NO_MAP(QuadConEQ, "acc:quadeq")
  STORE_CONSTRAINT_TYPE__NO_MAP(QuadConGE, "acc:quadge")

  /// Our own functional constraints: LFC, QFC
  STORE_CONSTRAINT_TYPE__WITH_MAP(
      LinearFunctionalConstraint, "acc:linfunccon")
  STORE_CONSTRAINT_TYPE__WITH_MAP(
      QuadraticFunctionalConstraint, "acc:quadfunccon")

  /// Flattened NL expressions
  STORE_CONSTRAINT_TYPE__WITH_MAP(MaxConstraint, "acc:max")
  STORE_CONSTRAINT_TYPE__WITH_MAP(MinConstraint, "acc:min")
  STORE_CONSTRAINT_TYPE__WITH_MAP(AbsConstraint, "acc:abs")
  STORE_CONSTRAINT_TYPE__WITH_MAP(AndConstraint,
                                  "acc:and acc:forall")
  STORE_CONSTRAINT_TYPE__WITH_MAP(OrConstraint,
                                  "acc:or acc:exists")

  STORE_CONSTRAINT_TYPE__WITH_MAP(CondLinConEQ, "acc:condlineq")
  STORE_CONSTRAINT_TYPE__WITH_MAP(CondLinConLE, "acc:condlinle")
  STORE_CONSTRAINT_TYPE__WITH_MAP(CondLinConLT, "acc:condlinlt")
  STORE_CONSTRAINT_TYPE__WITH_MAP(CondLinConGE, "acc:condlinge")
  STORE_CONSTRAINT_TYPE__WITH_MAP(CondLinConGT, "acc:condlingt")

  STORE_CONSTRAINT_TYPE__WITH_MAP(CondQuadConEQ, "acc:condquadeq")
  STORE_CONSTRAINT_TYPE__WITH_MAP(CondQuadConLE, "acc:condquadle")
  STORE_CONSTRAINT_TYPE__WITH_MAP(CondQuadConLT, "acc:condquadlt")
  STORE_CONSTRAINT_TYPE__WITH_MAP(CondQuadConGE, "acc:condquadge")
  STORE_CONSTRAINT_TYPE__WITH_MAP(CondQuadConGT, "acc:condquadgt")

  STORE_CONSTRAINT_TYPE__WITH_MAP(NotConstraint, "acc:not")
  STORE_CONSTRAINT_TYPE__WITH_MAP(DivConstraint, "acc:div")
  STORE_CONSTRAINT_TYPE__WITH_MAP(IfThenConstraint, "acc:ifthen")
  STORE_CONSTRAINT_TYPE__WITH_MAP(ImplicationConstraint, "acc:impl")
  STORE_CONSTRAINT_TYPE__WITH_MAP(AllDiffConstraint, "acc:alldiff")
  STORE_CONSTRAINT_TYPE__WITH_MAP(NumberofConstConstraint,
                                  "acc:numberofconst")
  STORE_CONSTRAINT_TYPE__WITH_MAP(NumberofVarConstraint,
                                  "acc:numberofvar")
  STORE_CONSTRAINT_TYPE__WITH_MAP(CountConstraint, "acc:count")

  STORE_CONSTRAINT_TYPE__WITH_MAP(ExpConstraint, "acc:exp")
  STORE_CONSTRAINT_TYPE__WITH_MAP(ExpAConstraint, "acc:expa acc:expA")
  STORE_CONSTRAINT_TYPE__WITH_MAP(LogConstraint, "acc:log")
  STORE_CONSTRAINT_TYPE__WITH_MAP(LogAConstraint, "acc:loga acc:logA")
  STORE_CONSTRAINT_TYPE__WITH_MAP(PowConstraint, "acc:pow")
  STORE_CONSTRAINT_TYPE__WITH_MAP(SinConstraint, "acc:sin")
  STORE_CONSTRAINT_TYPE__WITH_MAP(CosConstraint, "acc:cos")
  STORE_CONSTRAINT_TYPE__WITH_MAP(TanConstraint, "acc:tan")
  STORE_CONSTRAINT_TYPE__WITH_MAP(AsinConstraint, "acc:asin")
  STORE_CONSTRAINT_TYPE__WITH_MAP(AcosConstraint, "acc:acos")
  STORE_CONSTRAINT_TYPE__WITH_MAP(AtanConstraint, "acc:atan")
  STORE_CONSTRAINT_TYPE__WITH_MAP(SinhConstraint, "acc:sinh")
  STORE_CONSTRAINT_TYPE__WITH_MAP(CoshConstraint, "acc:cosh")
  STORE_CONSTRAINT_TYPE__WITH_MAP(TanhConstraint, "acc:tanh")
  STORE_CONSTRAINT_TYPE__WITH_MAP(AsinhConstraint, "acc:asinh")
  STORE_CONSTRAINT_TYPE__WITH_MAP(AcoshConstraint, "acc:acosh")
  STORE_CONSTRAINT_TYPE__WITH_MAP(AtanhConstraint, "acc:atanh")

  /// No maps for static constraints
  STORE_CONSTRAINT_TYPE__NO_MAP(
      IndicatorConstraintLinLE, "acc:indle acc:indlinle")
  STORE_CONSTRAINT_TYPE__NO_MAP(
      IndicatorConstraintLinEQ, "acc:indeq acc:indlineq")
  STORE_CONSTRAINT_TYPE__NO_MAP(
      IndicatorConstraintLinGE, "acc:indge acc:indlinge")
  STORE_CONSTRAINT_TYPE__NO_MAP(
      IndicatorConstraintQuadLE, "acc:indquadle")
  STORE_CONSTRAINT_TYPE__NO_MAP(
      IndicatorConstraintQuadEQ, "acc:indquadeq")
  STORE_CONSTRAINT_TYPE__NO_MAP(
      IndicatorConstraintQuadGE, "acc:indquadge")
  STORE_CONSTRAINT_TYPE__NO_MAP(PLConstraint,
                                "acc:pl acc:pwl acc:piecewise")
  STORE_CONSTRAINT_TYPE__NO_MAP(SOS1Constraint, "acc:sos1")
  STORE_CONSTRAINT_TYPE__NO_MAP(SOS2Constraint, "acc:sos2")
  STORE_CONSTRAINT_TYPE__NO_MAP(
      ComplementarityLinear, "acc:compl acc:compllin")
  STORE_CONSTRAINT_TYPE__NO_MAP(
      ComplementarityQuadratic, "acc:complquad")
  STORE_CONSTRAINT_TYPE__NO_MAP(
      QuadraticConeConstraint, "acc:quadcone")
  STORE_CONSTRAINT_TYPE__NO_MAP(
      RotatedQuadraticConeConstraint, "acc:rotatedquadcone")
  STORE_CONSTRAINT_TYPE__NO_MAP(
      PowerConeConstraint, "acc:powercone")
  STORE_CONSTRAINT_TYPE__NO_MAP(
      ExponentialConeConstraint, "acc:expcone")
  STORE_CONSTRAINT_TYPE__NO_MAP(
      GeometricConeConstraint, "acc:geomcone")


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
};


/// A 'final' flat converter in a CRTP hierarchy
template <template <typename, typename, typename> class FlatCvt,
          class Backend, class Model = FlatModel< > >
class FlatCvtImpl :
    public FlatCvt<FlatCvtImpl<FlatCvt, Backend, Model>, Backend, Model> {
public:
  using Base = FlatCvt<FlatCvtImpl<FlatCvt, Backend, Model>, Backend, Model>;
  FlatCvtImpl(Env& e) : Base(e) { }
};

} // namespace mp

#endif // CONVERTER_FLAT_H
