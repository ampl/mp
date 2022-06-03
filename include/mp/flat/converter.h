#ifndef CONVERTER_FLAT_H
#define CONVERTER_FLAT_H

#include <unordered_map>
#include <map>
#include <cstdio>
#include <cmath>
#include <utility>

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

namespace mp {

/// FlatConverter: preprocesses and manages flat constraints.
/// Such constraints might need to be converted to others, which is
/// handled by overloaded methods in derived classes.
/// @param Impl: the final CRTP class
/// @param ModelAPI: the solver's model API wrapper
/// @param FlatModel: internal representation of a flat model
template <class Impl, class ModelAPI,
          class FlatModel = BasicFlatModel< > >
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
  /// Fix a resulting variable of a logical expression as true
  /// and propagate positive ctx.
  /// TODO avoid creating resvar for root logical constraints
  void FixAsTrue(int resvar) {
    PropagateResultOfInitExpr(resvar, 1.0, 1.0, +Context());
  }


public:
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

  /// From am affine expression:
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


protected:
  //////////////////////////// CUSTOM CONSTRAINTS CONVERSION ////////////////////////////
  ///
  //////////////////////////// THE CONVERSION LOOP: BREADTH-FIRST ///////////////////////
  void ConvertItems() {
    try {
      MP_DISPATCH( ConvertAllConstraints() );
      // MP_DISPATCH( PreprocessIntermediate() );     // preprocess after each level
      MP_DISPATCH( ConvertMaps() );
      MP_DISPATCH( PreprocessFinal() );               // final prepro
    } catch (const ConstraintConversionFailure& cff) {
      MP_RAISE(cff.message());
    }
  }

  void ConvertAllConstraints() {
    GetModel(). ConvertAllConstraints(*this);
  }

  /// Default map conversions. Currently empty
  void ConvertMaps() { }

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
  /// Assume mixed context if not set in the constraint
  /// TODO Make sure context is always propagated for all constraints and objectives
  template <class Constraint>
  void RunConversion(const Constraint& con, int i) {
    if (con.HasContext())
      if (con.GetContext().IsNone())
        con.SetContext(Context::CTX_MIX);
    MP_DISPATCH(Convert(con, i));
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
          std::string("Convertion of constraint type '") +
            Constraint::GetTypeName() + "' not implemented");
  }

  //////////////////////////// SOME SPECIFIC CONSTRAINT CONVERTERS
  /// ///////////////////////////////////// ///////////////////////////

  /// If backend does not like LDC, we can redefine it
  void Convert(const LinearFunctionalConstraint& ldc) {
    this->AddConstraint(ldc.to_linear_constraint());
  }

  /// If backend does not like QDC, we can redefine it
  void Convert(const QuadraticFunctionalConstraint& qdc) {
    this->AddConstraint(qdc.to_quadratic_constraint());
  }


public:
  /// ADD CUSTOM CONSTRAINT
  ///
  /// Use only for non-mapped constraints. For functional constraints
  /// stored __WITH_MAP, use AssignResult(Var)2Args().
  /// TODO non-functional constraints __WITH_MAP.
  /// Takes ownership.
  /// @return Node reference for the stored constraint
  template <class Constraint>
  pre::NodeRange AddConstraint(Constraint con) {
    auto node_range =
        AddConstraintAndTryNoteResultVariable( std::move(con) );
    return node_range;
  }

  template <class Constraint>
  const Constraint& GetConstraint(int i) const {
    return GET_CONST_CONSTRAINT_KEEPER(Constraint).GetConstraint(i);
  }

protected:
  USE_BASE_MAP_FINDERS( BaseConverter )

  template <class Constraint>
  pre::NodeRange AddConstraintAndTryNoteResultVariable(Constraint&& con) {
    const auto resvar = con.GetResultVar();
    auto& ck = GET_CONSTRAINT_KEEPER( Constraint );
    auto i = ck.AddConstraint(std::move(con));
    ConInfo ci{&ck, i};
    if (resvar>=0)
      AddInitExpression(resvar, ci);
    /// Can also cache non-functional constraints
    if (! MP_DISPATCH( MapInsert(
                         MPD(template GetConstraint<Constraint>(i)), i ) ))
      MP_RAISE("Trying to MapInsert() duplicated constraint: " +
                             ck.GetDescription());
    return ck.GetValueNode().Select(i);
  }


public:
  void StartModelInput() { }

  void FinishModelInput() {
    MPD( ConvertModel() );
    if (relax())              // TODO value presolve link?
      GetModel().RelaxIntegrality();
    GetModel().PushModelTo(GetModelAPI());
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
    if (GetEnv().verbose_mode())
      GetEnv().PrintWarnings();
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

  std::vector<int> common_exprs_;               // variables equal to the result

public:

  //////////////////////////// CREATE OR FIND A FIXED VARIABLE //////////////////////////////
  pre::NodeRange MakeFixedVar(double value) {
    auto it = map_fixed_vars_.find(value);
    if (map_fixed_vars_.end()!=it)
      return GetVarValueNode().Select( it->second );
    auto v = MPD( DoAddVar(value, value) );
    map_fixed_vars_[value] = v;
    return GetVarValueNode().Select( v );
  }

  /// Create or find a fixed variable
  pre::NodeRange AddVar(double lb, double ub, var::Type type = var::CONTINUOUS) {
    if (lb!=ub)
      return DoAddVar(lb, ub, type);
    return MakeFixedVar(lb);
  }

  /// Add several variables
  /// @return value node range for them
  pre::NodeRange AddVars(const typename BaseFlatModel::VarBndVec& lbs,
               const typename BaseFlatModel::VarBndVec& ubs,
               const typename BaseFlatModel::VarTypeVec& types) {
    assert(0==BaseFlatModel::num_vars());                     // allow this only once
    BaseFlatModel::AddVars__basic(lbs, ubs, types);
    return GetVarValueNode().Add( lbs.size() );
  }

  /// Reuse Presolver's target nodes for all variables
  pre::ValueNode& GetVarValueNode()
  { return GetValuePresolver().GetTargetNodes().GetVarValues().MakeSingleKey(); }

  /// Constraint type's Value Node
  template <class Constraint>
  pre::ValueNode& GetValueNode(Constraint*)
  { return GET_CONSTRAINT_KEEPER(Constraint).GetValueNode(); }

public:
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
  /// Shortcut lb(var)
  void set_var_lb(int var, double lb) { this->GetModel().set_lb(var, lb); }
  /// Shortcut ub(var)
  void set_var_ub(int var, double ub) { this->GetModel().set_ub(var, ub); }

  /// Narrow variable domain range
  void NarrowVarBounds(int var, double lb, double ub) {
    auto& m = GetModel();
    m.set_lb(var, std::max(m.lb(var), lb));
    m.set_ub(var, std::min(m.ub(var), ub));
    if (m.lb(var)>m.ub(var))             // TODO write .sol, report .iis
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
    AffineExpr ae({{-1.0}, {bvar}}, 1.0); // TODO use map / FCC?
    return MP_DISPATCH( Convert2Var(std::move(ae)) );
  }

  /// Typedef ConInfo
  using ConInfo = AbstractConstraintLocation;

  /// Add variable. Type: var::CONTINUOUS by default
  pre::NodeRange DoAddVar(double lb=MinusInfty(), double ub=Infty(),
             var::Type type = var::CONTINUOUS) {
    int v = GetModel().AddVar__basic(lb, ub, type);
    return GetVarValueNode().Select( v );
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
  void AddInitExpression(int var, const ConInfo& vi) {
    var_info_.resize(std::max(var_info_.size(), (size_t)var+1));
    var_info_[var] = vi;
  }

public:
  /// Variable has an init expr?
  bool HasInitExpression(int var) const {
    return int(var_info_.size())>var && var_info_[var].HasId();
  }

  /// Get the init expr
  const ConInfo& GetInitExpression(int var) const {
    assert(HasInitExpression(var));
    return var_info_[var];
  }

  /// The internal flat model type
  using ModelType = FlatModel;
  /// The internal flat model object, const ref
  const ModelType& GetModel() const { return *this; }
  /// The internal flat model object, ref
  ModelType& GetModel() { return *this; }


  ///////////////////////////////////////////////////////////////////////
  /////////////////////// OPTIONS /////////////////////////
  ///
private:
  struct Options {
    int preprocessAnything_ = 1;
    int preprocessEqualityResultBounds_ = 1;
    int preprocessEqualityBvar_ = 1;
    int relax_ = 0;
  };
  Options options_;

protected:
  int relax() const { return options_.relax_; }

public:
  /// Init FlatConverter options
  void InitOptions() {
    InitOwnOptions();
    GetModelAPI().InitOptions();
  }

protected:
  const mp::OptionValueInfo values_relax_[2] = {
    {     "0", "No (default)", 0 },
    {     "1", "Yes: treat integer and binary variables as continuous.", 1}
  };

private:
  void InitOwnOptions() {
    GetEnv().AddOption("cvt:pre:all",
        "0/1*: Set to 0 to disable most presolve in the flat converter.",
        options_.preprocessAnything_, 0, 1);
    GetEnv().AddOption("cvt:pre:eqresult",
        "0/1*: Preprocess reified equality comparison's boolean result bounds.",
        options_.preprocessEqualityResultBounds_, 0, 1);
    GetEnv().AddOption("cvt:pre:eqbinary",
        "0/1*: Preprocess reified equality comparison with a binary variable.",
        options_.preprocessEqualityBvar_, 0, 1);
    GetEnv().AddOption("alg:relax relax",
        "0*/1: Whether to relax integrality of variables.",
        options_.relax_, 0, 1);
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


public:
  /// for tests. TODO make friends
  using ModelAPIType = ModelAPI;

  /// AddWarning. Strings should remain valid
  void AddWarning(const char* key, const char* msg) {
    GetEnv().AddWarning(key, msg);
  }


private:
  ModelAPIType modelapi_;      // We store modelapi in the converter for speed
  pre::ValuePresolver value_presolver_;
  pre::CopyLink copy_link_ { GetValuePresolver() };


protected:
  /////////////////////// CONSTRAINT KEEPERS /////////////////////////
  /// Constraint keepers and converters should be initialized after \a presolver_

  /// Define constraint keepers for all constraint types
  STORE_CONSTRAINT_TYPE__NO_MAP(LinConRange)
  STORE_CONSTRAINT_TYPE__NO_MAP(LinConLE)
  STORE_CONSTRAINT_TYPE__NO_MAP(LinConEQ)
  STORE_CONSTRAINT_TYPE__NO_MAP(LinConGE)

  STORE_CONSTRAINT_TYPE__NO_MAP(QuadConRange)
  STORE_CONSTRAINT_TYPE__NO_MAP(QuadConLE)
  STORE_CONSTRAINT_TYPE__NO_MAP(QuadConEQ)
  STORE_CONSTRAINT_TYPE__NO_MAP(QuadConGE)

  /// TODO Use FunctionalConstraintConverter with LFC, QFC
  STORE_CONSTRAINT_TYPE__WITH_MAP(LinearFunctionalConstraint)
  STORE_CONSTRAINT_TYPE__WITH_MAP(QuadraticFunctionalConstraint)

  /// With maps we store flattened NL expressions
  /// (functional constraints)
  STORE_CONSTRAINT_TYPE__WITH_MAP(MaxConstraint)
  STORE_CONSTRAINT_TYPE__WITH_MAP(MinConstraint)
  STORE_CONSTRAINT_TYPE__WITH_MAP(AbsConstraint)
  STORE_CONSTRAINT_TYPE__WITH_MAP(AndConstraint)
  STORE_CONSTRAINT_TYPE__WITH_MAP(OrConstraint)

  STORE_CONSTRAINT_TYPE__WITH_MAP(CondLinConEQ)
  STORE_CONSTRAINT_TYPE__WITH_MAP(CondLinConLE)
  STORE_CONSTRAINT_TYPE__WITH_MAP(CondLinConLT)
  STORE_CONSTRAINT_TYPE__WITH_MAP(CondLinConGE)
  STORE_CONSTRAINT_TYPE__WITH_MAP(CondLinConGT)

  STORE_CONSTRAINT_TYPE__WITH_MAP(CondQuadConEQ)
  STORE_CONSTRAINT_TYPE__WITH_MAP(CondQuadConLE)
  STORE_CONSTRAINT_TYPE__WITH_MAP(CondQuadConLT)
  STORE_CONSTRAINT_TYPE__WITH_MAP(CondQuadConGE)
  STORE_CONSTRAINT_TYPE__WITH_MAP(CondQuadConGT)

  STORE_CONSTRAINT_TYPE__WITH_MAP(NotConstraint)
  STORE_CONSTRAINT_TYPE__WITH_MAP(DivConstraint)
  STORE_CONSTRAINT_TYPE__WITH_MAP(IfThenConstraint)
  STORE_CONSTRAINT_TYPE__WITH_MAP(AllDiffConstraint)
  STORE_CONSTRAINT_TYPE__WITH_MAP(NumberofConstConstraint)
  STORE_CONSTRAINT_TYPE__WITH_MAP(NumberofVarConstraint)
  STORE_CONSTRAINT_TYPE__WITH_MAP(CountConstraint)

  STORE_CONSTRAINT_TYPE__WITH_MAP(ExpConstraint)
  STORE_CONSTRAINT_TYPE__WITH_MAP(ExpAConstraint)
  STORE_CONSTRAINT_TYPE__WITH_MAP(LogConstraint)
  STORE_CONSTRAINT_TYPE__WITH_MAP(LogAConstraint)
  STORE_CONSTRAINT_TYPE__WITH_MAP(PowConstraint)
  STORE_CONSTRAINT_TYPE__WITH_MAP(SinConstraint)
  STORE_CONSTRAINT_TYPE__WITH_MAP(CosConstraint)
  STORE_CONSTRAINT_TYPE__WITH_MAP(TanConstraint)

  /// No maps for static constraints
  STORE_CONSTRAINT_TYPE__NO_MAP(IndicatorConstraintLinLE)
  STORE_CONSTRAINT_TYPE__NO_MAP(IndicatorConstraintLinEQ)
  STORE_CONSTRAINT_TYPE__NO_MAP(IndicatorConstraintLinGE)
  STORE_CONSTRAINT_TYPE__NO_MAP(IndicatorConstraintQuadLE)
  STORE_CONSTRAINT_TYPE__NO_MAP(IndicatorConstraintQuadEQ)
  STORE_CONSTRAINT_TYPE__NO_MAP(IndicatorConstraintQuadGE)
  STORE_CONSTRAINT_TYPE__NO_MAP(PLConstraint)
  STORE_CONSTRAINT_TYPE__NO_MAP(SOS1Constraint)
  STORE_CONSTRAINT_TYPE__NO_MAP(SOS2Constraint)
  STORE_CONSTRAINT_TYPE__NO_MAP(ComplementarityLinear)
  STORE_CONSTRAINT_TYPE__NO_MAP(ComplementarityQuadratic)

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

  /// MapInsert.
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

public:
  /// Presolve link copying values between model items
  pre::CopyLink& GetCopyLink() { return copy_link_; }
};


/// A 'final' flat converter in a CRTP hierarchy
template <template <typename, typename, typename> class FlatCvt,
          class Backend, class Model = BasicFlatModel< > >
class FlatCvtImpl :
    public FlatCvt<FlatCvtImpl<FlatCvt, Backend, Model>, Backend, Model> {
public:
  using Base = FlatCvt<FlatCvtImpl<FlatCvt, Backend, Model>, Backend, Model>;
  FlatCvtImpl(Env& e) : Base(e) { }
};

} // namespace mp

#endif // CONVERTER_FLAT_H
