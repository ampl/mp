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
#include "mp/flat/convert_functional.h"
#include "mp/flat/constraint_keeper.h"
#include "mp/flat/constraints_std.h"
#include "mp/flat/model.h"
#include "mp/presolve.h"
#include "mp/flat/redef/std/range_con.h"

namespace mp {

/// FlatConverter: preprocesses and manages flat constraints.
/// Such constraints might need to be converted to others, which is
/// handled by overloaded methods in derived classes
/// @param Impl: the final CRTP class
/// @param Backend: the solver's model API wrapper
/// @param Model: internal representation of a flat model
template <class Impl, class Backend,
          class Model = BasicFlatModel< > >
class FlatConverter :
    public BasicFlatConverter,
    public Model,
    public EnvKeeper
{
public:
  static constexpr const char* name() { return "FlatConverter"; };

  using Var = typename Model::Var;
  static constexpr Var VoidVar() { return Model::VoidVar(); }

public:
  using VarArray = std::vector<int>;

protected:
  using ClassType = FlatConverter<Impl, Backend, Model>;
  using BaseConverter = BasicFlatConverter;
  using BaseFlatModel = Model;

public:
  static const char* GetConverterName() { return "FlatConverter"; }

  FlatConverter(Env& e) : EnvKeeper(e), backend_(e) { }


  //////////////////////////// CONVERTERS OF STANDRAD MP ITEMS //////////////////////////////
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////
public:
  /// Fix a resulting variable of a logical expression as true
  /// and propagate +Ctx
  /// TODO avoid creating resvar for root logical constraints
  void FixAsTrue(int resvar) {
    PropagateResultOfInitExpr(resvar, 1.0, 1.0, +Context());
  }

protected:
  void PropagateResultOfInitExpr(int var, double lb, double ub, Context ctx) {
    NarrowVarBounds(var, lb, ub);
    if (HasInitExpression(var)) {
      const auto& ckid = GetInitExpression(var);
      ckid.GetCK()->PropagateResult(*this, ckid.GetIndex(), lb, ub, ctx);
    }
  }

  void NarrowVarBounds(int var, double lb, double ub) {
    auto& m = GetModel();
    m.set_lb(var, std::max(m.lb(var), lb));
    m.set_ub(var, std::min(m.ub(var), ub));
    if (m.lb(var)>m.ub(var))             // TODO write .sol, report .iis
      throw std::logic_error("infeasibility: empty variable domain");
  }


public:
  //////////////////////////////////// VISITOR ADAPTERS /////////////////////////////////////////

  /// From an expression:
  /// Adds a result variable r and constraint r == expr
  int Convert2Var(QuadExp&& ee) {
    if (ee.is_variable())
      return ee.get_representing_variable();
    if (ee.is_constant())
      return MakeFixedVar(ee.constant_term());
    PreprocessInfoStd bnt = ComputeBoundsAndType(ee);
    auto r = MP_DISPATCH( AddVar(bnt.lb_, bnt.ub_, bnt.type_) );
    if (ee.is_affine())
      AddConstraint(LinearFunctionalConstraint(r, std::move(ee.GetAE())));
    else
      AddConstraint(QuadraticFunctionalConstraint(r, std::move(ee)));
    return r;
  }

  PreprocessInfoStd ComputeBoundsAndType(const QuadExp& ee) {
    auto bntAE = ComputeBoundsAndType(ee.GetAE());
    auto bntQT = ComputeBoundsAndType(ee.GetQT());
    return AddBoundsAndType(bntAE, bntQT);
  }

  PreprocessInfoStd ComputeBoundsAndType(const AffExp& ae) {
    PreprocessInfoStd result;
    result.lb_ = result.ub_ = ae.constant_term();    // TODO reuse bounds if supplied
    result.type_ = is_integer(result.lb_) ? var::INTEGER : var::CONTINUOUS;
    result.linexp_type_ = var::INTEGER;
    auto& model = MP_DISPATCH( GetModel() );
    for (auto i=ae.size(); i--; ) {
      auto v = ae.var(i);
      auto c = ae.coef(i);
      if (c >= 0.0) {
        result.lb_ += c * model.lb(v);
        result.ub_ += c * model.ub(v);
      } else {
        result.lb_ += c * model.ub(v);
        result.ub_ += c * model.lb(v);
      }
      if (var::INTEGER!=model.var_type(v) || !is_integer(c)) {
        result.type_=var::CONTINUOUS;
        result.linexp_type_=var::CONTINUOUS;
      }
    }
    return result;
  }

  PreprocessInfoStd ComputeBoundsAndType(const QuadTerms& qt) {
    PreprocessInfoStd result;
    result.lb_ = result.ub_ = 0.0;
    result.type_ = var::INTEGER;
    result.linexp_type_ = var::INTEGER;
    auto& model = MP_DISPATCH( GetModel() );
    for (auto i=qt.size(); i--; ) {
      auto coef = qt.coef(i);
      auto v1 = qt.var1(i);
      auto v2 = qt.var2(i);
      auto prodBnd = ProductBounds(v1, v2);
      if (coef >= 0.0) {
        result.lb_ += coef * prodBnd.first;
        result.ub_ += coef * prodBnd.second;
      } else {
        result.lb_ += coef * prodBnd.second;
        result.ub_ += coef * prodBnd.first;
      }
      if (var::INTEGER!=model.var_type(v1) ||
          var::INTEGER!=model.var_type(v2) ||
          !is_integer(coef)) {
        result.type_=var::CONTINUOUS;
        result.linexp_type_=var::CONTINUOUS;
      }
    }
    return result;
  }

  template <class Var>
  std::pair<double, double> ProductBounds(Var x, Var y) const {
    auto lx=lb(x), ly=lb(y), ux=ub(x), uy=ub(y);
    std::array<double, 4> pb{lx*ly, lx*uy, ux*ly, ux*uy};
    return {*std::min_element(pb.begin(), pb.end()), *std::max_element(pb.begin(), pb.end())};
  }

  PreprocessInfoStd AddBoundsAndType(const PreprocessInfoStd& bnt1,
                                     const PreprocessInfoStd& bnt2) {
    return {bnt1.lb()+bnt2.lb(), bnt1.ub()+bnt2.ub(),
      var::INTEGER==bnt1.type() && var::INTEGER==bnt2.type() ?
            var::INTEGER : var::CONTINUOUS};
  }

  /// Take FuncConstraint with arguments
  /// If the result of the function is known, return it
  /// Otherwise, create a result variable and add the constraint
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
      throw std::logic_error(cff.message());
    }
  }

  void ConvertAllConstraints() {
    GetModel(). ConvertAllConstraints(*this);
  }

  //////////////////////// WHOLE-MODEL PREPROCESSING /////////////////////////
  void PreprocessIntermediate() { }
  void ConvertMaps() { }
  void PreprocessFinal() { }


  //////////////////////////// CONSTRAINT PROPAGATORS ///////////////////////////////////

  /// Allow FCC to access Preprocess methods
  template <class Impl1, class Converter, class Constraint>
  friend class BasicFCC;

  template <class PreprocessInfo>
  void PreprocessConstraint(
      MinimumConstraint& c, PreprocessInfo& prepro) {
    auto& m = MP_DISPATCH( GetModel() );
    auto& args = c.GetArguments();
    prepro.narrow_result_bounds( m.lb_array(args),
                          m.ub_min_array(args) );
    prepro.set_result_type( m.common_type(args) );
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      MaximumConstraint& c, PreprocessInfo& prepro) {
    auto& m = MP_DISPATCH( GetModel() );
    auto& args = c.GetArguments();
    prepro.narrow_result_bounds( m.lb_max_array(args),
                          m.ub_array(args) );
    prepro.set_result_type( m.common_type(args) );
  }

  /// When the result variable is set,
  /// the constraint is skipped
  template <class PreprocessInfo>
  void PreprocessConstraint(
      AbsConstraint& c, PreprocessInfo& prepro) {
    const auto argvar = c.GetArguments()[0];
    const auto lb = this->lb(argvar),
        ub = this->ub(argvar);
    if (lb>=0.0) {  // When result var is set, \a c is skipped
      prepro.set_result_var(argvar);
      return;
    } else if (ub<=0.0) {
      auto res = AssignResult2Args(   // create newvar = -argvar
            LinearFunctionalConstraint({ {{-1.0}, {argvar}}, 0.0 }));
      prepro.set_result_var(res.get_var());
      return;
    }
    prepro.narrow_result_bounds(0.0, std::max(-lb, ub));
    prepro.set_result_type( var_type(argvar) );
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      EQ0Constraint& c, PreprocessInfo& prepro) {
    prepro.narrow_result_bounds(0.0, 1.0);
    prepro.set_result_type( var::INTEGER );
    if (0!=CanPreprocess( options_.preprocessEqualityResultBounds_ ))
      if (FixEqualityResult(c, prepro))
        return;
    PreprocessEqVarConst__unifyCoef(c);
    if (0!=CanPreprocess( options_.preprocessEqualityBvar_ ))
      if (ReuseEqualityBinaryVar(c, prepro))
        return;
  }

  template <class PreprocessInfo>
  bool FixEqualityResult(
      EQ0Constraint& c, PreprocessInfo& prepro) {
    AffExp& ae = c.GetArguments();
    if (ae.is_constant()) {               // const==0
      auto res = (double)int(0.0==ae.constant_term());
      /// TODO this depends on context???
      prepro.narrow_result_bounds(res, res);
      return true;
    }
    auto bndsNType = ComputeBoundsAndType(ae);
    if (bndsNType.lb() > 0.0 || bndsNType.ub() < 0.0) {
      /// TODO this depends on context???
      prepro.narrow_result_bounds(0.0, 0.0);
      return true;
    }
    if (bndsNType.lb()==0.0 && bndsNType.ub()==0.0) {
      /// TODO this depends on context???
      prepro.narrow_result_bounds(1.0, 1.0);
      return true;
    }
    if (var::INTEGER==bndsNType.linexp_type_ &&
        !is_integer(ae.constant_term())) {
      /// TODO this depends on context???
      prepro.narrow_result_bounds(0.0, 0.0);
      return true;
    }
    return false;
  }

  static void PreprocessEqVarConst__unifyCoef(EQ0Constraint& c) {
    AffExp& ae = c.GetArguments();
    if (1==ae.size()) {
      const double c = ae.coef(0);
      if (1.0!=c) {
        assert(0.0!=std::fabs(c));
        ae.constant_term(ae.constant_term() / c);
        ae.set_coef(0, 1.0);
      }
    }
  }

  template <class PreprocessInfo>
  bool ReuseEqualityBinaryVar(
      EQ0Constraint& c, PreprocessInfo& prepro) {
    auto& m = MP_DISPATCH( GetModel() );
    AffExp& ae = c.GetArguments();
    if (1==ae.size()) {                           // var==const
      assert( 1.0==ae.coef(0) );
      int var = ae.var(0);
      if (m.is_binary_var(var)) {            // See if this is binary var==const
        const double rhs = -ae.constant_term();
        if (1.0==rhs)
          prepro.set_result_var( var );
        else if (0.0==rhs)
          prepro.set_result_var( MakeComplementVar(var) );
        else
          prepro.narrow_result_bounds(0.0, 0.0);    // not 0/1 value, result false
        return true;
      }
    }
    return false;
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      LE0Constraint& , PreprocessInfo& prepro) {
    prepro.narrow_result_bounds(0.0, 1.0);
    prepro.set_result_type( var::INTEGER );
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      ConjunctionConstraint& , PreprocessInfo& prepro) {
    prepro.narrow_result_bounds(0.0, 1.0);
    prepro.set_result_type( var::INTEGER );
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      DisjunctionConstraint& , PreprocessInfo& prepro) {
    prepro.narrow_result_bounds(0.0, 1.0);
    prepro.set_result_type( var::INTEGER );
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      AllDiffConstraint& , PreprocessInfo& prepro) {
    prepro.narrow_result_bounds(0.0, 1.0);
    prepro.set_result_type( var::INTEGER );
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      NumberofConstConstraint& con, PreprocessInfo& prepro) {
    prepro.narrow_result_bounds(0.0, con.GetArguments().size());
    prepro.set_result_type( var::INTEGER );
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      NumberofVarConstraint& con, PreprocessInfo& prepro) {
    prepro.narrow_result_bounds(0.0, con.GetArguments().size());
    prepro.set_result_type( var::INTEGER );
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      CountConstraint& con, PreprocessInfo& prepro) {
    prepro.narrow_result_bounds(0.0, con.GetArguments().size());
    prepro.set_result_type( var::INTEGER );
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      NotConstraint& , PreprocessInfo& prepro) {
    prepro.narrow_result_bounds(0.0, 1.0);
    prepro.set_result_type( var::INTEGER );
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      IfThenConstraint& c, PreprocessInfo& prepro) {
    const auto& args = c.GetArguments();
    prepro.narrow_result_bounds(std::min(lb(args[1]), lb(args[2])),
        std::max(ub(args[1]), ub(args[2])));
    prepro.set_result_type( MP_DISPATCH(GetModel()).
                            common_type( { args[1], args[2] } ) );
  }

  ////////////////////// NONLINEAR FUNCTIONS //////////////////////
  template <class PreprocessInfo>
  void PreprocessConstraint(
      ExpConstraint& , PreprocessInfo& prepro) {
    prepro.narrow_result_bounds(0.0, this->Infty());
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      ExpAConstraint& , PreprocessInfo& prepro) {
    prepro.narrow_result_bounds(0.0, this->Infty());
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      LogConstraint& c, PreprocessInfo& ) {
    NarrowVarBounds(c.GetArguments()[0], 0.0, this->Infty());
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      LogAConstraint& c, PreprocessInfo& ) {
    NarrowVarBounds(c.GetArguments()[0], 0.0, this->Infty());
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      SinConstraint& , PreprocessInfo& prepro) {
    prepro.narrow_result_bounds(-1.0, 1.0);
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      CosConstraint& , PreprocessInfo& prepro) {
    prepro.narrow_result_bounds(-1.0, 1.0);
  }



  //////////////////////////// CUSTOM CONSTRAINTS //////////////////////
  ///
  //////////////////////////// SPECIFIC CONSTRAINT RESULT-TO-ARGUMENTS PROPAGATORS //////
  /// Currently we should propagate to all arguments, be it always the CTX_MIX.

  /// Allow ConstraintKeeper to PropagateResult()
  template <class , class , class >
  friend class ConstraintKeeper;

  /// By default, declare mixed context
  template <class Constraint>
  void PropagateResult(Constraint& con, double lb, double ub, Context ctx) {
    internal::Unused(&con, lb, ub, &ctx);
    con.SetContext(Context::CTX_MIX);
    for (const auto a: con.GetArguments())
      PropagateResultOfInitExpr(a, this->MinusInfty(), this->Infty(),
                                Context::CTX_MIX);
  }

  void PropagateResult(LinearFunctionalConstraint& con, double lb, double ub,
                       Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    for (auto v: con.GetAffineExpr().vars())
      PropagateResultOfInitExpr(v, this->MinusInfty(), this->Infty(),
                                Context::CTX_MIX);
  }

  void PropagateResult(QuadraticFunctionalConstraint& con, double lb, double ub,
                       Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    const auto& args = con.GetArguments();
    for (auto v: args.GetAE().vars())
      PropagateResultOfInitExpr(v, this->MinusInfty(), this->Infty(),
                                Context::CTX_MIX);
    const auto& qt = args.GetQT();
    for (auto i=qt.size(); i--; ) {
      PropagateResultOfInitExpr(qt.var1(i), this->MinusInfty(), this->Infty(),
                                Context::CTX_MIX);
      PropagateResultOfInitExpr(qt.var2(i), this->MinusInfty(), this->Infty(),
                                Context::CTX_MIX);
    }
  }

  void PropagateResult(QuadraticConstraint& con, double lb, double ub,
                       Context ctx) {
    internal::Unused(lb, ub, ctx);
    for (const auto v: con.vars())
      PropagateResultOfInitExpr(v, this->MinusInfty(), this->Infty(),
                                Context::CTX_MIX);
    for (const auto v: con.GetQPTerms().vars1())
      PropagateResultOfInitExpr(v, this->MinusInfty(), this->Infty(),
                                Context::CTX_MIX);
    for (const auto v: con.GetQPTerms().vars2())
      PropagateResultOfInitExpr(v, this->MinusInfty(), this->Infty(),
                                Context::CTX_MIX);
  }

  template <int sens>
  void PropagateResult(IndicatorConstraint< LinConRhs<sens> >& con,
                       double lb, double ub, Context ctx) {
    internal::Unused(lb, ub, ctx);
    PropagateResultOfInitExpr(con.get_binary_var(),
                              this->MinusInfty(), this->Infty(),
                              Context::CTX_MIX);
    for (const auto v: con.get_constraint().vars())
      PropagateResultOfInitExpr(v, this->MinusInfty(), this->Infty(),
                                Context::CTX_MIX);
  }

  template <int type>
  void PropagateResult(SOS_1or2_Constraint<type>& con, double lb, double ub, Context ) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    for (const auto v: con.get_vars())
      PropagateResultOfInitExpr(v, this->MinusInfty(), this->Infty(), Context::CTX_MIX);
  }

  void PropagateResult(NotConstraint& con, double lb, double ub, Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    for (const auto a: con.GetArguments())
      PropagateResultOfInitExpr(a, 1.0-ub, 1.0-lb, -ctx);
  }

  void PropagateResult(ConjunctionConstraint& con, double lb, double ub, Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    for (const auto a: con.GetArguments())
      PropagateResultOfInitExpr(a, lb, 1.0, +ctx);
  }

  void PropagateResult(DisjunctionConstraint& con, double lb, double ub, Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    for (const auto a: con.GetArguments())
      PropagateResultOfInitExpr(a, 0.0, ub, +ctx);
  }

  void PropagateResult(IfThenConstraint& con, double lb, double ub, Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    auto& args = con.GetArguments();
    PropagateResultOfInitExpr(args[0], 0.0, 1.0, Context::CTX_MIX);
    PropagateResultOfInitExpr(args[1], this->MinusInfty(), this->Infty(), +ctx);
    PropagateResultOfInitExpr(args[2], this->MinusInfty(), this->Infty(), -ctx);
  }

  void PropagateResult(AllDiffConstraint& con, double lb, double ub, Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    // TODO go into arguments
  }

  void PropagateResult(NumberofConstConstraint& con, double lb, double ub, Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    // TODO go into arguments
  }

  void PropagateResult(NumberofVarConstraint& con, double lb, double ub, Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    // TODO go into arguments
  }

  void PropagateResult(LE0Constraint& con, double lb, double ub, Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
  }

  void PropagateResult(EQ0Constraint& con, double lb, double ub, Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    auto& args = con.GetArguments().vars();
    for (auto v: args) {
      PropagateResultOfInitExpr(v, this->MinusInfty(), this->Infty(),
                                Context::CTX_MIX);
    }
  }


  //////////////////////////// CUSTOM CONSTRAINTS CONVERSION ////////////////////////////
  ///
  //////////////////////////// SPECIFIC CONSTRAINT CONVERTERS ///////////////////////////

  USE_BASE_CONSTRAINT_CONVERTERS(BasicFlatConverter);      // reuse default converters

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

  /// If backend does not like LDC, we can redefine it
  void Convert(const LinearFunctionalConstraint& ldc) {
    this->AddConstraint(ldc.to_linear_constraint());
  }

  /// If backend does not like QDC, we can redefine it
  void Convert(const QuadraticFunctionalConstraint& qdc) {
    this->AddConstraint(qdc.to_quadratic_constraint());
  }


public:
  //////////////////////// ADD CUSTOM CONSTRAINT ///////////////////////
  //////////////////////// Takes ownership /////////////////////////////
  template <class Constraint>
  pre::NodeRange AddConstraint(Constraint&& con) {
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
    ConstraintLocation<Constraint> cl{&ck, i};
    if (! MP_DISPATCH( MapInsert( cl ) ))
      throw std::logic_error("Trying to map_insert() duplicated constraint: " +
                             ck.GetDescription());
    return ck.GetValueNode().Select(i);
  }


public:
  void StartModelInput() { }

  void FinishModelInput() {
    MPD( ConvertModel() );
    if (relax())              // TODO bridge?
      GetModel().RelaxIntegrality();
    GetModel().PushModelTo(GetBackend());
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
  const Backend& GetBasicBackend() const { return backend_; }
  Backend& GetBasicBackend() { return backend_; }

  /// Expose Presolver
  const pre::Presolver& GetPresolver() const { return presolver_; }
  pre::Presolver& GetPresolver() { return presolver_; }

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
  { return GetPresolver().GetTargetNodes().GetVarValues().MakeSingleKey(); }

  /// Constraint type's Value Node
  template <class Constraint>
  pre::ValueNode& GetValueNode(Constraint*)
  { return GET_CONSTRAINT_KEEPER(Constraint).GetValueNode(); }

public:
  double lb(int var) const { return this->GetModel().lb(var); }
  double ub(int var) const { return this->GetModel().ub(var); }
  template <class VarArray>
  double lb_array(const VarArray& va) const { return this->GetModel().lb_array(va); }
  template <class VarArray>
  double ub_array(const VarArray& va) const { return this->GetModel().ub_array(va); }
  var::Type var_type(int var) const { return this->GetModel().var_type(var); }
  bool is_fixed(int var) const { return this->GetModel().is_fixed(var); }
  double fixed_value(int var) const { return this->GetModel().fixed_value(var); }

  int MakeComplementVar(int bvar) {
    if (! (lb(bvar)==0.0 && ub(bvar)==1.0) )
      throw std::logic_error("Asked to complement variable with bounds "
                             + std::to_string(lb(bvar)) + ".." + std::to_string(ub(bvar)));
    AffExp ae({{-1.0}, {bvar}}, 1.0); // TODO use map / FCC?
    return MP_DISPATCH( Convert2Var(std::move(ae)) );
  }

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
  bool HasInitExpression(int var) const {
    return int(var_info_.size())>var && var_info_[var].HasId();
  }

  const ConInfo& GetInitExpression(int var) const {
    assert(HasInitExpression(var));
    return var_info_[var];
  }


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
  void InitOptions() {
    InitOwnOptions();
    GetBackend().InitOptions();
  }

protected:
  const mp::OptionValueInfo values_relax_[2] = {
    {     "0", "No (default)", 0 },
    {     "1", "Yes: treat integer and binary variables as continuous.", 1}
  };

private:
  void InitOwnOptions() {
    GetEnv().AddOption("flat:pre:all",
        "0/1*: Set to 0 to disable all presolve in the flat converter.",
        options_.preprocessAnything_, 0, 1);
    GetEnv().AddOption("flat:pre:eqresult",
        "0/1*: Preprocess reified equality comparison's boolean result bounds.",
        options_.preprocessEqualityResultBounds_, 0, 1);
    GetEnv().AddOption("flat:pre:eqbinary",
        "0/1*: Preprocess reified equality comparison with a binary variable.",
        options_.preprocessEqualityBvar_, 0, 1);
    GetEnv().AddOption("alg:relax relax",
        "0*/1: Whether to relax integrality of variables.",
        options_.relax_, 0, 1);
  }

protected:
  bool CanPreprocess(int f) const {
    return 0!=options_.preprocessAnything_ && 0!=f;
  }

  using ModelType = Model;

public:
  /// for tests. TODO make friends
  using BackendType = Backend;

  /// AddWarning
  void AddWarning(const char* key, const char* msg) {
    GetEnv().AddWarning(key, msg);
  }


private:
  BackendType backend_;
  pre::Presolver presolver_;
  pre::CopyBridge copy_bridge_ { GetPresolver() };


protected:
  /// We store backend in the converter for speed
  const BackendType& GetBackend() const { return backend_; }
  BackendType& GetBackend() { return backend_; }

  /// The internal flat model
  const ModelType& GetModel() const { return *this; }
  ModelType& GetModel() { return *this; }

  /////////////////////// CONSTRAINT KEEPERS /////////////////////////
  /// Constraint keepers and converters should be initialized after \a presolver_

  /// Define constraint keepers for all constraint types
  STORE_CONSTRAINT_TYPE(RangeLinCon)
  STORE_CONSTRAINT_TYPE(LinConLE)
  STORE_CONSTRAINT_TYPE(LinConEQ)
  STORE_CONSTRAINT_TYPE(LinConGE)
  STORE_CONSTRAINT_TYPE(LinearFunctionalConstraint)

  STORE_CONSTRAINT_TYPE(QuadraticConstraint)
  STORE_CONSTRAINT_TYPE(QuadraticFunctionalConstraint)

  STORE_CONSTRAINT_TYPE(MaximumConstraint)
  STORE_CONSTRAINT_TYPE(MinimumConstraint)
  STORE_CONSTRAINT_TYPE(AbsConstraint)
  STORE_CONSTRAINT_TYPE(ConjunctionConstraint)
  STORE_CONSTRAINT_TYPE(DisjunctionConstraint)
  STORE_CONSTRAINT_TYPE(EQ0Constraint)
  STORE_CONSTRAINT_TYPE(LE0Constraint)
  STORE_CONSTRAINT_TYPE(NotConstraint)
  STORE_CONSTRAINT_TYPE(IfThenConstraint)
  STORE_CONSTRAINT_TYPE(AllDiffConstraint)
  STORE_CONSTRAINT_TYPE(NumberofConstConstraint)
  STORE_CONSTRAINT_TYPE(NumberofVarConstraint)
  STORE_CONSTRAINT_TYPE(CountConstraint)

  STORE_CONSTRAINT_TYPE(ExpConstraint)
  STORE_CONSTRAINT_TYPE(ExpAConstraint)
  STORE_CONSTRAINT_TYPE(LogConstraint)
  STORE_CONSTRAINT_TYPE(LogAConstraint)
  STORE_CONSTRAINT_TYPE(PowConstraint)
  STORE_CONSTRAINT_TYPE(SinConstraint)
  STORE_CONSTRAINT_TYPE(CosConstraint)
  STORE_CONSTRAINT_TYPE(TanConstraint)
  STORE_CONSTRAINT_TYPE(IndicatorConstraintLinLE)
  STORE_CONSTRAINT_TYPE(IndicatorConstraintLinEQ)
  STORE_CONSTRAINT_TYPE(PLConstraint)
  STORE_CONSTRAINT_TYPE(SOS1Constraint)
  STORE_CONSTRAINT_TYPE(SOS2Constraint)

  /////////////////////// CONSTRAINT CONVERTERS /////////////////////////
  /// Constraint keepers and converters should be initialized after \a presolver_

  /// Convert range constraints, if necessary
  INSTALL_ITEM_CONVERTER(RangeConstraintConverter)

public:
  /// Presolve bridge copying values between model items
  pre::CopyBridge& GetCopyBridge() { return copy_bridge_; }
};

/// A 'final' flat converter in a hierarchy
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
