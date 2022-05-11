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
#include "mp/flat/constraint_keeper.h"
#include "mp/flat/constraints_std.h"
#include "mp/presolve.h"
#include "mp/flat/redef/std/range_con.h"

namespace mp {

/// FlatConverter: preprocesses and manages flat constraints.
/// Such constraints might need to be converted to others, which is
/// handled by overloaded methods in derived classes
/// @param Impl: the final CRTP class
/// @param ModelAPI: the solver's model API wrapper
/// @param FlatModel: internal representation of a flat model
template <class Impl, class ModelAPI,
          class FlatModel = BasicFlatModel< > >
class FlatConverter :
    public BasicFlatConverter,
    public FlatModel,
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


  //////////////////////////// CONVERTERS OF STANDRAD MP ITEMS //////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////
public:
  /// Fix a resulting variable of a logical expression as true
  /// and propagate positive ctx
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


public:
  /// Narrow variable domain range
  void NarrowVarBounds(int var, double lb, double ub) {
    auto& m = GetModel();
    m.set_lb(var, std::max(m.lb(var), lb));
    m.set_ub(var, std::min(m.ub(var), ub));
    if (m.lb(var)>m.ub(var))             // TODO write .sol, report .iis
      MP_INFEAS("empty variable domain");
  }


public:
  //////////////////////////////////// VISITOR ADAPTERS /////////////////////////////////////////

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

  /// ComputeBoundsAndType(LinTerms)
  PreprocessInfoStd ComputeBoundsAndType(const LinTerms& lt) {
    PreprocessInfoStd result;
    result.lb_ = result.ub_ = 0.0;    // TODO reuse bounds if supplied
    result.type_ = var::INTEGER;
    auto& model = MP_DISPATCH( GetModel() );
    for (auto i=lt.size(); i--; ) {
      auto v = lt.var(i);
      auto c = lt.coef(i);
      if (c >= 0.0) {
        result.lb_ += c * model.lb(v);
        result.ub_ += c * model.ub(v);
      } else {
        result.lb_ += c * model.ub(v);
        result.ub_ += c * model.lb(v);
      }
      if (var::INTEGER!=model.var_type(v) || !is_integer(c)) {
        result.type_=var::CONTINUOUS;
      }
    }
    return result;
  }

  /// ComputeBoundsAndType(AlgebraicExpr<>)
  template <class Body>
  PreprocessInfoStd ComputeBoundsAndType(
      const AlgebraicExpression<Body>& ae) {
    PreprocessInfoStd result = ComputeBoundsAndType(ae.GetBody());
    result.lb_ += ae.constant_term();    // TODO reuse bounds if supplied
    result.ub_ += ae.constant_term();
    if (!is_integer(ae.constant_term()))
      result.type_ = var::CONTINUOUS;
    return result;
  }

  /// ComputeBoundsAndType(QuadTerms)
  PreprocessInfoStd ComputeBoundsAndType(const QuadTerms& qt) {
    PreprocessInfoStd result;
    result.lb_ = result.ub_ = 0.0;
    result.type_ = var::INTEGER;
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
      }
    }
    return result;
  }

  /// ComputeBoundsAndType(QuadAndLinearTerms)
  PreprocessInfoStd ComputeBoundsAndType(const QuadAndLinTerms& qlt) {
    auto bntLT = ComputeBoundsAndType(qlt.GetLinTerms());
    auto bntQT = ComputeBoundsAndType(qlt.GetQPTerms());
    return AddBoundsAndType(bntLT, bntQT);
  }

  /// Product bounds
  template <class Var>
  std::pair<double, double> ProductBounds(Var x, Var y) const {
    auto lx=lb(x), ly=lb(y), ux=ub(x), uy=ub(y);
    std::array<double, 4> pb{lx*ly, lx*uy, ux*ly, ux*uy};
    return { *std::min_element(pb.begin(), pb.end()),
          *std::max_element(pb.begin(), pb.end()) };
  }

  /// Add / merge bounds and type
  PreprocessInfoStd AddBoundsAndType(const PreprocessInfoStd& bnt1,
                                     const PreprocessInfoStd& bnt2) {
    return {bnt1.lb()+bnt2.lb(), bnt1.ub()+bnt2.ub(),
      var::INTEGER==bnt1.type() && var::INTEGER==bnt2.type() ?
            var::INTEGER : var::CONTINUOUS};
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

  /// PreprocessConstraint
  /// TODO rename 'PropagateUp' and define together with the constraint

  template <class PreprocessInfo>
  void PreprocessConstraint(
      LinearFunctionalConstraint& c, PreprocessInfo& prepro) {
    auto pre = ComputeBoundsAndType(c.GetAffineExpr());
    prepro.narrow_result_bounds( pre.lb(), pre.ub() );
    prepro.set_result_type( pre.type() );
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      QuadraticFunctionalConstraint& c, PreprocessInfo& prepro) {
    auto pre = ComputeBoundsAndType(c.GetQuadExpr());
    prepro.narrow_result_bounds( pre.lb(), pre.ub() );
    prepro.set_result_type( pre.type() );
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      PowConstraint& c, PreprocessInfo& prepro) {
    auto& m = MP_DISPATCH( GetModel() );
    auto arg = c.GetArguments()[0];
    auto prm = c.GetParameters()[0];
    auto lb = std::pow(m.lb(arg), prm),
        ub = std::pow(m.ub(arg), prm);
    prepro.narrow_result_bounds( std::min(lb, ub),
                          std::max(lb, ub) );
    prepro.set_result_type( m.var_type(arg) );
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      TanConstraint& , PreprocessInfo& ) {
    // TODO improve
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      MinConstraint& c, PreprocessInfo& prepro) {
    auto& m = MP_DISPATCH( GetModel() );
    auto& args = c.GetArguments();
    prepro.narrow_result_bounds( m.lb_array(args),
                          m.ub_min_array(args) );
    prepro.set_result_type( m.common_type(args) );
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      MaxConstraint& c, PreprocessInfo& prepro) {
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

  /// Preprocess CondLinConEQ
  template <class PreprocessInfo>
  void PreprocessConstraint(
      CondLinConEQ& c, PreprocessInfo& prepro) {
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

  /// Preprocess CondQuadConEQ
  template <class PreprocessInfo>
  void PreprocessConstraint(
      CondQuadConEQ& c, PreprocessInfo& prepro) {
    prepro.narrow_result_bounds(0.0, 1.0);
    prepro.set_result_type( var::INTEGER );
    if (0!=CanPreprocess( options_.preprocessEqualityResultBounds_ ))
      if (FixEqualityResult(c, prepro))
        return;
  }

  /// Try and fix conditional equality result
  /// @return true if success
  template <class PreprocessInfo, class CondAlgCon>
  bool FixEqualityResult(
      CondAlgCon& c, PreprocessInfo& prepro) {
    const auto& con = c.GetConstraint();
    const auto& body = con.GetBody();
    const auto rhs = con.rhs();
    // TODO expr is empty. Possible?
    auto bndsNType = ComputeBoundsAndType(body);
    if (bndsNType.lb() > rhs || bndsNType.ub() < rhs) {
      /// TODO this depends on context???
      prepro.narrow_result_bounds(0.0, 0.0);
      return true;
    }
    if (bndsNType.lb()==rhs && bndsNType.ub()==rhs) {
      /// TODO this depends on context???
      prepro.narrow_result_bounds(1.0, 1.0);
      return true;
    }
    if (var::INTEGER==bndsNType.type_ &&
        !is_integer(con.rhs())) {
      /// TODO this depends on context???
      prepro.narrow_result_bounds(0.0, 0.0);
      return true;
    }
    return false;
  }

  /// Normalize conditional equality coef * var == const
  static void PreprocessEqVarConst__unifyCoef(CondLinConEQ& c) {
    auto& con = c.GetConstraint();
    auto& body = con.GetBody();
    if (1==body.size()) {
      const double coef = body.coef(0);
      if (1.0!=coef) {
        assert(0.0!=std::fabs(coef));
        con.set_rhs(con.rhs() / coef);
        body.set_coef(0, 1.0);
      }
    }
  }

  /// Simplify conditional equality bin_var==0/1
  /// by reusing bin_var or its complement
  template <class PreprocessInfo>
  bool ReuseEqualityBinaryVar(
      CondLinConEQ& c, PreprocessInfo& prepro) {
    auto& m = MP_DISPATCH( GetModel() );
    const auto& con = c.GetConstraint();
    const auto& body = con.GetBody();
    if (1==body.size()) {                           // var==const
      assert( 1.0==body.coef(0) );                  // is normalized
      int var = body.var(0);
      if (m.is_binary_var(var)) {            // See if this is binary var==const
        const double rhs = con.rhs();
        if (1.0==rhs)
          prepro.set_result_var( var );      // TODO need to know var is a result var?
        else if (0.0==std::fabs(rhs))
          prepro.set_result_var( MakeComplementVar(var) );
        else
          prepro.narrow_result_bounds(0.0, 0.0);    // not 0/1 value, result false
        return true;
      }
    }
    return false;
  }

  /// Preprocess other conditional comparisons TODO
  template <class PreprocessInfo, class Body, int kind>
  void PreprocessConstraint(
      ConditionalConstraint<
        AlgebraicConstraint< Body, AlgConRhs<kind> > >& ,
      PreprocessInfo& prepro) {
    prepro.narrow_result_bounds(0.0, 1.0);
    prepro.set_result_type( var::INTEGER );
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      AndConstraint& , PreprocessInfo& prepro) {
    prepro.narrow_result_bounds(0.0, 1.0);
    prepro.set_result_type( var::INTEGER );
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      OrConstraint& , PreprocessInfo& prepro) {
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
    prepro.narrow_result_bounds(0.0,     // size()-1: 1st arg is the ref var
                                con.GetArguments().size()-1);
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
      DivConstraint& c, PreprocessInfo& prepro) {
    auto& m = MPD( GetModel() );
    auto v1 = c.GetArguments()[0], v2 = c.GetArguments()[1];
    const auto l1=m.lb(v1), u1=m.ub(v1), l2=m.lb(v2), u2=m.ub(v2);
    if (l1 > this->PracticallyMinusInfty() &&
        u1 < this->PracticallyInfty() &&
        l2 > this->PracticallyMinusInfty() &&
        u2 < this->PracticallyInfty() &&
        l2 * u2 > 0.0) {
      auto l0 = std::numeric_limits<double>::max();
      auto u0 = std::numeric_limits<double>::min();
      {
        l0 = std::min(l0, l1 / l2);
        l0 = std::min(l0, l1 / u2);
        l0 = std::min(l0, u1 / l2);
        l0 = std::min(l0, u1 / u2);
        u0 = std::max(u0, l1 / l2);
        u0 = std::max(u0, l1 / u2);
        u0 = std::max(u0, u1 / l2);
        u0 = std::max(u0, u1 / u2);
      }
      prepro.narrow_result_bounds( l0, u0 );
    }
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

  template <class PreprocessInfo>
  void PreprocessConstraint(
      PLConstraint& , PreprocessInfo& ) {
  }



  //////////////////////////// CUSTOM CONSTRAINTS //////////////////////
  ///
  //////////////////////////// SPECIFIC CONSTRAINT RESULT-TO-ARGUMENTS PROPAGATORS //////
  /// Currently we should propagate to all arguments, be it always the CTX_MIX.

  /// Allow ConstraintKeeper to PropagateResult(), use GetBackend() etc
  template <class , class , class >
  friend class ConstraintKeeper;

public:
  /// By default, set mixed context for argument variables
  template <class Constraint>
  void PropagateResult(Constraint& con, double lb, double ub, Context ctx) {
    internal::Unused(&con, lb, ub, ctx);
    con.SetContext(ctx);
    PropagateResult2Args(con.GetArguments(),
                         this->MinusInfty(), this->Infty(), Context::CTX_MIX);
  }

  void PropagateResult(LinearFunctionalConstraint& con, double lb, double ub,
                       Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    PropagateResult2LinTerms(con.GetAffineExpr(),
                             this->MinusInfty(), this->Infty(), +ctx);
  }

  void PropagateResult(QuadraticFunctionalConstraint& con, double lb, double ub,
                       Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    const auto& args = con.GetArguments();
    PropagateResult2LinTerms(args.GetLinTerms(),
                             this->MinusInfty(), this->Infty(), +ctx);
    PropagateResult2QuadTerms(args.GetQPTerms(),
                              this->MinusInfty(), this->Infty(), +ctx);
  }

  void PropagateResult(QuadConRange& con, double lb, double ub,
                       Context ctx) {
    internal::Unused(lb, ub, ctx);
    PropagateResult2LinTerms(con.GetLinTerms(), // TODO sense dep. on bounds
                             this->MinusInfty(), this->Infty(), ctx);
    PropagateResult2QuadTerms(con.GetQPTerms(), // TODO bounds?
                              this->MinusInfty(), this->Infty(), ctx);
  }

  template <class Body, int sens>
  void PropagateResult(IndicatorConstraint<
                         AlgebraicConstraint< Body, AlgConRhs<sens> > >& con,
                       double lb, double ub, Context ctx) {
    internal::Unused(lb, ub, ctx);
    PropagateResultOfInitExpr(con.get_binary_var(),
                              this->MinusInfty(), this->Infty(),
                              1==con.get_binary_value() ?  // b==1 means b in CTX_NEG
                                Context::CTX_NEG : Context::CTX_POS);
    PropagateResult2Args(con.get_constraint().GetBody(),   // Assume Con::BodyType is handled
                             this->MinusInfty(), this->Infty(),
                             0==sens ? Context::CTX_MIX :
                                       0<sens ? +ctx : -ctx);
  }

  template <int type>
  void PropagateResult(SOS_1or2_Constraint<type>& con, double lb, double ub, Context ) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    PropagateResult2Vars(con.get_vars(),
                         this->MinusInfty(), this->Infty(), Context::CTX_MIX);
  }

  void PropagateResult(ComplementarityLinear& con, double lb, double ub,
                       Context ctx) {
    internal::Unused(lb, ub, ctx);
    PropagateResult2LinTerms(con.GetExpression().GetLinTerms(),
                         lb, ub, Context::CTX_MIX);
    PropagateResultOfInitExpr(con.GetVariable(), lb, ub, Context::CTX_MIX);
  }

  void PropagateResult(ComplementarityQuadratic& con, double lb, double ub,
                       Context ctx) {
    internal::Unused(lb, ub, ctx);
    PropagateResult2LinTerms(con.GetExpression().GetLinTerms(),
                    lb, ub, Context::CTX_MIX);
    PropagateResult2QuadTerms(con.GetExpression().GetQPTerms(),
                              lb, ub, Context::CTX_MIX);
    PropagateResultOfInitExpr(con.GetVariable(), lb, ub, Context::CTX_MIX);
  }

  void PropagateResult(NotConstraint& con, double lb, double ub, Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    PropagateResultOfInitExpr(con.GetArguments()[0], 1.0-ub, 1.0-lb, -ctx);
  }

  void PropagateResult(AndConstraint& con, double lb, double ub, Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    PropagateResult2Vars(con.GetArguments(), lb, 1.0, +ctx);
  }

  void PropagateResult(OrConstraint& con, double lb, double ub, Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    PropagateResult2Vars(con.GetArguments(), 0.0, ub, +ctx);
  }

  void PropagateResult(IfThenConstraint& con, double lb, double ub, Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    auto& args = con.GetArguments();
    /// TODO consider bounds for then/else for the context:
    PropagateResultOfInitExpr(args[0], 0.0, 1.0, Context::CTX_MIX);
    PropagateResultOfInitExpr(args[1], this->MinusInfty(), this->Infty(), +ctx);
    PropagateResultOfInitExpr(args[2], this->MinusInfty(), this->Infty(), -ctx);
  }

  void PropagateResult(AllDiffConstraint& con, double lb, double ub, Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    PropagateResult2Vars(con.GetArguments(), this->MinusInfty(), this->Infty(),
                         Context::CTX_MIX);
  }

  void PropagateResult(NumberofConstConstraint& con, double lb, double ub, Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    PropagateResult2Vars(con.GetArguments(), this->MinusInfty(), this->Infty(),
                         Context::CTX_MIX);
  }

  void PropagateResult(NumberofVarConstraint& con, double lb, double ub, Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    PropagateResult2Vars(con.GetArguments(), this->MinusInfty(), this->Infty(),
                         Context::CTX_MIX);
  }

  void PropagateResult(CondLinConEQ& con, double lb, double ub, Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    PropagateResult2LinTerms(con.GetConstraint().GetBody(),
                         lb, ub, Context::CTX_MIX);
  }

  void PropagateResult(CondQuadConEQ& con, double lb, double ub, Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    PropagateResult2QuadAndLinTerms(con.GetConstraint().GetBody(),
                         lb, ub, Context::CTX_MIX);
  }

  template <class Body, int kind>
  void PropagateResult(
      ConditionalConstraint<
        AlgebraicConstraint< Body, AlgConRhs<kind> > >& con,
      double lb, double ub, Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    PropagateResult2Args(con.GetConstraint().GetBody(), lb, ub,
                             kind>0 ? ctx : -ctx);
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
      PropagateResultOfInitExpr(lint.var(i),      /// TODO bounds as well
                                this->MinusInfty(), this->Infty(),
                                (lint.coef(i)>=0.0) ? +ctx : -ctx);
    }
  }

  /// Propagate given bounds & context into a vector of variables
  /// @param lb, ub: bounds for each variable
  template <class Vec>
  void PropagateResult2Vars(const Vec& vars, double lb, double ub, Context ctx) {
    for (auto v: vars) {
      PropagateResultOfInitExpr(v, lb, ub, ctx);
    }
  }

  /// Propagate result into QuadTerms
  void PropagateResult2QuadTerms(const QuadTerms& quadt, double , double , Context ) {
    for (auto i=quadt.size(); i--; ) {             /// TODO context for special cases
      PropagateResultOfInitExpr(quadt.var1(i),     /// TODO bounds as well
                                this->MinusInfty(), this->Infty(), Context::CTX_MIX);
      PropagateResultOfInitExpr(quadt.var2(i),
                                this->MinusInfty(), this->Infty(), Context::CTX_MIX);
    }
  }


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
    if (relax())              // TODO bridge?
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
  /// var_type()
  var::Type var_type(int var) const { return this->GetModel().var_type(var); }
  /// is_fixed()
  bool is_fixed(int var) const { return this->GetModel().is_fixed(var); }
  /// fixed_value()
  double fixed_value(int var) const { return this->GetModel().fixed_value(var); }

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
        "0/1*: Set to 0 to disable all presolve in the flat converter.",
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

protected:
  bool CanPreprocess(int f) const {
    return 0!=options_.preprocessAnything_ && 0!=f;
  }

  using ModelType = FlatModel;

public:
  /// for tests. TODO make friends
  using ModelAPIType = ModelAPI;

  /// AddWarning. Strings should remain valid
  void AddWarning(const char* key, const char* msg) {
    GetEnv().AddWarning(key, msg);
  }


private:
  ModelAPIType modelapi_;      // We store modelapi in the converter for speed
  pre::Presolver presolver_;
  pre::CopyBridge copy_bridge_ { GetPresolver() };


protected:
  /// The internal flat model
  const ModelType& GetModel() const { return *this; }
  ModelType& GetModel() { return *this; }

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
  /// Presolve bridge copying values between model items
  pre::CopyBridge& GetCopyBridge() { return copy_bridge_; }
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
