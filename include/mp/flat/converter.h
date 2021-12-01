#ifndef CONVERTER_FLAT_H
#define CONVERTER_FLAT_H

#include <unordered_map>
#include <map>
#include <cstdio>
#include <cmath>

#include "mp/format.h"
#include "mp/flat/converter_proxy_base.h"
#include "mp/flat/backend_base.h"
#include "mp/flat/convert_functional.h"
#include "mp/flat/constraint_keeper.h"
#include "mp/flat/model.h"
#include "mp/flat/std_constr.h"
#include "mp/flat/preprocess.h"

namespace mp {

/// BasicMPFlatConverter: it "flattens" most expressions
/// by replacing them by a result variable and constraints.
/// Such constraints might need to be converted to others, which is
/// handled by overloaded methods in derived classes
template <class Impl, class Backend,
          class Model = BasicFlatModel< > >
class FlatConverter :
    public FlatConverterProxy,
    public BasicFlatConverter,
    public Model
{
public:
  static constexpr const char* name() { return "Flat Converter"; };

  using Var = typename Model::Var;
  static constexpr Var VoidVar() { return Model::VoidVar(); }

public:
  using VarArray = std::vector<int>;
  template <class Constraint>
    using ConstraintKeeperType = ConstraintKeeper<Impl, Backend, Constraint>;

protected:
  using ClassType = FlatConverter<Impl, Backend, Model>;
  using BaseConverter = BasicFlatConverter;
  using BaseFlatModel = Model;

public:
  static const char* GetConverterName() { return "BasicMPFlatConverter"; }

  FlatConverter()
  { GetBackend().ProvideFlatConverterProxyObject(this); }


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
    if (HasInitExpression(var))
      GetInitExpression(var)->PropagateResult(*this, lb, ub, ctx);
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
  int Convert2Var(QuadExpr&& ee) {
    if (ee.is_variable())
      return ee.get_representing_variable();
    if (ee.is_constant())
      return MakeFixedVar(ee.constant_term());
    PreprocessInfoStd bnt = ComputeBoundsAndType(ee);
    auto r = MP_DISPATCH( AddVar(bnt.lb_, bnt.ub_, bnt.type_) );
    if (ee.is_affine())
      AddConstraint(LinearDefiningConstraint(r, std::move(ee.GetAE())));
    else
      AddConstraint(QuadraticDefiningConstraint(r, std::move(ee)));
    return r;
  }

  PreprocessInfoStd ComputeBoundsAndType(const QuadExpr& ee) {
    auto bntAE = ComputeBoundsAndType(ee.GetAE());
    auto bntQT = ComputeBoundsAndType(ee.GetQT());
    return AddBoundsAndType(bntAE, bntQT);
  }

  PreprocessInfoStd ComputeBoundsAndType(const AffineExpr& ae) {
    PreprocessInfoStd result;
    result.lb_ = result.ub_ = ae.constant_term();    // TODO reuse bounds if supplied
    result.type_ = is_integer(result.lb_) ? var::INTEGER : var::CONTINUOUS;
    result.linexp_type_ = var::INTEGER;
    auto& model = MP_DISPATCH( GetModel() );
    for (const auto& term: ae) {
      auto v = term.var_index();
      if (term.coef() >= 0.0) {
        result.lb_ += term.coef() * model.lb(v);
        result.ub_ += term.coef() * model.ub(v);
      } else {
        result.lb_ += term.coef() * model.ub(v);
        result.ub_ += term.coef() * model.lb(v);
      }
      if (var::INTEGER!=model.var_type(v) || !is_integer(term.coef())) {
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
    for (int i=0; i<qt.num_terms(); ++i) {
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
      for (int endConstraintsThisLoop = 0, endPrevious = 0;
           (endConstraintsThisLoop = this->GetModel().num_custom_cons()) > endPrevious;
           endPrevious = endConstraintsThisLoop
           ) {
        MP_DISPATCH( PreprocessIntermediate() );                        // preprocess before each level
        ConvertExtraItemsInRange(endPrevious, endConstraintsThisLoop);
      }
      MP_DISPATCH( ConvertMaps(); );
      MP_DISPATCH( PreprocessFinal() );                                 // final prepro
    } catch (const ConstraintConversionFailure& cff) {
      throw std::logic_error(cff.message());
    }
  }

  void ConvertExtraItemsInRange(int first, int after_last) {
    for (; first<after_last; ++first) {
      auto* pConstraint = this->GetModel().custom_con(first);
      if (!pConstraint->IsRemoved()) {
        const auto acceptanceLevel =
            pConstraint->BackendAcceptance(this->GetBackend());
        if (NotAccepted == acceptanceLevel) {
          pConstraint->ConvertWith(*this);
          pConstraint->Remove();
        }
        else if (AcceptedButNotRecommended == acceptanceLevel) {
          try {
            pConstraint->ConvertWith(*this);
            pConstraint->Remove();
          } catch (const ConstraintConversionFailure& ccf) {
            printf( fmt::format(      // TODO use Env
                           "WARNING: {}. Will pass the constraint "
                           "to the backend {}. Continuing\n",
                           ccf.message(),
                           MP_DISPATCH( GetBackend() ).GetBackendName() ).c_str() );
          }
        }
      }
    }
  }


  //////////////////////// WHOLE-MODEL PREPROCESSING /////////////////////////
  void PreprocessIntermediate() { }
  void ConvertMaps() { }
  void PreprocessFinal() { }



  //////////////////////////// CUSTOM CONSTRAINTS ////////////////////////////
  ///
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

  template <class PreprocessInfo>
  void PreprocessConstraint(
      AbsConstraint& c, PreprocessInfo& prepro) {
    const auto argvar = c.GetArguments()[0];
    const auto lb = this->lb(argvar),
        ub = this->ub(argvar);
    if (lb>=0.0) {
      prepro.set_result_var(argvar);
      return;
    } else if (ub<=0.0) {
      auto res = AssignResult2Args(   // create newvar = -argvar
            LinearDefiningConstraint({ {-1.0}, {argvar}, 0.0 }));
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
    AffineExpr& ae = c.GetArguments();
    if (ae.is_constant()) {                  // const==0
      auto res = (double)int(0.0==ae.constant_term());
      prepro.narrow_result_bounds(res, res);
      return true;
    }
    auto bndsNType = ComputeBoundsAndType(ae);
    if (bndsNType.lb() > 0.0 || bndsNType.ub() < 0.0) {
      prepro.narrow_result_bounds(0.0, 0.0);
      return true;
    }
    if (bndsNType.lb()==0.0 && bndsNType.ub()==0.0) {
      prepro.narrow_result_bounds(1.0, 1.0);
      return true;
    }
    if (var::INTEGER==bndsNType.linexp_type_ &&
        !is_integer(ae.constant_term())) {
      prepro.narrow_result_bounds(0.0, 0.0);
      return true;
    }
    return false;
  }

  static void PreprocessEqVarConst__unifyCoef(EQ0Constraint& c) {
    AffineExpr& ae = c.GetArguments();
    if (1==ae.num_terms()) {
      const double coef = ae.coef(0);
      if (1.0!=coef) {
        assert(0.0!=std::fabs(coef));
        ae.constant_term(ae.constant_term() / coef);
        ae.set_coef(0, 1.0);
      }
    }
  }

  template <class PreprocessInfo>
  bool ReuseEqualityBinaryVar(
      EQ0Constraint& c, PreprocessInfo& prepro) {
    auto& m = MP_DISPATCH( GetModel() );
    AffineExpr& ae = c.GetArguments();
    if (1==ae.num_terms()) {                           // var==const
      assert( 1.0==ae.coef(0) );
      int var = ae.var_index(0);
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
    internal::Unused(con, lb, ub, ctx);
    con.SetContext(Context::CTX_MIX);
    for (const auto a: con.GetArguments())
      PropagateResultOfInitExpr(a, this->MinusInfty(), this->Infty(), Context::CTX_MIX);
  }

  void PropagateResult(LinearDefiningConstraint& con, double lb, double ub, Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    for (const auto& term: con.GetAffineExpr())
      PropagateResultOfInitExpr(term.var_index(), this->MinusInfty(), this->Infty(), Context::CTX_MIX);
  }

  void PropagateResult(QuadraticDefiningConstraint& con, double lb, double ub, Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
    const auto& args = con.GetArguments();
    for (const auto& term: args.GetAE())
      PropagateResultOfInitExpr(term.var_index(), this->MinusInfty(), this->Infty(), Context::CTX_MIX);
    const auto& qt = args.GetQT();
    for (int i=0; i<qt.num_terms(); ++i) {
      PropagateResultOfInitExpr(qt.var1(i), this->MinusInfty(), this->Infty(), Context::CTX_MIX);
      PropagateResultOfInitExpr(qt.var2(i), this->MinusInfty(), this->Infty(), Context::CTX_MIX);
    }
  }

  void PropagateResult(LinearConstraint& con, double lb, double ub, Context ctx) {
    internal::Unused(lb, ub, ctx);
    for (const auto v: con.vars())
      PropagateResultOfInitExpr(v, this->MinusInfty(), this->Infty(), Context::CTX_MIX);
  }

  void PropagateResult(QuadraticConstraint& con, double lb, double ub, Context ctx) {
    internal::Unused(lb, ub, ctx);
    for (const auto v: con.vars())
      PropagateResultOfInitExpr(v, this->MinusInfty(), this->Infty(), Context::CTX_MIX);
    for (const auto v: con.GetQPTerms().vars1())
      PropagateResultOfInitExpr(v, this->MinusInfty(), this->Infty(), Context::CTX_MIX);
    for (const auto v: con.GetQPTerms().vars2())
      PropagateResultOfInitExpr(v, this->MinusInfty(), this->Infty(), Context::CTX_MIX);
  }

  template <int sens>
  void PropagateResult(IndicatorConstraintLin<sens>& con, double lb, double ub, Context ctx) {
    internal::Unused(lb, ub, ctx);
    PropagateResultOfInitExpr(con.get_binary_var(),
                              this->MinusInfty(), this->Infty(), Context::CTX_MIX);
    for (const auto v: con.get_lin_vars())
      PropagateResultOfInitExpr(v, this->MinusInfty(), this->Infty(), Context::CTX_MIX);
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

  void PropagateResult(LE0Constraint& con, double lb, double ub, Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
  }

  void PropagateResult(EQ0Constraint& con, double lb, double ub, Context ctx) {
    MPD( NarrowVarBounds(con.GetResultVar(), lb, ub) );
    con.AddContext(ctx);
  }


  //////////////////////////// CUSTOM CONSTRAINTS CONVERSION ////////////////////////////
  ///
  //////////////////////////// SPECIFIC CONSTRAINT CONVERTERS ///////////////////////////

  USE_BASE_CONSTRAINT_CONVERTERS(BasicFlatConverter);      // reuse default converters

public: // for ConstraintKeeper
  /// Assume mixed context if not set in the constraint
  /// TODO Make sure context is always propagated for all constraints and objectives
  template <class Constraint>
  void RunConversion(const Constraint& con) {
    if (con.HasContext())
      if (con.GetContext().IsNone())
        con.SetContext(Context::CTX_MIX);
    MP_DISPATCH(Convert(con));
  }


  /// If backend does not like LDC, we can redefine it
  void Convert(const LinearDefiningConstraint& ldc) {
    this->AddConstraint(ldc.to_linear_constraint());
  }

  /// If backend does not like QDC, we can redefine it
  void Convert(const QuadraticDefiningConstraint& qdc) {
    this->AddConstraint(qdc.to_quadratic_constraint());
  }

public:
  //////////////////////// ADD CUSTOM CONSTRAINT ///////////////////////
  //////////////////////// Takes ownership /////////////////////////////
  template <class Constraint>
  void AddConstraint(Constraint&& con) {
    const auto pck = makeConstraintKeeper<Impl, Constraint>(std::forward<Constraint>(con));
    AddConstraintAndTryNoteResultVariable(pck);
//    orig_lin_constr_.push_back(n_alg_constr_);
//    orig_qp_constr_.push_back(n_alg_constr_);
//    ++n_alg_constr_;
  }

protected:
  template <class ConstraintKeeper>
  void AddConstraintAndTryNoteResultVariable(ConstraintKeeper* pbc) {
    MP_DISPATCH( GetModel() ).AddConstraint(pbc);
    const auto resvar = pbc->GetResultVar();
    if (resvar>=0)
      AddInitExpression(resvar, pbc);
    if (! MP_DISPATCH( MapInsert(pbc) ))
      throw std::logic_error("Trying to map_insert() duplicated constraint: " +
                             pbc->GetDescription());
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
    MP_DISPATCH( PrepareConversion() );
    MPD( ConvertItems() );
  }

  void PrepareConversion() {
    MP_DISPATCH( MemorizeModelSize() );
  }

  /// TODO wrong, Visit() may have produced aux constraints
  void MemorizeModelSize() {
    n_constr_orig_ = GetModel().num_custom_cons();
    for (size_t i=0; i<n_constr_orig_; ++i) {
      auto cc = GetModel().custom_con(i)->ConstraintClass();
      if (1==cc) {
        orig_lin_constr_.push_back(i);
      } else if (2==cc) {
        orig_qp_constr_.push_back(i);
      }
    }
  }

  /// Elementary constraint classes
  /// TODO generalize
  static constexpr int ConstraintClass(BasicConstraint* )
  { return 0; }
  static constexpr int ConstraintClass(LinearConstraint* )
  { return 1; }
  static constexpr int ConstraintClass(LinearDefiningConstraint* )
  { return 1; }
  static constexpr int ConstraintClass(QuadraticConstraint* )
  { return 2; }
  static constexpr int ConstraintClass(QuadraticDefiningConstraint* )
  { return 2; }

  //////////////////////////// UTILITIES /////////////////////////////////
  ///
public:
  /// Expose abstract Backend
  const BasicBackend& GetBasicBackend() const { return backend_; }
  BasicBackend& GetBasicBackend() { return backend_; }

private:
  std::unordered_map<double, int> map_fixed_vars_;

  std::vector<int> common_exprs_;               // variables equal to the result

public:

  //////////////////////////// CREATE OR FIND A FIXED VARIABLE //////////////////////////////
  int MakeFixedVar(double value) {
    auto it = map_fixed_vars_.find(value);
    if (map_fixed_vars_.end()!=it)
      return it->second;
    auto v = MPD( DoAddVar(value, value) );
    map_fixed_vars_[value] = v;
    return v;
  }

  /// Create or find a fixed variable
  int AddVar(double lb, double ub, var::Type type = var::CONTINUOUS) {
    if (lb!=ub)
      return DoAddVar(lb, ub, type);
    return MakeFixedVar(lb);
  }

  using BaseFlatModel::AddVars;

protected:
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
    AffineExpr ae({-1.0}, {bvar}, 1.0); // TODO use map / FCC?
    return MP_DISPATCH( Convert2Var(std::move(ae)) );
  }

  struct VarInfo {
    BasicConstraintKeeper *pInitExpr=nullptr;
  };

  /// Add variable. Type: var::CONTINUOUS by default
  int DoAddVar(double lb=MinusInfty(), double ub=Infty(),
             var::Type type = var::CONTINUOUS) {
    return GetModel().AddVar(lb, ub, type);
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
  std::vector<VarInfo> var_info_;

protected:
  void AddInitExpression(int var, BasicConstraintKeeper* pie) {
    var_info_.resize(std::max(var_info_.size(), (size_t)var+1));
    var_info_[var].pInitExpr = pie;
  }

  bool HasInitExpression(int var) const {
    return int(var_info_.size())>var &&
        nullptr!=var_info_[var].pInitExpr;
  }

  BasicConstraintKeeper* GetInitExpression(int var) {
    assert(HasInitExpression(var));
    return var_info_[var].pInitExpr;
  }

protected:
  /// Convenience methods
  /// Gurobi reports duals etc separately for linear and QCP constraints
  /// Implement FlatConverterProxy's interface

  /////////////////////////////////////////////////////////////////////////
  /// PRESOLVE ///
  /////////////////////////////////////////////////////////////////////////

  /// From original NL model's suffix or duals
  std::vector<double>
  ExtractLinConValues(ArrayRef<double> allval) override
  { return ExtractLinConValues<>(allval); }

  std::vector<int>
  ExtractLinConValues(ArrayRef<int> allval) override
  { return ExtractLinConValues<>(allval); }

  std::vector<double> ExtractQCValues(ArrayRef<double> allval) override
  { return ExtractQCValues<>(allval); }

  /////////////////////////////////////////////////////////////////////////
  /// POSTSOLVE ///
  /////////////////////////////////////////////////////////////////////////

  /// To original NL model's indexing
  std::vector<double> MakeConstrValuesFromLPAndQCP(
        ArrayRef<double> pi, ArrayRef<double> qcpi) override
  { return MakeConstrValuesFromLPAndQCP<>(pi, qcpi); }
  std::vector<int> MakeConstrValuesFromLPAndQCP(
        ArrayRef<int> pi, ArrayRef<int> qcpi) override
  { return MakeConstrValuesFromLPAndQCP<>(pi, qcpi); }

  /////////////////////////////////////////////////////////////////////////
  /// IMPLEMENTATIONS OF PRE- / POSTSOLVE ///
  /////////////////////////////////////////////////////////////////////////

  template <class Vec>
  auto MakeConstrValuesFromLPAndQCP(
        const Vec& pi, const Vec& qcpi) ->
          std::vector<
            typename std::decay< decltype(*pi.data()) >::type>
  {
    std::vector<typename std::decay< decltype(*pi.data()) >::type>
        result(NumValuedOrigConstr());
    FillByIndex(result, pi, GetIndexesOfValuedLinearConstraints());
    FillByIndex(result, qcpi, GetIndexesOfValuedQPConstraints());
    return result;
  }
  template <class Vec>
  auto ExtractLinConValues(const Vec& allval) ->
  std::vector<typename std::decay< decltype(*allval.data()) >::type>
  {
    return ExtractByIndex(allval, GetIndexesOfValuedLinearConstraints());
  }
  template <class Vec>
  auto ExtractQCValues(const Vec& allval) ->
  std::vector<typename std::decay< decltype(*allval.data()) >::type>
  {
    return ExtractByIndex(allval, GetIndexesOfValuedQPConstraints());
  }

  template <class Vec1, class Vec2, class Vec3>
  static void FillByIndex(Vec1& dest, const Vec2& src, const Vec3& idx) {
    for (size_t ii=0; ii<src.size() && ii<idx.size(); ++ii)
      dest.at(idx[ii]) = src[ii];
  }
  template <class Vec2, class Vec3>
  static auto ExtractByIndex(const Vec2& src, const Vec3& idx) ->
  std::vector<typename std::decay< decltype(*src.data()) >::type>
  {
    std::vector<typename std::decay< decltype(*src.data()) >::type> rsl(idx.size());
    for (size_t ii=0; ii<idx.size(); ++ii)
      if (idx[ii] < src.size())
        rsl[ii] = src[idx[ii]];
    return rsl;
  }

  /// Number of original NL algebraic constraints relevant for suffixes
  size_t NumValuedOrigConstr() const override { return n_constr_orig_; }
  /// Gurobi handles linear constraints as a separate class.
  /// AMPL provides suffixes for all constraints together.
  /// The method returns the indexes of linear constraints
  /// which are relevant for suffix / dual value exchange,
  /// in the overall constraints list.
  const std::vector<size_t>& GetIndexesOfValuedLinearConstraints() const
  { return orig_lin_constr_; }
  /// Same for QP constraints
  const std::vector<size_t>& GetIndexesOfValuedQPConstraints() const
  { return orig_qp_constr_; }


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
  /// Chance for the Backend to init solver environment, etc
  void InitOptionParsing() { GetBackend().InitOptionParsing(); }
  /// Chance to consider options immediately (open cloud, etc)
  void FinishOptionParsing() { GetBackend().FinishOptionParsing(); }

  template <class OptionManager>
  void InitOptions(OptionManager& opt) {
    MPD( GetBackend() ).InitMetaInfoAndOptions();
    InitOwnOptions(opt);
  }

protected:
  const mp::OptionValueInfo values_relax_[2] = {
    {     "0", "No (default)", 0 },
    {     "1", "Yes: treat integer and binary variables as continuous.", 1}
  };

private:
  template <class OptionManager>
  void InitOwnOptions(OptionManager& opt) {
    opt.AddOption("flat:pre:all",
        "0/1*: Set to 0 to disable all presolve in the flat converter.",
        options_.preprocessAnything_, 0, 1);
    opt.AddOption("flat:pre:eqresult",
        "0/1*: Preprocess reified equality comparison's boolean result bounds.",
        options_.preprocessEqualityResultBounds_, 0, 1);
    opt.AddOption("flat:pre:eqbinary",
        "0/1*: Preprocess reified equality comparison with a binary variable.",
        options_.preprocessEqualityBvar_, 0, 1);
    opt.AddOption("alg:relax relax",
        "0*/1: Whether to relax integrality of variables.",
        options_.relax_, 0, 1);
  }

protected:
  bool CanPreprocess(int f) const {
    return 0!=options_.preprocessAnything_ && 0!=f;
  }

  using ModelType = Model;

public:  // for tests. TODO make friends
  using BackendType = Backend;

private:
  BackendType backend_;

  /// TODO replace by pre- / postsolve
  /// Indices of NL original linear constr in the total constr ordering
  std::vector<size_t> orig_lin_constr_;
  /// Indices of NL original QP constr in the total constr ordering
  std::vector<size_t> orig_qp_constr_;
  size_t n_constr_orig_=0;

protected:
  /// We store backend in the converter for speed
  const BackendType& GetBackend() const { return backend_; }
  BackendType& GetBackend() { return backend_; }

  /// The internal flat model
  const ModelType& GetModel() const { return *this; }
  ModelType& GetModel() { return *this; }

public:
  /// TODO use universal Env instead
  /// (which can well use these "MPUtils",
  /// but ideally an appr new base class of Solver)
  using MPUtils = typename BackendType::MPUtils;
  const MPUtils& GetMPUtils() const { return GetBackend().GetMPUtils(); }
  MPUtils& GetMPUtils() { return GetBackend().GetMPUtils(); }
};

/// A 'final' converter in a hierarchy
template <template <typename, typename, typename> class Converter,
          class Backend, class Model = BasicFlatModel< > >
class ConverterImpl :
    public Converter<ConverterImpl<Converter, Backend, Model>, Backend, Model> { };

template <template <typename, typename, typename> class Converter,
          class Backend, class Model = BasicFlatModel< > >
using Interface = ConverterImpl<Converter, Backend, Model>;

} // namespace mp

#endif // CONVERTER_FLAT_H