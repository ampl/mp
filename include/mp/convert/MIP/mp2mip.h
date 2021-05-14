#ifndef MP2MIP_H
#define MP2MIP_H

#include "mp/convert/converter_flat.h"

namespace mp {

/// MPToMIPConverter: one of the converters requiring a "minimal" output interface
template <class Impl, class Backend,
          class Model = BasicModel< > >
class MPToMIPConverter
    : public BasicMPFlatConverter<Impl, Backend, Model>
{
public:
  using BaseConverter = BasicMPFlatConverter<Impl, Backend, Model>;
  template <class Constraint>
    using ConstraintKeeperType = typename
      BaseConverter::template ConstraintKeeperType<Constraint>;

public:
  static const char* GetConverterName() { return "MPToMIPConverter"; }
  MPToMIPConverter() {
    InitOptions();
  }

  ///////////////////// SPECIALIZED CONSTRAINT CONVERTERS //////////////////
  USE_BASE_CONSTRAINT_CONVERTERS( BaseConverter )           // reuse default ones

  template <int sense, class MinOrMaxConstraint>
  void ConvertMinOrMax(const MinOrMaxConstraint& mc) {
    const auto& args = mc.GetArguments();
    const std::size_t nargs = args.size();
    const auto flags = this->AddVars(nargs, 0.0, 1.0, var::Type::INTEGER);   // binary flags
    MP_DISPATCH( AddConstraint(LinearConstraint(std::vector<double>(nargs, 1.0),    // sum of the flags >= 1
                        flags, 1.0, this->Infty())) );
    const auto resvar = mc.GetResultVar();
    for (int i=0; i<nargs; ++i) {
      MP_DISPATCH( AddConstraint(LinearConstraint({1.0*sense, -1.0*sense},
                          {args[i], resvar}, this->MinusInfty(), 0.0)) );
      MP_DISPATCH( AddConstraint(IndicatorConstraintLinLE{flags[i], 1,
                          {1.0*sense, -1.0*sense}, {resvar, args[i]}, 0.0}) );
    }
  }

  void Convert(const MaximumConstraint& mc) {
    ConvertMinOrMax<1>(mc);
  }

  void Convert(const MinimumConstraint& mc) {
    ConvertMinOrMax<-1>(mc);
  }

  void Convert(const AbsConstraint& ac) {
    const int arg = ac.GetArguments()[0];
    const int res = ac.GetResultVar();
    this->AddConstraint(LinearConstraint({1.0, 1.0}, {res, arg}, 0.0, this->Infty()));
    this->AddConstraint(LinearConstraint({1.0, -1.0}, {res, arg}, 0.0, this->Infty()));
    const int flag = this->AddVar(0.0, 1.0, var::INTEGER);
    this->AddConstraint(IndicatorConstraintLinLE(flag, 1, {1.0, 1.0}, {res, arg}, 0.0));
    this->AddConstraint(IndicatorConstraintLinLE(flag, 0, {1.0, -1.0}, {res, arg}, 0.0));
  }

  void Convert(const NotConstraint& nc) {
    MP_DISPATCH( AddConstraint(LinearDefiningConstraint(
      nc.GetResultVar(), {{-1.0}, {nc.GetArguments()[0]}, 1.0})) );
  }

  void Convert(const LE0Constraint& le0c) {
    auto& m = this->GetModel();
    if (m.is_fixed(le0c.GetResultVar()))
      throw std::logic_error("LE0Constraint: result fixed, not implemented");
    assert(!le0c.GetContext().IsNone());
    if (le0c.GetContext().HasPositive())
      ConvertImplied(le0c);
    if (le0c.GetContext().HasNegative())
      ConvertReverseImplied(le0c);
  }

  void ConvertImplied(const LE0Constraint& le0c) {
    auto& m = this->GetModel();
    const auto& ae = le0c.GetArguments();
    if (ae.is_constant()) {
      if (ae.constant_term() > 0.0)
        m.narrow_var_bounds(le0c.GetResultVar(), 0.0, 0.0);
    } else {
      LinearExprUnzipper le(ae);
      MP_DISPATCH( AddConstraint(IndicatorConstraintLinLE(
                   le0c.GetResultVar(), 1,
                   le.coefs(), le.var_indexes(), -ae.constant_term())) );
    }
  }

  void ConvertReverseImplied(const LE0Constraint& le0c) {
    auto& m = this->GetModel();
    auto ae = le0c.GetArguments();
    ae.Negate();
    if (ae.is_constant()) {
      if (ae.constant_term() >= 0.0)
        m.narrow_var_bounds(le0c.GetResultVar(), 1.0, 1.0);
    } else {
      auto bNt = MP_DISPATCH( ComputeBoundsAndType(ae) );
      double cmpEps = var::INTEGER==bNt.get_result_type() ? 1.0 : 1e-6;
      double d = ae.constant_term() + cmpEps;
      LinearExprUnzipper le(ae);
      MP_DISPATCH( AddConstraint(IndicatorConstraintLinLE(
                   le0c.GetResultVar(), 0,
                   le.coefs(), le.var_indexes(), -d)) );
    }
  }

  double ComparisonEpsilon(int var) {
    return (MP_DISPATCH( GetModel() ).is_integer_var(var)) ? 1.0 : 1e-6; // TODO param
  }

  void Convert(const EQ0Constraint& eq0c) {
    auto& m = this->GetModel();
    if (m.is_fixed(eq0c.GetResultVar()))
      throw std::logic_error("EQ0Constraint: result fixed, not implemented");
    /// TODO assert(!eq0c.GetContext().IsNone());
    if (1<eq0c.GetArguments().num_terms()) {
      throw std::logic_error("MIP conversion of linexp==0 comparison not implemented");
    }
    /// Stop here, for var==const rest done by postprocessing
  }

  void Convert(const IndicatorConstraintLinLE& indc) {
    auto binvar=indc.get_binary_var();
    if (indc.is_binary_value_1())                  /// If binval==1, complement the variable
      binvar = this->MakeComplementVar(binvar);
    /// Convert indc's linear inequality to 'cmpvar<=0'
    /// Could use the full inequality instead of the new var
    int cmpvar = MP_DISPATCH( Convert2Var(indc.to_lhs_affine_expr()) );
    if (this->ub(cmpvar) >= this->Infty())
      throw ConstraintConversionFailure("Cannot convert indicator constraint with variable " +
                             std::to_string(cmpvar) + " having infinite upper bound."
                             " Define finite upper bound or use solver built-in indicator");
    MP_DISPATCH( AddConstraint(LinearConstraint(          /// Big-M constraint cmpvar <= ub(cmpvar)*binvar
        {1.0, -this->ub(cmpvar)}, {cmpvar, binvar}, this->MinusInfty(), 0.0)) );
  }

  void Convert(const DisjunctionConstraint& disj) {
    assert(!disj.GetContext().IsNone());
    if (disj.GetContext().HasPositive())
      ConvertImplied(disj);
    if (disj.GetContext().HasNegative())
      ConvertReverseImplied(disj);
  }

  void ConvertImplied(const DisjunctionConstraint& disj) {
    const auto& args = disj.GetArguments();
    auto flags = args;
    flags.push_back(disj.GetResultVar());
    std::vector<double> ones(args.size(), 1.0);
    ones.push_back(-1.0);
    MP_DISPATCH( AddConstraint(LinearConstraint(ones, flags, 0.0, MP_DISPATCH( Infty() ))) );
  }

  void ConvertReverseImplied(const DisjunctionConstraint& disj) {
    std::array<double, 2> coefs{1.0, -1.0};
    std::array<int, 2> vars{-1, disj.GetResultVar()};
    for (auto arg: disj.GetArguments()) {
      vars[0] = arg;
      MP_DISPATCH( AddConstraint(LinearConstraint(coefs, vars, MP_DISPATCH( MinusInfty() ), 0.0)) );
    }
  }

  void Convert(const IfThenConstraint& itc) {
    assert(!itc.GetContext().IsNone());
    const auto& args = itc.GetArguments();
    if (!this->is_fixed(args[1]) || !this->is_fixed(args[2]))
      throw std::logic_error("MP2MIP: IfThen with variable then/else arguments not implemented");
    else
      ConvertIfThen_constantThenElse(itc);
  }

  void ConvertIfThen_constantThenElse(const IfThenConstraint& itc) {
    const auto& args = itc.GetArguments();
    assert((this->is_fixed(args[1]) && this->is_fixed(args[2])));
    const double const1 = this->fixed_value(args[1]);
    const double const2 = this->fixed_value(args[2]);
    this->AddConstraint( LinearDefiningConstraint(
                           itc.GetResultVar(),
    { {const1-const2}, {args[0]}, const2 } ) );
  }

  //////////////////// ALLDIFF ///////////////////////
  void Convert(const AllDiffConstraint& alld) {
    if (!this->is_fixed(alld.GetResultVar()) ||
        1.0!=this->fixed_value(alld.GetResultVar()))
      throw std::logic_error("MP2MIP: only static alldifferent implemented.");
    else
      Convert_staticAllDiff(alld);
  }

  void Convert_staticAllDiff(const AllDiffConstraint& alld) {
    assert (this->is_fixed(alld.GetResultVar()) &&
        1.0==this->fixed_value(alld.GetResultVar()));
    const auto& args = alld.GetArguments();
    const auto lba_dbl = this->lb_array(args);
    const auto uba_dbl = this->ub_array(args);
    if (lba_dbl==this->MinusInfty() || uba_dbl==this->Infty())
      throw std::logic_error("MP2MIP: AllDiff on unbounded variables not implemented");
    if (lba_dbl<std::numeric_limits<int>::min() || uba_dbl>std::numeric_limits<int>::max())
      throw std::logic_error("MP2MIP: AllDiff on variables with domain out of integer range not implemented");
    const int lba = (int)std::round(lba_dbl);
    const int uba = (int)std::round(uba_dbl);
    std::vector<double> coefs(args.size(), 1.0);
    std::vector<int> flags(args.size());                // unary encoding flags
    for (int v=lba; v!=uba+1; ++v) {                    // for each value in the domain union
      for (size_t ivar = 0; ivar < args.size(); ++ivar) {
        flags[ivar] = this->AssignResultToArguments(
              EQ0Constraint( { {1.0}, {args[ivar]}, -double(v) } ) ).
            get_representing_variable();
      }
      this->AddConstraint( LinearConstraint(
          coefs, flags, this->MinusInfty(), 1.0 ) );
    }
  }

  ///////////////////////////////////////////////////////////////////////
  /////////////////////////// MAPS //////////////////////////
  ///
private:
  /// For a single variable, map its equality comparisons
  using SingleVarEqConstMap = std::unordered_map<double,
                     const ConstraintKeeperType<EQ0Constraint>*>;
  /// A map keeping such maps for certain variables
  using VarsEqConstMap = std::unordered_map<int, SingleVarEqConstMap>;

  VarsEqConstMap map_vars_eq_const_;

public:
  ///////////////////////////////////////////////////////////////////////
  ///
  USE_BASE_MAP_FINDERS( BaseConverter )

  const BasicConstraintKeeper* MapFind(const EQ0Constraint& eq0c) const {
    const auto isVCC = IsVarConstCmp( eq0c );
    if (isVCC.first) {                    // only var==const comparisons
      return MapFind__VarConstCmp(isVCC.second.first, isVCC.second.second);
    }
    return nullptr;
  }

  const BasicConstraintKeeper* MapFind__VarConstCmp(int var, double val) const {
    auto itVar = map_vars_eq_const_.find(var);
    if (map_vars_eq_const_.end() != itVar) {
      auto itCmp = itVar->second.find( val );
      if (itVar->second.end() != itCmp)
        return itCmp->second;
    }
    return nullptr;
  }

  bool MapInsert(const ConstraintKeeperType<EQ0Constraint>* pck) {
    const auto isVCC = IsVarConstCmp( pck->GetConstraint() );
    if (isVCC.first) {                    // only var==const comparisons
      return MapInsert__VarConstCmp(isVCC.second.first, isVCC.second.second, pck);
    }
    return true;
  }

  bool MapInsert__VarConstCmp(int var, double val,
                              const ConstraintKeeperType<EQ0Constraint>* pck) {
    auto result = map_vars_eq_const_[var].
        insert( std::make_pair( val, pck ) );
    return result.second;
  }

  using VarConstCmp = std::pair<int, double>;
  static std::pair<bool, VarConstCmp> IsVarConstCmp(const EQ0Constraint& cons) {
    const AffineExpr& args = cons.GetArguments();
    if (1==args.num_terms()) {
      assert(1.0==args.coef(0));
      return { true, { args.var_index(0), -args.constant_term() } };
    }
    return { false, {} };
  }

  void ConvertMaps() {
    MP_DISPATCH( ConvertEqVarConstMaps() );
  }

  void ConvertEqVarConstMaps() {
    for (const auto& m: map_vars_eq_const_) {
      ConvertEqVarConstMap(m.first, m.second);
    }
  }

  void ConvertEqVarConstMap(int var, const SingleVarEqConstMap& map) {
    CreateUnaryEncoding(var, map);
  }

  void CreateUnaryEncoding(int var,  const SingleVarEqConstMap& map) {
    const Model& model = MP_DISPATCH( GetModel() );
    if (!model.is_integer_var(var))
      throw std::logic_error("MP2MIP: Equality-comparing non-integer variables not implemented");
    const auto lb_dbl = this->lb(var);
    const auto ub_dbl = this->ub(var);
    if (lb_dbl==this->MinusInfty() || ub_dbl==this->Infty())
      throw std::logic_error("MP2MIP: Equality-comparing unbounded variables not implemented");
    if (lb_dbl<std::numeric_limits<int>::min() || ub_dbl>std::numeric_limits<int>::max())
      throw std::logic_error("MP2MIP: Equality-comparing variables with domain out of integer range not implemented");
    const int lb = (int)std::round(lb_dbl);
    const int ub = (int)std::round(ub_dbl);
    std::vector<int> unaryEncVars(ub-lb+1);
    int nTaken=0;
    for (int v=lb; v!=ub+1; ++v) {
      auto itV = map.find(v);
      if (map.end() != itV) {
        ++nTaken;
        unaryEncVars[v-lb] = itV->second->GetResultVar();
      } else {
        unaryEncVars[v-lb] = this->AddVar(0.0, 1.0, var::INTEGER);
      }
    }
    assert(map.size()==(size_t)nTaken);
    std::vector<double> coefs(ub_dbl-lb_dbl+1, 1.0);
    this->AddConstraint(LinearConstraint(coefs, unaryEncVars, 1.0, 1.0));
    unaryEncVars.push_back(var);
    for (int v=lb; v!=ub+1; ++v) {
      coefs[v-lb] = v;
    }
    coefs.push_back(-1.0);
    this->AddConstraint(LinearConstraint(coefs, unaryEncVars, 0.0, 0.0));
  }

  ///////////////////////////////////////////////////////////////////////
  /////////////////////// OPTIONS /////////////////////////
  ///
private:
  struct Options {
  };
  Options options_;

  void InitOptions() {
    this->add_to_long_name(" with MP-to-MIP Converter Layer");
    this->add_to_version("\nMP-to-MIP Converter Layer for AMPL");
    this->add_to_option_header(
          "\n"
          "Including MP-to-MIP Converter Layer Options\n"
          "-------------------------------------------\n"
          );
  }

};

} // namespace mp

#endif // MP2MIP_H
