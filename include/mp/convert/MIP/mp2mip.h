#ifndef MP2MIP_H
#define MP2MIP_H

#include "mp/convert/converter_flat.h"

namespace mp {

/// MPToMIPConverter: one of the converters requiring a "minimal" output interface
template <class Impl, class Backend,
          class Model = BasicModel<std::allocator<char> > >
class MPToMIPConverter
    : public BasicMPFlatConverter<Impl, Backend, Model>
{
public:
  MPToMIPConverter() {
    InitOptions();
  }

  using BaseConverter = BasicMPFlatConverter<Impl, Backend, Model>;

  ///////////////////// SPECIALIZED CONSTRAINT CONVERTERS //////////////////
  USE_BASE_CONSTRAINT_CONVERTERS( BaseConverter )           // reuse default ones

  template <int sense, class MinOrMaxConstraint>
  void ConvertMinOrMax(const MinOrMaxConstraint& mc) {
    const auto& args = mc.GetArguments();
    const int nargs = args.size();
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

  void Convert(const NotConstraint& nc) {
    MP_DISPATCH( AddConstraint(LinearDefiningConstraint(
      nc.GetResultVar(), {{-1.0}, {nc.GetArguments()[0]}, 1.0})) );
  }

  void Convert(const LE0Constraint& le0c) {
    auto& m = this->GetModel();
    if (m.is_fixed(le0c.GetResultVar()))
      throw std::logic_error("LEConstraint: result fixed, not implemented");
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
      throw std::logic_error("LEConstraint: result fixed, not implemented");
    assert(!eq0c.GetContext().IsNone());
    /// Stop here, rest done by postprocessing
  }

  void Convert(const IndicatorConstraintLinLE& indc) {
    auto binvar=indc.get_binary_var();
    if (indc.is_binary_value_1())                  /// If binval==1, complement the variable
      binvar = this->MakeComplementVar(binvar);
    /// Convert indc's linear inequality to 'cmpvar<=0'
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


  ///////////////////////////////////////////////////////////////////////
  /////////////////////////// MAPS //////////////////////////
  ///
private:
  /// For a single variable, map its equality comparisons
  using SingleVarEqConstMap = std::unordered_map<double,
                     const ConstraintKeeper<Impl, Backend, EQ0Constraint>*>;
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
      auto itVar = map_vars_eq_const_.find(isVCC.second.first);
      if (map_vars_eq_const_.end() != itVar) {
        auto itCmp = itVar->second.find( isVCC.second.second );
        if (itVar->second.end() != itCmp)
          return itCmp->second;
      }
    }
    return nullptr;
  }

  bool MapInsert(const ConstraintKeeper<Impl, Backend, EQ0Constraint>* pck) {
    const auto isVCC = IsVarConstCmp( pck->GetConstraint() );
    if (isVCC.first) {                    // only var==const comparisons
      auto result = map_vars_eq_const_[isVCC.second.first].
          insert( std::make_pair( isVCC.second.second, pck ) );
      return result.second;
    }
    return true;
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
      throw std::logic_error("Equality-comparing non-integer variables not implemented");
    const auto lb_dbl = this->lb(var);
    const auto ub_dbl = this->ub(var);
    if (lb_dbl==this->MinusInfty() || ub_dbl==this->Infty())
      throw std::logic_error("Equality-comparing unbounded variables not implemented");
    if (lb_dbl<std::numeric_limits<int>::min() || ub_dbl>std::numeric_limits<int>::max())
      throw std::logic_error("Equality-comparing variables with domain out of integer range not implemented");
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
