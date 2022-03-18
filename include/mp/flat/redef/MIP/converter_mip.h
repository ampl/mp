#ifndef MP2MIP_H
#define MP2MIP_H

#include "mp/flat/converter.h"
#include "mp/flat/redef/MIP/redefs_mip_std.h"

namespace mp {

/// MIPFlatConverter: converts flattened expressions for MIP
template <class Impl, class Backend,
          class Model = BasicFlatModel< > >
class MIPFlatConverter
    : public FlatConverter<Impl, Backend, Model>
{
public:
  /// Class name for diagnostics
  static constexpr const char* GetConverterName() { return "MIPFlatConverter"; }

  /// BaseConverter typedef
  using BaseConverter = FlatConverter<Impl, Backend, Model>;

  /// Constructor
  MIPFlatConverter(Env& e) : BaseConverter(e) {  }


  ///////////////////// SPECIALIZED CONSTRAINT CONVERTERS //////////////////
  /// Reuse default (empty) ones
  USE_BASE_CONSTRAINT_CONVERTERS( BaseConverter );


  /// Strict comparison tolerance.
  /// Need a big eps to avoid misinterpretation,
  /// at least the solver's feasibility tolerance
  double ComparisonEps(int var) const {
    return MPCD(is_var_integer(var)) ? 1.0 : cmpEpsContinuous();
  }
  /// Strict comparison tolerance
  double ComparisonEps(var::Type vartype) const {
    return var::INTEGER==vartype ? 1.0 : cmpEpsContinuous();
  }


  ///////////////////////////////////////////////////////////////////////
  /////////////////////////// MAPS //////////////////////////
  ///
private:
  /// For a single variable, map its equality comparisons
  /// for the comparison value (double), map the EQ0Constraint index
  using SingleVarEqConstMap = std::unordered_map<double, int>;
  /// A map keeping such maps for certain variables
  using VarsEqConstMap = std::unordered_map<int, SingleVarEqConstMap>;

  VarsEqConstMap map_vars_eq_const_;


public:
  //////////////////////////////// CONSTRAINT MAPS ///////////////////////////////////////
  ///
  USE_BASE_MAP_FINDERS( BaseConverter )

  AbstractConstraintLocation MapFind(const EQ0Constraint& eq0c) {
    const auto isVCC = IsVarConstCmp( eq0c );
    if (isVCC.first) {                    // only var==const comparisons
      return MapFind__VarConstCmp(isVCC.second.first, isVCC.second.second);
    }
    return { };
  }

  AbstractConstraintLocation MapFind__VarConstCmp(int var, double val) {
    auto itVar = map_vars_eq_const_.find(var);
    if (map_vars_eq_const_.end() != itVar) {
      auto itCmp = itVar->second.find( val );
      if (itVar->second.end() != itCmp)
        /// Make sure we store the comparisons in EQ0Con's
        return {&GET_CONSTRAINT_KEEPER(EQ0Constraint), itCmp->second};
    }
    return { };
  }

  bool MapInsert(ConstraintLocation<EQ0Constraint> ck) {
    const auto isVCC = IsVarConstCmp( ck.GetConstraint() );
    if (isVCC.first) {                    // only var==const comparisons
      return MapInsert__VarConstCmp(isVCC.second.first, isVCC.second.second, ck);
    }
    return true;
  }

  bool MapInsert__VarConstCmp(int var, double val,
                              ConstraintLocation<EQ0Constraint> ck) {
    auto result = map_vars_eq_const_[var].
        insert( std::make_pair( val, ck.GetIndex() ) );
    return result.second;
  }

  /// Result of IsVarConstCmp(): var, const
  using VarConstCmp = std::pair<int, double>;

  /// Check if the eq0c is a var==const
  static std::pair<bool, VarConstCmp> IsVarConstCmp(const EQ0Constraint& con) {
    const AffExp& args = con.GetArguments();
    if (1==args.size()) {
      assert(1.0==args.coef(0));
      return { true, { args.var(0), -args.constant_term() } };
    }
    return { false, {} };
  }

  void ConvertMaps() {
    MP_DISPATCH( ConvertEqVarConstMaps() );
  }

  void ConvertEqVarConstMaps() {
    for (const auto& m: map_vars_eq_const_) {
      if (IfUseEqualityEncodingForVar(m.first))
        ConvertEqVarConstMap(m.first, m.second);
    } // Otherwise, indicators / big-Ms should have been applied
  }

  bool IfUseEqualityEncodingForVar(int var) const {
    if (!MPCD(is_var_integer(var)))
      return false;
    const auto lb_dbl = this->lb(var);
    const auto ub_dbl = this->ub(var);
    return (lb_dbl>std::numeric_limits<int>::min() && // TODO enable outside
            ub_dbl<std::numeric_limits<int>::max()) &&
        (ub_dbl-lb_dbl <= 100000);      // TODO param
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
    for (int v=lb; v!=ub+1; ++v) { // TODO run thru the map first
      auto itV = map.find(v);
      if (map.end() != itV) {
        ++nTaken;
        unaryEncVars[v-lb] =
          GET_CONSTRAINT_KEEPER(EQ0Constraint).GetResultVar(itV->second);
      } else {
        unaryEncVars[v-lb] = this->AddVar(0.0, 1.0, var::INTEGER);
      }
    }
    assert(map.size()==(size_t)nTaken);
    std::vector<double> coefs(ub_dbl-lb_dbl+1, 1.0);
    this->AddConstraint(LinConEQ({coefs, unaryEncVars}, 1.0));
    unaryEncVars.push_back(var);
    for (int v=lb; v!=ub+1; ++v) {
      coefs[v-lb] = v;
    }
    coefs.push_back(-1.0);
    this->AddConstraint(LinConEQ({coefs, unaryEncVars}, 0.0));
  }

public:
  /////////////////////// CONSTRAINT CONVERTERS /////////////////////////

  /// Abs
  INSTALL_ITEM_CONVERTER(AbsConverter_MIP)
  /// And
  INSTALL_ITEM_CONVERTER(AndConverter_MIP)
  /// AllDiff
  INSTALL_ITEM_CONVERTER(AllDiffConverter_MIP)
  /// Count
  INSTALL_ITEM_CONVERTER(CountConverter_MIP)
  /// EQ0
  INSTALL_ITEM_CONVERTER(EQ0Converter_MIP)
  /// IfThenElse
  INSTALL_ITEM_CONVERTER(IfThenElseConverter_MIP)
  /// ImplLE0
  INSTALL_ITEM_CONVERTER(IndicatorLinLEConverter_MIP)
  /// ImplEQ0
  INSTALL_ITEM_CONVERTER(IndicatorLinEQConverter_MIP)
  /// LE0
  INSTALL_ITEM_CONVERTER(LE0Converter_MIP)
  /// LT0
  INSTALL_ITEM_CONVERTER(LT0Converter_MIP)
  /// Min
  INSTALL_ITEM_CONVERTER(MinConverter_MIP)
  /// Max
  INSTALL_ITEM_CONVERTER(MaxConverter_MIP)
  /// Not
  INSTALL_ITEM_CONVERTER(NotConverter_MIP)
  /// NumberofConst
  INSTALL_ITEM_CONVERTER(NumberofConstConverter_MIP)
  /// NumberofVar
  INSTALL_ITEM_CONVERTER(NumberofVarConverter_MIP)
  /// Or
  INSTALL_ITEM_CONVERTER(OrConverter_MIP)


  ///////////////////////////////////////////////////////////////////////
  /////////////////////// OPTIONS /////////////////////////
  ///
public:
  void InitOptions() {
    BaseConverter::InitOptions();
    InitOwnOptions();
  }

protected:
  double cmpEpsContinuous() const { return options_.cmpEps_; }

private:
  struct Options {
    double cmpEps_ = 1e-3;
  };
  Options options_;

  void InitOwnOptions() {
    this->GetEnv().AddOption("cvt:mip:eps",
                       "Tolerance for strict comparison of continuous variables for MIP. "
                       "Ensure larger than the solver's feasibility tolerance.",
                       options_.cmpEps_, 0.0, 1e100);
  }
};

} // namespace mp

#endif // MP2MIP_H
