#ifndef MP2MIP_H
#define MP2MIP_H

#include "mp/flat/converter.h"
#include "mp/flat/redef/MIP/redefs_mip_std.h"
#include "mp/flat/redef/encodings.h"

namespace mp {

/// MIPFlatConverter: converts flattened expressions for MIP
template <class Impl, class ModelAPI,
          class Model = BasicFlatModel< > >
class MIPFlatConverter
    : public FlatConverter<Impl, ModelAPI, Model>
{
public:
  /// Class name for diagnostics
  static constexpr const char* GetTypeName() { return "MIPFlatConverter"; }

  /// BaseConverter typedef
  using BaseConverter = FlatConverter<Impl, ModelAPI, Model>;

  /// Constructor
  MIPFlatConverter(Env& e) : BaseConverter(e) {  }


  ///////////////////// MIP CONSTRAINT CONVERTERS ///////////////////////
  /// Reuse default (empty) ones
  USE_BASE_CONSTRAINT_CONVERTERS( BaseConverter );

  /// Abs
  INSTALL_ITEM_CONVERTER(AbsConverter_MIP)
  /// And
  INSTALL_ITEM_CONVERTER(AndConverter_MIP)
  /// AllDiff
  INSTALL_ITEM_CONVERTER(AllDiffConverter_MIP)
  /// Complementarity linear
  INSTALL_ITEM_CONVERTER(ComplCvtLin_MIP)
  /// Complementarity quadratic
  INSTALL_ITEM_CONVERTER(ComplCvtQuad_MIP)
  /// Count
  INSTALL_ITEM_CONVERTER(CountConverter_MIP)
  /// Div
  INSTALL_ITEM_CONVERTER(DivConverter_MIP)
  /// CondLinEQ
  INSTALL_ITEM_CONVERTER(CondLinEQConverter_MIP)
  /// CondQuadEQ
  INSTALL_ITEM_CONVERTER(CondQuadEQConverter_MIP)
  /// CondLinLE
  INSTALL_ITEM_CONVERTER(CondLinLEConverter_MIP)
  /// CondQuadLE
  INSTALL_ITEM_CONVERTER(CondQuadLEConverter_MIP)
  /// CondLinLT
  INSTALL_ITEM_CONVERTER(CondLinLTConverter_MIP)
  /// CondQuadLT
  INSTALL_ITEM_CONVERTER(CondQuadLTConverter_MIP)
  /// CondLinGE
  INSTALL_ITEM_CONVERTER(CondLinGEConverter_MIP)
  /// CondQuadGE
  INSTALL_ITEM_CONVERTER(CondQuadGEConverter_MIP)
  /// CondLinGT
  INSTALL_ITEM_CONVERTER(CondLinGTConverter_MIP)
  /// CondQuadGT
  INSTALL_ITEM_CONVERTER(CondQuadGTConverter_MIP)
  /// IfThenElse
  INSTALL_ITEM_CONVERTER(IfThenElseConverter_MIP)
  /// ImplLE
  INSTALL_ITEM_CONVERTER(IndicatorLinLEConverter_MIP)
  /// ImplEQ
  INSTALL_ITEM_CONVERTER(IndicatorLinEQConverter_MIP)
  /// ImplGE
  INSTALL_ITEM_CONVERTER(IndicatorLinGEConverter_MIP)
  /// ImplQuadLE
  INSTALL_ITEM_CONVERTER(IndicatorQuadLEConverter_MIP)
  /// ImplQuadEQ
  INSTALL_ITEM_CONVERTER(IndicatorQuadEQConverter_MIP)
  /// ImplQuadGE
  INSTALL_ITEM_CONVERTER(IndicatorQuadGEConverter_MIP)
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
  /// PLConstraint
  INSTALL_ITEM_CONVERTER(PLConverter_MIP)
  /// PowConstExp
  INSTALL_ITEM_CONVERTER(PowConstExponentConverter_MIP)
  /// SOS2
  INSTALL_ITEM_CONVERTER(SOS2Converter_MIP)

  /// Smooth nonlinear

  /// Exp
  INSTALL_ITEM_CONVERTER(FuncConConverter_MIP_Exp)
  /// ExpA
  INSTALL_ITEM_CONVERTER(FuncConConverter_MIP_ExpA)
  /// Log
  INSTALL_ITEM_CONVERTER(FuncConConverter_MIP_Log)
  /// LogA
  INSTALL_ITEM_CONVERTER(FuncConConverter_MIP_LogA)
  /// Sin
  INSTALL_ITEM_CONVERTER(FuncConConverter_MIP_Sin)
  /// Cos
  INSTALL_ITEM_CONVERTER(FuncConConverter_MIP_Cos)
  /// Tan
  INSTALL_ITEM_CONVERTER(FuncConConverter_MIP_Tan)
  /// ASin
  INSTALL_ITEM_CONVERTER(FuncConConverter_MIP_Asin)
  /// ACos
  INSTALL_ITEM_CONVERTER(FuncConConverter_MIP_Acos)
  /// ATan
  INSTALL_ITEM_CONVERTER(FuncConConverter_MIP_Atan)
  /// Sinh
  INSTALL_ITEM_CONVERTER(FuncConConverter_MIP_Sinh)
  /// Cosh
  INSTALL_ITEM_CONVERTER(FuncConConverter_MIP_Cosh)
  /// Tanh
  INSTALL_ITEM_CONVERTER(FuncConConverter_MIP_Tanh)
  /// ASinh
  INSTALL_ITEM_CONVERTER(FuncConConverter_MIP_Asinh)
  /// ACosh
  INSTALL_ITEM_CONVERTER(FuncConConverter_MIP_Acosh)
  /// ATanh
  INSTALL_ITEM_CONVERTER(FuncConConverter_MIP_Atanh)


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
  /////////////////////////// MIP CONSTRAINT MAPS ///////////////////////
public:
  USE_BASE_MAP_FINDERS( BaseConverter )

  /// Specialize MapFind for CondLinConEQ.
  /// We distinguish the case var==const
  int MapFind(const CondLinConEQ& eq0c) {
    const auto isVCC = IsVarConstCmp( eq0c );
    if (isVCC.first)                    // only var==const comparisons
      return MapFind__VarConstCmp(isVCC.second.first, isVCC.second.second);
    return MPD( MapFind__Impl(eq0c) );
  }

  /// Specialize MapInsert for CondLiConEQ
  bool MapInsert(const CondLinConEQ& eq0c, int i) {
    const auto isVCC = IsVarConstCmp( eq0c );
    if (isVCC.first)                    // only var==const comparisons
      return MapInsert__VarConstCmp(isVCC.second.first, isVCC.second.second, i);
    return MPD( MapInsert__Impl(eq0c, i) );
  }

  /// Convert MIP-specific maps
  void ConvertMaps() {
    BaseConverter::ConvertMaps();
    MP_DISPATCH( ConvertEqVarConstMaps() );
  }

  /// Whether to use equality encoding for \a var
  bool IfUseEqualityEncodingForVar(int var) const {
    if (!MPCD(is_var_integer(var)))
      return false;
    const auto lb_dbl = this->lb(var);
    const auto ub_dbl = this->ub(var);
    return (lb_dbl>std::numeric_limits<int>::min() &&
            ub_dbl<std::numeric_limits<int>::max()) &&
        (ub_dbl-lb_dbl <= 100000);
  }

  /// Obtain extended column \a k of ZZI encoding C^r
  std::vector<double>
  GetZZIExtendedColumn(int r, int k, int v0, int v1) {
    return GetExtendedColumn(*p_zzi_, r, k, v0, v1);
  }


private:
  /// For a single variable, map its equality comparisons
  /// for the comparison value (double), map the EQ0Constraint index
  using SingleVarEqConstMap = std::unordered_map<double, int>;
  /// A map keeping such maps for certain variables
  using VarsEqConstMap = std::unordered_map<int, SingleVarEqConstMap>;

  VarsEqConstMap map_vars_eq_const_;

  P_ZZI_Encoding p_zzi_ { MakeZZIEncoding() };


protected:
  int MapFind__VarConstCmp(int var, double val) {
    auto itVar = map_vars_eq_const_.find(var);
    if (map_vars_eq_const_.end() != itVar) {
      auto itCmp = itVar->second.find( val );
      if (itVar->second.end() != itCmp)
        /// Make sure we store the comparisons in CondConEQ's
        return itCmp->second;
    }
    return -1;
  }

  bool MapInsert__VarConstCmp(int var, double val, int i) {
    auto result = map_vars_eq_const_[var].
        insert( { val, i } );
    return result.second;
  }

  /// Result of IsVarConstCmp(): var, const
  using VarConstCmp = std::pair<int, double>;

  /// Check if \a con is a conditional var==const
  static std::pair<bool, VarConstCmp>
  IsVarConstCmp(const CondLinConEQ& con) {
    const auto& linEQ = con.GetConstraint();
    if (1==linEQ.size()) {
      assert(1.0==linEQ.coef(0));
      return { true, { linEQ.var(0), linEQ.rhs() } };
    }
    return { false, {} };
  }

  void ConvertEqVarConstMaps() {
    for (const auto& m: map_vars_eq_const_) {
      if (IfUseEqualityEncodingForVar(m.first))
        ConvertEqVarConstMap(m.first, m.second);
    } // Otherwise, indicators / big-Ms should have been applied.
  }   // But why the map then?

  void ConvertEqVarConstMap(int var, const SingleVarEqConstMap& map) {
    CreateUnaryEncoding(var, map);
  }

  void CreateUnaryEncoding(int var, const SingleVarEqConstMap& map) {
    const Model& model = MP_DISPATCH( GetModel() );
    if (!model.is_integer_var(var))
      MP_RAISE("MP2MIP: Equality encoding: comparing non-integer variables not implemented");
    const auto lb_dbl = this->lb(var);
    const auto ub_dbl = this->ub(var);
    if (lb_dbl==this->MinusInfty() || ub_dbl==this->Infty())
      MP_RAISE("MP2MIP: Equality-comparing unbounded variables not implemented");
    if (lb_dbl<std::numeric_limits<int>::min() || ub_dbl>std::numeric_limits<int>::max())
      MP_RAISE("MP2MIP: Equality-comparing variables with domain out of integer range not implemented");
    const int lb = (int)std::round(lb_dbl);
    const int ub = (int)std::round(ub_dbl);
    std::vector<int> unaryEncVars(ub-lb+1);
    int nTaken=0;
    for (int v=lb; v!=ub+1; ++v) {  // Running thru the map faster?
      auto itV = map.find(v);
      if (map.end() != itV) {
        ++nTaken;
        unaryEncVars[v-lb] =
          GET_CONSTRAINT_KEEPER(CondLinConEQ).GetResultVar(itV->second);
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

  ///////////////////////////////////////////////////////////////////////
  /////////////////////// OPTIONS /////////////////////////           ///
public:
  /// Init MIPFlatConverter options
  void InitOptions() {
    BaseConverter::InitOptions();
    InitOwnOptions();
  }

protected:
  double cmpEpsContinuous() const { return options_.cmpEps_; }

private:
  struct Options {
    double cmpEps_ = 1e-4;
  };
  Options options_;

  void InitOwnOptions() {
    this->GetEnv().AddOption("cvt:mip:eps cvt:cmp:eps",
                       "Tolerance for strict comparison of continuous variables for MIP. "
                       "Ensure larger than the solver's feasibility tolerance.",
                       options_.cmpEps_, 0.0, 1e100);
  }
};

} // namespace mp

#endif // MP2MIP_H
