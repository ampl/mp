#ifndef MP2MIP_H
#define MP2MIP_H

#include "mp/flat/converter.h"

#include "mp/flat/redef/MIP/alldiff.h"

namespace mp {

/// MIPFlatConverter: converts flattened expressions for MIP
template <class Impl, class Backend,
          class Model = BasicFlatModel< > >
class MIPFlatConverter
    : public FlatConverter<Impl, Backend, Model>
{
public:
  static constexpr const char* name() { return "MIPFlatConverter"; };

public:
  using BaseConverter = FlatConverter<Impl, Backend, Model>;

public:
  static const char* GetConverterName() { return "MIPFlatConverter"; }
  MIPFlatConverter(Env& e) : BaseConverter(e) {  }

  ///////////////////// SPECIALIZED CONSTRAINT CONVERTERS //////////////////
  USE_BASE_CONSTRAINT_CONVERTERS( BaseConverter );        ///< reuse default ones

  template <int sense, class MinOrMaxConstraint>
  void ConvertMinOrMax(const MinOrMaxConstraint& mc) {
    const auto& args = mc.GetArguments();
    const std::size_t nargs = args.size();
    const auto flags =
        this->AddVars_returnIds(nargs, 0.0, 1.0, var::Type::INTEGER);   // binary flags
    MP_DISPATCH( AddConstraint(LinConGE(
                                 { std::vector<double>(nargs, 1.0),  // sum of the flags >= 1
                                   flags }, 1.0)) );
    const auto resvar = mc.GetResultVar();
    for (size_t i=0; i<nargs; ++i) {
      MP_DISPATCH( AddConstraint(LinConLE(
                         { {1.0*sense, -1.0*sense},
                           {args[i], resvar} }, 0.0)) );
      MP_DISPATCH( AddConstraint(
                     IndicatorConstraintLinLE{flags[i], 1,
                         { { {1.0*sense, -1.0*sense}, {resvar, args[i]} }, 0.0 }}) );
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
    this->AddConstraint(LinConGE({{1.0, 1.0}, {res, arg}}, {0.0}));
    this->AddConstraint(LinConGE({{1.0, -1.0}, {res, arg}}, {0.0}));
    const int flag = this->AddVar(0.0, 1.0, var::INTEGER);
    this->AddConstraint(
          IndicatorConstraintLinLE(flag, 1, {{{1.0, 1.0}, {res, arg}}, 0.0}));
    this->AddConstraint(
          IndicatorConstraintLinLE(flag, 0, {{{1.0, -1.0}, {res, arg}}, 0.0}));
  }

  void Convert(const NotConstraint& nc) {
    MP_DISPATCH( AddConstraint(LinearFunctionalConstraint(
      nc.GetResultVar(), {{{-1.0}, {nc.GetArguments()[0]}}, 1.0})) );
  }

  void Convert(const LE0Constraint& le0c) {
    assert(!le0c.GetContext().IsNone());
    if (le0c.GetContext().HasPositive())
      ConvertImplied(le0c);
    if (le0c.GetContext().HasNegative())
      ConvertReverseImplied(le0c);
  }

  void ConvertImplied(const LE0Constraint& le0c) {
    const auto& ae = le0c.GetArguments();
    if (ae.is_constant()) {
      if (ae.constant_term() > 0.0)
        MPD( NarrowVarBounds(le0c.GetResultVar(), 0.0, 0.0) );
    } else {
      if (MPD(is_fixed(le0c.GetResultVar()))) {
        if (MPD(fixed_value(le0c.GetResultVar()))) {      // fixed to 1
          MP_DISPATCH( AddConstraint( ExtractConstraint(le0c) ) );
        }
      } else {
        MP_DISPATCH( AddConstraint(IndicatorConstraintLinLE(
                                     le0c.GetResultVar(), 1,
                                     ExtractConstraint(le0c))) );
      }
    }
  }

  /// b==0 --> c'x > d
  void ConvertReverseImplied(const LE0Constraint& le0c) {
    auto ae = le0c.GetArguments();
    if (ae.is_constant()) {
      if (ae.constant_term() <= 0.0)
        MPD( NarrowVarBounds(le0c.GetResultVar(), 1.0, 1.0) );
    } else {
      if (MPD(is_fixed(le0c.GetResultVar()))) {
        if (!MPD(fixed_value(le0c.GetResultVar()))) {      // fixed to 0
          MP_DISPATCH( AddConstraint(LinConGE(
                                       LinTerms(ae),
                                       -ae.constant_term()+1)) );
        }
      } else {
        ae.negate();
        auto bNt = MP_DISPATCH( ComputeBoundsAndType(ae) );
        double cmpEps = var::INTEGER==bNt.get_result_type() ? 1.0 : 1e-6;
        double d = ae.constant_term() + cmpEps;
        MP_DISPATCH( AddConstraint(IndicatorConstraintLinLE(
                                     le0c.GetResultVar(), 0,
                                     { LinTerms(ae), -d })) );
      }
    }
  }

  double ComparisonEps(int var) const {
    return MPCD(is_var_integer(var)) ? 1.0 : 1e-6; // TODO param
  }
  double ComparisonEps(var::Type vartype) const {
    return var::INTEGER==vartype ? 1.0 : 1e-6; // TODO param
  }

  void Convert(const EQ0Constraint& eq0c) {
    assert(!eq0c.GetContext().IsNone());
    const auto& args = eq0c.GetArguments();
    if (1<args.size() ||
        !IfUseEqualityEncodingForVar(
          args.var(0))) {
      auto ctx = eq0c.GetContext();
      if (ctx.HasPositive())
        ConvertImplied(eq0c);
      if (ctx.HasNegative())
        ConvertReverseImplied(eq0c);
    } // else, using unary encoding whose flags are,
  }   // in the fixed case, fixed by PropagateResult()

  /// resvar==1 => c'x==d
  void ConvertImplied(const EQ0Constraint& eq0c) {
    const auto& ae = eq0c.GetArguments();
    if (ae.is_constant()) {         // TODO consider resvar+context
      if (std::fabs(ae.constant_term()) != 0.0)
        MPD( NarrowVarBounds(eq0c.GetResultVar(), 0.0, 0.0) );
    } else {
      if (MPD(is_fixed(eq0c.GetResultVar()))) {
        if (MPD(fixed_value(eq0c.GetResultVar()))) {     // fixed to 1
          MP_DISPATCH( AddConstraint( ExtractConstraint(eq0c) ) );
        } // else, skip
      }
      MP_DISPATCH( AddConstraint(IndicatorConstraintLinEQ(
                   eq0c.GetResultVar(), 1,
                   ExtractConstraint(eq0c))) );
    }
  }

  /// resvar==0 ==> c'x!=d
  void ConvertReverseImplied(const EQ0Constraint& eq0c) {
    auto ae = eq0c.GetArguments();
    if (ae.is_constant()) {
      if (std::fabs(ae.constant_term()) != 0.0)
        MPD( NarrowVarBounds(eq0c.GetResultVar(), 1.0, 1.0) );
      // TODO use resvar + context
    } else if ( !MPD(is_fixed(eq0c.GetResultVar())) ||       // not fixed, or
                !MPD(fixed_value(eq0c.GetResultVar())) )     // fixed to 0
    { // TODO We are in MIP so doing algebra, not DisjunctiveConstr. Why?
      // Well in party1.mod, although this results in more fixed variables,
      // Gurobi 9.5 runs 31s vs 91s.
      auto newvars = MPD( AddVars_returnIds(2, 0.0, 1.0, var::INTEGER) );
      newvars.push_back( eq0c.GetResultVar() );
      MPD( AddConstraint( LinConGE(   // b1+b2+resvar >= 1
                            {{1.0, 1.0, 1.0}, newvars},
                                         1.0 ) ) );
      auto bNt = MP_DISPATCH( ComputeBoundsAndType(ae) );
      double cmpEps = MPD( ComparisonEps( bNt.get_result_type() ) );
      {
        MP_DISPATCH( AddConstraint(IndicatorConstraintLinLE(
                                     newvars[0], 1,
                                     { {ae.coefs(), ae.vars()},
                                       -ae.constant_term() - cmpEps })) );
      }
      ae.negate();
      MP_DISPATCH( AddConstraint(IndicatorConstraintLinLE(
                                   newvars[1], 1,
                                   { {ae.coefs(), ae.vars()},
                                     -ae.constant_term() - cmpEps })) );
    } // else, skip
  }



  void Convert(const IndicatorConstraintLinLE& indc) {
    auto binvar=indc.get_binary_var();
    auto ae = indc.to_lhs_expr();
    auto bnds = MPD( ComputeBoundsAndType(ae) );
    ConvertImplicationLE(binvar, indc.get_binary_value(), bnds, std::move(ae));
  }

  /// (b==val ==> ae<=0)
  void ConvertImplicationLE(int b, int val,
                   const PreprocessInfoStd& bnds, AffExp ae) {
    /// TODO fail if lb>0 +report .iis if requested
    /// TODO skip if ub<0
    if (bnds.ub() >= this->PracticallyInfty())
      throw ConstraintConversionFailure( "IndicatorInfBound",
          "The redefinition of a (possibly auxiliary) indicator constraint failed"
          " so it had to be passed to the solver."
          " Provide tight bounds on variables entering logical expressions, "
          "or set acc:ind_le=2");
    if (val)            // left condition is b==1
      ae += {{bnds.ub(), b}, -bnds.ub()};
    else
      ae += {{-bnds.ub(), b}, 0.0};
    MP_DISPATCH( AddConstraint(LinConLE(     /// Big-M constraint
        (LinTerms&&)ae, -ae.constant_term() )) );
  }

  /// b==val ==> c'x==d
  void Convert(const IndicatorConstraintLinEQ& indc) {
    auto binvar=indc.get_binary_var();
    auto ae = indc.to_lhs_expr();
    auto bnds = MPD( ComputeBoundsAndType(ae) );
    ConvertImplicationLE(binvar, indc.get_binary_value(), bnds, ae);
    ae.negate();
    bnds.NegateBounds();
    ConvertImplicationLE(binvar, indc.get_binary_value(), bnds, std::move(ae));
  }

  void Convert(const ConjunctionConstraint& conj) {
    assert(!conj.GetContext().IsNone());
    if (conj.GetContext().HasPositive())
      ConvertImplied(conj);
    if (conj.GetContext().HasNegative())
      ConvertReverseImplied(conj);
  }

  void ConvertReverseImplied(const ConjunctionConstraint& conj) {
    const auto& args = conj.GetArguments();
    auto flags = args;
    flags.push_back(conj.GetResultVar());
    std::vector<double> ones(args.size(), 1.0); // res+n-1 >= sum(args) in CTX-
    ones.push_back(-1.0);
    MP_DISPATCH( AddConstraint(
                   LinConLE({ones, flags},
                               {(double)args.size()-1} )) );
  }

  void ConvertImplied(const ConjunctionConstraint& conj) {
    std::array<double, 2> coefs{-1.0, 1.0};
    std::array<int, 2> vars{-1, conj.GetResultVar()};
    for (auto arg: conj.GetArguments()) {       // res <= arg[i] in CTX+
      vars[0] = arg;
      MP_DISPATCH( AddConstraint(
                     LinConLE({coefs, vars}, {0.0} )) );
    }
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
    std::vector<double> ones(args.size(), 1.0);  // res <= sum(args) in CTX+
    ones.push_back(-1.0);
    MP_DISPATCH( AddConstraint(
                   LinConGE({ones, flags}, {0.0} )) );
  }

  void ConvertReverseImplied(const DisjunctionConstraint& disj) {
    std::array<double, 2> coefs{1.0, -1.0};
    std::array<int, 2> vars{-1, disj.GetResultVar()};
    for (auto arg: disj.GetArguments()) {        // res >= arg[i] in CTX-
      vars[0] = arg;
      MP_DISPATCH( AddConstraint(
                     LinConLE({coefs, vars}, {0.0} )) );
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
    this->AddConstraint( LinearFunctionalConstraint(
                           itc.GetResultVar(),
    { {{const1-const2}, {args[0]}}, const2 } ) );
  }

  //////////////////// NUMBEROF CONST ///////////////////////
  void Convert(const NumberofConstConstraint& nocc) {
    const auto& args = nocc.GetArguments();
    const double k = nocc.GetParameters()[0];
    std::vector<double> coefs(args.size()+1, 1.0);
    std::vector<int> flags(args.size()+1, nocc.GetResultVar());
    for (size_t ivar = 0; ivar < args.size(); ++ivar) {
      flags[ivar] = this->AssignResultVar2Args(   // flag = (args[i]==k)
            EQ0Constraint( { {{1.0}, {args[ivar]}}, -k } ) );
    }
    coefs.back() = -1.0;
    this->AddConstraint( LinConEQ( {coefs, flags}, {0.0} ) );
  }

  //////////////////// NUMBEROF VAR ///////////////////////
  /// Very basic, could be improved
  void Convert(const NumberofVarConstraint& novc) {
    const auto& args = novc.GetArguments();
    std::vector<double> coefs(args.size(), 1.0);
    coefs.front() = -1.0;
    std::vector<int> flags(args.size(), novc.GetResultVar());
    for (size_t ivar = 1; ivar < args.size(); ++ivar) {
      flags[ivar] = this->AssignResultVar2Args(   // flag = (args[i]==args[0])
            EQ0Constraint(
                 { { {1.0, -1.0}, {args[ivar], args[0]} }, 0.0 } ) );
    }
    this->AddConstraint( LinConEQ( {coefs, flags}, {0.0} ) );
  }

  void Convert(const CountConstraint& cc) {
    const auto& args = cc.GetArguments();
    std::vector<double> coefs(args.size()+1, 1.0);
    coefs.back() = -1.0;
    std::vector<int> flags(args.size()+1, cc.GetResultVar());
    for (size_t ivar = 0; ivar < args.size(); ++ivar) {
      flags[ivar] = args[ivar];
      /// Force booleanize
      /// Either reify !=0 or constrain to 0..1? TODO param? TODO warning?
      if (!MPD( is_binary_var(args[ivar]) )) {
        auto feq0 = this->AssignResultVar2Args(   // feq0 = (args[i]==0)
            EQ0Constraint( { {{1.0}, {args[ivar]}}, 0.0 } ) );
        flags[ivar] = this->AssignResultVar2Args(   // flag = (args[i]!=0)
            NotConstraint( {feq0} ));
      }
    }
    this->AddConstraint( LinConEQ( {coefs, flags}, 0.0 ) );
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

  /// AllDiff
  INSTALL_ITEM_CONVERTER(AllDiffConverter_MIP)


  ///////////////////////////////////////////////////////////////////////
  /////////////////////// OPTIONS /////////////////////////
  ///
public:
  void InitOptions() {
    BaseConverter::InitOptions();
  }

};

} // namespace mp

#endif // MP2MIP_H
