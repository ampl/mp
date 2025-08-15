#ifndef CONSTR_PREPRO_H
#define CONSTR_PREPRO_H

/*
 * Preprocess flat constraints before adding.
 *
 * Possible tasks:
 * 1. Simplify constraints.
 * 2. Replace a functional constraint by a different one,
 *    via returning its result variable from another
 *    (see, e.g., conditional inequalities, pow).
 * 3. Narrow domain of the result variable,
 *    so that reformulations can be tight
 *    (variable bounds, integrality.)
 *
 * Possible design for future:
 * A. Standardize fail procedure for infeasibility
 *    (which should be optional,
 *     by default we might leave it to solver?)
 * B. Make prepro code reentrable for presolve.
 * C. Consider existing result bounds.
 *    Example:
 *
 *      max(x, y) <= 5;
 *
 * After functional constraints, this file currently contains
 * preprocessing for static constraints.
 */

#include <cmath>
#include <algorithm>

#include "mp/flat/preprocess.h"
#include "mp/flat/constr_std.h"

namespace mp {

/// A mix-in base class
/// providing preprocessors of flat constraints.
/// Currently used before adding a constraint
/// (if not simplified to nothing).
template <class Impl>
class ConstraintPreprocessors {
public:

  /// Preprocess LFC
  template <class PreprocessInfo>
  void PreprocessConstraint(
      LinearFunctionalConstraint& c, PreprocessInfo& prepro) {
    auto pre = MPD( ComputeBoundsAndType(c.GetAffineExpr()) );
    prepro.narrow_result_bounds( pre.lb(), pre.ub() );
    prepro.set_result_type( pre.type() );
  }

  /// Preprocess QFC
  template <class PreprocessInfo>
  void PreprocessConstraint(
      QuadraticFunctionalConstraint& c, PreprocessInfo& prepro) {
    auto pre = MPD( ComputeBoundsAndType(c.GetQuadExpr()) );
    prepro.narrow_result_bounds( pre.lb(), pre.ub() );
    prepro.set_result_type( pre.type() );
  }

  /// Preprocess Pow
  template <class PreprocessInfo>
  void PreprocessConstraint(
      PowConstraint& c, PreprocessInfo& prepro) {
    auto x = c.GetArguments()[0];
    auto lbx = MPCD( lb(x) );
    if (lbx >= 0.0) {
      auto y = c.GetArguments()[1];
      auto ubx = MPCD( ub(x) );
      auto lby = MPCD( lb(y) );
      auto uby = MPCD( ub(y) );
      if (lbx >= 1e-10 || lby >= 1e-10) {
        std::array<double, 4> rb
            {std::pow(lbx, lby), std::pow(lbx, uby),
                                 std::pow(ubx, lby), std::pow(ubx, uby)};
        prepro.narrow_result_bounds(
            *std::min_element(rb.begin(), rb.end()),
            *std::max_element(rb.begin(), rb.end()) );
      } else {
        std::array<double, 4> rb
            {std::pow(1e-10, lby), std::pow(1e-10, uby),
             std::pow(ubx, lby), std::pow(ubx, uby)};
        prepro.narrow_result_bounds(
            *std::min_element(rb.begin(), rb.end()),
            *std::max_element(rb.begin(), rb.end()) );
      }
    }
  }

  /// Preprocess PowConstExp
  template <class PreprocessInfo>
  void PreprocessConstraint(
      PowConstExpConstraint& c, PreprocessInfo& prepro) {
    auto pwr = c.GetParameters()[0];
    if (0.0==std::fabs(pwr)) {              // decidable case
      prepro.narrow_result_bounds(1.0, 1.0);
      return;
    }
    auto arg = c.GetArguments()[0];
    if (1.0==pwr) {                         // decidable case
      prepro.set_result_var(arg);
      return;
    }
    auto& m = MP_DISPATCH( GetModel() );
    bool lbx_neg = m.lb(arg)<0.0;
    bool ubx_pos = m.ub(arg)>0.0;
    bool pow_int = MPD( is_integer_value(pwr) );
    if ((!pow_int && lbx_neg)          // a is fractional, lbx<0
        || (pwr<0 && lbx_neg)) {       // a<0, lbx<=0
      // We _COULD_ PL approximate when pwr<0, pow_int, lbx<0.
      // But we leave it here (don't even narrow result).
      // Gurobi 10 does not handle a<0 && lbx<0.
    } else {
      auto lbr = std::pow(m.lb(arg), pwr),
          ubr = std::pow(m.ub(arg), pwr);
      if (pow_int && pwr>=0.0) {
        // result integer if x integer, a>=0
        prepro.set_result_type( m.var_type(arg) );
      }
      if (MPD( is_integer_value(pwr / 2.0) )) {  // a even
        if (lbx_neg && ubx_pos) {
          ubr = std::max(lbr, ubr);
          lbr = 0.0;
        }
      }  // else, monotone
      prepro.narrow_result_bounds( std::min(lbr, ubr),
                                   std::max(lbr, ubr) );
    }
  }

  /// Preprocess Min
  template <class PreprocessInfo>
  void PreprocessConstraint(
      MinConstraint& c, PreprocessInfo& prepro) {
    auto& m = MP_DISPATCH( GetModel() );
    auto& args = c.GetArguments();
    if (MPCD( IfPreproUnnest() & 8 ))
      IntegrateNested(c);
    if (MPCD( IfPreproSortUnify() ))
      SortAndUnify(args);
    prepro.narrow_result_bounds( m.lb_array(args),
                          m.ub_min_array(args) );
    prepro.set_result_type( m.common_type(args) );
  }

  /// Preprocess Max
  template <class PreprocessInfo>
  void PreprocessConstraint(
      MaxConstraint& c, PreprocessInfo& prepro) {
    auto& m = MP_DISPATCH( GetModel() );
    auto& args = c.GetArguments();
    if (MPCD( IfPreproUnnest() & 8 ))
      IntegrateNested(c);
    if (MPCD( IfPreproSortUnify() ))
      SortAndUnify(args);
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
    const auto lb = MPD( lb(argvar) ),
        ub = MPD( ub(argvar) );
    if (lb>=0.0) {  // When result var is set, \a c is skipped
      prepro.set_result_var(argvar);
      return;
    } else if (ub<=0.0) {
      auto res = MPD( AssignResult2Args(   // create newvar = -argvar
            LinearFunctionalConstraint({ {{-1.0}, {argvar}}, 0.0 })) );
      prepro.set_result_var(res.get_var());
      return;
    }
    prepro.narrow_result_bounds(0.0, std::max(-lb, ub));
    prepro.set_result_type( MPD( var_type(argvar) ) );
  }

  /// Preprocess CondLinConEQ
  template <class PreprocessInfo>
  void PreprocessConstraint(
      CondLinConEQ& c, PreprocessInfo& prepro) {
    if (CheckEmptySubCon(c, prepro))
      return;
    prepro.narrow_result_bounds(0.0, 1.0);
    prepro.set_result_type( var::INTEGER );
    if (!IsNormalized(c))
      c.GetConstraint().negate();   // for equality
    if (0!=MPD( IfPreproEqResBounds() ))
      if (FixEqualityResult(c, prepro))
        return;
    PreprocessEqVarConst__unifyCoef(c);
    if (0!=MPD( IfPreproEqBinVar() ))
      if (ReuseEqualityBinaryVar(c, prepro))
        return;
  }

  /// See if the argument of a conditional
  /// algebraic constraint is normalized
  template <class Body, int kind>
  bool IsNormalized(
      ConditionalConstraint<
        AlgebraicConstraint< Body, AlgConRhs<kind> > >& cc) {
    auto& arg = cc.GetConstraint();
    return arg.is_normalized();
  }

  /// Preprocess CondQuadConEQ
  template <class PreprocessInfo>
  void PreprocessConstraint(
      CondQuadConEQ& c, PreprocessInfo& prepro) {
    if (CheckEmptySubCon(c, prepro))
      return;
    prepro.narrow_result_bounds(0.0, 1.0);
    prepro.set_result_type( var::INTEGER );
    if (0!=MPD( IfPreproEqResBounds() ))
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
    auto bndsNType = MPD( ComputeBoundsAndType(body) );
    if (bndsNType.lb() > rhs || bndsNType.ub() < rhs) {
      prepro.narrow_result_bounds(0.0, 0.0);
      return true;
    }
    if (bndsNType.lb()==rhs && bndsNType.ub()==rhs) {
      prepro.narrow_result_bounds(1.0, 1.0);
      return true;
    }
    if (var::INTEGER==bndsNType.type_ &&
        !is_integer(con.rhs())) {
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
          prepro.set_result_var( var );
        else if (0.0==std::fabs(rhs))
          prepro.set_result_var( MPD( MakeComplementVar(var) ) );
        else
          prepro.narrow_result_bounds(0.0, 0.0);    // not 0/1 value, result false
        return true;
      }
    }
    return false;
  }

  /// (Non)strict inequalities
  template <class PreprocessInfo, class Body, int kind>
  void PreprocessConstraint(
      ConditionalConstraint<
        AlgebraicConstraint< Body, AlgConRhs<kind> > >& cc,
      PreprocessInfo& prepro) {
    if (CheckEmptySubCon(cc, prepro))
      return;
    prepro.narrow_result_bounds(0.0, 1.0);
    prepro.set_result_type( var::INTEGER );
    assert(kind);
    auto& algc = cc.GetArguments();
    if (!IsNormalized(cc)) {
      auto arg1 = algc; // Add the normalized one instead
      arg1.negate();    // Negate the terms and sense
      prepro.set_result_var(
            MPD( AssignResultVar2Args(
                   ConditionalConstraint<
                     AlgebraicConstraint< Body, AlgConRhs<
                   -kind> > > { {
                   std::move(arg1.GetBody()), arg1.rhs()
                 } } ) ));
      return;
    }
    // See if we need to round the constant term
    auto rhs = algc.rhs();
    auto bnt = MPD(
          ComputeBoundsAndType(algc.GetBody()) );
    if (MPCD( IfPreproIneqRHS() )
        && var::INTEGER == bnt.get_result_type()
        && std::floor(rhs) != std::ceil(rhs)) {  // rhs not int
      if (1==kind)  // algc is >=
        algc.set_rhs( std::ceil(rhs) );
      else if (-1==kind)
        algc.set_rhs( std::floor(rhs) );
      else if (2==kind)  // algc is >
        algc.set_rhs( std::floor(rhs) ); // > floor(rhs)
      else {
        assert(-2==kind);
        algc.set_rhs( std::ceil(rhs) );
      }
    }
    // Decidable cases
    if (MPCD( IfPreproIneqResBounds() )) {
      if ((1==kind && bnt.lb()>=rhs)
          || (-1==kind && bnt.ub()<=rhs)
          || (2==kind && bnt.lb()>rhs)
          || (-2==kind && bnt.ub()<rhs)) {
        prepro.narrow_result_bounds(1.0, 1.0);   // TRUE
      } else {
        if ((1==kind && bnt.ub()<rhs)
            || (-1==kind && bnt.lb()>rhs)
            || (2==kind && bnt.ub()<=rhs)
            || (-2==kind && bnt.lb()>=rhs)) {
          prepro.narrow_result_bounds(0.0, 0.0);   // FALSE
        }
      }
    }
  }

  /// Dave experiments with logic presolve.
  ///
  /// @return true if constraint is empty
  /// and we should return the fixed result.
  template <class SubCon, class PreprocessInfo>
  bool CheckEmptySubCon(
      const ConditionalConstraint<SubCon>& cc,
      PreprocessInfo& prepro) {
    if (cc.GetConstraint().empty()) {
#ifndef NDEBUG
      MPD(AddWarning("empty_cmp",
                     "Empty comparison in a logical constraint\n  of type '"
                     + std::string(cc.GetTypeName())
                     + "'.\n  Please contact AMPL support."));
#endif
      auto res = ComputeValue(cc, std::vector<double>{});
      prepro.narrow_result_bounds(res, res);
      return true;
    }
    return false;
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      AndConstraint& con, PreprocessInfo& prepro) {
    prepro.narrow_result_bounds(0.0, 1.0);
    prepro.set_result_type( var::INTEGER );
    // Remove fixed variables for XPRESS (solvers/#61).
    auto n01 = count_fixed_01(con.GetArguments());
    if (n01.first) {               // AND = false
      prepro.narrow_result_bounds(0.0, 0.0);
      return;
    }
    if ((int)con.GetArguments().size() == n01.second) {
      prepro.narrow_result_bounds(1.0, 1.0);
      return;
    }
    if (n01.second) {
      std::vector<int> arg1;
      arg1.reserve(con.GetArguments().size() - n01.second);
      for (auto x: con.GetArguments()) {
        if (MPCD( lb(x) ) <= 0.0)    // not fixed
          arg1.push_back(x);
      }
      con.GetArguments() = arg1;
    }
    if (MPCD( IfPreproUnnest() & 1 ))
      IntegrateNested(con);            // flatten nested
    if (MPCD( IfPreproSortUnify() ))
      SortAndUnify(con.GetArguments());
    if (con.GetArguments().empty())
      prepro.narrow_result_bounds(1.0, 1.0);  // empty conjunction
    if (1 == con.GetArguments().size())
      prepro.set_result_var(con.GetArguments()[0]);
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      OrConstraint& con, PreprocessInfo& prepro) {
    prepro.narrow_result_bounds(0.0, 1.0);
    prepro.set_result_type( var::INTEGER );
    // Remove fixed variables for XPRESS (solvers/#61).
    auto n01 = count_fixed_01(con.GetArguments());
    if (n01.second) {               // OR = true
      prepro.narrow_result_bounds(1.0, 1.0);
      return;
    }
    if ((int)con.GetArguments().size() == n01.first) {
      prepro.narrow_result_bounds(0.0, 0.0);
      return;
    }
    if (n01.first) {
      std::vector<int> arg1;
      arg1.reserve(con.GetArguments().size() - n01.first);
      for (auto x: con.GetArguments()) {
        if (MPCD( ub(x) ) >= 1.0)    // not fixed
          arg1.push_back(x);
      }
      con.GetArguments() = arg1;
    }
    if (MPCD( IfPreproUnnest() & 1 ))
      IntegrateNested(con);            // flatten nested
    if (MPCD( IfPreproSortUnify() ))
      SortAndUnify(con.GetArguments());
    if (con.GetArguments().empty())
      prepro.narrow_result_bounds(0.0, 0.0);  // empty disjunction
    if (1 == con.GetArguments().size())
      prepro.set_result_var(con.GetArguments()[0]);
  }

  /// Count N fixed binary vars
  template <class Vec>
  std::pair<int, int> count_fixed_01(const Vec& vec) const {
    std::pair<int, int> result {0, 0};
    for (auto x: vec) {
      assert(MPCD( is_binary_var(x) ));
      if (MPCD( ub(x) ) <= 0.0)
        ++ result.first;
      if (MPCD( lb(x) ) >= 1.0)
        ++ result.second;
    }
    assert(result.first + result.second <= (int)vec.size());
    return result;
  }

  /// Integrate nested constraints (And/Or)
  template <class Con>
  void IntegrateNested(Con& con) {
    bool fChanges = false;
    VarArray args_new;
    args_new.reserve(con.GetArguments().size());
    for (auto v: con.GetArguments()) {
      if (auto pNested = MPD( template
          GetInitExpressionOfType<Con>(v) )) {
        fChanges = true;
        const auto& args2 = pNested->GetArguments();
        args_new.insert(args_new.end(),
                        args2.begin(), args2.end());
        // Not any more #201 #266: now usage++ by FixAsTrue()
        // MPD( DecrementVarUsage(v) );
      } else {
        args_new.push_back(v);
      }
    }
    if (fChanges)
      con.GetArguments() = args_new;
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      EquivalenceConstraint& , PreprocessInfo& prepro) {
    prepro.narrow_result_bounds(0.0, 1.0);
    prepro.set_result_type( var::INTEGER );
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      AllDiffConstraint& c, PreprocessInfo& prepro) {
    if (MPCD( IfPreproSortUnify() ))
      Sort(c.GetArguments());        // @todo warn if duplicates
    prepro.narrow_result_bounds(0.0, 1.0);
    prepro.set_result_type( var::INTEGER );
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      NumberofConstConstraint& con, PreprocessInfo& prepro) {
    if (MPCD( IfPreproSortUnify() ))
      Sort(con.GetArguments());
    prepro.narrow_result_bounds(0.0, (double)con.GetArguments().size());
    prepro.set_result_type( var::INTEGER );
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      NumberofVarConstraint& con, PreprocessInfo& prepro) {
    if (MPCD( IfPreproSortUnify() ))
      Sort(con.GetArguments());
    prepro.narrow_result_bounds(0.0,     // size()-1: 1st arg is the ref var
                                (double)con.GetArguments().size()-1);
    prepro.set_result_type( var::INTEGER );
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      CountConstraint& con, PreprocessInfo& prepro) {
    if (MPCD( IfPreproSortUnify() ))
      Sort(con.GetArguments());
    prepro.narrow_result_bounds(0.0, (double)con.GetArguments().size());
    prepro.set_result_type( var::INTEGER );
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      NotConstraint& nc, PreprocessInfo& prepro) {
    auto nfixed = count_fixed_01(nc.GetArguments());
    if (nfixed.first) {
      prepro.narrow_result_bounds(1.0, 1.0);
    } else if (nfixed.second) {
      prepro.narrow_result_bounds(0.0, 0.0);
    } else {
      prepro.narrow_result_bounds(0.0, 1.0);
      prepro.set_result_type( var::INTEGER );
    }
  }

  /// Preprocess Div
  template <class PreprocessInfo>
  void PreprocessConstraint(
      DivConstraint& c, PreprocessInfo& prepro) {
    auto& m = MPD( GetModel() );
    auto v1 = c.GetArguments()[0], v2 = c.GetArguments()[1];
    const auto l1=m.lb(v1), u1=m.ub(v1), l2=m.lb(v2), u2=m.ub(v2);
		if (l1 > MPD( PracticallyMinusInf() ) &&
				u1 < MPD( PracticallyInf() ) &&
				l2 > MPD( PracticallyMinusInf() ) &&
				u2 < MPD( PracticallyInf() ) &&
        l2 * u2 > 0.0) {
      auto l0 = std::numeric_limits<double>::max();
      auto u0 = std::numeric_limits<double>::lowest();
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

  /// Preprocess IfThen
  template <class PreprocessInfo>
  void PreprocessConstraint(
      IfThenConstraint& c, PreprocessInfo& prepro) {
    const auto& args = c.GetArguments();
    prepro.narrow_result_bounds(
          std::min(MPD( lb(args[1]) ), MPD( lb(args[2]) )),
        std::max(MPD( ub(args[1]) ), MPD( ub(args[2]) )));
    prepro.set_result_type( MP_DISPATCH(GetModel()).
                            common_type( { args[1], args[2] } ) );
  }

  /// Preprocess IfThen
  template <class PreprocessInfo>
  void PreprocessConstraint(
      ImplicationConstraint& , PreprocessInfo& prepro) {
    prepro.narrow_result_bounds(0.0, 1.0); // no prepro yet
    prepro.set_result_type( var::INTEGER );
  }

  ////////////////////// NONLINEAR FUNCTIONS //////////////////////
  template <class PreprocessInfo>
  void PreprocessConstraint(
      ExpConstraint& , PreprocessInfo& prepro) {
    prepro.narrow_result_bounds(0.0, MPD( Infty() ));
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      ExpAConstraint& , PreprocessInfo& prepro) {
    prepro.narrow_result_bounds(0.0, MPD( Infty() ));
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      LogConstraint& c, PreprocessInfo& ) {
    if (MPCD( IfBoundLogArg() )) {
      auto x = c.GetArguments()[0];  // if no positive lb,
      MPD( NarrowVarBounds(x, 0.0, MPD( Infty() )) );
    }
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      LogAConstraint& c, PreprocessInfo& ) {
    if (MPCD( IfBoundLogArg() ))
      MPD( NarrowVarBounds(
          c.GetArguments()[0], 0.0, MPD( Infty() )) );
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

  /// Preprocess Tan
  template <class PreprocessInfo>
  void PreprocessConstraint(
      TanConstraint& , PreprocessInfo& ) {
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      AsinConstraint& , PreprocessInfo& prepro) {
    prepro.narrow_result_bounds(-MPD(Pi())/2, MPD(Pi())/2);
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      AcosConstraint& , PreprocessInfo& prepro) {
    prepro.narrow_result_bounds(0.0, MPD(Pi()));
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      AtanConstraint& , PreprocessInfo& prepro) {
    prepro.narrow_result_bounds(-MPD(Pi())/2, MPD(Pi())/2);
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      SinhConstraint& , PreprocessInfo& ) {
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      CoshConstraint& , PreprocessInfo& prepro) {
    prepro.narrow_result_bounds(1.0, MPD(Infty()));
  }

  /// Preprocess Tan
  template <class PreprocessInfo>
  void PreprocessConstraint(
      TanhConstraint& , PreprocessInfo& prepro) {
    prepro.narrow_result_bounds(-1.0, 1.0);
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      AsinhConstraint& , PreprocessInfo& ) {
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      AcoshConstraint& , PreprocessInfo& prepro) {
    prepro.narrow_result_bounds(0.0, MPD(Infty()));
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      AtanhConstraint& , PreprocessInfo& ) {
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      PLConstraint& , PreprocessInfo& ) {
  }


  /// Static constraints.
  /// @return true iff the constraints has been presolved
  ///   into something different, no need to keep it.
  /// @note could accept non-const ref and modify the argument.
  template <class Con>
  bool PreprocessStaticConstraint(const Con& )
  { return false; }

  /// Preprocess Indicator
  /// @note Necessary for XPRESS 9.4.2
  template <class SubCon>
  bool PreprocessStaticConstraint(
      const IndicatorConstraint<SubCon>& indc) {
    if (MPCD( is_fixed(indc.get_binary_var()) )) {
      int fixed_var_val = MPCD( fixed_value(indc.get_binary_var()) );
      if (fixed_var_val == indc.get_binary_value()) {
        MPD( AddConstraint(indc.get_constraint()) );
      } else {
        // forget
      }
      return true;
    }
    return false;
  }

};

} // namespace mp

#endif // CONSTR_PREPRO_H
