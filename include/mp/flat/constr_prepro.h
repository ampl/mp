#ifndef CONSTR_PREPRO_H
#define CONSTR_PREPRO_H

/**
 * Preprocess flat constraints before adding.
 *
 * Possible tasks:
 * 1. Simplify constraints
 * 2. Replace a functional constraint by a different one,
 *    via returning its result variable from another
 *    (see conditional inequalities).
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
    auto lb = std::pow(m.lb(arg), pwr),
        ub = std::pow(m.ub(arg), pwr);
		if (MPD( is_integer_value(pwr) )) {
			if (pwr>=0.0) {
				// result integer if arg is and integer, >=0 exponent
				prepro.set_result_type( m.var_type(arg) );
				if (MPD( is_integer_value(pwr / 2.0) )) {  // exponent is even, >=0
					bool lb_neg = m.lb(arg)<0.0;
					bool ub_pos = m.ub(arg)>0.0;
					if (lb_neg && ub_pos) {
						ub = std::max(lb, ub); lb = 0.0;
					} else if (lb_neg) {
						std::swap(lb, ub);
					}
				}
			}
		} else {                      // fractional power
			if (lb<0.0)
				lb = 0.0;
		}
    prepro.narrow_result_bounds( std::min(lb, ub),
                          std::max(lb, ub) );
  }

  /// Preprocess Min
  template <class PreprocessInfo>
  void PreprocessConstraint(
      MinConstraint& c, PreprocessInfo& prepro) {
    auto& m = MP_DISPATCH( GetModel() );
    auto& args = c.GetArguments();
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
    prepro.narrow_result_bounds(0.0, 1.0);
    prepro.set_result_type( var::INTEGER );
    assert(kind);
    auto& algc = cc.GetArguments();
    if (!IsNormalized(cc)) {
      auto arg1 = algc;
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
    auto bnt_body = MPD(
          ComputeBoundsAndType(algc.GetBody()) );
    if (var::INTEGER == bnt_body.get_result_type()
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
    if (con.GetArguments().size() == n01.second) {
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
    if (con.GetArguments().size() == n01.first) {
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
    return result;
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
    auto x = c.GetArguments()[0];  // if no positive lb,
    if (MPD(lb(x))<=0.0)           // set (fixed) lb=1e-6
      MPD( NarrowVarBounds(x, 1e-6, MPD( Infty() )) );
  }

  template <class PreprocessInfo>
  void PreprocessConstraint(
      LogAConstraint& c, PreprocessInfo& ) {
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
    prepro.narrow_result_bounds(-MPD(Pi())/2, MPD(Pi()));
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


};

} // namespace mp

#endif // CONSTR_PREPRO_H
