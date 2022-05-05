#ifndef CONSTRAINTS_GENERAL_H
#define CONSTRAINTS_GENERAL_H

/**
  * Static general constraints
  */

#include <vector>
#include <string>

#include "mp/flat/constraints_algebraic.h"


namespace mp {

////////////////////////////////////////////////////////////////////////
/// Indicator: b==bv -> [constraint]
template <class Con>
class IndicatorConstraint: public BasicConstraint {
public:
  static const std::string& GetTypeName() {
    static std::string name
      { "IndicatorConstraint[" + Con::GetTypeName() + ']' };
    return name;
  }

  /// Constructor
  IndicatorConstraint(int b, int bv, Con con) noexcept :
    b_(b), bv_(bv), con_(std::move(con)) { assert(check()); }
  bool check() const { return (b_>=0) && (bv_==0 || bv_==1); }

  /// Getters
  int get_binary_var() const { return b_; }
  int get_binary_value() const { return bv_; }
  bool is_binary_value_1() const { return 1==get_binary_value(); }
  const Con& get_constraint() const { return con_; }


private:
  const int b_=-1;                            // the indicator variable
  const int bv_=1;                            // the value, 0/1
  const Con con_;
};


/// Typedef indicator<LinConLE>
using IndicatorConstraintLinLE = IndicatorConstraint<LinConLE>;

/// Typedef indicator<LinConEQ>
using IndicatorConstraintLinEQ = IndicatorConstraint<LinConEQ>;

/// Typedef indicator<QuadConLE>
using IndicatorConstraintQuadLE = IndicatorConstraint<QuadConLE>;

/// Typedef indicator<QuadConEQ>
using IndicatorConstraintQuadEQ = IndicatorConstraint<QuadConEQ>;


////////////////////////////////////////////////////////////////////////
/// SOS1, SOS2
template <int type>
class SOS_1or2_Constraint: public BasicConstraint {
  static constexpr const char* name1_ = "SOS1Constraint";
  static constexpr const char* name2_ = "SOS2Constraint";

  const std::vector<int> v_;
  const std::vector<double> w_;
public:
  static const char* GetTypeName()
  { return 1==type ? name1_ : name2_; }

  int get_sos_type() const { return type; }
  int size() const { return (int)v_.size(); }
  const std::vector<int>& get_vars() const { return v_; }
  const std::vector<double>& get_weights() const { return w_; }

  /// Constructor
  template <class VV=std::vector<int>, class WV=std::vector<double> >
  SOS_1or2_Constraint(VV&& v, WV&& w) noexcept :
    v_(std::forward<VV>(v)), w_(std::forward<WV>(w))
  { assert(check()); }
  bool check() const { return type>=1 && type<=2 &&
                     v_.size()==w_.size(); }
};


/// Typedef SOS1Constraint
using SOS1Constraint = SOS_1or2_Constraint<1>;

/// Typedef SOS2Constraint
using SOS2Constraint = SOS_1or2_Constraint<2>;


////////////////////////////////////////////////////////////////////////
/// Complementarity constraint.
/// <Expr> complements a variable.
/// @param Expr: an affine or quadratic functional expression
template <class Expr>
class ComplementarityConstraint : public BasicConstraint {
public:
  /// The expression type
  using ExprType = Expr;

  /// Name
  static const std::string& GetTypeName() {
    static std::string name
      { std::string("ComplementarityConstraint[") +
          Expr::GetTypeName() + ']' };
    return name;
  }

  /// Constructor
  ComplementarityConstraint(ExprType expr, int var) :
    compl_expr_(std::move(expr)), compl_var_(var) { }

  /// Get constraint
  const ExprType& GetExpression() const { return compl_expr_; }

  /// Get variable
  int GetVariable() const { return compl_var_; }


private:
  ExprType compl_expr_;
  int compl_var_;
};


/// Typedef ComplementarityLinRange
using ComplementarityLinear = ComplementarityConstraint<AffineExpr>;

/// Typedef ComplementarityQuadRange
using ComplementarityQuadratic = ComplementarityConstraint<QuadraticExpr>;


} // namespace mp

#endif // CONSTRAINTS_GENERAL_H
