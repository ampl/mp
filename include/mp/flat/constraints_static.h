#ifndef CONSTRAINTS_STATIC_H
#define CONSTRAINTS_STATIC_H

/**
  * Static (non-functional) flat constraints
  */

#include <vector>
#include <algorithm>
#include <map>
#include <string>
#include <cmath>

#include "mp/error.h"
#include "mp/flat/constraint_base.h"
#include "mp/flat/expr_quadratic.h"


namespace mp {

////////////////////////////////////////////////////////////////////////
/// Generic algebraic constraint
/// @param Body: linear or linear and higher-order terms
/// @param RhsOrRange: rhs or range
template <class Body, class RhsOrRange>
class AlgebraicConstraint :
    public BasicConstraint, public Body, public RhsOrRange {
public:
  static const std::string& GetTypeName() {
    static std::string name {
      std::string("AlgebraicConstraint:") +
      Body::GetTypeName() + "::" +
      RhsOrRange::GetTypeName() };
    return name;
  }

  /// BodyType
  using BodyType = Body;

  /// RhsOrRangeType
  using RhsOrRangeType = RhsOrRange;

  /// Constructor.
  /// By default (\a fSort = true), it sorts terms.
  /// Pass \a fSort = false to skip if you complement the terms list
  /// but do sorting later.
  /// @param le: linear / linear + higher-order terms
  /// @param rr: rhs or range
  AlgebraicConstraint(Body le, RhsOrRange rr, bool fSort=true)
    : Body(std::move(le)), RhsOrRange(std::move(rr))
  { if (fSort) sort_terms(); }

  /// Body: linear or linear + higher-order terms, const
  const Body& GetBody() const { return (const Body&)(*this); }

  /// Body: linear or linear + higher-order terms
  Body& GetBody() { return (Body&)(*this); }

  /// Range or RHS. Used for hash<>
  RhsOrRange GetRhsOrRange() const { return (RhsOrRange)(*this); }

  /// Synonym, For PropagateResult()
  const Body& GetArguments() const { return GetBody(); }

  /// If no variable terms in the body
  bool empty() const { return Body::empty(); }

  /// Compute lower slack
  double ComputeLowerSlack(ArrayRef<double> x) const {
    return Body::ComputeValue(x) - RhsOrRange::lb();
  }

  /// Sorting and merging terms, some solvers require
  void sort_terms() { Body::sort_terms(); }

  /// Negate
  void negate() { Body::negate(); RhsOrRange::negate(); }

  /// Testing API
  bool operator==(const AlgebraicConstraint& lc) const {
    return Body::equals(lc) && RhsOrRange::equals(lc);
  }
};


/// Algebraic constraint range (template parameter)
class AlgConRange {
public:
  /// Class name
  static constexpr const char* GetTypeName() { return "Range"; }
  /// Constructor
  AlgConRange(double l, double u) : lb_(l), ub_(u) { }
  /// range lb()
  double lb() const { return lb_; }
  /// range ub()
  double ub() const { return ub_; }
  /// negate
  void negate() { auto tmp=ub_; ub_=-lb_; lb_=-tmp; }
  /// operator==
  bool equals(const AlgConRange& r) const
  { return lb()==r.lb() && ub()==r.ub(); }
private:
  double lb_, ub_;
};


/// Algebraic constraint right-hand side (template parameter).
/// Kind: -2/-1/0/1 for <= / == / >=
template <int kind_>
class AlgConRhs {
  static constexpr const char* kind_str_[] =
  { "LT", "LE", "EQ", "GE", "GT" };
public:
  /// name
  static std::string GetTypeName()
  { return std::string("Rhs") + kind_str_[kind_+2]; }
  /// Constructor
  AlgConRhs(double r) : rhs_(r) { }
  /// Kind
  int kind() const { return kind_; }
  /// rhs()
  double rhs() const { return rhs_; }
  /// lb(): this is a specialization of the range constraint
  double lb() const {
    return kind_<0 ? -INFINITY : rhs();
  }
  /// ub(): this is a specialization of the range constraint
  double ub() const {
    return kind_>0 ? INFINITY : rhs();
  }
  /// Set rhs
  void set_rhs(double v) { rhs_ = v; }
  /// negate
  void negate() { rhs_ = -rhs_; }
  /// operator==
  bool equals(const AlgConRhs& r) const
  { return rhs()==r.rhs(); }
private:
  double rhs_;
};


////////////////////////////////////////////////////////////////////////
/// Linear range constraint
using LinConRange = AlgebraicConstraint<LinTerms, AlgConRange>;

/// Convenience typedef
template <int sens>
using LinConRhs = AlgebraicConstraint< LinTerms, AlgConRhs<sens> >;
/// Linear constraint c'x <  d
using LinConLT = LinConRhs<-2>;
/// Linear constraint c'x <= d
using LinConLE = LinConRhs<-1>;
/// Linear constraint c'x == d
using LinConEQ = LinConRhs< 0>;
/// Linear constraint c'x >= d
using LinConGE = LinConRhs< 1>;

/// Where applicable, produces expr
/// so that the constraint is equivalent to expr<=>0.0
template <int sens>
AffExp ToLhsExpr(
    const LinConRhs<sens>& lc) {
  return { LinTerms(lc), -lc.rhs() };
}


////////////////////////////////////////////////////////////////////////
/// Quadratic range constraint
using QuadConRange =
    AlgebraicConstraint<QuadAndLinTerms, AlgConRange>;

/// Convenience typedef
template <int sens>
using QuadConRhs =
    AlgebraicConstraint< QuadAndLinTerms, AlgConRhs<sens> >;
/// Quadratic constraint c'x+x'Qx <  d
using QuadConLT = QuadConRhs<-2>;
/// Quadratic constraint c'x+x'Qx <= d
using QuadConLE = QuadConRhs<-1>;
/// Quadratic constraint c'x+x`Qx == d
using QuadConEQ = QuadConRhs< 0>;
/// Quadratic constraint c'x+x`Qx >= d
using QuadConGE = QuadConRhs< 1>;


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

  /// Getters
  int get_binary_var() const { return b_; }
  int get_binary_value() const { return bv_; }
  bool is_binary_value_1() const { return 1==get_binary_value(); }
  const Con& get_constraint() const { return con_; }

  /// Where applicable, produces expr
  /// so that the constraint is equivalent to expr<=>0.0
  auto to_lhs_expr() const ->decltype (ToLhsExpr(get_constraint())) {
    return ToLhsExpr(get_constraint());
  }

  /// Constructor
  IndicatorConstraint(int b, int bv, Con con) noexcept :
    b_(b), bv_(bv), con_(std::move(con)) { assert(check()); }
  bool check() const { return (b_>=0) && (bv_==0 || bv_==1); }


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
using ComplementarityLinear = ComplementarityConstraint<AffExp>;

/// Typedef ComplementarityQuadRange
using ComplementarityQuadratic = ComplementarityConstraint<QuadExp>;

} // namespace mp

#endif // CONSTRAINTS_STATIC_H
