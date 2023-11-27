#ifndef CONSTRAINTS_GENERAL_H
#define CONSTRAINTS_GENERAL_H

/**
  * Static general constraints
  */

#include <vector>
#include <string>
#include <map>

#include "mp/error.h"
#include "mp/flat/constr_algebraic.h"


namespace mp {

////////////////////////////////////////////////////////////////////////
/// Indicator: b==bv -> [constraint]
template <class Con>
class IndicatorConstraint: public BasicConstraint {
public:
  /// Constraint type name
  static const std::string& GetTypeName() {
    static std::string name
      { "IndicatorConstraint[" + Con::GetTypeName() + ']' };
    return name;
  }

  /// Is logical?
  static bool IsLogical() { return true; }

  /// Constructor
  IndicatorConstraint(int b, int bv, Con con) noexcept :
    b_(b), bv_(bv), con_(std::move(con)) { assert(check()); }
  bool check() const { return (b_>=0) && (bv_==0 || bv_==1); }

  /// Getters
  int get_binary_var() const { return b_; }
  int get_binary_value() const { return bv_; }
  bool is_binary_value_1() const
  { return 1==get_binary_value(); }
  const Con& get_constraint() const { return con_; }

  /// Compute violation
  template <class VarInfo>
  Violation ComputeViolation(const VarInfo& x) const {
    assert(b_<(int)x.size());
    auto xb = x[b_];
    if (std::round(xb) == bv_)      // Implication needed
      return con_.ComputeViolation(x);
    return {0.0, 0.0};
  }


private:
  const int b_=-1;                  // the indicator variable
  const int bv_=1;                  // the value, 0/1
  const Con con_;
};


/// Typedef indicator<LinConLE>
using IndicatorConstraintLinLE = IndicatorConstraint<LinConLE>;

/// Typedef indicator<LinConEQ>
using IndicatorConstraintLinEQ = IndicatorConstraint<LinConEQ>;

/// Typedef indicator<LinConGE>
using IndicatorConstraintLinGE = IndicatorConstraint<LinConGE>;

/// Typedef indicator<QuadConLE>
using IndicatorConstraintQuadLE = IndicatorConstraint<QuadConLE>;

/// Typedef indicator<QuadConEQ>
using IndicatorConstraintQuadEQ = IndicatorConstraint<QuadConEQ>;

/// Typedef indicator<QuadConGE>
using IndicatorConstraintQuadGE = IndicatorConstraint<QuadConGE>;


////////////////////////////////////////////////////////////////////////
/// SOS1, SOS2

/// SOS constraint extra info, supplied for better conversion
struct SOSExtraInfo {
  /// Bounds on the sum of variables
  struct Bounds {
    Bounds(double lb=-1e100, double ub=1e100) :
      lb_(lb), ub_(ub) { }
    double lb_, ub_;
    bool operator==(const Bounds& b) const
    { return lb_==b.lb_ && ub_==b.ub_; }
  } bounds_;
  /// Constructor
  SOSExtraInfo(Bounds b={}) : bounds_(b) { }
  /// operator==
  bool operator==(const SOSExtraInfo& ei) const
  { return bounds_==ei.bounds_; }
};

/// SOS1, SOS2 constraints
template <int type>
class SOS_1or2_Constraint: public BasicConstraint {
  static constexpr const char* name1_ = "SOS1Constraint";
  static constexpr const char* name2_ = "SOS2Constraint";

  std::vector<int> v_;
  std::vector<double> w_;
  const SOSExtraInfo extra_info_;
public:
  /// Constraint type name
  static const char* GetTypeName()
  { return 1==type ? name1_ : name2_; }

  /// Is logical?
  static bool IsLogical() { return true; }

  int get_sos_type() const { return type; }
  int size() const { return (int)v_.size(); }
  /// Returns vector of variables, sorted by weights
  const std::vector<int>& get_vars() const { return v_; }
  /// Returns weights, sorted
  const std::vector<double>& get_weights() const { return w_; }
  /// SOS2 extra info
  const SOSExtraInfo& get_extra_info() const { return extra_info_; }
  /// Sum of vars range, extra info if supplied
  SOSExtraInfo::Bounds get_sum_of_vars_range() const
  { return extra_info_.bounds_; }

  /// Constructor
  template <class VV=std::vector<int>, class WV=std::vector<double> >
  SOS_1or2_Constraint(VV v, WV w,
                      SOSExtraInfo ei, std::string nm={}) :
    v_(v), w_(w), extra_info_(ei) {
    SetName(std::move(nm));
    sort();
    assert(check());
  }
  /// Constructor
  template <class VV=std::vector<int>, class WV=std::vector<double> >
  SOS_1or2_Constraint(VV v, WV w, std::string nm={}) :
    v_(v), w_(w) {
    SetName(std::move(nm));
    sort();
    assert(check());
  }
  /// Check data
  bool check() const { return type>=1 && type<=2 &&
                     v_.size()==w_.size(); }

  /// Compute violation
  template <class VarInfo>
  Violation ComputeViolation(const VarInfo& x) const {
    if (1==type)
      return ComputeViolationSOS1(x);
    return ComputeViolationSOS2(x);
  }

  /// Compute violation
  template <class VarInfo>
  Violation ComputeViolationSOS1(const VarInfo& x) const {
    int nnz=0;
    for (int i=(int)v_.size(); i--; ) {
      if (x.is_nonzero(v_[i]))
        ++nnz;
    }
    return {(double)std::max(0, nnz-1), 0.0};
  }

  /// Compute violation
  template <class VarInfo>
  Violation ComputeViolationSOS2(const VarInfo& x) const {
    int npos=0, pos1=-1, pos2=-1;
    int posDist=1;     // should be 1
    for (int i=(int)v_.size(); i--; ) {
      if (x.is_positive(v_[i])) {
        ++npos;
        if (pos1<0)
          pos1=i;
        else if (pos2<0) {
          pos2=i;
          posDist = pos1-pos2;
        }
      }
    }
    return {double(std::max(0, npos-2)
        + std::abs(1 - posDist)), 0.0};
  }


protected:
  /// Sort by weights
  void sort() {
    std::map<double, int> sorted;
    for (auto i=v_.size(); i--; )
      MP_ASSERT_ALWAYS(
            sorted.insert({ w_.at(i), v_.at(i) }).second,
            "SOS2: weights not unique");
    v_.clear();
    w_.clear();
    for (const auto& wv: sorted) {
      v_.push_back(wv.second);
      w_.push_back(wv.first);
    }
  }
};


/// Typedef SOS1Constraint
using SOS1Constraint = SOS_1or2_Constraint<1>;

/// Typedef SOS2Constraint
using SOS2Constraint = SOS_1or2_Constraint<2>;


////////////////////////////////////////////////////////////////////////
/// Complementarity constraint.
/// \a Expr complements a variable.
/// @param Expr: an affine or quadratic functional expression
template <class Expr>
class ComplementarityConstraint : public BasicConstraint {
public:
  /// The expression type
  using ExprType = Expr;

  /// Constraint type name
  static const std::string& GetTypeName() {
    static std::string name
      { std::string("ComplementarityConstraint[") +
          Expr::GetTypeName() + ']' };
    return name;
  }

  /// Is logical?
  static bool IsLogical() { return true; }

  /// Constructor
  ComplementarityConstraint(ExprType expr, int var) :
    compl_expr_(std::move(expr)), compl_var_(var) { }

  /// Get constraint
  const ExprType& GetExpression() const { return compl_expr_; }

  /// Get variable
  int GetVariable() const { return compl_var_; }

  /// Compute violation
  template <class VarInfo>
  Violation ComputeViolation(const VarInfo& x) const {
    auto ve = compl_expr_.ComputeValue(x);
    if (x.is_at_lb(compl_var_))
      return {-ve, 0.0};
    else if (x.is_at_ub(compl_var_))
      return {ve, 0.0};
    return {std::fabs(ve), 0.0};
  }


private:
  ExprType compl_expr_;
  int compl_var_;
};


/// Typedef ComplementarityLinRange
using ComplementarityLinear = ComplementarityConstraint<AffineExpr>;

/// Typedef ComplementarityQuadRange
using ComplementarityQuadratic = ComplementarityConstraint<QuadraticExpr>;

/// Quadratic cone
DEF_STATIC_CONSTR_WITH_PRM( QuadraticConeConstraint, VarArray, DblParamArray,
                            "Quadratic cone p1*x1 >= sqrt((p2*x2)^2 + ...)),"
                            " with factors p1..pn");
/// Rotated quadratic cone
DEF_STATIC_CONSTR_WITH_PRM( RotatedQuadraticConeConstraint, VarArray, DblParamArray,
                            "Rotated quadratic cone "
                            "2 * p1*x1*p2*x2 >= (p3*x3)^2 + ...),"
                            " x1, x2 >= 0, with factors p1..pn");
/// Exponential cone
DEF_STATIC_CONSTR_WITH_PRM( ExponentialConeConstraint, VarArray3, DblParamArray3,
                            "Exponential cone ax >= by exp(cz / (by)),"
                            " where ax, by >= 0, with factors a,b,c");
/// Power cone
DEF_STATIC_CONSTR_WITH_PRM( PowerConeConstraint, VarArray, DblParamArray,
                            "Power cone with factors ");
/// Geometric cone
DEF_STATIC_CONSTR_WITH_PRM( GeometricConeConstraint, VarArray, DblParamArray,
                            "Geometric with factors ");

} // namespace mp

#endif // CONSTRAINTS_GENERAL_H
