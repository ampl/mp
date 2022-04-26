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
#include "mp/arrayref.h"
#include "mp/flat/constraint_base.h"
#include "mp/flat/expr_quadratic.h"


namespace mp {

////////////////////////////////////////////////////////////////////////
/// Generic algebraic constraint
/// @param Body: linear / higher-order terms
/// @param RhsOrRange: rhs or range
template <class Body, class RhsOrRange>
class AlgebraicConstraint :
    public BasicConstraint, public Body, public RhsOrRange {
public:
  static const std::string& GetTypeName() {
    static std::string name {
      "AlgebraicConstraint:" +
      Body::GetTypeName() + "::" +
      RhsOrRange::GetTypeName() };
    return name;
  }

  /// Constructor.
  /// By default (\a fSort = true), it sorts terms.
  /// Pass \a fSort = false to skip if you complement the terms list
  /// but do sorting later.
  /// @param le: linear / linear + higher-order terms
  /// @param rr: rhs or range
  AlgebraicConstraint(Body le, RhsOrRange rr, bool fSort=true)
    : Body(std::move(le)), RhsOrRange(std::move(rr))
  { if (fSort) sort_terms(); }

  /// Body: linear or linear + higher-order terms
  const Body& GetBody() const { return (const Body&)(*this); }

  /// Synonym, For PropagateResult()
  const Body& GetArguments() const { return GetBody(); }

  /// Compute lower slack
  double ComputeLowerSlack(ArrayRef<double> x) const {
    double s=0.0;
    for (size_t i=LinTerms::coefs().size(); i--; )
      s += LinTerms::coefs()[i] * x[LinTerms::vars()[i]];
    return s - RhsOrRange::lb();
  }

  /// Sorting and merging terms, some solvers require
  void sort_terms() { Body::sort_terms(); }

  /// Testing API
  bool operator==(const AlgebraicConstraint& lc) const {
    return Body::equals(lc) && RhsOrRange::equals(lc);
  }
};


/// Algebraic constraint range (template parameter)
class AlgConRange {
public:
  /// Class name
  static std::string GetTypeName() { return "Range"; }
  /// Constructor
  AlgConRange(double l, double u) : lb_(l), ub_(u) { }
  /// range lb()
  double lb() const { return lb_; }
  /// range ub()
  double ub() const { return ub_; }
  /// operator==
  bool equals(const AlgConRange& r) const
  { return lb()==r.lb() && ub()==r.ub(); }
private:
  double lb_, ub_;
};


/// Algebraic constraint right-hand side (template parameter).
/// Kind: -1/0/1 for <= / == / >=
template <int kind_>
class AlgConRhs {
  static constexpr const char* kind_str_[] =
  { "LE", "EQ", "GE" };
public:
  /// name
  static std::string GetTypeName()
  { return std::string(kind_str_[kind_+1]) + "Rhs"; }
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
  /// operator==
  bool equals(const AlgConRhs& r) const
  { return rhs()==r.rhs(); }
private:
  double rhs_;
};

/// Range linear constraint
using LinConRange = AlgebraicConstraint<LinTerms, AlgConRange>;
/// Convenience typedef
template <int sens>
using LinConRhs = AlgebraicConstraint< LinTerms, AlgConRhs<sens> >;
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
/// Standard quadratic constraint
/// TODO make range/rhs versions
class QuadConRange : public LinConRange {
  QuadTerms qt_;
public:
  static std::string GetTypeName() { return "QuadraticConstraint"; }

  /// Construct from a linear constraint and QP terms.
  /// Always sort terms.
  QuadConRange(LinConRange&& lc, QuadTerms&& qt) :
    LinConRange(std::move(lc)), qt_(std::move(qt)) {
    sort_qp_terms();         // LinearConstr sorts them itself
  }

  /// Constructor for testing.
  /// Always sort terms.
  QuadConRange(std::initializer_list<std::pair<double, int>> lin_terms,
                      std::initializer_list<std::tuple<double, int, int>> quad_terms,
                      double lb, double ub) :
    LinConRange({}, {lb, ub}), qt_(quad_terms) {
    LinConRange::add_terms(lin_terms);
    sort_qp_terms();
  }

  const QuadTerms& GetQPTerms() const { return qt_; }

  // To enable comparison. Also eliminates zeros
  void sort_terms() {
    LinConRange::sort_terms();
    sort_qp_terms();
  }

  void sort_qp_terms() {
    qt_.sort_terms();
  }

  /// Testing API
  bool operator==(const QuadConRange& qc) const {
    return LinConRange::operator==(qc) && qt_==qc.qt_;
  }
};


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

////////////////////////////////////////////////////////////////////////
/// AMPL represents PWL by a list of slopes
/// and breakpoints between them, assuming (X0,Y0) is on the line
class PLSlopes {
  const std::vector<double> breakpoints_, slopes_;
  const double X0_, Y0_;                        // some point on the PWL
public:
  template <class Vec>
  PLSlopes(Vec&& bp, Vec&& sl, double x, double y) noexcept :
    breakpoints_(std::forward<Vec>(bp)), slopes_(std::forward<Vec>(sl)),
    X0_(x), Y0_(y) { assert(check()); }
  const std::vector<double>& GetBP() const { return breakpoints_; }
  const std::vector<double>& GetSlopes() const { return slopes_; }
  double GetX0() const { return X0_; }
  double GetY0() const { return Y0_; }
  int GetNBP() const { return GetBP().size(); }
  int GetNSlopes() const { return GetSlopes().size(); }
  bool check() const { return GetNBP()>0 && GetNSlopes()==GetNBP()+1; }
};

/// Representing a PWL by points
struct PLPoints {
  std::vector<double> x_, y_;
  PLPoints(const PLSlopes& pls);
};

DEF_NUMERIC_FUNC_CONSTR_WITH_PRM( PLConstraint,
                  VarArray1, PLSlopes, "r = piecewise_linear(x)");


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
/// <RangeCon> complements a variable.
/// @param RangeCon: a linear or quadratic range constraint
template <class RangeCon>
class ComplementarityConstraint : public BasicConstraint {
public:
  /// The algebraic constraint
  using ConType = RangeCon;

  /// Name
  static const std::string& GetTypeName() {
    static std::string name
      { "ComplementarityConstraint[" + RangeCon::GetTypeName() + ']' };
    return name;
  }

  /// Constructor
  ComplementarityConstraint(ConType con, int var) :
    compl_con_(std::move(con)), compl_var_(var) { }

  /// Get constraint
  const ConType& GetConstraint() const { return compl_con_; }

  /// Get variable
  int GetVariable() const { return compl_var_; }

private:
  ConType compl_con_;
  int compl_var_;
};

/// Typedef ComplementarityLinRange
using ComplementarityLinRange = ComplementarityConstraint<LinConRange>;

/// Typedef ComplementarityQuadRange
using ComplementarityQuadRange = ComplementarityConstraint<QuadConRange>;

} // namespace mp

#endif // CONSTRAINTS_STATIC_H
