#ifndef STD_CONSTR_H
#define STD_CONSTR_H

#include <vector>
#include <algorithm>
#include <map>
#include <string>
#include <limits>

#include "mp/arrayref.h"
#include "mp/flat/basic_constr.h"
#include "mp/flat/quad_expr.h"

namespace mp {

////////////////////////////////////////////////////////////////////////
/// Generic linear constraint
template <class LinTerms, class RhsOrRange>
class LinearConstraint :
    public BasicConstraint, public LinTerms, public RhsOrRange {
public:
  static const std::string& GetConstraintName() {
    static std::string name { "LinCon" + RhsOrRange::name() };
    return name;
  }
  LinearConstraint(LinTerms le, RhsOrRange rr) noexcept
    : LinTerms(std::move(le)), RhsOrRange(std::move(rr))
  { /* preprocess(); */ }

  /// For PropagateResult()
  const std::vector<int>& GetArguments() const
  { return LinTerms::vars(); }

  /// Compute lower slack
  double ComputeLowerSlack(ArrayRef<double> x) const {
    double s=0.0;
    for (size_t i=LinTerms::coefs().size(); i--; )
      s += LinTerms::coefs()[i] * x[LinTerms::vars()[i]];
    return s - RhsOrRange::lb();
  }

  /// Testing API
  bool operator==(const LinearConstraint& lc) const {
    return LinTerms::coefs()==lc.coefs() && LinTerms::vars()==lc.vars() &&
        RhsOrRange::equals(lc);
  }
};

class AlgConRange {
public:
  /// name
  static std::string name() { return "Range"; }
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

/// Kind: -1/0/1 for <= / == / >=
template <int kind_>
class AlgConRhs {
  static constexpr const char* kind_str_[] =
  { "LE", "EQ", "GE" };
public:
  /// name
  static std::string name()
  { return std::string(kind_str_[kind_+1]) + "Rhs"; }
  /// Constructor
  AlgConRhs(double r) : rhs_(r) { }
  /// Kind
  int kind() const { return kind_; }
  /// rhs()
  double rhs() const { return rhs_; }
  /// lb(): this is a specialization of the range constraint
  double lb() const {
    return kind_<0 ? -std::numeric_limits<double>::infinity() : rhs();
  }
  /// ub(): this is a specialization of the range constraint
  double ub() const {
    return kind_>0 ? std::numeric_limits<double>::infinity() : rhs();
  }
  /// operator==
  bool equals(const AlgConRhs& r) const
  { return rhs()==r.rhs(); }
private:
  double rhs_;
};

/// Range linear constraint
using RangeLinCon = LinearConstraint<LinTerms, AlgConRange>;
/// Convenience typedef
template <int sens>
using LinConRhs = LinearConstraint< LinTerms, AlgConRhs<sens> >;
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
class QuadraticConstraint : public RangeLinCon {
  QuadTerms qt_;
public:
  static const char* GetConstraintName() { return "QuadraticConstraint"; }
  /// Construct from a linear constraint and QP terms
  QuadraticConstraint(RangeLinCon&& lc, QuadTerms&& qt) noexcept :
    RangeLinCon(std::move(lc)), qt_(std::move(qt)) {
    sort_qp_terms(); // LinearConstr sorts them itself
  }
  /// Constructor for testing
  QuadraticConstraint(std::initializer_list<std::pair<double, int>> lin_exp,
                      std::initializer_list<std::tuple<double, int, int>> quad_terms,
                      double lb, double ub) :
    RangeLinCon({}, {lb, ub}), qt_(quad_terms) {
    RangeLinCon::add_terms(lin_exp);
    sort_qp_terms();
  }

  const QuadTerms& GetQPTerms() const { return qt_; }

  // To enable comparison. Also eliminates zeros
  void sort_terms() {
    RangeLinCon::sort_terms();
    sort_qp_terms();
  }

  void sort_qp_terms() {
    qt_.sort_terms();
  }

  /// Testing API
  bool operator==(const QuadraticConstraint& qc) const {
    return RangeLinCon::operator==(qc) && qt_==qc.qt_;
  }
};


////////////////////////////////////////////////////////////////////////
/// Linear Defining Constraint: r = affine_expr
class LinearDefiningConstraint :
    public DefiningConstraint {
  AffExp affine_expr_;
public:
  static const char* GetConstraintName() { return "LinearDefiningConstraint"; }
  using Arguments = AffExp;
  using DefiningConstraint::GetResultVar;
  /// A constructor ignoring result variable: use AssignResultToArguments() then
  LinearDefiningConstraint(AffExp&& ae) noexcept :
    affine_expr_(std::move(ae)) {  // TODO sort+merge elements?
  }
  LinearDefiningConstraint(int r, AffExp&& ae) noexcept :
    DefiningConstraint(r), affine_expr_(std::move(ae)) {
    /// TODO sort+merge elements
  }
  const AffExp& GetAffineExpr() const { return affine_expr_; }
  const Arguments& GetArguments() const { return GetAffineExpr(); }
  LinConEQ to_linear_constraint() const {
    const auto& ae = GetAffineExpr();
    auto le = ae.get_lin_exp();
    le.add_term(-1.0, DefiningConstraint::GetResultVar());
    LinConEQ lc { le, -ae.constant_term() };
    return lc;
  }
};

////////////////////////////////////////////////////////////////////////
/// Quadratic Defining Constraint: r = quad_expr
class QuadraticDefiningConstraint :
    public DefiningConstraint {
  QuadExp quad_expr_;
public:
  static const char* GetConstraintName() { return "QuadraticDefiningConstraint"; }
  using Arguments = QuadExp;
  using DefiningConstraint::GetResultVar;
  /// A constructor ignoring result variable: use AssignResultToArguments() then
  QuadraticDefiningConstraint(QuadExp&& qe) noexcept :
    quad_expr_(std::move(qe)) {
    /// TODO sort+merge elements
  }
  QuadraticDefiningConstraint(int r, QuadExp&& qe) noexcept :
    DefiningConstraint(r), quad_expr_(std::move(qe)) {
    /// TODO sort+merge elements
  }
  const QuadExp& GetQuadExpr() const { return quad_expr_; }
  const Arguments& GetArguments() const { return GetQuadExpr(); }
  QuadraticConstraint to_quadratic_constraint() const {
    const auto& qe = GetQuadExpr();
    const auto& ae = qe.GetAE();
    auto le = ae.get_lin_exp();
    le.add_term(-1.0, DefiningConstraint::GetResultVar());
    auto qt = qe.GetQT();
    return {RangeLinCon(std::move(le),
                             {-ae.constant_term(), -ae.constant_term()}),
            std::move(qt)};
  }
};

////////////////////////////////////////////////////////////////////////
DEFINE_CUSTOM_DEFINING_CONSTRAINT( MaximumConstraint, VarArray,
                                   "r = max(v1, v2, ..., vn)");

////////////////////////////////////////////////////////////////////////
DEFINE_CUSTOM_DEFINING_CONSTRAINT( MinimumConstraint, VarArray,
                                   "r = min(v1, v2, ..., vn)");

////////////////////////////////////////////////////////////////////////
DEFINE_CUSTOM_DEFINING_CONSTRAINT( AbsConstraint, VarArray1,
                                   "r = abs(v)");

////////////////////////////////////////////////////////////////////////
DEFINE_CUSTOM_DEFINING_CONSTRAINT( ConjunctionConstraint, VarArray,
                                   "r = forall({vi})");

////////////////////////////////////////////////////////////////////////
DEFINE_CUSTOM_DEFINING_CONSTRAINT( DisjunctionConstraint, VarArray,
                                   "r = exists({vi})");

////////////////////////////////////////////////////////////////////////
/// Storing AffExp instead of LinConEQ because big-M is straightforwardly
/// computed for (aff_exp) <= 0:
/// b -> ae<=0 is linearized as ae <= ub(ae)*(1-b).
/// If we stored LinConEQ:
/// b -> lin_exp<=d would be linearized as
/// le <= d + (ub(le)-d)*(1-b)  <==>
/// le <= d + ub_le - d - ub_le*b + d*b  <==>
/// le <= ub_le + (d-ub_le)*b.  Not too complex.
/// Keep it with AffineExpr, indicators need that
/// and we don't want quadratics with big-M's?
/// Or, add QuadraticEq0Constraint?
/// TODO Use LinConEq / QuadConEQ ?
/// TODO Have the actual conditional constraint as a template parameter?
/// TODO Diff to Indicator?
DEFINE_CUSTOM_DEFINING_CONSTRAINT( EQ0Constraint, AffExp,
                                   "r = (expr == 0)");

/// Extract underlying constraint.
/// Can be done more general if using something like
/// ConditionalConstraint<> instead if EQ0C / LE0C
inline LinConEQ ExtractConstraint(const EQ0Constraint& eq0c) {
  const auto& ae=eq0c.GetArguments();
  return { (LinTerms)ae, -ae.constant_term() };
}

/// Not using: var1 != var2.
/// Represented by Not { Eq0Constraint... }
////////////////////////////////////////////////////////////////////////
DEFINE_CUSTOM_DEFINING_CONSTRAINT( NEConstraint__unused, VarArray2,
                                   "r = (v1 != v2)");

////////////////////////////////////////////////////////////////////////
////// Keep it with AffineExpr, indicators need that
/// and we don't want quadratics with big-M's?
DEFINE_CUSTOM_DEFINING_CONSTRAINT( LE0Constraint, AffExp,
                                   "r = (expr <= 0)");

/// Extract underlying constraint.
/// Can be done more general if using something like
/// ConditionalConstraint<> instead if EQ0C / LE0C
inline LinConLE ExtractConstraint(const LE0Constraint& le0c) {
  const auto& ae=le0c.GetArguments();
  return { (LinTerms)ae, -ae.constant_term() };
}

////////////////////////////////////////////////////////////////////////
DEFINE_CUSTOM_DEFINING_CONSTRAINT( NotConstraint, VarArray1,
                                   "r = !v");

////////////////////////////////////////////////////////////////////////
DEFINE_CUSTOM_DEFINING_CONSTRAINT( IfThenConstraint, VarArrayN<3>,
                                  "if (cond) then (expr1) else (expr2)");

////////////////////////////////////////////////////////////////////////
DEFINE_CUSTOM_DEFINING_CONSTRAINT( AllDiffConstraint, VarArray,
                                  "alldiff({})");

////////////////////////////////////////////////////////////////////////
DEFINE_CUSTOM_DEFINING_CONSTRAINT_WITH_PARAMS( NumberofConstConstraint,
                                  VarArray, DblParamArray1,
                                  "numberof_const(k (=x0), {x1...xn})");

////////////////////////////////////////////////////////////////////////
DEFINE_CUSTOM_DEFINING_CONSTRAINT( NumberofVarConstraint, VarArray,
                                  "numberof_var(x0, {x1...xn})");

////////////////////////////////////////////////////////////////////////
DEFINE_CUSTOM_DEFINING_CONSTRAINT( CountConstraint, VarArray,
                                   "count({x0...xn})");


//////////////////// NONLINEAR FUNCTIONS //////////////////////
////////////////////////////////////////////////////////////////////////
DEFINE_CUSTOM_DEFINING_CONSTRAINT( ExpConstraint, VarArray1,
                                   "r = exp(v)");

////////////////////////////////////////////////////////////////////////
DEFINE_CUSTOM_DEFINING_CONSTRAINT_WITH_PARAMS( ExpAConstraint,
                  VarArray1, DblParamArray1, "r = a**v");

////////////////////////////////////////////////////////////////////////
DEFINE_CUSTOM_DEFINING_CONSTRAINT( LogConstraint, VarArray1,
                                   "r = log(v)");

////////////////////////////////////////////////////////////////////////
DEFINE_CUSTOM_DEFINING_CONSTRAINT_WITH_PARAMS( LogAConstraint,
                  VarArray1, DblParamArray1, "r = log(v)/log(a)");

////////////////////////////////////////////////////////////////////////
DEFINE_CUSTOM_DEFINING_CONSTRAINT_WITH_PARAMS( PowConstraint,
                  VarArray1, DblParamArray1, "r = v ** a");

////////////////////////////////////////////////////////////////////////
DEFINE_CUSTOM_DEFINING_CONSTRAINT( SinConstraint, VarArray1,
                                   "r = sin(v)");

////////////////////////////////////////////////////////////////////////
DEFINE_CUSTOM_DEFINING_CONSTRAINT( CosConstraint, VarArray1,
                                   "r = cos(v)");

////////////////////////////////////////////////////////////////////////
DEFINE_CUSTOM_DEFINING_CONSTRAINT( TanConstraint, VarArray1,
                                   "r = tan(v)");


/// Where applicable, produces expr
/// so that the constraint is equivalent to expr<=>0.0
template <class Constraint>
int ToLhsExpr(const Constraint& ) {
  throw std::runtime_error(
        std::string("Cannot produce a lhs expr from ") +
          typeid (Constraint).name());
  return 0;
}

////////////////////////////////////////////////////////////////////////
/// Indicator: b==bv -> [constraint]
template <class Con>
class IndicatorConstraint: public BasicConstraint {
public:
  static const std::string& GetConstraintName() {
    static std::string name
      { "IndicatorConstraint[" + Con::name() + ']' };
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

using IndicatorConstraintLinLE = IndicatorConstraint<LinConLE>;
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

DEFINE_CUSTOM_DEFINING_CONSTRAINT_WITH_PARAMS( PLConstraint,
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
  static const char* GetConstraintName()
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

using SOS1Constraint = SOS_1or2_Constraint<1>;
using SOS2Constraint = SOS_1or2_Constraint<2>;

} // namespace mp

#endif // STD_CONSTR_H
