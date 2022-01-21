#ifndef STD_CONSTR_H
#define STD_CONSTR_H

#include <vector>
#include <algorithm>
#include <map>
//#include <iostream>

#include "mp/flat/basic_constr.h"
#include "mp/flat/quad_expr.h"

namespace mp {

////////////////////////////////////////////////////////////////////////
/// Standard linear constraint

class LinearConstraint : public BasicConstraint {
  std::vector<double> coefs_;
  std::vector<int> vars_;
public:
  static const std::string& name() {
    static std::string name { "LinCon" };
    return name;
  }
  template <class CV=std::vector<double>, class VV=std::vector<int> >
  LinearConstraint(CV&& c, VV&& v) noexcept : coefs_(std::forward<CV>(c)), vars_(std::forward<VV>(v))
  { assert(coefs_.size()==vars_.size()); preprocess(); }
  template <size_t N>
  LinearConstraint(std::array<double, N>& c, std::array<int, N>& v) :
     coefs_(c.begin(), c.end()), vars_(v.begin(), v.end())
  { assert(coefs_.size()==vars_.size()); preprocess(); }
  LinearConstraint(std::initializer_list<std::pair<double, int>> lin_exp)  {
    coefs_.reserve(lin_exp.size());
    vars_.reserve(lin_exp.size());
    for (const auto& term: lin_exp) {
      coefs_.push_back(term.first);
      vars_.push_back(term.second);
    }
    preprocess();
  }
  int nnz() const { return (int)coefs_.size(); }
  const double* pcoefs() const { return coefs_.data(); }
  const int* pvars() const { return vars_.data(); }
  const std::vector<double>& coefs() const { return coefs_; }
  const std::vector<int>& vars() const { return vars_; }
  /// For PropagateResult()
  const std::vector<int>& GetArguments() const { return vars_; }

  void AddTerm(double c, int v) { coefs_.push_back(c); vars_.push_back(v); }

  void preprocess() { sort_terms(); /*eliminate_zeros();*/ }
  /// So Gurobi does not complain about 0 coefs and duplicate variables
  void eliminate_zeros() {
    auto ic1 = std::find_if(coefs_.begin(), coefs_.end(),
                            [](double c){return 0.0==std::fabs(c);});
    if (coefs_.end()!=ic1) {
      auto iv1 = vars_.begin() + (ic1-coefs_.begin());
      auto iv = iv1;
      auto ic = ic1;
      for( ; ++iv, ++ic != coefs_.end(); )
        if (0.0!=std::fabs(*ic)) {
          *ic1++ = *ic;
          *iv1++ = *iv;
        }
      coefs_.resize(ic1-coefs_.begin());
      vars_.resize(ic1-coefs_.begin());
      assert(coefs_.size()==vars_.size());
    }
  }
  /// Add same variables, eliminate 0's
  void sort_terms() {
    std::map<int, double> var_coef_map;
    for (int i=0; i<nnz(); ++i)
      if (0.0!=std::fabs(coefs_[i]))
        var_coef_map[vars_[i]] += coefs_[i];
    if (true) {                                  // would check size if using hash map
      coefs_.clear();
      vars_.clear();
      for (const auto& vc: var_coef_map) {
        if (0.0!=std::fabs(vc.second)) {         // Need tolerance?
          coefs_.push_back(vc.second);
          vars_.push_back(vc.first);
        }
      }
    }
  }
  
  /// Testing API
  bool operator==(const LinearConstraint& lc) const {
      return this->equal_to(lc);
  }
protected:
    virtual bool equal_to(LinearConstraint const& lc) const {
        return coefs_ == lc.coefs_ && vars_ == lc.vars_;
   }
   //TODO
//  void print(std::ostream& os) const {
//    os << lb_ << " <= ";
//    for (int i=0; i<nnz(); ++i) {
//      os << coefs_[i] << "*[" << vars_[i] << ']';
//      if (i<nnz()-1)
//        os << " + ";
//    }
//    os << " <= " << ub_;
//  }
};


class RangeLinCon : public LinearConstraint {
public:
    template <class CV = std::vector<double>, class VV = std::vector<int> >
    RangeLinCon(CV&& c, VV&& v, std::initializer_list<double> r) noexcept
        : LinearConstraint(c, v) {
        // Of course here i wouldn't use std::initializer_list but simply
        // RangeLinCon(CV&& c, VV&& v, double lb, double ub) noexcept
        // This was just to have the code compiling without refactoring all the calls to these constructors
        auto it = r.begin();
        lb_ = *it++;
        ub_ = *it;
    }
    template <size_t N>
    RangeLinCon(std::array<double, N>& c, std::array<int, N>& v, std::initializer_list<double> r)
        : LinearConstraint(c, v) {
        auto it = r.begin();
        lb_ = *it++;
        ub_ = *it;
    }
    RangeLinCon(std::initializer_list<std::pair<double, int>> lin_exp,
        std::initializer_list<double> r) : LinearConstraint(lin_exp) {
        auto it = r.begin();
        lb_ = *it++;
        ub_ = *it;
    }

    /// name
    static std::string name() { return  LinearConstraint::name() + "Range"; }
    /// Constructor
    /// range lb()
    double lb() const { return lb_; }
    /// range ub()
    double ub() const { return ub_; }
    /// operator==
    bool equals(const RangeLinCon& r) const
    {
        return equal_to(r) && lb() == r.lb() && ub() == r.ub();
    }
private:
    double lb_ = std::numeric_limits<double>::min(), ub_ = std::numeric_limits<double>::max();
};
/// Kind: -1/0/1 for <= / == / >=
template <int kind_> class RHSLinCon : public LinearConstraint {
    static constexpr const char* kind_str_[] =
    { "LE", "EQ", "GE" };
public:
    template <class CV = std::vector<double>, class VV = std::vector<int> >
    RHSLinCon(CV&& c, VV&& v, double rhs) noexcept
        : LinearConstraint(c, v), rhs_(rhs) { }
    template <size_t N>
    RHSLinCon(std::array<double, N>& c, std::array<int, N>& v, double rhs)
        : LinearConstraint(c, v), rhs_(rhs) {}
    RHSLinCon(std::initializer_list<std::pair<double, int>> lin_exp,
        double rhs) : LinearConstraint(lin_exp), rhs_(rhs) {}
    /// name
    static std::string name()
    {
        return LinearConstraint::name() + std::string(kind_str_[kind_ + 1]) + "Rhs";
    }
    /// Kind
    int kind() const { return kind_; }
    /// rhs()
    double rhs() const { return rhs_; }
    /// operator==
    bool equals(const RHSLinCon& r) const
    {
        return equal_to(r) && rhs() == r.rhs();
    }
private:
    double rhs_;
};

using LinConLE = RHSLinCon< -1>;
using LinConEQ = RHSLinCon< 0>;
using LinConGE = RHSLinCon<  1>;

////////////////////////////////////////////////////////////////////////
/// Standard quadratic constraint
class QuadraticConstraint : public RangeLinCon {
  QuadTerms qt_;
public:
  static const char* GetConstraintName() { return "QuadraticConstraint"; }
  QuadraticConstraint(RangeLinCon&& lc, QuadTerms&& qt) noexcept :
    RangeLinCon(std::move(lc)), qt_(std::move(qt)) {
    sort_qp_terms(); // LinearConstr sorts them itself
  }
  QuadraticConstraint(std::initializer_list<std::pair<double, int>> lin_exp,
                      std::initializer_list<std::tuple<double, int, int>> quad_terms,
                      double lb, double ub) :
      RangeLinCon(lin_exp, { lb, ub }), qt_(quad_terms) {
    sort_qp_terms();
  }

  const QuadTerms& GetQPTerms() const { return qt_; }

  // To enable comparison. Also eliminates zeros
  void sort_terms() {
    RangeLinCon::sort_terms();
    sort_qp_terms();
  }

  void sort_qp_terms() {
    qt_.SortTerms();
  }

  /// Testing API
  bool operator==(const QuadraticConstraint& qc) const {
    return RangeLinCon::operator==(qc) && qt_==qc.qt_;
  }
//  void print(std::ostream& os) const {
//    os << lb() << " <= ";
//    for (int i=0; i<nnz(); ++i) {
//      os << coefs()[i] << "*[" << vars()[i] << ']';
//      if (i<nnz()-1)
//        os << " + ";
//    }
//    for (int i=0; i<qt_.num_terms(); ++i) {
//      os << " + "
//         << qt_.coef(i) << "*[" << qt_.var1(i) << "]*[" << qt_.var2(i) << "]";
//    }
//    os << " <= " << ub();
//  }
};

////////////////////////////////////////////////////////////////////////
/// Converting linear expr to 2 vectors.
struct LinearExprUnzipper {
  std::vector<double> c_;
  std::vector<int> v_;
  LinearExprUnzipper() { }
  LinearExprUnzipper(const LinearExpr& e) { AddTerms(e); }
  std::size_t num_terms() const { return c_.size(); }
  const std::vector<double>& coefs() const { return c_; }
  const std::vector<int>& var_indexes() const { return v_; }
  void Reserve(std::size_t s) { c_.reserve(s); v_.reserve(s); }
  void AddTerm(int v, double c) { c_.push_back(c); v_.push_back(v); }
  void AddTerms(const LinearExpr& e) {
    Reserve(num_terms() + e.num_terms());
    for (LinearExpr::const_iterator it=e.begin(); it!=e.end(); ++it) {
      AddTerm(it->var_index(), it->coef());
    }
  }
};

////////////////////////////////////////////////////////////////////////
/// Linear Defining Constraint: r = affine_expr
class LinearDefiningConstraint :
    public DefiningConstraint {
  AffineExpr affine_expr_;
public:
  static const char* GetConstraintName() { return "LinearDefiningConstraint"; }
  using Arguments = AffineExpr;
  using DefiningConstraint::GetResultVar;
  /// A constructor ignoring result variable: use AssignResultToArguments() then
  LinearDefiningConstraint(AffineExpr&& ae) noexcept :
    affine_expr_(std::move(ae)) {
    /// TODO sort+merge elements
  }
  LinearDefiningConstraint(int r, AffineExpr&& ae) noexcept :
    DefiningConstraint(r), affine_expr_(std::move(ae)) {
    /// TODO sort+merge elements
  }
  const AffineExpr& GetAffineExpr() const { return affine_expr_; }
  const Arguments& GetArguments() const { return GetAffineExpr(); }
  RangeLinCon to_linear_constraint() const {
      const auto& ae = GetAffineExpr();
      LinearExprUnzipper aeu(ae);
      aeu.AddTerm(DefiningConstraint::GetResultVar(), -1.0);
      return RangeLinCon(std::move(aeu.c_), std::move(aeu.v_),
          { -ae.constant_term(), -ae.constant_term() });
  }
};

////////////////////////////////////////////////////////////////////////
/// Quadratic Defining Constraint: r = quad_expr
class QuadraticDefiningConstraint :
    public DefiningConstraint {
  QuadExpr quad_expr_;
public:
  static const char* GetConstraintName() { return "QuadraticDefiningConstraint"; }
  using Arguments = QuadExpr;
  using DefiningConstraint::GetResultVar;
  /// A constructor ignoring result variable: use AssignResultToArguments() then
  QuadraticDefiningConstraint(QuadExpr&& qe) noexcept :
    quad_expr_(std::move(qe)) {
    /// TODO sort+merge elements
  }
  QuadraticDefiningConstraint(int r, QuadExpr&& qe) noexcept :
    DefiningConstraint(r), quad_expr_(std::move(qe)) {
    /// TODO sort+merge elements
  }
  const QuadExpr& GetQuadExpr() const { return quad_expr_; }
  const Arguments& GetArguments() const { return GetQuadExpr(); }
  QuadraticConstraint to_quadratic_constraint() const {
    const auto& qe = GetQuadExpr();
    const auto& ae = qe.GetAE();
    LinearExprUnzipper aeu(ae);
    aeu.AddTerm(DefiningConstraint::GetResultVar(), -1.0);
    auto qt = qe.GetQT();
    return {RangeLinCon(std::move(aeu.c_), std::move(aeu.v_),
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
/// \brief DEFINE_CUSTOM_DEFINING_CONSTRAINT
/// Keep it with AffineExpr, indicators need that
/// and we don't want quadratics with big-M's?
/// Or, add QuadraticEq0Constraint?
DEFINE_CUSTOM_DEFINING_CONSTRAINT( EQ0Constraint, AffineExpr,
                                   "r = (expr == 0)");

/// Not using: var1 != var2.
/// Represented by Not { Eq0Constraint... }
////////////////////////////////////////////////////////////////////////
DEFINE_CUSTOM_DEFINING_CONSTRAINT( NEConstraint__unused, VarArray2,
                                   "r = (v1 != v2)");

////////////////////////////////////////////////////////////////////////
/// \brief DEFINE_CUSTOM_DEFINING_CONSTRAINT
////// Keep it with AffineExpr, indicators need that
/// and we don't want quadratics with big-M's?
DEFINE_CUSTOM_DEFINING_CONSTRAINT( LE0Constraint, AffineExpr,
                                   "r = (expr <= 0)");

////////////////////////////////////////////////////////////////////////
DEFINE_CUSTOM_DEFINING_CONSTRAINT( NotConstraint, VarArray1,
                                   "r = !v");

////////////////////////////////////////////////////////////////////////
DEFINE_CUSTOM_DEFINING_CONSTRAINT( IfThenConstraint, VarArrayN<3>,
                                  "if (cond) then (expr1) else (expr2)");

////////////////////////////////////////////////////////////////////////
DEFINE_CUSTOM_DEFINING_CONSTRAINT( AllDiffConstraint, VarArray,
                                  "alldiff({})");



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



////////////////////////////////////////////////////////////////////////
/// Indicator: b==bv -> c'x ? rhs
template <int sens>
class IndicatorConstraintLin: public BasicConstraint {
public:
  static const char* GetConstraintName() { return "IndicatorConstraint"; }
  const int b_=-1;                            // the indicator variable
  const int bv_=1;                            // the value, 0/1
  const std::vector<double> c_;
  const std::vector<int> v_;
  const double rhs_;
  /// Getters
  int get_binary_var() const { return b_; }
  int get_binary_value() const { return bv_; }
  bool is_binary_value_1() const { return 1==get_binary_value(); }
  const std::vector<double>& get_lin_coefs() const { return c_; }
  const std::vector<int>& get_lin_vars() const { return v_; }
  double get_lin_rhs() const { return rhs_; }
  constexpr int get_sense() const { return sens; }
  /// Produces affine expr ae so that the inequality is equivalent to ae<=0.0
  AffineExpr to_lhs_affine_expr() const {
    return {get_lin_coefs(), get_lin_vars(), -get_lin_rhs()};
  }
  /// Constructor
  template <class CV=std::vector<double>, class VV=std::vector<int> >
  IndicatorConstraintLin(int b, int bv,
                           CV&& c, VV&& v,
                           double rhs) noexcept :
    b_(b), bv_(bv), c_(std::forward<CV>(c)), v_(std::forward<VV>(v)), rhs_(rhs)
  { assert(check()); }
  bool check() const { return (b_>=0) && (bv_==0 || bv_==1); }
};

using IndicatorConstraintLinLE = IndicatorConstraintLin<-1>;
using IndicatorConstraintLinEQ = IndicatorConstraintLin<0>;

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
