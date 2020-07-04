#ifndef STD_CONSTR_H
#define STD_CONSTR_H

#include <vector>
#include <algorithm>

#include "mp/convert/basic_constr.h"
#include "mp/convert/affine_expr.h"

namespace mp {

////////////////////////////////////////////////////////////////////////
/// Standard linear constraint
class LinearConstraint : public BasicConstraint {
  std::vector<double> coefs_;
  std::vector<int> vars_;
  double lb_, ub_;
public:
  static const char* GetConstraintName() { return "LinearConstraint"; }
  template <class CV=std::vector<double>, class VV=std::vector<int> >
  LinearConstraint(CV&& c, VV&& v, double l, double u)
    : coefs_(std::forward<CV>(c)), vars_(std::forward<VV>(v)),
      lb_(l), ub_(u) { assert(coefs_.size()==vars_.size()); preprocess(); }
  template <size_t N>
  LinearConstraint(std::array<double, N>& c, std::array<int, N>& v, double l, double u)
    : coefs_(c.begin(), c.end()), vars_(v.begin(), v.end()),
      lb_(l), ub_(u) { assert(coefs_.size()==vars_.size()); preprocess(); }
  int nnz() const { return (int)coefs_.size(); }
  const double* coefs() const { return coefs_.data(); }
  const int* vars() const { return vars_.data(); }
  double lb() const { return lb_; }
  double ub() const { return ub_; }

  void preprocess() { eliminate_zeros(); }          // TODO check duplicates?
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
};

////////////////////////////////////////////////////////////////////////
/// Converting linear expr to 2 vectors.
struct LinearExprUnzipper {
  std::vector<double> c_;
  std::vector<int> v_;
  LinearExprUnzipper() { }
  LinearExprUnzipper(const LinearExpr& e) {
    Reserve(e.num_terms());
    for (LinearExpr::const_iterator it=e.begin(); it!=e.end(); ++it) {
      AddTerm(it->var_index(), it->coef());
    }
  }
  int num_terms() const { return c_.size(); }
  const std::vector<double>& coefs() const { return c_; }
  const std::vector<int>& var_indexes() const { return v_; }
  void Reserve(size_t s) { c_.reserve(s); v_.reserve(s); }
  void AddTerm(int v, double c) { c_.push_back(c); v_.push_back(v); }
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
  LinearDefiningConstraint(int r, AffineExpr&& ae) :
    DefiningConstraint(r), affine_expr_(std::move(ae)) {
    /// TODO sort elements
  }
  const AffineExpr& GetAffineExpr() const { return affine_expr_; }
  LinearConstraint to_linear_constraint() const {
    const auto& ae = GetAffineExpr();
    LinearExprUnzipper aeu(ae);
    aeu.AddTerm(DefiningConstraint::GetResultVar(), -1.0);
    return LinearConstraint(std::move(aeu.c_), std::move(aeu.v_),
                            -ae.constant_term(), -ae.constant_term());
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
DEFINE_CUSTOM_DEFINING_CONSTRAINT( EQ0Constraint, AffineExpr,
                                   "r = (expr == 0)");

////////////////////////////////////////////////////////////////////////
DEFINE_CUSTOM_DEFINING_CONSTRAINT( NEConstraint__unused, VarArray2,
                                   "r = (v1 != v2)");

////////////////////////////////////////////////////////////////////////
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
/// Indicator: b==bv -> c'x <= rhs
class IndicatorConstraintLinLE: public BasicConstraint {
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
  /// Produces affine expr ae so that the inequality is equivalent to ae<=0.0
  AffineExpr to_lhs_affine_expr() const {
    return {get_lin_coefs(), get_lin_vars(), -get_lin_rhs()};
  }
  /// Constructor
  template <class CV=std::vector<double>, class VV=std::vector<int> >
  IndicatorConstraintLinLE(int b, int bv,
                           CV&& c, VV&& v,
                           double rhs) :
    b_(b), bv_(bv), c_(std::forward<CV>(c)), v_(std::forward<VV>(v)), rhs_(rhs)
  { assert(check()); }
  bool check() const { return (b_>=0) && (bv_==0 || bv_==1); }
};

} // namespace mp

#endif // STD_CONSTR_H
