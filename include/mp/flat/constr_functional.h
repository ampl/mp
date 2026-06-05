#ifndef CONSTRAINTS_FUNCTIONAL_H
#define CONSTRAINTS_FUNCTIONAL_H

/*
  * Functional flat constraints
  */

#include "mp/error.h"
#include "mp/flat/constr_static.h"


namespace mp {


////////////////////////////////////////////////////////////////////////
DEF_NUMERIC_FUNC_CONSTR( Max, VarArray,
                                   "r = max(v1, v2, ..., vn)");

////////////////////////////////////////////////////////////////////////
DEF_NUMERIC_FUNC_CONSTR( Min, VarArray,
                                   "r = min(v1, v2, ..., vn)");

////////////////////////////////////////////////////////////////////////
DEF_NUMERIC_FUNC_CONSTR( Abs, VarArray1,
                                   "r = abs(v)");

////////////////////////////////////////////////////////////////////////
DEF_LOGICAL_FUNC_CONSTR( And, VarArray,
                                   "r = forall({vi})");

////////////////////////////////////////////////////////////////////////
DEF_LOGICAL_FUNC_CONSTR( Or, VarArray,
                                   "r = exists({vi})");

////////////////////////////////////////////////////////////////////////
/// NL equivalence constraint: r <==> (expr(var1) <==> expr(var2)).
/// We have to do this for NL output because flat constraints
/// represent equivalence algebraically.
DEF_LOGICAL_FUNC_CONSTR( Equivalence, VarArray2,
                        "r = (v1 <==> v2)");


////////////////////////////////////////////////////////////////////////
DEF_LOGICAL_FUNC_CONSTR( Not, VarArray1,
                                   "r = !v");

////////////////////////////////////////////////////////////////////////
/// \brief DivConstraint and DivExpression
DEF_NUMERIC_FUNC_CONSTR( Div, VarArray2,
                                  "r = v1 / v2 and v2!=0");

////////////////////////////////////////////////////////////////////////
DEF_NUMERIC_FUNC_CONSTR( IfThen, VarArrayN<3>,
                         "Expr-valued: if (cond) then (expr1) else (expr2)");

////////////////////////////////////////////////////////////////////////
DEF_LOGICAL_FUNC_CONSTR( Implication, VarArrayN<3>,
                         "Logic-valued: if (cond) then (con1) else (con2)");

////////////////////////////////////////////////////////////////////////
DEF_LOGICAL_FUNC_CONSTR( AllDiff, VarArray,
                                  "alldiff({})");

////////////////////////////////////////////////////////////////////////
DEF_NUMERIC_FUNC_CONSTR_WITH_PRM( NumberofConst,
                                  VarArray, DblParamArray1,
                                  "numberof_const(k, {x0...xn})");

////////////////////////////////////////////////////////////////////////
DEF_NUMERIC_FUNC_CONSTR( NumberofVar, VarArray,
                                  "numberof_var(x0, {x1...xn})");

////////////////////////////////////////////////////////////////////////
DEF_NUMERIC_FUNC_CONSTR( Count, VarArray,
                                   "count({x0...xn})");


//////////////////// NONLINEAR FUNCTIONS //////////////////////
////////////////////////////////////////////////////////////////////////
DEF_NUMERIC_FUNC_CONSTR( Exp, VarArray1,
                                   "r = exp(v)");

////////////////////////////////////////////////////////////////////////
DEF_NUMERIC_FUNC_CONSTR_WITH_PRM( ExpA,
                  VarArray1, DblParamArray1, "r = a**v");

////////////////////////////////////////////////////////////////////////
DEF_NUMERIC_FUNC_CONSTR( Log, VarArray1,
                                   "r = log(v)");

////////////////////////////////////////////////////////////////////////
DEF_NUMERIC_FUNC_CONSTR_WITH_PRM( LogA,
                  VarArray1, DblParamArray1, "r = log(v)/log(a)");

////////////////////////////////////////////////////////////////////////
DEF_NUMERIC_FUNC_CONSTR_WITH_PRM( PowConstExp,
                                 VarArray1, DblParamArray1,
                                 "r = v ** a (a is constant)");

////////////////////////////////////////////////////////////////////////
DEF_NUMERIC_FUNC_CONSTR_WITH_PRM( SignpowConstExp,
                                 VarArray1, DblParamArray1,
                                 "r = sign(v) * abs(v) ^ a (a is constant)");

////////////////////////////////////////////////////////////////////////
DEF_NUMERIC_FUNC_CONSTR_WITH_PRM( Logistic,
                                 VarArray1, ParamArray0,
                                 "r = logistic_fn(v) := 1 / (1+e^(-v))");

////////////////////////////////////////////////////////////////////////
DEF_NUMERIC_FUNC_CONSTR_WITH_PRM( Pow,
                                 VarArray2, ParamArray0,
                                 "r = x ** y (both variable)");

////////////////////////////////////////////////////////////////////////
DEF_NUMERIC_FUNC_CONSTR( Sin, VarArray1,
                                   "r = sin(v)");

////////////////////////////////////////////////////////////////////////
DEF_NUMERIC_FUNC_CONSTR( Cos, VarArray1,
                                   "r = cos(v)");

////////////////////////////////////////////////////////////////////////
DEF_NUMERIC_FUNC_CONSTR( Tan, VarArray1,
                                   "r = tan(v)");

////////////////////////////////////////////////////////////////////////
DEF_NUMERIC_FUNC_CONSTR( Asin, VarArray1,
                                   "r = asin(v)");

////////////////////////////////////////////////////////////////////////
DEF_NUMERIC_FUNC_CONSTR( Acos, VarArray1,
                                   "r = acos(v)");

////////////////////////////////////////////////////////////////////////
DEF_NUMERIC_FUNC_CONSTR( Atan, VarArray1,
                                   "r = atan(v)");

////////////////////////////////////////////////////////////////////////
DEF_NUMERIC_FUNC_CONSTR( Sinh, VarArray1,
                                   "r = sinh(v)");

////////////////////////////////////////////////////////////////////////
DEF_NUMERIC_FUNC_CONSTR( Cosh, VarArray1,
                                   "r = cosh(v)");

////////////////////////////////////////////////////////////////////////
DEF_NUMERIC_FUNC_CONSTR( Tanh, VarArray1,
                                   "r = tanh(v)");

////////////////////////////////////////////////////////////////////////
DEF_NUMERIC_FUNC_CONSTR( Asinh, VarArray1,
                                   "r = asinh(v)");

////////////////////////////////////////////////////////////////////////
DEF_NUMERIC_FUNC_CONSTR( Acosh, VarArray1,
                                   "r = acosh(v)");

////////////////////////////////////////////////////////////////////////
DEF_NUMERIC_FUNC_CONSTR( Atanh, VarArray1,
                                   "r = atanh(v)");

////////////////////////////////////////////////////////////////////////
DEF_NUMERIC_FUNC_CONSTR_WITH_PRM( Call, VarArray, IntParamArray1,
                                 "r = call[function[i]](args)");

////////////////////////////////////////////////////////////////////////
DEF_NUMERIC_FUNC_CONSTR_WITH_PRM( SDPDotProd, VarArray0, IntParamArray,
                                 "r = sum {i} SDPDotProduct(A[p[i*2]], X[p[i*2+1]]), "
                                 "with A[p[i*2]] coef matrix and "
                                 "X[p[i*2+1]] SDP variable. The parameter vector p "
                                 "has consecutive pairs of indexes "
                                 "of A's and X's.");


/// Not using: var1 != var2.
/// Represented by Not { Eq0Constraint... }
////////////////////////////////////////////////////////////////////////
// DEF_LOGICAL_FUNC_CONSTR( NE__unused, VarArray2,
//                                   "r = (v1 != v2)");


////////////////////////////////////////////////////////////////////////
DEF_CONDITIONAL_CONSTRAINT_WRAPPER(CondLinLT, LinConLT);
/// Typedef CondLinConLT for legacy
using CondLinConLT = CondLinLTConstraint;
/// Typedef CondLTExpression: its LHS is an expression
using CondLTExpression = CondLinLTExpression;
////////////////////////////////////////////////////////////////////////
DEF_CONDITIONAL_CONSTRAINT_WRAPPER(CondLinLE, LinConLE);
using CondLinConLE = CondLinLEConstraint;
using CondLEExpression = CondLinLEExpression;
////////////////////////////////////////////////////////////////////////
DEF_CONDITIONAL_CONSTRAINT_WRAPPER(CondLinEQ, LinConEQ);
using CondLinConEQ = CondLinEQConstraint;
using CondEQExpression = CondLinEQExpression;
////////////////////////////////////////////////////////////////////////
DEF_CONDITIONAL_CONSTRAINT_WRAPPER(CondLinGE, LinConGE);
using CondLinConGE = CondLinGEConstraint;
using CondGEExpression = CondLinGEExpression;
////////////////////////////////////////////////////////////////////////
DEF_CONDITIONAL_CONSTRAINT_WRAPPER(CondLinGT, LinConGT);
using CondLinConGT = CondLinGTConstraint;
using CondGTExpression = CondLinGTExpression;

////////////////////////////////////////////////////////////////////////
DEF_CONDITIONAL_CONSTRAINT_WRAPPER(CondQuadLT, QuadConLT);
using CondQuadConLT = CondQuadLTConstraint;
////////////////////////////////////////////////////////////////////////
DEF_CONDITIONAL_CONSTRAINT_WRAPPER(CondQuadLE, QuadConLE);
using CondQuadConLE = CondQuadLEConstraint;
////////////////////////////////////////////////////////////////////////
DEF_CONDITIONAL_CONSTRAINT_WRAPPER(CondQuadEQ, QuadConEQ);
using CondQuadConEQ = CondQuadEQConstraint;
////////////////////////////////////////////////////////////////////////
DEF_CONDITIONAL_CONSTRAINT_WRAPPER(CondQuadGE, QuadConGE);
using CondQuadConGE = CondQuadGEConstraint;
////////////////////////////////////////////////////////////////////////
DEF_CONDITIONAL_CONSTRAINT_WRAPPER(CondQuadGT, QuadConGT);
using CondQuadConGT = CondQuadGTConstraint;


////////////////////////////////////////////////////////////////////////
/// Linear Functional Constraint: r = affine_expr
class LinearFunctionalConstraint :
    public FunctionalConstraint,
    public NumericFunctionalConstraintTraits {
  AffineExpr affine_expr_;
public:
  /// Constraint type name
  static const char* GetTypeName()
  { return "LinearFunctionalConstraint"; }
  /// Typedef Arguments
  using Arguments = AffineExpr;
  /// using GetResultVar()
  using FunctionalConstraint::GetResultVar;
  /// A constructor ignoring result variable:
  /// use AssignResultToArguments() then.
  /// Not sorting+merging
  LinearFunctionalConstraint(AffineExpr&& ae) noexcept :
    affine_expr_(std::move(ae)) {
  }
  /// Constructor with result variable known.
  /// Not sorting+merging
  LinearFunctionalConstraint(int r, AffineExpr&& ae) noexcept :
    FunctionalConstraint(r), affine_expr_(std::move(ae)) { }
  /// Get the affine expr
  const AffineExpr& GetAffineExpr() const { return affine_expr_; }
  /// Get the arguments (affine expr)
  const Arguments& GetArguments() const { return GetAffineExpr(); }
  /// Modify constant term
  void add_to_constant(double d) { affine_expr_.add_to_constant(d); }
  /// Produce corresp linear constraint
  LinConEQ to_linear_constraint() const {
    const auto& ae = GetAffineExpr();
    auto le = ae.GetLinTerms();
    le.add_term(-1.0, FunctionalConstraint::GetResultVar());
    return { std::move(le), -ae.constant_term() };
  }
};


/// Typedef NLAffineExpression.
/// Like AffineExpr, but arguments are expressions.
using NLAffineExpression = ExprWrapper<LinearFunctionalConstraint>;

/// Write LFC without name.
template <class Writer>
inline void WriteModelItem(Writer& wrt,
                    const LinearFunctionalConstraint& lfc,
                    ItemNamer& vnam) {
  wrt << vnam.at(lfc.GetResultVar()) << " == ";
  WriteModelItem(wrt, lfc.GetArguments(), vnam);
}

/// Write a LinFuncCon
inline void WriteJSON(JSONW jw,
                      const LinearFunctionalConstraint& lfc) {
  jw["res_var"] = lfc.GetResultVar();
  WriteJSON(jw["expr"], lfc.GetAffineExpr());
}


////////////////////////////////////////////////////////////////////////
/// Quadratic Functional Constraint: r = quad_expr
class QuadraticFunctionalConstraint :
    public FunctionalConstraint,
    public NumericFunctionalConstraintTraits {
  QuadraticExpr quad_expr_;
public:
  /// Constraint type name
  static const char* GetTypeName() { return "QuadraticFunctionalConstraint"; }
  /// Typedef Arguments
  using Arguments = QuadraticExpr;
  /// using GetResultVar()
  using FunctionalConstraint::GetResultVar;

  /// A constructor ignoring result variable:
  /// use AssignResultToArguments() then.
  /// Not sorting+merging
  QuadraticFunctionalConstraint(QuadraticExpr&& qe) noexcept :
    quad_expr_(std::move(qe)) { }

  /// Constructor: known result var + body.
  /// Not sorting+merging
  QuadraticFunctionalConstraint(int r, QuadraticExpr&& qe) noexcept :
    FunctionalConstraint(r), quad_expr_(std::move(qe)) { }

  /// Getter: quad expr
  const QuadraticExpr& GetQuadExpr() const { return quad_expr_; }
  /// GetArguments(): get quad expr
  const Arguments& GetArguments() const { return GetQuadExpr(); }

  /// Get constant
  double constant_term() const { return quad_expr_.constant_term(); }
  /// Get linear part
  const LinTerms& GetLinTerms() const { return quad_expr_.GetLinTerms(); }
  /// Get linear part
  const QuadTerms& GetQPTerms() const { return quad_expr_.GetQPTerms(); }
  /// Modify constant term
  void add_to_constant(double d) { quad_expr_.add_to_constant(d); }

  /// add respective static constraint to a converter.
  /// Use >=< depending on context.
  template <class Converter>
  Context AddQuadraticConstraint(Converter& cvt) const {
    auto le = GetQuadExpr().GetLinTerms();
    le.add_term(-1.0, FunctionalConstraint::GetResultVar());
    auto qt = GetQuadExpr().GetQPTerms();
    if (GetContext().IsMixed())
      cvt.AddConstraint( QuadConEQ{ { std::move(le), std::move(qt) },
      -GetQuadExpr().constant_term() } );
    else if (GetContext().HasPositive())
      cvt.AddConstraint( QuadConGE{ { std::move(le), std::move(qt) },
      -GetQuadExpr().constant_term() } );
    else if (GetContext().HasNegative())
      cvt.AddConstraint( QuadConLE{ { std::move(le), std::move(qt) },
      -GetQuadExpr().constant_term() } );
    else
      MP_RAISE("QuadraticFuncCon: no context");
    return GetContext();
  }
};


/// Typedef NLQuadExpression.
/// Like QuadExpr, but arguments are expressions.
using NLQuadExpression = ExprWrapper<QuadraticFunctionalConstraint>;


/// Function definition
class FuncDef {
public:
  /// Construct
  FuncDef(std::string_view nm, int na) : name_(nm), n_args_(na) { }
  /// Name
  const char* Name() const { return name_.c_str(); }
  /// N args
  int NumArgs() const { return n_args_; }
  /// Type: always numeric
  int Type() const { return 0; }
private:
  std::string name_;
  int n_args_;
};


/// Shortcut to make an LFC from linear body and constant
inline LinearFunctionalConstraint
MakeFunctionalConstraint(AffineExpr ae) {
  return {std::move(ae)};
}

/// Shortcut to make a QFC from quadr body and constant
inline QuadraticFunctionalConstraint
MakeFunctionalConstraint(QuadraticExpr qe) {
  return {std::move(qe)};
}

/// Write LFC without name.
template <class Writer>
inline void WriteModelItem(Writer& wrt,
                    const QuadraticFunctionalConstraint& qfc,
                    ItemNamer& vnam) {
  wrt << vnam.at(qfc.GetResultVar()) << " == ";
  WriteModelItem(wrt, qfc.GetArguments(), vnam);
}

/// Write a QuadrFuncCon
inline void WriteJSON(JSONW jw,
                      const QuadraticFunctionalConstraint& qfc) {
  jw["res_var"] = qfc.GetResultVar();
  WriteJSON(jw["expr"], qfc.GetQuadExpr());
}


struct PLPoints;

////////////////////////////////////////////////////////////////////////
/// AMPL represents PWL by a list of slopes
/// and breakpoints between them, assuming (X0,Y0) is on the line
class PLSlopes {
public:
  /// Default constructor
  PLSlopes() { }
  /// Construct from breakpoints, slopes, and a sample point
  template <class Vec>
  PLSlopes(Vec&& bp, Vec&& sl, double x, double y) noexcept :
    breakpoints_(std::forward<Vec>(bp)), slopes_(std::forward<Vec>(sl)),
    X0_(x), Y0_(y) { assert(check()); }
	/// Construct from PLPoints
	PLSlopes(const PLPoints& );
	/// Check if have information
	bool empty() const { assert(!GetNBP() || check()); return !GetNBP(); }
	/// Get breakpoints
  const std::vector<double>& GetBP() const { return breakpoints_; }
  /// Get slopes
  const std::vector<double>& GetSlopes() const { return slopes_; }
  /// Get sample point's X
  double GetX0() const { return X0_; }
	/// Get sample point's Y
  double GetY0() const { return Y0_; }
  /// Get number of bp
  int GetNBP() const { return (int)GetBP().size(); }
  /// Get number slopes
  int GetNSlopes() const { return (int)GetSlopes().size(); }
  /// Validate
  bool check() const { return GetNBP()>0 && GetNSlopes()==GetNBP()+1; }

private:
  std::vector<double> breakpoints_, slopes_;
  double X0_, Y0_;                             // some point on the PWL
};


/// Representing a PWL by points
struct PLPoints {
  /// Check if have information
  bool empty() const { return x_.empty(); }
  /// size()
	int size() const
	{ assert(x_.size()==y_.size()); return (int)x_.size(); }
  /// The x, y coordinates of the PL function
  std::vector<double> x_, y_;
  /// Default construct
  PLPoints() { }
  /// Construct from 2 vectors
  PLPoints(std::vector<double> x, std::vector<double> y) :
    x_{std::move(x)}, y_{std::move(y)} { }
  /// Construct from PLSlopes
  PLPoints(const PLSlopes& pls);
  /// Add point
  void AddPoint(double x, double y) {
    if (!empty())
      assert(x > x_.back());
    if (empty() || x>x_.back()+1e-4) {  // skip near points for Gurobi
      if (size()>=2 &&                  // simple check: 3rd equal y
         y_[size()-1]==y && y_[size()-2]==y) {
        x_.back() = x;                  // update last x
      } else {
        x_.push_back(x);
        y_.push_back(y);
      }
    }
  }
  /// Clear out
  void clear() {
    x_.clear();
    y_.clear();
  }
  /// Get preslope
  double PreSlope() const {
    return (x_.size()<=1 || x_[0]>=x_[1]) ?
        0.0 :
        (y_[1]-y_[0]) / (x_[1]-x_[0]);
  }
  /// Get postslope
  double PostSlope() const {
    auto i1=x_.size()-1, i0=i1-1;
    return (x_.size()<=1 || x_[i0]>=x_[i1]) ?
        0.0 :
        (y_[i1]-y_[i0]) / (x_[i1]-x_[i0]);
  }

  /// operator==
  bool operator==(const PLPoints& plp) const
  { return x_==plp.x_ && y_==plp.y_; }
};


/// Parameters of a PL constraint: keeps any or both of
/// PLSlopes and PLPoints
class PLConParams {
public:
  /// Construct from PLSlopes
  PLConParams(PLSlopes pls) : pls_(std::move(pls)) { }
  /// Construct from PLPoints
  PLConParams(PLPoints plp) : plp_(std::move(plp)) { }
  /// Produce PLPoints, either stored or from PLSlopes
  const PLPoints& GetPLPoints() const {
    if (plp_.empty())
      plp_ = pls_;
    return plp_;
  }
	/// Produce PLSlopes, either stored or from PLPoints
	const PLSlopes& GetPLSlopes() const {
		if (pls_.empty())
			pls_ = plp_;
		return pls_;
	}
	/// operator==
  bool operator==(const PLConParams& plcp) const
  { return GetPLPoints() == plcp.GetPLPoints(); }

private:
	mutable PLSlopes pls_;
  mutable PLPoints plp_;
};


/// Is PL convex?
bool IsConvex(const PLSlopes& pls);
/// Is PL concave?
bool IsConcave(const PLSlopes& pls);

/// Is PL convex?
inline bool IsConvex(const PLConParams& plp)
{ return IsConvex(plp.GetPLSlopes()); }
/// Is PL concave?
inline bool IsConcave(const PLConParams& plp)
{ return IsConcave(plp.GetPLSlopes()); }

/// Define PLConstraint
DEF_NUMERIC_FUNC_CONSTR_WITH_PRM( PL,
                  VarArray1, PLConParams, "r = piecewise_linear(x)");

/// Write flat expr/obj/con parameters:
/// PLConParams
template <class Writer>
void WriteModelItemParameters(
    Writer& wrt, const PLConParams& plp) {
  wrt << "plpoints[";
  const auto& p = plp.GetPLPoints();
  for (size_t i=0; i<p.x_.size(); ++i) {
    if (i)
      wrt << ", ";
    wrt << '(' << p.x_[i] << ", " << p.y_[i] << ')';
  }
  wrt << ']';
}

/// Write PLConParams
inline void WriteJSON(JSONW jw,
                      const PLConParams& plp) {
  const auto& p = plp.GetPLPoints();
  jw["pl_x"] = p.x_;
  jw["pl_y"] = p.y_;
}

} // namespace mp

#endif // CONSTRAINTS_FUNCTIONAL_H
