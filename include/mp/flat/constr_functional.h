#ifndef CONSTRAINTS_FUNCTIONAL_H
#define CONSTRAINTS_FUNCTIONAL_H

/**
  * Functional flat constraints
  */

#include "mp/error.h"
#include "mp/flat/constr_static.h"


namespace mp {


////////////////////////////////////////////////////////////////////////
DEF_NUMERIC_FUNC_CONSTR( MaxConstraint, VarArray,
                                   "r = max(v1, v2, ..., vn)");

////////////////////////////////////////////////////////////////////////
DEF_NUMERIC_FUNC_CONSTR( MinConstraint, VarArray,
                                   "r = min(v1, v2, ..., vn)");

////////////////////////////////////////////////////////////////////////
DEF_NUMERIC_FUNC_CONSTR( AbsConstraint, VarArray1,
                                   "r = abs(v)");

////////////////////////////////////////////////////////////////////////
DEF_LOGICAL_FUNC_CONSTR( AndConstraint, VarArray,
                                   "r = forall({vi})");

////////////////////////////////////////////////////////////////////////
DEF_LOGICAL_FUNC_CONSTR( OrConstraint, VarArray,
                                   "r = exists({vi})");


////////////////////////////////////////////////////////////////////////
DEF_LOGICAL_FUNC_CONSTR( NotConstraint, VarArray1,
                                   "r = !v");

////////////////////////////////////////////////////////////////////////
/// \brief DivConstraint
DEF_NUMERIC_FUNC_CONSTR( DivConstraint, VarArray2,
                                  "r = v1 / v2 and v2!=0");

////////////////////////////////////////////////////////////////////////
DEF_LOGICAL_FUNC_CONSTR( IfThenConstraint, VarArrayN<3>,
                                  "if (cond) then (expr1) else (expr2)");

////////////////////////////////////////////////////////////////////////
DEF_LOGICAL_FUNC_CONSTR( AllDiffConstraint, VarArray,
                                  "alldiff({})");

////////////////////////////////////////////////////////////////////////
DEF_NUMERIC_FUNC_CONSTR_WITH_PRM( NumberofConstConstraint,
                                  VarArray, DblParamArray1,
                                  "numberof_const(k (=x0), {x1...xn})");

////////////////////////////////////////////////////////////////////////
DEF_NUMERIC_FUNC_CONSTR( NumberofVarConstraint, VarArray,
                                  "numberof_var(x0, {x1...xn})");

////////////////////////////////////////////////////////////////////////
DEF_NUMERIC_FUNC_CONSTR( CountConstraint, VarArray,
                                   "count({x0...xn})");


//////////////////// NONLINEAR FUNCTIONS //////////////////////
////////////////////////////////////////////////////////////////////////
DEF_NUMERIC_FUNC_CONSTR( ExpConstraint, VarArray1,
                                   "r = exp(v)");

////////////////////////////////////////////////////////////////////////
DEF_NUMERIC_FUNC_CONSTR_WITH_PRM( ExpAConstraint,
                  VarArray1, DblParamArray1, "r = a**v");

////////////////////////////////////////////////////////////////////////
DEF_NUMERIC_FUNC_CONSTR( LogConstraint, VarArray1,
                                   "r = log(v)");

////////////////////////////////////////////////////////////////////////
DEF_NUMERIC_FUNC_CONSTR_WITH_PRM( LogAConstraint,
                  VarArray1, DblParamArray1, "r = log(v)/log(a)");

////////////////////////////////////////////////////////////////////////
DEF_NUMERIC_FUNC_CONSTR_WITH_PRM( PowConstraint,
                  VarArray1, DblParamArray1, "r = v ** a");

////////////////////////////////////////////////////////////////////////
DEF_NUMERIC_FUNC_CONSTR( SinConstraint, VarArray1,
                                   "r = sin(v)");

////////////////////////////////////////////////////////////////////////
DEF_NUMERIC_FUNC_CONSTR( CosConstraint, VarArray1,
                                   "r = cos(v)");

////////////////////////////////////////////////////////////////////////
DEF_NUMERIC_FUNC_CONSTR( TanConstraint, VarArray1,
                                   "r = tan(v)");


/// Not using: var1 != var2.
/// Represented by Not { Eq0Constraint... }
////////////////////////////////////////////////////////////////////////
// DEF_LOGICAL_FUNC_CONSTR( NEConstraint__unused, VarArray2,
//                                   "r = (v1 != v2)");


////////////////////////////////////////////////////////////////////////
DEF_CONDITIONAL_CONSTRAINT_WRAPPER(CondLinConLT, LinConLT);

////////////////////////////////////////////////////////////////////////
DEF_CONDITIONAL_CONSTRAINT_WRAPPER(CondLinConLE, LinConLE);

////////////////////////////////////////////////////////////////////////
DEF_CONDITIONAL_CONSTRAINT_WRAPPER(CondLinConEQ, LinConEQ);

////////////////////////////////////////////////////////////////////////
DEF_CONDITIONAL_CONSTRAINT_WRAPPER(CondLinConGE, LinConGE);

////////////////////////////////////////////////////////////////////////
DEF_CONDITIONAL_CONSTRAINT_WRAPPER(CondLinConGT, LinConGT);


////////////////////////////////////////////////////////////////////////
DEF_CONDITIONAL_CONSTRAINT_WRAPPER(CondQuadConLT, QuadConLT);

////////////////////////////////////////////////////////////////////////
DEF_CONDITIONAL_CONSTRAINT_WRAPPER(CondQuadConLE, QuadConLE);

////////////////////////////////////////////////////////////////////////
DEF_CONDITIONAL_CONSTRAINT_WRAPPER(CondQuadConEQ, QuadConEQ);

////////////////////////////////////////////////////////////////////////
DEF_CONDITIONAL_CONSTRAINT_WRAPPER(CondQuadConGE, QuadConGE);

////////////////////////////////////////////////////////////////////////
DEF_CONDITIONAL_CONSTRAINT_WRAPPER(CondQuadConGT, QuadConGT);



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
  /// Produce corresp linear constraint
  LinConEQ to_linear_constraint() const {
    const auto& ae = GetAffineExpr();
    auto le = ae.GetLinTerms();
    le.add_term(-1.0, FunctionalConstraint::GetResultVar());
    return { std::move(le), -ae.constant_term() };
  }
};

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

  /// add respective static constraint to a converter.
  /// Use >=< depending on context.
  template <class Converter>
  void AddQuadraticConstraint(Converter& cvt) const {
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
  }
};


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
  /// Get breakpoints
  const std::vector<double>& GetBP() const { return breakpoints_; }
  /// Get slopes
  const std::vector<double>& GetSlopes() const { return slopes_; }
  /// Get sample point's X
  double GetX0() const { return X0_; }
  /// Get sample poont's Y
  double GetY0() const { return Y0_; }
  /// Get number of bp
  int GetNBP() const { return GetBP().size(); }
  /// Get number slopes
  int GetNSlopes() const { return GetSlopes().size(); }
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
  /// The x, y coordinates of the PL function
  std::vector<double> x_, y_;
  /// Default construct
  PLPoints() { }
  /// Construct from PLSlopes
  PLPoints(const PLSlopes& pls);
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
  operator PLPoints() const { return plp_.empty() ? pls_ : plp_; }

private:
  PLSlopes pls_;
  PLPoints plp_;
};

DEF_NUMERIC_FUNC_CONSTR_WITH_PRM( PLConstraint,
                  VarArray1, PLConParams, "r = piecewise_linear(x)");

} // namespace mp

#endif // CONSTRAINTS_FUNCTIONAL_H
