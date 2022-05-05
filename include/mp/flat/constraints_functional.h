#ifndef CONSTRAINTS_FUNCTIONAL_H
#define CONSTRAINTS_FUNCTIONAL_H

/**
  * Functional flat constraints
  */

#include "mp/flat/constraints_static.h"


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
/// TODO Keep the full expression of v1/2,
/// or use it in MIP to multiply out
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


////////////////////////////////////////////////////////////////////////
/// OLD: Storing AffExp instead of LinConEQ because big-M is straightforwardly
/// computed for (aff_exp) <= 0:
/// b -> ae<=0 is linearized as ae <= ub(ae)*(1-b) <==> le-d <= (ub(le)-d)*(1-b)
/// <==> le <= d + ub(le) - d + (d-ub(le))*b
/// If we stored LinConEQ:
/// b -> lin_exp<=d would be linearized as
/// le <= d + (ub(le)-d)*(1-b)  <==>
/// le <= d + ub_le - d - ub_le*b + d*b  <==>
/// le <= ub_le + (d-ub_le)*b.  Not too complex.
/// Keep it with AffineExpr, indicators need that
/// and we don't want quadratics with big-M's?
/// TODO Diff to Indicator?

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
DEF_CONDITIONAL_CONSTRAINT_WRAPPER(CondQuadConLT, QuadConLT);

////////////////////////////////////////////////////////////////////////
DEF_CONDITIONAL_CONSTRAINT_WRAPPER(CondQuadConLE, QuadConLE);

////////////////////////////////////////////////////////////////////////
DEF_CONDITIONAL_CONSTRAINT_WRAPPER(CondQuadConEQ, QuadConEQ);



/// TODO use macros for FLC / FQC too?

////////////////////////////////////////////////////////////////////////
/// Linear Functional Constraint: r = affine_expr
class LinearFunctionalConstraint :
    public FunctionalConstraint {
  AffineExpr affine_expr_;
public:
  static const char* GetName() { return "LinearFunctionalConstraint"; }
  using Arguments = AffineExpr;
  using FunctionalConstraint::GetResultVar;
  /// A constructor ignoring result variable: use AssignResultToArguments() then
  LinearFunctionalConstraint(AffineExpr&& ae) noexcept :
    affine_expr_(std::move(ae)) {  // TODO sort+merge elements?
  }
  LinearFunctionalConstraint(int r, AffineExpr&& ae) noexcept :
    FunctionalConstraint(r), affine_expr_(std::move(ae)) {
    /// TODO sort+merge elements
  }
  const AffineExpr& GetAffineExpr() const { return affine_expr_; }
  const Arguments& GetArguments() const { return GetAffineExpr(); }
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
    public FunctionalConstraint {
  QuadraticExpr quad_expr_;
public:
  static const char* GetName() { return "QuadraticFunctionalConstraint"; }
  using Arguments = QuadraticExpr;
  using FunctionalConstraint::GetResultVar;

  /// A constructor ignoring result variable: use AssignResultToArguments() then
  QuadraticFunctionalConstraint(QuadraticExpr&& qe) noexcept :
    quad_expr_(std::move(qe)) {
    /// TODO sort+merge elements
  }

  /// Constructor: result var + body
  QuadraticFunctionalConstraint(int r, QuadraticExpr&& qe) noexcept :
    FunctionalConstraint(r), quad_expr_(std::move(qe)) {
    /// TODO sort+merge elements
  }

  /// Getters
  const QuadraticExpr& GetQuadExpr() const { return quad_expr_; }
  const Arguments& GetArguments() const { return GetQuadExpr(); }

  /// produce respective static constraint
  QuadConEQ to_quadratic_constraint() const {
    auto le = GetQuadExpr().GetLinTerms();
    le.add_term(-1.0, FunctionalConstraint::GetResultVar());
    auto qt = GetQuadExpr().GetQPTerms();
    return { { std::move(le), std::move(qt) },
      -GetQuadExpr().constant_term() };
  }
};


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

} // namespace mp

#endif // CONSTRAINTS_FUNCTIONAL_H
