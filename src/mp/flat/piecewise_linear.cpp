/*
 MP PL facilities implementation.

 Copyright (C) 2022 AMPL Optimization Inc

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization Inc disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Gleb Belov
 */
#include <cmath>
#include <array>
#include <vector>
#include <set>
#include <string>
#include <algorithm>

#include "mp/flat/constr_std.h"
#include "mp/flat/redef/MIP/core/lin_approx_core.h"

namespace mp {

/// Declare function-specific approximator
template <class Con>
class PLApproximator;

/// Define generic approximation call
template <class FuncCon>
void PLApproximate(const FuncCon& con,
                             PLApproxParams& laPrm) {
  PLApproximator<FuncCon> pla{con, laPrm};
  pla.Run();
}


/// Interface class for univariate function's
/// piecewise-linear smooth approximator:
/// domain, periodicity, derivatives.
/// The domain (or the period, if relevant)
/// should be split into subintervals where both
/// the function and the 1st drivative are monotone.
/// @param FuncCon: the functional constraint
template <class FuncCon>
class BasicPLApproximator {
public:
  /// Constructor
  BasicPLApproximator(const FuncCon& con, PLApproxParams& laPrm) :
    con_(con), laPrm_(laPrm) { }
  /// Virtual destructor
  virtual ~BasicPLApproximator() { }
  /// Constraint type name
  virtual const char* GetConTypeName() const
  { return con_.GetTypeName(); }
  /// Run the approximation
  void Run();


  ////////////////////// USER API ///////////////////////////
  /// Typically, this is enough to define an approximator ///
  ///////////////////////////////////////////////////////////
protected:
  /// Largest accepted argument domain.
  /// If the argument domain is not contained, fail.
  /// Use if the approximation is artificially restricted,
  /// e.g., only positive bases with negative powers
  virtual Range GetLargestAcceptedArgumentRange() const
  { return {-1e100, 1e100}; }
  /// Ample, but realistic, function graph domain
  /// (values should not overflow when computing (pre-)image).
  /// Depends on the solver's integrality and feasibility
  /// tolerances. For values involved in indicator-like
  /// constraints, should be below 1e6, sometimes 1e5
  virtual FuncGraphDomain GetFuncGraphDomain() const = 0;
  /// Monotone? Then we can clip image range.
  /// This should be superseded by walking thru subintervals.
  virtual bool IsMonotone() const { return false; }
  /// Normally periodic?
  virtual bool IsPeriodic() const { return false; }
  /// The base period interval, if relevant.
  /// For example, [-pi/2, pi/2] for tan
  virtual Range GetDefaultPeriod() const
  { return {-1e100, 1e100}; }
  /// Breakpoint list, defines subintervals: [bp0, bp1], [bp1, bp2], ...
  /// Include domain/period margins,
  /// either for the global domain,
  /// or, if IsPeriodic(), in the GetDefaultPeriod().
  /// Should have realistic margins and
  /// internal breakpoints for (local) extrema and
  /// deflection points.
  /// For example, [ -pi/2+eps, 0.0, pi/2-eps ] for tan
  using BreakpointList = std::vector<double>;
  /// Default breakpoints for the function.
  /// The default implementation returns the domain margins.
  virtual BreakpointList GetDefaultBreakpoints() const {
    auto grDom = GetFuncGraphDomain();
    return {grDom.lbx, grDom.ubx};
  }
  /// Evaluate function at \a x
  virtual double eval(double x) const = 0;
  /// Evaluate inv(f) at \a y. Consider current subinterval
  virtual double inverse(double y) const = 0;
  /// Evaluate f' at \a x
  virtual double eval_1st(double x) const = 0;
  /// Evaluate inv(f') at \a y. Consider current subinterval
  virtual double inverse_1st(double y) const = 0;
  /// Evaluate f'' at \a x
  virtual double eval_2nd(double x) const = 0;


  ////////////////// INTERNAL API /////////////////////

  virtual double GetEps() const { return 1e-6; }

  /// Constraint parameters, such as the logarithm base
  const typename FuncCon::Parameters& GetConParams() const {
    return con_.GetParameters();
  }
  /// Function graph domain:
  /// 1. Check if x outside of the largest accepted range
  /// 2. Clip to the user-provided domain, and
  /// 3. to its function values.
  /// The latter only if the function is monotone
  virtual void ClipFuncGraphDomain() {
    auto rngAcc = GetLargestAcceptedArgumentRange();
    auto& lbx = laPrm_.grDom.lbx;
    auto& ubx = laPrm_.grDom.ubx;
    MP_ASSERT_ALWAYS(lbx >= rngAcc.lb && rngAcc.ub >= ubx,
                     fmt::format(
                       "PLApproximator<{}>: argument range "
                       "[{}, {}] outside of the accepted [{}, {}]",
                       GetConTypeName(),
                       lbx, ubx, rngAcc.lb, rngAcc.ub));
    laPrm_.grDom.intersect(GetFuncGraphDomain());
    if (IsMonotone())     // TODO in the reduced domain only
      ClipWithFunctionValues(laPrm_.grDom);
    lbx_ = lbx;
    ubx_ = ubx;
    laPrm_.grDomOut = laPrm_.grDom;        // communicate new domain
  }
  /// Check domain, throw if infeas, or return trivial PL
  /// @throw if infeasible
  /// @return false iff trivial
  /// @param result: the trivial PL
  virtual bool CheckDomainReturnFalseIfTrivial(PLPoints& result);

  /// Consider and initialize periodic approximation
  /// @return Whether decided to go periodic
  virtual bool InitPeriodic();
  /// Initialize non-periodic approximation
  virtual void InitNonPeriodic();
  /// initialize loop through subintervals
  virtual void InitSubintervalLoop();
  /// init next subinterval
  virtual bool NextSubinterval();
  /// Approximate in the chosen subinterval
  virtual void ApproximateSubinterval();
  /// ClipGraphFinal.
  /// Is necessary beacuse initial domain clipping is
  /// only done for monotone functions.
  virtual void ClipGraphFinal();
  /// See if we can use argument's integrality
  virtual void ConsiderIntegrality();

  /// Provide initial step length to the right from
  /// the current point, trying have the error nearly ok
  virtual double ComputeInitialStepLength(double x0);
  /// Increase step length while error small enough
  /// @param dx0: in-out parameter, step length
  virtual void IncreaseStepWhileErrorSmallEnough(
      double x0, double f0, double& dx0);
  /// Decrease step length while error too high
  /// @param dx0: in-out parameter, step length
  virtual void DecreaseStepWhileErrorTooBig(
      double x0, double f0, double& dx0);
  /// Compare error on the given linear segment
  /// to the relevant upper bound.
  /// The error can be relative or absolute.
  /// @return -1, 0, 1 for <, ==, >
  virtual int CompareError(
      double x0, double y0, double x1, double y1);

  /// Maximal absolute error (along Y axis) to the linear segment
  /// (x0, y0) -> (x1, y1).
  /// Standard implementation checks in the ends
  /// and in the middle value point
  virtual double maxErrorAbs(
      double x0, double y0, double x1, double y1);
  /// Maximal relative error (along Y axis) of the linear segment
  /// (x0, y0) -> (x1, y1) to the function.
  /// For function values between +-1, the error bound is taken
  /// absolute
  virtual double maxErrorRelAbove1(
      double x0, double y0, double x1, double y1);
  /// Clip the graph domain, considering the actual function.
  virtual void ClipWithFunctionValues(FuncGraphDomain& grDom) const;

  /// Call inverse() and check result to be in lb_sub()..ub_sub()
  virtual double inverse_with_check(double y) const;
  /// Call inverse_1st and check result to be in lb_sub()..ub_sub()
  virtual double inverse_1st_with_check(double y) const;

  /// Function interval to approximate: lb
  double lbx() const { return lbx_; }
  /// Function interval to approximate: ub
  double ubx() const { return ubx_; }

  /// Current subinterval index
  int GetSubIntvIndex() const { return iSubIntv_; }

  /// Subinterval to approximate: lb
  double lb_sub() const { return breakpoints_.at(iSubIntv_); }
  /// Subinterval to approximate: ub
  double ub_sub() const { return breakpoints_.at(iSubIntv_+1); }

  /// Reference to the output PL, const
  const PLPoints& GetPL() const { return laPrm_.plPoints; }
  /// Reference to the output PL
  PLPoints& GetPL() { return laPrm_.plPoints; }

  /// Remainder of division, non-negative
  static double nonneg_fmod(double num, double denom) {
    return num - std::floor(num / denom) * denom;
  }


private:
  const FuncCon& con_;
  PLApproxParams& laPrm_;
  double lbx_ = -1e100, ubx_=1e100;
  int iSubIntv_ { -100 };        // chosen subinterval
  BreakpointList breakpoints_;   // subintervals
};


/// Define value of \a pi
constexpr double pi = 3.14159265358979323846;


template <class FuncCon>
void BasicPLApproximator<FuncCon>::Run() {
   ClipFuncGraphDomain();
   if (CheckDomainReturnFalseIfTrivial(GetPL())) {
     if (!InitPeriodic())
       InitNonPeriodic();
     InitSubintervalLoop();
     do {
       ApproximateSubinterval();
     } while (NextSubinterval());
     ClipGraphFinal();
     ConsiderIntegrality();
   }
 }

template <class FuncCon>
bool BasicPLApproximator<FuncCon>::InitPeriodic() {
  if (IsPeriodic()) {

    breakpoints_ = GetDefaultBreakpoints();

    auto per = GetDefaultPeriod();

    laPrm_.fUsePeriod = true;
    laPrm_.periodLength = per.ub - per.lb;
    laPrm_.periodRemainderRange =
      {breakpoints_.front(), breakpoints_.back()};
    laPrm_.periodicFactorRange = {
      std::floor( (lbx()-per.lb) / laPrm_.periodLength ),
      std::ceil( (ubx()-per.lb) / laPrm_.periodLength )
    };

    return true;
  }
  return false;
}

template <class FuncCon>
void BasicPLApproximator<FuncCon>::InitNonPeriodic() {
  // If using breakpoints, need to narrow / exclude subinterval
  // according to lbx() / ubx()
  laPrm_.fUsePeriod = false;
  auto bp_default = GetDefaultBreakpoints();
  std::set<float> bpl_set(bp_default.begin(), bp_default.end());
  auto it = bpl_set.insert(lbx()).first;
  bpl_set.erase(bpl_set.begin(), it);      // remove points before lbx()
  it = bpl_set.insert(ubx()).first;
  bpl_set.erase(++it, bpl_set.end());      // remove points after ubx()
  breakpoints_.assign(bpl_set.begin(), bpl_set.end());
  assert(breakpoints_.size() >= 2);
}

template <class FuncCon>
void BasicPLApproximator<FuncCon>::InitSubintervalLoop() {
  iSubIntv_ = 0;
  assert(breakpoints_.size() >= 2);
  auto x0 = lb_sub();                   // current left point
  auto f0 = eval(x0);
  GetPL().AddPoint(x0, f0);
}

template <class FuncCon>
bool BasicPLApproximator<FuncCon>::NextSubinterval() {
  return ++iSubIntv_+1 < (int)breakpoints_.size();
}

/// Approximate chosen subinterval
template <class FuncCon>
void BasicPLApproximator<FuncCon>::ApproximateSubinterval() {
  assert(GetPL().size());
  auto x0 = GetPL().x_.back();
  auto f0 = GetPL().y_.back();    // can be != f(x0)
  /// Simple: breakpoints on the function
  do {
    auto dx0 = ComputeInitialStepLength(x0);
    assert(x0+dx0 <= ub_sub());
    IncreaseStepWhileErrorSmallEnough(x0, f0, dx0);
    DecreaseStepWhileErrorTooBig(x0, f0, dx0);
    x0 += dx0;
    if (ub_sub()-x0 < 1e-6)        // constant param
      x0 = ub_sub();
    GetPL().AddPoint(x0, f0 = eval(x0));
  } while (x0 < ub_sub());
}

template <class FuncCon>
void BasicPLApproximator<FuncCon>::ClipGraphFinal() { }

template <class FuncCon>
void BasicPLApproximator<FuncCon>::ConsiderIntegrality() {
  /// Simple case when N breakpoints >= N integer points
  if (laPrm_.is_x_int && !laPrm_.fUsePeriod) {
    auto x0=std::ceil(laPrm_.grDomOut.lbx);
    auto xN=std::floor(laPrm_.grDomOut.ubx);
    auto N = int(xN - x0 + 1);
    if (N <= laPrm_.plPoints.size()) {
      laPrm_.plPoints.clear();
      for (int k=0; k<N; ++k)       // use double + int
        laPrm_.plPoints.AddPoint(x0+k, eval(x0+k));
    }
  }
}

template <class FuncCon>
bool BasicPLApproximator<FuncCon>::CheckDomainReturnFalseIfTrivial(
    PLPoints& result) {
  if (lbx() > ubx()+1e-6)
    MP_INFEAS(fmt::format("PLApprox {}: "
                          "empty argument domain [{}, {}]",
                          GetConTypeName(), lbx(), ubx()));
  /// Domain ~ single point
  if (lbx() > ubx()-1e-6) {
    result = { {(lbx()+ubx())/2.0}, {eval((lbx()+ubx())/2.0)} };
    return false;
  }
  return true;
}

template <class FuncCon>
double BasicPLApproximator<FuncCon>::ComputeInitialStepLength(
    double x0) {
  /// Quadratic approximation from current left point
  auto f2 = eval_2nd(x0);             // f''(x0)
  if (std::fabs(f2)<1e-100)           // cannot do,
    return (ub_sub()-x0)/100;         // return interval/100
  /// Initial step
  auto dx0 = std::sqrt(
        std::fabs(laPrm_.ubErr * 8.0 / 3.0 / f2) );
  if (x0+dx0 > ub_sub())
    dx0 = ub_sub()-x0;
  if (dx0<1e-10)
    dx0 = (ub_sub()-x0) / 100.0;
  return dx0;
}

template <class FuncCon>
void BasicPLApproximator<FuncCon>::IncreaseStepWhileErrorSmallEnough(
    double x0, double f0, double& dx0) {
  double f1;
  while ( (f1=eval(x0+dx0))==f0 ||           // same value
          0 > CompareError(x0, f0, x0+dx0, f1) ) {
    dx0 *= 1.2;
    if (x0+dx0 > ub_sub()) {
      dx0 = ub_sub()-x0;
      break;
    }
  }
}

template <class FuncCon>
void BasicPLApproximator<FuncCon>::DecreaseStepWhileErrorTooBig(
    double x0, double f0, double& dx0) {
  double f1;
  while ( (f1=eval(x0+dx0))!=f0 &&          // same value
          0 < CompareError(x0, f0, x0+dx0, f1) ) {
    dx0 *= (1.0 / 1.1);
  }
}

template <class FuncCon>
int BasicPLApproximator<FuncCon>::CompareError(
    double x0, double y0, double x1, double y1) {
//  auto err = maxErrorAbs(x0, y0, x1, y1);
  auto err = maxErrorRelAbove1(x0, y0, x1, y1);
  auto ub = laPrm_.ubErr;
  if (err<ub)
    return -1;
  if (err>ub)
    return 1;
  return 0;
}

template <class FuncCon>
double BasicPLApproximator<FuncCon>::maxErrorAbs(
    double x0, double y0, double x1, double y1) {
  MP_ASSERT_ALWAYS(x1>x0,
                   "PLApprox maxErrAbs(): degenerate segment, x0>=x1");
  auto f0 = eval(x0);
  auto f1 = eval(x1);
  auto errMax = std::max( std::fabs(f0-y0), std::fabs(f1-y1) );
  auto xMid = inverse_1st_with_check( (y1-y0)/(x1-x0) );
  auto fMid = eval(xMid);
  auto yMid = y0 + (y1-y0) * (xMid-x0) / (x1-x0);
  errMax = std::max( errMax, std::fabs(fMid-yMid) );
  return errMax;
}

template <class FuncCon>
double BasicPLApproximator<FuncCon>::maxErrorRelAbove1(
    double x0, double y0, double x1, double y1) {
  MP_ASSERT_ALWAYS(x1>x0,
                   fmt::format(
                     "PLApprox maxErrRel(): degenerate segment,"
                     " x0>=x1: {}, {}", x0, x1));
  MP_ASSERT_ALWAYS(laPrm_.ubErr>0.0,
                   "PLApprox maxErrRel(): ubErr<=0");
  std::vector< std::array<double, 2> > points;
  auto f0 = eval(x0);
  auto f1 = eval(x1);
  points.push_back({f0, y0});
  points.push_back({f1, y1});
  auto slope = (y1-y0)/(x1-x0);        // segment slope
  auto xMid = inverse_1st_with_check( slope );    // middle-value point
  points.push_back( { eval(xMid), y0 + (xMid-x0) * slope } );
  /// Assuming f' is monotone on the subinterval,
  /// check points of tilted slope assuming they are present.
  /// An alternative to evaluating f' on the endpoints
  /// would be to have a domain/subintervals for f'.
  auto fp0 = eval_1st(x0);
  auto fp1 = eval_1st(x1);
  if (fp0>fp1)
    std::swap(fp0, fp1);
  {
    auto slopeTiltedAway = slope / (1+laPrm_.ubErr);
    if (slopeTiltedAway>=fp0 && slopeTiltedAway<=fp1) {
      auto xMidUp = inverse_1st_with_check( slopeTiltedAway );
      points.push_back( { eval(xMidUp), y0 + (xMidUp-x0) * slope } );
    }
  }
  if (1.0!=laPrm_.ubErr) {
    auto slopeTiltedTo = slope / (1-laPrm_.ubErr);
    if (slopeTiltedTo>=fp0 && slopeTiltedTo<=fp1) {
      auto xMidDn = inverse_1st_with_check( slopeTiltedTo );
      points.push_back( { eval(xMidDn), y0 + (xMidDn-x0) * slope } );
    }
  }
  if (f0<1.0 && f1>1.0) {
    auto x_preim_1 = inverse_with_check(1.0);
    MP_ASSERT_ALWAYS(x0<x_preim_1 && x1>x_preim_1,
                     "PLApprox maxErrRel(): preim(1.0) outside");
    points.push_back( { 1.0, y0 + (x_preim_1-x0) * slope } );
  }
  if (f0<-1.0 && f1>-1.0) {
    auto x_preim_1 = inverse_with_check(-1.0);
    MP_ASSERT_ALWAYS(x0<x_preim_1 && x1>x_preim_1,
                     "PLApprox maxErrRel(): preim(-1.0) outside");
    points.push_back( { -1.0, y0 + (x_preim_1-x0) * slope } );
  }
  /// Check rel / abs errors in these points
  double errMax=0.0;
  for (const auto& pt: points) {
    auto f=pt[0];
    auto y=pt[1];
    auto err = (f>=-1.0 && f<=1.0) ?
          std::fabs(f-y) :                    // abs error in -1..1
          std::fabs(f-y) / std::fabs(f);      // rel error otherwise
    errMax = std::max(errMax, err);
  }
  return errMax;
}

template <class FuncCon>
double BasicPLApproximator<FuncCon>::
inverse_with_check(double y) const {
  auto x = inverse(y);
  assert(x >= lb_sub()-1e-6 && x <= ub_sub()+1e-6);
  return x;
}

template <class FuncCon>
double BasicPLApproximator<FuncCon>::
inverse_1st_with_check(double y) const {
  auto x = inverse_1st(y);
  assert(x >= lb_sub()-1e-6 && x <= ub_sub()+1e-6);
  return x;
}


////////////////// SPECIALIZED PLApproximators<> ////////////////////

/// PLApproximator<ExpConstraint>
template <>
class PLApproximator<ExpConstraint> :
    public BasicPLApproximator<ExpConstraint> {
public:
  PLApproximator(const ExpConstraint& con, PLApproxParams& p) :
    BasicPLApproximator<ExpConstraint>(con, p) { }
  FuncGraphDomain GetFuncGraphDomain() const override
  { return { -230.2585, 230.2585, 1e-100, 1e100 }; }
  bool IsMonotone() const override { return true; }
  double eval(double x) const override { return std::exp(x); }
  double inverse(double y) const override { return std::log(y); }
  double eval_1st(double x) const override { return std::exp(x); }
  double inverse_1st(double y) const override { return std::log(y); }
  double eval_2nd(double x) const override { return std::exp(x); }
};

/// Instantiate PLApproximate<ExpConstraint>
template
void PLApproximate<ExpConstraint>(
    const ExpConstraint& con, PLApproxParams& laPrm);


/// PLApproximator<LogConstraint>
template <>
class PLApproximator<LogConstraint> :
    public BasicPLApproximator<LogConstraint> {
public:
  PLApproximator(const LogConstraint& con, PLApproxParams& p) :
    BasicPLApproximator<LogConstraint>(con, p) { }
  FuncGraphDomain GetFuncGraphDomain() const override
  { return { 1e-6, 1e100, -230.2585, 230.2585 }; }
  bool IsMonotone() const override { return true; }
  double eval(double x) const override { return std::log(x); }
  double inverse(double y) const override { return std::exp(y); }
  double eval_1st(double x) const override { return 1.0/x; }
  double inverse_1st(double y) const override { return 1.0/y; }
  double eval_2nd(double x) const override { return -1.0/(x*x); }
};

/// Instantiate PLApproximate<LogConstraint>
template
void PLApproximate<LogConstraint>(
    const LogConstraint& con, PLApproxParams& laPrm);


/// PLApproximator<ExpAConstraint>
template <>
class PLApproximator<ExpAConstraint> :
    public BasicPLApproximator<ExpAConstraint> {
public:
  PLApproximator(const ExpAConstraint& con, PLApproxParams& p) :
    BasicPLApproximator<ExpAConstraint>(con, p),
    A_(GetConParams()[0]), logA_(std::log(A_)) { }
  FuncGraphDomain GetFuncGraphDomain() const override
  { return { -1e100, 1e100, 1e-100, 1e100 }; }
  bool IsMonotone() const override { return true; }
  double eval(double x) const override { return std::pow(A_, x); }
  double inverse(double y) const override { return std::log(y)/logA_; }
  double eval_1st(double x) const override { return std::pow(A_, x)*logA_; }
  double inverse_1st(double y) const override { return std::log(y/logA_)/logA_; }
  double eval_2nd(double x) const override { return std::pow(A_, x)*logA_*logA_; }


private:
  double A_;
  double logA_;
};

/// Instantiate PLApproximate<ExpAConstraint>
template
void PLApproximate<ExpAConstraint>(
    const ExpAConstraint& con, PLApproxParams& laPrm);


/// PLApproximator<LogAConstraint>
template <>
class PLApproximator<LogAConstraint> :
    public BasicPLApproximator<LogAConstraint> {
public:
  PLApproximator(const LogAConstraint& con, PLApproxParams& p) :
    BasicPLApproximator<LogAConstraint>(con, p),
    A_(GetConParams()[0]), logA_(std::log(A_)) { }
  FuncGraphDomain GetFuncGraphDomain() const override
  { return { 1e-6, 1e100, -1e100, 1e100 }; }
  bool IsMonotone() const override { return true; }
  double eval(double x) const override
  { return std::log(x) / logA_; }
  double inverse(double y) const override { return std::pow(A_, y); }
  double eval_1st(double x) const override { return 1.0/(x*logA_); }
  double inverse_1st(double y) const override { return 1.0/(y*logA_); }
  double eval_2nd(double x) const override { return -1.0/(x*x*logA_); }


private:
  double A_;
  double logA_;
};

/// Instantiate PLApproximate<LogConstraint>
template
void PLApproximate<LogAConstraint>(
    const LogAConstraint& con, PLApproxParams& laPrm);


/// PLApproximator<PowConstraint>
template <>
class PLApproximator<PowConstraint> :
    public BasicPLApproximator<PowConstraint> {
public:
  PLApproximator(const PowConstraint& con, PLApproxParams& p) :
    BasicPLApproximator<PowConstraint>(con, p) {
    InitPowDomain();
  }
  /// Constraint name. This should normally go via some
  /// printing facilities
  const char* GetConTypeName() const override {
    static std::string nm {
      "PowConstraint ^ " + std::to_string(GetConParams()[0])
    };
    return nm.c_str();
  }
  Range GetLargestAcceptedArgumentRange() const override
  { return rngAccepted_; }
  FuncGraphDomain GetFuncGraphDomain() const override
  { return fgd_; }
  /// Return false in monotone, otherwise fails in func value clipping
  bool IsMonotone() const override { return false; }
  BreakpointList GetDefaultBreakpoints() const override
  { return bpl_; }

  double eval(double x) const override
  { return std::pow(x, GetConParams()[0]); }
  double inverse(double y) const override {
    auto x00 = std::pow(std::fabs(y), 1.0 / GetConParams()[0]);
    if (lb_sub()<0.0)           // subinterval before 0
      return -x00;
    return x00;
  }
  double eval_1st(double x) const override {
    return GetConParams()[0] * std::pow(x, GetConParams()[0]-1);
  }
  double inverse_1st(double y) const override {
    auto x00 = std::pow(
          std::fabs(y / GetConParams()[0]),
        1.0 / (GetConParams()[0]-1) );
    if (lb_sub()<0.0)           // subinterval before 0
      return -x00;
    return x00;
  }
  double eval_2nd(double x) const override {
    return GetConParams()[0] *
        (GetConParams()[0]-1.0) *
        std::pow(x, GetConParams()[0]-2.0);
  }


protected:
  /// Consider cases
  void InitPowDomain() {
    auto pwr = GetConParams()[0];       // the exponent
    auto ubx = std::min(1e5, std::pow(1e5, 1/pwr));
    bool fPowInt = std::floor(pwr)==pwr;
    if (pwr<0.0) {                      // a<0
      rngAccepted_ = {0.0, 1e100};      // don' bother with x<0
      ubx = std::min(1e5, std::pow(1e-5, 1/pwr));  // smaller for HiGHS
      fgd_ = {1e-3, ubx, 0.0, 1e5};
      fMonotone = true;
      bpl_ = {{1e-3, 1e5}};
    }
    else if (!fPowInt) {                // a fractional
      rngAccepted_ = {0.0, 1e100};      // don' bother with x<0
      fgd_ = {0.0, ubx, 0.0, 1e5};
      fMonotone = true;
      bpl_ = {{0.0, 1e5}};
    } else {
      fgd_ = {-ubx, ubx, -1e5, 1e5};
    }
  }


private:
  Range rngAccepted_ {-1e100, 1e100};
  FuncGraphDomain fgd_ {-1e6, 1e6, -1e6, 1e6};
  bool fMonotone {false};
  BreakpointList bpl_{{-1e6, 0.0, 1e6}};
};

/// Instantiate PLApproximate<PowConstraint>
template
void PLApproximate<PowConstraint>(
    const PowConstraint& con, PLApproxParams& laPrm);


/// PLApproximator<SinConstraint>
template <>
class PLApproximator<SinConstraint> :
    public BasicPLApproximator<SinConstraint> {
public:
  PLApproximator(const SinConstraint& con, PLApproxParams& p) :
    BasicPLApproximator<SinConstraint>(con, p) { }
  FuncGraphDomain GetFuncGraphDomain() const override
  { return { -1e100, 1e100, -1.0, 1.0 }; }
  bool IsMonotone() const override { return false; }
  bool IsPeriodic() const override { return true; }
  Range GetDefaultPeriod() const override { return {-pi/2.0, pi*1.5}; }
  BreakpointList GetDefaultBreakpoints() const override
  { return {-pi/2, 0.0, pi/2, pi, pi*1.5}; }

  double eval(double x) const override
  { return std::sin(x); }
  /// Distinguish subinterval 0..3 (coord plane quarter):
  double inverse(double y) const override {
    assert(std::fabs(y) <= 1.0);
    if (GetSubIntvIndex()<2) {
      return std::asin(y);
    }
    return pi - std::asin(y);
  }
  double eval_1st(double x) const override { return std::cos(x); }
  /// Distinguish subinterval 0..3 (coord plane quarter):
  double inverse_1st(double y) const override {
    assert(std::fabs(y) <= 1.0);
    if (GetSubIntvIndex()<1) {
      return -std::acos(y);
    }
    if (GetSubIntvIndex()>2) {
      return pi*2 - std::acos(y);
    }
    return std::acos(y);
  }
  double eval_2nd(double x) const override { return -std::sin(x); }
};

/// Instantiate PLApproximate<LogConstraint>
template
void PLApproximate<SinConstraint>(
    const SinConstraint& con, PLApproxParams& laPrm);


/// PLApproximator<CosConstraint>
template <>
class PLApproximator<CosConstraint> :
    public BasicPLApproximator<CosConstraint> {
public:
  PLApproximator(const CosConstraint& con, PLApproxParams& p) :
    BasicPLApproximator<CosConstraint>(con, p) { }
  FuncGraphDomain GetFuncGraphDomain() const override
  { return { -1e100, 1e100, -1.0, 1.0 }; }
  bool IsMonotone() const override { return false; }
  bool IsPeriodic() const override { return true; }
  Range GetDefaultPeriod() const override { return {-pi/2.0, pi*1.5}; }
  BreakpointList GetDefaultBreakpoints() const override
  { return {-pi/2, 0.0, pi/2, pi, pi*1.5}; }

  double eval(double x) const override
  { return std::cos(x); }
  /// Distinguish subinterval 0..3 (coord plane quarter):
  double inverse(double y) const override {
    assert(std::fabs(y) <= 1.0);
    if (GetSubIntvIndex()<1) {
      return -std::acos(y);
    }
    if (GetSubIntvIndex()>2) {
      return pi*2 - std::acos(y);
    }
    return std::acos(y);
  }
  double eval_1st(double x) const override { return -std::sin(x); }
  /// Distinguish subinterval 0..3 (coord plane quarter):
  double inverse_1st(double y) const override {
    assert(std::fabs(y) <= 1.0);
    if (GetSubIntvIndex()<2) {
      return std::asin(-y);
    }
    return pi - std::asin(-y);
  }
  double eval_2nd(double x) const override { return -std::cos(x); }
};

/// Instantiate PLApproximate<CosConstraint>
template
void PLApproximate<CosConstraint>(
    const CosConstraint& con, PLApproxParams& laPrm);


/// PLApproximator<TanConstraint>
template <>
class PLApproximator<TanConstraint> :
    public BasicPLApproximator<TanConstraint> {
public:
  PLApproximator(const TanConstraint& con, PLApproxParams& p) :
    BasicPLApproximator<TanConstraint>(con, p) { }
  FuncGraphDomain GetFuncGraphDomain() const override
  { return { -1e100, 1e100, -1e100, 1e100 }; }
  bool IsMonotone() const override { return false; }
  bool IsPeriodic() const override { return true; }
  Range GetDefaultPeriod() const override { return {-pi/2.0, pi/2.0}; }
  BreakpointList GetDefaultBreakpoints() const override
  { return {-pi/2+1.4e-3, 0.0, pi/2-1.4e-3}; }

  double eval(double x) const override
  { return std::tan(x); }
  double inverse(double y) const override {
    return std::atan(y);
  }
  double eval_1st(double x) const override {
    return 1.0 / std::pow( std::cos(x), 2.0 );
  }
  /// Distinguish subinterval 0..3 (coord plane quarter):
  double inverse_1st(double y) const override {
    assert(std::fabs(y) >= 1.0);
    if (0==GetSubIntvIndex()) {                            // x<0
      return -std::acos( std::sqrt(1.0/y) );
    }
    return std::acos( std::sqrt(1.0/y) );
  }
  double eval_2nd(double x) const override {
    return 2 * std::tan(x) / std::pow( std::cos(x), 2.0 );
  }
};

/// Instantiate PLApproximate<TanConstraint>
template
void PLApproximate<TanConstraint>(
    const TanConstraint& con, PLApproxParams& laPrm);


/// PLApproximator<AsinConstraint>
template <>
class PLApproximator<AsinConstraint> :
    public BasicPLApproximator<AsinConstraint> {
public:
  PLApproximator(const AsinConstraint& con, PLApproxParams& p) :
    BasicPLApproximator<AsinConstraint>(con, p) { }
  FuncGraphDomain GetFuncGraphDomain() const override
  { return { -1.0, 1.0, -pi/2, pi/2 }; }
  BreakpointList GetDefaultBreakpoints() const override
  { return {-1.0, 0.0, 1.0}; }
  double eval(double x) const override { return std::asin(x); }
  double inverse(double y) const override { return std::sin(y); }
  double eval_1st(double x) const override { return std::pow(1-x*x, -1/2.0); }
  double inverse_1st(double y) const override {
    if (lb_sub() >= 0.0)
      return std::sqrt(1.0 - 1.0/y/y);
    return -std::sqrt(1.0 - 1.0/y/y);
  }
  double eval_2nd(double x) const override { return x*std::pow(1-x*x, -3/2.0); }
};

/// Instantiate PLApproximate<AsinConstraint>
template
void PLApproximate<AsinConstraint>(
    const AsinConstraint& con, PLApproxParams& laPrm);


/// PLApproximator<AcosConstraint>
template <>
class PLApproximator<AcosConstraint> :
    public BasicPLApproximator<AcosConstraint> {
public:
  PLApproximator(const AcosConstraint& con, PLApproxParams& p) :
    BasicPLApproximator<AcosConstraint>(con, p) { }
  FuncGraphDomain GetFuncGraphDomain() const override
  { return { -1.0, 1.0, 0.0, pi }; }
  BreakpointList GetDefaultBreakpoints() const override
  { return {-1.0, 0.0, 1.0}; }
  double eval(double x) const override { return std::acos(x); }
  double inverse(double y) const override { return std::cos(y); }
  double eval_1st(double x) const override { return -std::pow(1-x*x, -1/2.0); }
  double inverse_1st(double y) const override {
    if (lb_sub() >= 0.0)
      return std::sqrt(1.0 - 1.0/y/y);
    return -std::sqrt(1.0 - 1.0/y/y);
  }
  double eval_2nd(double x) const override { return -x*std::pow(1-x*x, -3/2.0); }
};

/// Instantiate PLApproximate<AcosConstraint>
template
void PLApproximate<AcosConstraint>(
    const AcosConstraint& con, PLApproxParams& laPrm);


/// PLApproximator<AtanConstraint>
template <>
class PLApproximator<AtanConstraint> :
    public BasicPLApproximator<AtanConstraint> {
public:
  PLApproximator(const AtanConstraint& con, PLApproxParams& p) :
    BasicPLApproximator<AtanConstraint>(con, p) { }
  FuncGraphDomain GetFuncGraphDomain() const override
  { return { -1e5, 1e5, -pi/2, pi/2 }; }      // HiGHS cannot do +-1e6, Sept 2022
  BreakpointList GetDefaultBreakpoints() const override
  { return {-1e5, 0.0, 1e5}; }
  double eval(double x) const override { return std::atan(x); }
  double inverse(double y) const override { return std::tan(y); }
  double eval_1st(double x) const override { return 1.0 / (1+x*x); }
  double inverse_1st(double y) const override {
    if (lb_sub() >= 0.0)
      return std::sqrt(1.0/y - 1.0);
    return -std::sqrt(1.0/y - 1.0);
  }
  double eval_2nd(double x) const override { return -2*x*std::pow(1+x*x, -2.0); }
};

/// Instantiate PLApproximate<AtanConstraint>
template
void PLApproximate<AtanConstraint>(
    const AtanConstraint& con, PLApproxParams& laPrm);


/// PLApproximator<SinhConstraint>
template <>
class PLApproximator<SinhConstraint> :
    public BasicPLApproximator<SinhConstraint> {
public:
  PLApproximator(const SinhConstraint& con, PLApproxParams& p) :
    BasicPLApproximator<SinhConstraint>(con, p) { }
  FuncGraphDomain GetFuncGraphDomain() const override
  { return { -14, 14, -1e100, 1e100 }; }
  BreakpointList GetDefaultBreakpoints() const override
  { return {-1e100, 0.0, 1e100}; }

  bool IsMonotone() const override { return true; }

  double eval(double x) const override
  { return std::sinh(x); }
  double inverse(double y) const override
  { return std::asinh(y); }
  double eval_1st(double x) const override { return std::cosh(x); }
  double inverse_1st(double y) const override
  { return lb_sub()>=0.0 ? std::acosh(y) : -std::acosh(y); }
  double eval_2nd(double x) const override { return std::sinh(x); }
};

/// Instantiate PLApproximate<SinhConstraint>
template
void PLApproximate<SinhConstraint>(
    const SinhConstraint& con, PLApproxParams& laPrm);


/// PLApproximator<CoshConstraint>
template <>
class PLApproximator<CoshConstraint> :
    public BasicPLApproximator<CoshConstraint> {
public:
  PLApproximator(const CoshConstraint& con, PLApproxParams& p) :
    BasicPLApproximator<CoshConstraint>(con, p) { }
  FuncGraphDomain GetFuncGraphDomain() const override
  { return { -14, 14, 1.0, 1e100 }; }
  BreakpointList GetDefaultBreakpoints() const override
  { return {-1e100, 0.0, 1e100}; }

  double eval(double x) const override
  { return std::cosh(x); }
  /// Distinguish subinterval 0..3 (coord plane quarter):
  double inverse(double y) const override
  { return lb_sub()>=0.0 ? std::acosh(y) : -std::acosh(y); }
  double eval_1st(double x) const override { return std::sinh(x); }
  double inverse_1st(double y) const override
  { return std::asinh(y); }
  double eval_2nd(double x) const override { return std::cosh(x); }
};

/// Instantiate PLApproximate<CoshConstraint>
template
void PLApproximate<CoshConstraint>(
    const CoshConstraint& con, PLApproxParams& laPrm);


/// PLApproximator<TanhConstraint>
template <>
class PLApproximator<TanhConstraint> :
    public BasicPLApproximator<TanhConstraint> {
public:
  PLApproximator(const TanhConstraint& con, PLApproxParams& p) :
    BasicPLApproximator<TanhConstraint>(con, p) { }
  FuncGraphDomain GetFuncGraphDomain() const override
  { return { -1e4, 1e4, -1, 1 }; }
  BreakpointList GetDefaultBreakpoints() const override
  { return {-1e100, 0.0, 1e100}; }

  bool IsMonotone() const override { return true; }

  double eval(double x) const override
  { return std::tanh(x); }
  double inverse(double y) const override
  { return std::atanh(y); }
  double eval_1st(double x) const override
  { return std::pow( std::cosh(x), -2.0 ); }
  /// Distinguish subinterval 0..3 (coord plane quarter):
  double inverse_1st(double y) const override {
    assert(y >= 0.0);
    if (lb_sub()>=0.0) {                            // x<0
      return std::acosh( std::sqrt(1.0/y) );
    }
    return -std::acosh( std::sqrt(1.0/y) );
  }
  double eval_2nd(double x) const override {
    return -2 * std::tanh(x) * (1.0-std::pow( std::tanh(x), 2.0 ));
  }
};

/// Instantiate PLApproximate<TanhConstraint>
template
void PLApproximate<TanhConstraint>(
    const TanhConstraint& con, PLApproxParams& laPrm);


/// PLApproximator<AsinhConstraint>
template <>
class PLApproximator<AsinhConstraint> :
    public BasicPLApproximator<AsinhConstraint> {
public:
  PLApproximator(const AsinhConstraint& con, PLApproxParams& p) :
    BasicPLApproximator<AsinhConstraint>(con, p) { }
  FuncGraphDomain GetFuncGraphDomain() const override
  { return { -1e6, 1e6, -14.5087, 14.5087 }; }
  BreakpointList GetDefaultBreakpoints() const override
  { return {-1e100, 0.0, 1e100}; }
  double eval(double x) const override { return std::asinh(x); }
  double inverse(double y) const override { return std::sinh(y); }
  double eval_1st(double x) const override { return std::pow(1+x*x, -1/2.0); }
  double inverse_1st(double y) const override {
    if (lb_sub() >= 0.0)
      return std::sqrt(1.0/y/y - 1.0);
    return -std::sqrt(1.0/y/y - 1.0);
  }
  double eval_2nd(double x) const override
  { return -x*std::pow(1+x*x, -3/2.0); }
};

/// Instantiate PLApproximate<AsinhConstraint>
template
void PLApproximate<AsinhConstraint>(
    const AsinhConstraint& con, PLApproxParams& laPrm);


/// PLApproximator<AcoshConstraint>
template <>
class PLApproximator<AcoshConstraint> :
    public BasicPLApproximator<AcoshConstraint> {
public:
  PLApproximator(const AcoshConstraint& con, PLApproxParams& p) :
    BasicPLApproximator<AcoshConstraint>(con, p) { }
  FuncGraphDomain GetFuncGraphDomain() const override
  { return { 1.0, 1e6, 0.0, 14.5087 }; }
  double eval(double x) const override { return std::acosh(x); }
  double inverse(double y) const override { return std::cosh(y); }
  double eval_1st(double x) const override { return std::pow(x*x-1, -1/2.0); }
  double inverse_1st(double y) const override
  { return std::sqrt(1.0 + 1.0/y/y); }
  double eval_2nd(double x) const override { return -x*std::pow(x*x-1, -3/2.0); }
};

/// Instantiate PLApproximate<AcoshConstraint>
template
void PLApproximate<AcoshConstraint>(
    const AcoshConstraint& con, PLApproxParams& laPrm);


/// PLApproximator<AtanhConstraint>
template <>
class PLApproximator<AtanhConstraint> :
    public BasicPLApproximator<AtanhConstraint> {
public:
  PLApproximator(const AtanhConstraint& con, PLApproxParams& p) :
    BasicPLApproximator<AtanhConstraint>(con, p) { }
  FuncGraphDomain GetFuncGraphDomain() const override
  { return { -0.999, 0.999, -1e100, 1e100 }; }  // need tighter bounds for highs
  BreakpointList GetDefaultBreakpoints() const override
  { return {-1, 0.0, 1}; }
  double eval(double x) const override { return std::atanh(x); }
  double inverse(double y) const override { return std::tanh(y); }
  double eval_1st(double x) const override { return 1.0 / (1-x*x); }
  double inverse_1st(double y) const override {
    if (lb_sub() >= 0.0)
      return std::sqrt(1.0 - 1.0/y);
    return -std::sqrt(1.0 - 1.0/y);
  }
  double eval_2nd(double x) const override { return 2*x*std::pow(1-x*x, -2.0); }
};

/// Instantiate PLApproximate<AtanhConstraint>
template
void PLApproximate<AtanhConstraint>(
    const AtanhConstraint& con, PLApproxParams& laPrm);




///////////////////////////////////////////////////////////////////
/////////////////////////// More technical code ///////////////////
void FuncGraphDomain::intersect(const FuncGraphDomain &grDom) {
  lbx = std::max(lbx, grDom.lbx);
  ubx = std::min(ubx, grDom.ubx);
  lby = std::max(lby, grDom.lby);
  uby = std::min(uby, grDom.uby);
}


/// Basic version, assume fn monotone
template <class FuncCon>
void BasicPLApproximator<FuncCon>::ClipWithFunctionValues(
    FuncGraphDomain &grDom) const {
  auto imlbx = eval(grDom.lbx);
  auto imubx = eval(grDom.ubx);
  auto prelby = inverse(grDom.lby);
  auto preuby = inverse(grDom.uby);
  grDom.lbx = std::max( grDom.lbx, std::min(prelby, preuby) );
  grDom.ubx = std::min( grDom.ubx, std::max(prelby, preuby) );
  grDom.lby = std::max( grDom.lby, std::min(imlbx, imubx) );
  grDom.uby = std::min( grDom.uby, std::max(imlbx, imubx) );
}


PLPoints::PLPoints(const PLSlopes &pls) {
  constexpr auto eps = 1.0;           // works?
  const auto& bp = pls.GetBP();
  const auto& sl = pls.GetSlopes();
  const auto nsl = sl.size();
  const auto X0 = pls.GetX0(), Y0 = pls.GetY0();
  x_.resize(nsl+1);
  y_.resize(nsl+1);
  /// Copy and add dummy breakpoints on both ends
  std::copy(bp.begin(), bp.end(), x_.begin()+1);
  x_[0] = x_[1] - eps;
  x_[nsl] = x_[nsl-1] + eps;
  y_[0] = 0.0;                        // initialize leftmost point
  /// Lift the line by this
  double deltaH {};
  if (x_[0] > X0)                     // if left x > reference point
    deltaH = sl[0] * (x_[0]-X0) + Y0;
  for (size_t i = 0; i < nsl; ++i) {
    assert( x_[i+1] > x_[i] );
    y_[i+1] = y_[i] + sl[i] * (x_[i+1]-x_[i]);
    if (x_[i]<=X0 && (x_[i+1]>=X0 || i==nsl-1))
      deltaH = Y0 - (y_[i] + sl[i]*(-x_[i]-X0));
  }
  for (size_t i = 0; i <= nsl; ++i) {
    y_[i] += deltaH;
  }
}

}  // namepsace mp
