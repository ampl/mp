/*
 A mathematical optimization solver.

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

#include "mp/flat/constr_std.h"
#include "mp/flat/redef/MIP/core/lin_approx_core.h"

namespace mp {

/// Declare function-specific approximator
template <class Con>
class PLApproximator;

/// Define generic approximation call
template <class FuncCon>
PLPoints PLApproximate(const FuncCon& con,
                             LinApproxParams& laPrm) {
  PLApproximator<FuncCon> pla{con, laPrm};
  return pla.Run();
}


/// Interface class for univariate function's
/// piecewise-linear smooth approximator:
/// domain, periodicity, derivatives.
/// The basic version assumes monotonicity
/// of the function and the 1st derivative.
template <class FuncCon>
class BasicPLApproximator {
public:
  /// Constructor
  BasicPLApproximator(const FuncCon& con, LinApproxParams& laPrm) :
    con_(con), laPrm_(laPrm) { }
  /// Virtual destructor
  virtual ~BasicPLApproximator() { }
  /// Constraint type name
  virtual const char* GetConTypeName() const
  { return con_.GetTypeName(); }
  /// Run the approximation
  PLPoints Run();


protected:
  /// Ample, but realistic, function graph domain
  /// (values should not overflow when computing (pre-)image)
  virtual FuncGraphDomain GetFuncGraphDomain() const = 0;
  /// Function graph domain: clip to the user-provided domain
  virtual void ClipFuncGraphDomain() {
    laPrm_.grDom.intersect(GetFuncGraphDomain());
    ClipWithFunctionValues(laPrm_.grDom);
    lbx_ = laPrm_.grDom.lbx;
    ubx_ = laPrm_.grDom.ubx;
  }
  /// Check domain, throw if infeas, or return trivial PL
  /// @throw if infeasible
  /// @return false iff trivial
  /// @param result: the trivial PL
  bool CheckDomainReturnFalseIfTrivial(PLPoints& result);
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

  /// Evaluate function at \a x
  virtual double eval(double x) const = 0;
  /// Evaluate inv(f) at \a y
  virtual double inverse(double y) const = 0;
  /// Evaluate f' at \a x
  virtual double eval_1st(double x) const = 0;
  /// Evaluate inv(f') at \a y
  virtual double inverse_1st(double y) const = 0;
  /// Evaluate f'' at \a x
  virtual double eval_2nd(double x) const = 0;

  /// Maximal absolute error (along Y axis) to the linear segment
  /// (x0, y0) - (x1, y1).
  /// Standard implementation checks in the ends
  /// and in the middle value point
  virtual double maxErrorAbs(
      double x0, double y0, double x1, double y1);
  /// Clip the graph domain, considering the actual function.
  virtual void ClipWithFunctionValues(FuncGraphDomain& grDom) const;

  /// Function interval to approximate: lb
  double lbx() const { return lbx_; }
  /// Function interval to approximate: ub
  double ubx() const { return ubx_; }

private:
  const FuncCon& con_;
  LinApproxParams& laPrm_;
  double lbx_ = -1e100, ubx_=1e100;
};

/// Approximation execution
template <class FuncCon>
PLPoints BasicPLApproximator<FuncCon>::Run() {
  PLPoints result;
  ClipFuncGraphDomain();
  if (CheckDomainReturnFalseIfTrivial(result)) {
    auto x0 = lbx();                      // current left point
    auto f0 = eval(x0);
    result.AddPoint(x0, f0);
    /// Simple: breakpoints on the function
    do {
      auto dx0 = ComputeInitialStepLength(x0);
      assert(x0+dx0 <= ubx());
      IncreaseStepWhileErrorSmallEnough(x0, f0, dx0);
      DecreaseStepWhileErrorTooBig(x0, f0, dx0);
      x0 += dx0;
      result.AddPoint(x0, f0 = eval(x0));
    } while (x0 < ubx());
  }
  return result;
}

template <class FuncCon>
bool BasicPLApproximator<FuncCon>::CheckDomainReturnFalseIfTrivial(
    PLPoints& result) {
  if (lbx() > ubx()+1e-6)
    MP_INFEAS(fmt::format("LinApprox {}: "
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
  MP_ASSERT_ALWAYS(std::fabs(f2)>1e-100,
                   fmt::format("LinApprox {}: f''(x0)==0",
                               GetConTypeName()));
  /// Initial step
  auto dx0 = std::sqrt(
        std::fabs(laPrm_.ubErrAbs * 8.0 / 3.0 / f2) );
  if (x0+dx0 > ubx())
    dx0 = ubx()-x0;
  return dx0;
}

template <class FuncCon>
void BasicPLApproximator<FuncCon>::IncreaseStepWhileErrorSmallEnough(
    double x0, double f0, double& dx0) {
  while ( maxErrorAbs(x0, f0, x0+dx0, eval(x0+dx0))
          < laPrm_.ubErrAbs ) {
    dx0 *= 1.2;
    if (x0+dx0 > ubx()) {
      dx0 = ubx()-x0;
      break;
    }
  }
}

template <class FuncCon>
void BasicPLApproximator<FuncCon>::DecreaseStepWhileErrorTooBig(
    double x0, double f0, double& dx0) {
  while ( maxErrorAbs(x0, f0, x0+dx0, eval(x0+dx0))
          > laPrm_.ubErrAbs ) {
    dx0 *= (1.0 / 1.1);
  }
}

template <class FuncCon>
double BasicPLApproximator<FuncCon>::maxErrorAbs(
    double x0, double y0, double x1, double y1) {
  auto f0 = eval(x0);
  auto f1 = eval(x1);
  auto errMax = std::max( std::fabs(f0-y0), std::fabs(f1-y1) );
  MP_ASSERT_ALWAYS(x1>x0,
                   "LinApprox maxErrAbs(): degenerate segment, x0>=x1");
  auto xMid = inverse_1st( (y1-y0)/(x1-x0) );
  auto fMid = eval(xMid);
  auto yMid = y0 + (y1-y0) * (xMid-x0) / (x1-x0);
  errMax = std::max( errMax, std::fabs(fMid-yMid) );
  return errMax;
}


////////////////// SPECIALIZED PLApproximators<> ////////////////////

/// PLApproximator<ExpConstraint>
template <>
class PLApproximator<ExpConstraint> :
    public BasicPLApproximator<ExpConstraint> {
public:
  PLApproximator(const ExpConstraint& con, LinApproxParams& p) :
    BasicPLApproximator<ExpConstraint>(con, p) { }
  FuncGraphDomain GetFuncGraphDomain() const override
  { return { -230.2585, 230.2585, 1e-100, 1e100 }; }
  double eval(double x) const override { return std::exp(x); }
  double inverse(double y) const override { return std::log(y); }
  double eval_1st(double x) const override { return std::exp(x); }
  double inverse_1st(double y) const override { return std::log(y); }
  double eval_2nd(double x) const override { return std::exp(x); }
};

/// Instantiate PLApproximate<ExpConstraint>
template
PLPoints PLApproximate<ExpConstraint>(
    const ExpConstraint& con, LinApproxParams& laPrm);


/// PLApproximator<LogConstraint>
template <>
class PLApproximator<LogConstraint> :
    public BasicPLApproximator<LogConstraint> {
public:
  PLApproximator(const LogConstraint& con, LinApproxParams& p) :
    BasicPLApproximator<LogConstraint>(con, p) { }
  FuncGraphDomain GetFuncGraphDomain() const override
  { return { 1e-5, 1e100, -230.2585, 230.2585 }; }
  double eval(double x) const override { return std::log(x); }
  double inverse(double y) const override { return std::exp(y); }
  double eval_1st(double x) const override { return 1.0/x; }
  double inverse_1st(double y) const override { return 1.0/y; }
  double eval_2nd(double x) const override { return -1.0/(x*x); }
};

/// Instantiate PLApproximate<LogConstraint>
template
PLPoints PLApproximate<LogConstraint>(
    const LogConstraint& con, LinApproxParams& laPrm);


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
