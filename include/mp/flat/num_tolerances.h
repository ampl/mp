#ifndef NUM_TOLERANCES_H
#define NUM_TOLERANCES_H

/**
 * Tolerances of the flat converter's own numeric decisions:
 * integrality (cvt:pre:inteps), floating-point margins of computed
 * values (cvt:pre:ulps) and comparisons (absolute cvt:pre:feastol).
 * The converter sets them once, after option parsing; the model
 * (rounding of integer bounds), the preprocessors and the MIP
 * reformulations read them from one place.
 */

#include <cmath>
#include <algorithm>

namespace mp {

/// Tolerances of the converter's numeric decisions
class NumericTolerances {
public:
  /// Largest integrality tolerance: below 0.5, so that a value
  /// with fractional part 0.5 is never taken as an integer
  static constexpr double MaxIntEps() { return 0.49999999999999994; }
  /// Default integrality tolerance (cvt:pre:inteps)
  static constexpr double DefaultIntEps() { return 1e-6; }
  /// Default floating-point margin in ulp (cvt:pre:ulps)
  static constexpr int DefaultUlps() { return 16; }
  /// Default comparison tolerance (cvt:pre:feastol)
  static constexpr double DefaultCmpTol() { return 1e-6; }

  /// Set all values.
  /// @param inteps: cvt:pre:inteps; negative: only the margin \a ulps.
  /// @param ulps: cvt:pre:ulps; 0: no floating-point margin.
  /// @param cmptol: absolute comparison tolerance (cvt:pre:feastol).
  void Set(double inteps, int ulps, double cmptol) {
    inteps_ = inteps;
    ulps_ = ulps;
    cmptol_ = cmptol;
  }

  /// Integrality tolerance as set (negative: exact up to FPMargin())
  double inteps() const { return inteps_; }
  /// Floating-point margin in units in the last place
  int ulps() const { return ulps_; }
  /// Absolute comparison tolerance
  double cmp_tol() const { return cmptol_; }

  /// The spacing of doubles at |b|: the size of a rounding error
  static double Ulp(double b) {
    const double a = std::fabs(b);
    return std::nextafter(a, INFINITY) - a;
  }

  /// Floating-point margin of a value computed by the converter
  /// whose terms have magnitude up to \a a:
  /// cvt:pre:ulps units in the last place of \a a.
  /// 0 for infinite \a a.
  double FPMargin(double a) const {
    return std::isfinite(a) ? ulps_ * Ulp(a) : 0.0;
  }

  /// Whether a value near an integer can be taken as that integer
  /// (false with cvt:pre:inteps < 0 and cvt:pre:ulps = 0:
  /// integrality is then exact and nothing is rounded)
  bool CanRoundToInteger() const { return inteps_ >= 0.0 || ulps_ > 0; }

  /// Integrality tolerance at \a b: max(inteps, FPMargin(b)),
  /// only FPMargin(b) with inteps < 0, at most MaxIntEps()
  double IntTol(double b) const {
    return std::min(MaxIntEps(), std::max(inteps_, FPMargin(b)));
  }

  /// Is finite \a v within IntTol(v) of an integer?
  /// Exact integers always are.
  bool IsNearInteger(double v) const {
    if (!std::isfinite(v))
      return false;
    if (std::trunc(v) == v)
      return true;
    return CanRoundToInteger()
        && std::fabs(v - std::round(v)) <= IntTol(v);
  }

  /// The integer near \a v if IsNearInteger(v), otherwise \a v
  double SnapToInteger(double v) const {
    return IsNearInteger(v) ? std::round(v) : v;
  }

  /// Margin d by which the computed bounds of a comparison's body may
  /// cross its right-hand side \a c without deciding the comparison:
  /// d = max(cmp_tol(), FPMargin(c)). cmp_tol() is the comparison
  /// tolerance in the units of the body; FPMargin(c), from the compared
  /// constant only, absorbs rounding errors of computed bounds near c,
  /// e.g., 1.9999999999999998 for 2.
  /// @return d, or a negative value if the comparison must not be
  ///   decided: FPMargin(c) exceeds a positive cmp_tol(), so that
  ///   rounding errors at the magnitude of c exceed the tolerance
  ///   (left to the solver)
  double DecisionMargin(double c) const {
    const double m = FPMargin(c);
    if (cmptol_ > 0.0 && m > cmptol_)
      return -1.0;
    return std::max(cmptol_, m);
  }

private:
  double inteps_ {DefaultIntEps()};
  int ulps_ {DefaultUlps()};
  double cmptol_ {DefaultCmpTol()};
};

} // namespace mp

#endif // NUM_TOLERANCES_H
