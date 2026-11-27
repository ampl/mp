#ifndef CPSATCOMMON_H
#define CPSATCOMMON_H

#include <string>
#include <vector>
#include <memory>
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>

#include "mp/backend-to-model-api.h"
#include "mp/flat/num_tolerances.h"

#include "ortools/sat/cp_model.h"
#include "ortools/sat/cp_model.pb.h"
#include "ortools/sat/sat_parameters.pb.h"
#include "ortools/sat/model.h"

#include "mp/format.h"

namespace mp {

/// Absolute tolerance for taking a value as the nearest integer:
/// coefficients, bounds and fixed values computed in floating point
/// (0.1*3*10 is 3.0000000000000004) are meant to be integers.
constexpr double kCpsatIntTol = 1e-9;

/// For large magnitudes, the rounding error of a few operations, in units
/// in the last place (ulp) of the value. 4 ulp keeps 0.07*1e15 and rejects
/// real fractions such as 1e12 + 0.001, 1e13 + 1/64 or 2^48 + 0.5 (8 ulp).
constexpr double kCpsatIntUlps = 4.0;

/// Upper limit of the tolerance: a value ending in .5 is never taken as an
/// integer, and a bound is never moved outward by 0.25 or more. From 2^51
/// (spacing 0.5) only exact integers pass; from 2^52 every double is one.
constexpr double kCpsatIntMaxTol = 0.25;

/// Is \a v within \a abs_tol of an integer?
inline bool IsInteger(double v, double abs_tol = kCpsatIntTol) {
  return std::fabs(v - std::round(v)) <= abs_tol;
}

/// Integrality tolerance for a finite \a v:
/// max(1e-9, min(4 ulp(v), 0.25)).
inline double IntegerTolerance(double v) {
  const double a = std::fabs(v);
  const double ulp =
      std::nextafter(a, std::numeric_limits<double>::infinity()) - a;
  return std::max(kCpsatIntTol,
                  std::min(kCpsatIntUlps * ulp, kCpsatIntMaxTol));
}

/// Is \a v an integer value, up to IntegerTolerance()? The one rule for
/// coefficients, row and variable bounds and fixed values.
inline bool IsIntegerValue(double v) {
  return std::isfinite(v) && IsInteger(v, IntegerTolerance(v));
}

/// Convert an integer value to int64_t; a value within the tolerance
/// of an integer is taken as that integer.
/// @param what: description of the value, for the error message.
/// @throw std::runtime_error if \a v is not an integer value
/// (IsIntegerValue()) or is out of range. The message shows all digits
/// of \a v, so that a value such as 2.9999999 does not print as "3".
inline int64_t ToInt64(double v, const char* what) {
  if (!IsIntegerValue(v))
    throw std::runtime_error(fmt::format(
        "CP-SAT requires integer {}, got {:.17g}", what, v));
  const double r = std::round(v);
  if (std::fabs(r) >= 4.0e18)
    throw std::runtime_error(fmt::format(
        "CP-SAT: {} out of range: {:.17g}", what, v));
  return static_cast<int64_t>(r);
}

/// Tolerance for taking a bound \a b of an integer-valued quantity
/// as the nearest integer: IntegerTolerance(b), or the converter's
/// integrality tolerance \a inteps (cvt:pre:inteps) if larger,
/// below 0.5. Negative \a inteps: IntegerTolerance(b) only.
inline double BoundTolerance(double b, double inteps) {
  const double t = IntegerTolerance(b);
  return inteps > t ? std::min(inteps, NumericTolerances::MaxIntEps()) : t;
}

/// Round a bound of an integer-valued quantity to an integer value:
/// to the nearest one if within BoundTolerance(), otherwise inward
/// (a lower bound up, an upper bound down). Infinite bounds are kept.
/// @param is_lower: whether \a b is a lower bound.
/// @param inteps: cvt:pre:inteps, negative: not used.
inline double RoundBoundInward(double b, bool is_lower,
                               double inteps = -1.0) {
  if (std::isfinite(b) && IsInteger(b, BoundTolerance(b, inteps)))
    return std::round(b);
  return is_lower ? std::ceil(b) : std::floor(b);
}

/// Convert an integer-valued bound to int64_t, clamping to the int64 range.
/// A bound at or beyond the range (also an infinite one) cannot bind,
/// or cannot be satisfied, either way: CP-SAT expressions stay within int64.
/// Not for NaN.
inline int64_t ClampToInt64(double v) {
  constexpr double kLimit = 9223372036854775808.0;   // 2^63, exact in double
  if (v >= kLimit)
    return std::numeric_limits<int64_t>::max();
  if (v <= -kLimit)
    return std::numeric_limits<int64_t>::min();
  return static_cast<int64_t>(v);
}

/// Default magnitude used for infinite variable bounds.
/// CP-SAT needs finite int64 domains for its overflow checks to succeed.
constexpr double kCpsatDefaultVarBound = 1e12;

/// Largest magnitude of a variable bound given to CP-SAT. Its model
/// checker requires domains within [-kint64max/2, kint64max/2] (about
/// 4.6e18); 4e18 is below that and exact in double.
constexpr double kCpsatMaxVarBound = 4e18;

/// Data shared by the Backend and the ModelAPI
struct CpsatData {
  /// The model under construction
  operations_research::sat::CpModelBuilder builder;
  /// One CP-SAT variable per AMPL (flat) variable
  std::vector<operations_research::sat::IntVar> vars;
  /// Bounds of the variables, as given to CP-SAT
  std::vector<int64_t> var_lb, var_ub;
  /// Value of a variable fixed at a non-integer value (MP's constants,
  /// e.g., a fractional objective constant), else NaN. Such a variable
  /// enters only the objective's offset; CP-SAT has a placeholder.
  std::vector<double> var_const;
  /// Solver parameters
  operations_research::sat::SatParameters params;
  /// Replaces infinite variable bounds (cvt:mip_max_bound)
  double var_bound = kCpsatDefaultVarBound;
  /// Integrality tolerance of the converter (cvt:pre:inteps) for
  /// row and variable bounds and numberof values; negative: not used
  double inteps = -1.0;
  /// A variable with an infinite bound replaced by cvt:mip_max_bound
  struct VarBoundSet {
    int var;               ///< index of the variable
    std::string name;      ///< its name, if known
    bool lb, ub;           ///< which bounds were set
  };
  /// The variables whose bounds cvt:mip_max_bound set, in order
  std::vector<VarBoundSet> varbound_set;
  /// Result of the last solve
  operations_research::sat::CpSolverResponse response;
  /// Set only by the interrupt (signal) handler. An atomic store is
  /// async-signal-safe, as MP requires of handlers. While solving, a watcher
  /// thread stops the CP-SAT search when it is set (see cpsatbackend.cc).
  /// Reset when the handler is installed, see CpsatBackend::SetInterrupter().
  std::atomic<bool> interrupted{false};
};


/// Information shared by both
/// `CpsatBackend` and `CpsatModelAPI`
struct CpsatCommonInfo {
  CpsatData* data() const { return data_; }
  void set_data(CpsatData* d) { data_ = d; }
private:
  CpsatData* data_ = nullptr;
};


/// Common API for Cpsat classes
class CpsatCommon :
    public Backend2ModelAPIConnector<CpsatCommonInfo> {
public:
  /// These methods access CP-SAT parameters (fields of SatParameters).
  /// Used by AddSolverOption()
  void GetSolverOption(const char* key, int& value) const;
  void SetSolverOption(const char* key, int value);
  void GetSolverOption(const char* key, double& value) const;
  void SetSolverOption(const char* key, double value);
  void GetSolverOption(const char* key, std::string& value) const;
  void SetSolverOption(const char* key, const std::string& value);

  static constexpr double Infinity() { return INFINITY; }
  static constexpr double MinusInfinity() { return -INFINITY; }

protected:
  int NumVars() const;
  int NumObjs() const { return 1; }

  /// The CP-SAT domain [lb, ub] of a variable with bounds \a lb, \a ub.
  /// Finite bounds are kept, rounded inward (lb up, ub down; within
  /// cvt:pre:inteps to the nearest integer, RoundBoundInward()); an infinite
  /// bound becomes -+var_bound, or var_bound beyond the other bound when
  /// that bound lies beyond var_bound itself. The result is clamped to
  /// +-kCpsatMaxVarBound, the domain limit of CP-SAT.
  /// @throw std::runtime_error if a finite bound excludes every
  /// value within that limit (lb above it, or ub below its negative).
  std::pair<int64_t, int64_t> VarDomainToInt64(double lb, double ub) const;
};

} // namespace mp

#endif // CPSATCOMMON_H
