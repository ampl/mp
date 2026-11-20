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

#include "ortools/sat/cp_model.h"
#include "ortools/sat/cp_model.pb.h"
#include "ortools/sat/sat_parameters.pb.h"
#include "ortools/sat/model.h"

#include "mp/format.h"

namespace mp {

/// Absolute tolerance for taking a value as the nearest integer:
/// coefficients, bounds, fixed values and exponents computed in floating
/// point (0.1*3*10 is 3.0000000000000004) are meant to be integers.
constexpr double kCpsatIntTol = 1e-9;

/// The same for large magnitudes, where rounding errors exceed kCpsatIntTol:
/// a few units of floating-point rounding, relative to the value. Kept
/// that small so that a real fraction (1000000000.5) is not taken as noise.
constexpr double kCpsatIntRelTol = 8 * std::numeric_limits<double>::epsilon();

/// Is \a v within \a abs_tol of an integer?
inline bool IsInteger(double v, double abs_tol = kCpsatIntTol) {
  return std::fabs(v - std::round(v)) <= abs_tol;
}

/// Is \a v an integer value, up to the tolerance (kCpsatIntTol, or
/// kCpsatIntRelTol relative to \a v when that is larger)? The one rule for
/// coefficients, row and variable bounds, fixed values and exponents.
inline bool IsIntegerValue(double v) {
  return std::isfinite(v) && IsInteger(
        v, std::max(kCpsatIntTol, kCpsatIntRelTol * std::fabs(v)));
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

/// Round a bound of an integer-valued quantity to an integer value:
/// to the nearest one if within the tolerance, otherwise inward
/// (a lower bound up, an upper bound down). Infinite bounds are kept.
/// @param is_lower: whether \a b is a lower bound.
inline double RoundBoundInward(double b, bool is_lower) {
  if (IsIntegerValue(b))
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
  /// Magnitude for infinite bounds
  double var_bound = kCpsatDefaultVarBound;
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

  /// Convert a lower / upper bound to an int64, rounding inwards
  /// and clamping to +-var_bound
  int64_t LowerBoundToInt64(double lb) const;
  int64_t UpperBoundToInt64(double ub) const;
};

} // namespace mp

#endif // CPSATCOMMON_H
