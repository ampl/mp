#ifndef CPSATCOMMON_H
#define CPSATCOMMON_H

#include <string>
#include <vector>
#include <memory>
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <stdexcept>

#include "mp/backend-to-model-api.h"

#include "ortools/sat/cp_model.h"
#include "ortools/sat/cp_model.pb.h"
#include "ortools/sat/sat_parameters.pb.h"
#include "ortools/sat/model.h"

#include "mp/format.h"

namespace mp {

/// Tolerance for considering a coefficient or a value integer.
/// Relative for large values, see IsIntegerValue().
constexpr double kCpsatIntTol = 1e-9;

/// Is \a v within \a abs_tol of an integer?
inline bool IsInteger(double v, double abs_tol = kCpsatIntTol) {
  return std::fabs(v - std::round(v)) <= abs_tol;
}

/// Is \a v an integer value, with a tolerance relative to its magnitude
/// (at least kCpsatIntTol)? Used for coefficients and bounds.
inline bool IsIntegerValue(double v) {
  return IsInteger(v, kCpsatIntTol * std::max(1.0, std::fabs(v)));
}

/// Convert an integer value to int64_t.
/// @param what: description of the value, for the error message.
/// @throw std::runtime_error if \a v is not an integer value
/// (IsIntegerValue()) or is out of range.
inline int64_t ToInt64(double v, const char* what) {
  if (!IsIntegerValue(v))
    throw std::runtime_error(fmt::format(
        "CP-SAT requires integer {}, got {}", what, v));
  const double r = std::round(v);
  if (std::fabs(r) >= 4.0e18)
    throw std::runtime_error(fmt::format(
        "CP-SAT: {} out of range: {}", what, v));
  return (int64_t)r;
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
