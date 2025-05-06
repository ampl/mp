#ifndef CUOPTLPCOMMON_H
#define CUOPTLPCOMMON_H

#include <string>

#include "mp/backend-to-model-api.h"

//extern "C" {
// TODO Typically import here the solver's C API headers

#include <cuopt/linear_programming/cuopt_c.h>

// Here we use a cplus plus stub instead
   #include "cuoptlp-solvermodel.h"
//}

#include "mp/format.h"

namespace mp {

struct ProblemData {
  cuOptOptimizationProblem problem;
  cuOptSolverSettings settings;
  cuOptSolution solution;


  std::vector<cuopt_float_t> lower_bounds;
  std::vector<cuopt_float_t> upper_bounds;
  std::vector<char> variable_types;

  std::vector<cuopt_float_t> constraint_matrix_coefficients;
  std::vector<cuopt_int_t> constraint_matrix_row_offsets;
  std::vector<cuopt_int_t> constraint_matrix_column_indices;

  cuopt_int_t num_constraints;
  cuopt_int_t num_variables;
  cuopt_int_t nnz;

  std::vector<char> constraint_sense;
  std::vector<cuopt_float_t> rhs;

  std::vector<cuopt_float_t> objective_coefficients;
  cuopt_int_t objective_sense;
  cuopt_float_t objective_offset;
};


/// Information shared by both
/// `CuoptlpBackend` and `CuoptlpModelAPI`
struct CuoptlpCommonInfo {

  // TODO provide accessors to the solver's in-memory model/environment
  //cuoptlp_env* env() const { return env_; }
  ProblemData* lp() const { return lp_; }

  // TODO provide accessors to the solver's in-memory model/environment
  //void set_env(cuoptlp_env* e) { env_ = e; }
  void set_lp(ProblemData* lp) { lp_ = lp; }

  ProblemData* lp_ = nullptr;

private:
  // TODO provide accessors to the solver's in-memory model/environment
  //cuoptlp_env*      env_ = NULL;
};


/// Common API for Cuoptlp classes
class CuoptlpCommon :
    public Backend2ModelAPIConnector<CuoptlpCommonInfo> {
public:
  /// These methods access Cuoptlp options. Used by AddSolverOption()
  void GetSolverOption(const char* key, int& value) const;
  void SetSolverOption(const char* key, int value);
  void GetSolverOption(const char* key, double& value) const;
  void SetSolverOption(const char* key, double value);
  void GetSolverOption(const char* key, std::string& value) const;
  void SetSolverOption(const char* key, const std::string& value);

  /// TODO Typically solvers define their own infinity; use them here
  static constexpr double Infinity() { return INFINITY;  }
  static constexpr double MinusInfinity() { return -INFINITY; }

protected:

  int NumLinCons() const;
  int NumVars() const;
  int NumObjs() const;
  int NumQPCons() const;


protected:
  // TODO if desirable, provide function to create the solver's environment
  // with own license
  // int (*createEnv) (solver_env**) = nullptr;

};


/// Convenience macro
// TODO This macro is useful to automatically throw an error if a function in the
// solver API does not return a valid errorcode. In this mock driver, we define it
// ourselves, normally this constant would be defined in the solver's API.
#define CUOPT_CCALL( call ) do { if (int e = (call) != CUOPT_SUCCESS) \
  throw std::runtime_error( \
    fmt::format("  Call failed: '{}' with code {}", #call, e )); } while (0)

} // namespace mp

#endif // CUOPTLPCOMMON_H
