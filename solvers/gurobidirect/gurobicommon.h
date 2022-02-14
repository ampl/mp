#ifndef GUROBICOMMON_H
#define GUROBICOMMON_H

#include <string>
#include <vector>

/// Common stuff for GurobiBackend and GurobiModelAPI
extern "C" {
  #include "gurobi_c.h"
}

#include "mp/env.h"
#include "mp/arrayref.h"
#include "mp/format.h"
#include "mp/error.h"
#include "mp/flat/model_api_base.h"
#include "mp/model-mgr-base.h"
#include "mp/presolve-base.h"

namespace mp {

/// Common ancestor for GurobiBackend and GurobiModelAPI
class GurobiCommon : public BasicFlatModelAPI {
public:
  ////////////////////////////////////////////////////////////
  //////////////////////// Metadata //////////////////////////
  ////////////////////////////////////////////////////////////

  /// +inf
  static constexpr double Infinity() { return GRB_INFINITY; }
  /// -inf
  static constexpr double MinusInfinity() { return -GRB_INFINITY; }

  /// Gurobi separates constraint classes
  int NumLinCons() const;
  int NumQPCons() const;
  int NumSOSCons() const;
  int NumGenCons() const;
  int NumVars() const;
  int NumObjs() const;
  int ModelSense() const;

  /// Public option API.
  /// These methods access Gurobi options. Used by AddSolverOption()
  void GetSolverOption(const char* key, int& value) const;
  void SetSolverOption(const char* key, int value);
  void GetSolverOption(const char* key, double& value) const;
  void SetSolverOption(const char* key, double value);
  void GetSolverOption(const char* key, std::string& value) const;
  void SetSolverOption(const char* key, const std::string& value);

  /// REMEMBER Gurobi does not update attributes before calling optimize() etc
  /// Scalar attributes. If (flag), set *flag <-> success,
  /// otherwise fail on error
  int GrbGetIntAttr(const char* attr_id, bool *flag=nullptr) const;
  /// If (flag), set *flag <-> success, otherwise fail on error
  double GrbGetDblAttr(const char* attr_id, bool *flag=nullptr) const;
  /// Vector attributes. Return empty vector on failure
  std::vector<int> GrbGetIntAttrArray(const char* attr_id,
    std::size_t size, std::size_t offset=0) const;
  std::vector<double> GrbGetDblAttrArray(const char* attr_id,
    std::size_t size, std::size_t offset=0) const;
  std::vector<int> GrbGetIntAttrArray(GRBmodel* mdl, const char* attr_id,
    std::size_t size, std::size_t offset=0) const;
  std::vector<double> GrbGetDblAttrArray(GRBmodel* mdl, const char* attr_id,
    std::size_t size, std::size_t offset=0) const;

  /// varcon: 0 - vars, 1 - constraints
  std::vector<double> GrbGetDblAttrArray_VarCon(
      const char* attr, int varcon) const;
  /// varcon: 0 - vars, 1 - constraints
  std::vector<double> GrbGetDblAttrArray_VarCon(GRBmodel* mdl,
      const char* attr, int varcon) const;

  /// Set attributes.
  /// Return false on failure
  void GrbSetIntAttr(const char* attr_id, int val);
  void GrbSetDblAttr(const char* attr_id, double val);
  ///  Silently ignore empty vector arguments.
  void GrbSetIntAttrArray(const char* attr_id,
                               ArrayRef<int> values, std::size_t start=0);
  void GrbSetDblAttrArray(const char* attr_id,
                               ArrayRef<double> values, std::size_t start=0);
  ///  Silently ignore empty vector arguments.
  void GrbSetIntAttrList(const char* attr_id,
                         const std::vector<int>& idx, const std::vector<int>& val);
  void GrbSetDblAttrList(const char* attr_id,
                         const std::vector<int>& idx, const std::vector<double>& val);

  //////////// Wrappers for Get/SetSolverOption(). Assume model_ is set
  int GrbGetIntParam(const char* key) const;
  double GrbGetDblParam(const char* key) const;
  std::string GrbGetStrParam(const char* key) const;
  void GrbSetIntParam(const char* key, int value);
  void GrbSetDblParam(const char* key, double value);
  void GrbSetStrParam(const char* key, const std::string& value);

  GurobiCommon *other_gurobi() { return other_; }
  void set_other_gurobi(GurobiCommon* o) { other_ = o; }


protected:
  void copy_handlers_from_other_gurobi();
  void copy_handlers_to_other_gurobi();

  GRBenv *env() const { return env_; }
  GRBenv *&env_ref() { return env_; }
  void set_env(GRBenv* e) { env_ = e; }
  GRBmodel *model() const { return model_; }
  GRBmodel *&model_ref() { return model_; }
  void set_model(GRBmodel* m) { model_ = m; }

private:
  GRBenv *env_ = nullptr;
  GRBmodel *model_ = nullptr;
  GurobiCommon *other_ = nullptr;
};


/// Data from creating a GurobiModelAPI
struct GurobiModelAPIData {
  std::unique_ptr<BasicModelManager> pmm_;
  pre::BasicPresolver* pre_;
};

/// Create Gurobi Model Manager
/// @param gc the Gurobi Backend
GurobiModelAPIData
CreateGurobiModelMgr(GurobiCommon& gc, Env& e);

} // namespace mp

#define GRB_CALL_MSG( call, msg ) do { if (int e=call) MP_RAISE( \
    fmt::format( \
      "Call failed: '{}' with code {},\n" \
      "Gurobi message: {}, hint: {}", #call, e, \
           GRBgeterrormsg(env()), msg ) \
  ); } while (0)
#define GRB_CALL( call ) do { if (int e=call) MP_RAISE( \
    fmt::format( \
      "Call failed: '{}' with code {}, Gurobi message: {}", #call, \
        e, GRBgeterrormsg(env()) ) \
  ); } while (0)

#endif // GUROBICOMMON_H
