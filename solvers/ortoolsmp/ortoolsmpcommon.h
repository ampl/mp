#ifndef ORTOOLSCOMMON_H
#define ORTOOLSCOMMON_H

#include <string>
#include <limits>
#include <memory> // for auto_ptr

#include "mp/backend-to-model-api.h"

#include "ortools/linear_solver/linear_solver.h"
namespace orr = operations_research;

/// The below would go into actual ...common.h:

#include "mp/format.h"

namespace mp {
/// Information shared by both
/// `OrtoolsBackend` and `OrtoolsModelAPI`
struct OrtoolsCommonInfo {
  orr::MPSolver* lp() const { return lp_; }
  void set_lp(orr::MPSolver* lp) { lp_ = lp; }

 private:
  orr::MPSolver* lp_ = NULL;
};


/// Common API for Ortools classes
class OrtoolsCommon :
    public Backend2ModelAPIConnector<OrtoolsCommonInfo> {
public:
  /// These methods access Ortools options. Used by AddSolverOption()
  void GetSolverOption(const char* key, int& value) const;
  void SetSolverOption(const char* key, int value);
  void GetSolverOption(const char* key, double& value) const;
  void SetSolverOption(const char* key, double value);
  void GetSolverOption(const char* key, std::string& value) const;
  void SetSolverOption(const char* key, const std::string& value);

  double Infinity() { return lp()->infinity(); }
  double MinusInfinity() { return -Infinity(); }

protected:
  std::map<std::string, int> paramNames_ = {
      {"RELATIVE_MIP_GAP",
       orr::MPSolverParameters::RELATIVE_MIP_GAP},
      {"PRIMAL_TOLERANCE",
       orr::MPSolverParameters::PRIMAL_TOLERANCE},
      {"DUAL_TOLERANCE",
       orr::MPSolverParameters::DUAL_TOLERANCE},

      {"PRESOLVE", orr::MPSolverParameters::PRESOLVE},
      {"LP_ALGORITHM", orr::MPSolverParameters::LP_ALGORITHM},
      {"INCREMENTALITY",
       orr::MPSolverParameters::INCREMENTALITY},
      {"SCALING", orr::MPSolverParameters::SCALING}};

  orr::MPSolverParameters params_;

  int getIntAttr(int name) const;
  double getDblAttr(const char* name) const;

  int NumLinCons() const;
  int NumVars() const;
  int NumObjs() const;
  int NumQPCons() const;
  int NumSOSCons() const;
  int NumIndicatorCons() const;

  class OrtoolsOptionsManager {

  public:
    struct Param {
      std::string id_;
      const std::string& id() const { return id_; }
      Param() : id_() {}
      Param(const std::string& id) : id_(id) {}
      virtual std::string value() { throw std::runtime_error("Not implemented in base class."); }
      virtual std::string toString() { return fmt::format("{}:{}", id(), value()); }
    };

    struct IntParam : public Param {
      IntParam(const std::string& id, int value) : Param(id), value_(value) { }
      int value_;
      std::string value() {
        return fmt::format("{}", value_);
      }
    };
    struct DoubleParam : public Param {
      DoubleParam(const std::string& id, double value) : Param(id), value_(value) { }
      double value_;
      std::string value() {
        return fmt::format("{}", value_);
      }
    };
    struct StringParam : public Param {
      StringParam(const std::string& id, const std::string& value) : Param(id), value_(value) { }
      std::string value_;
      std::string value() {
        return value_;
      }
    };


    void set(const std::string& id, int value) {
      params_[id]=std::make_unique<IntParam>(id, value);
    }
    void set(const std::string& id, double value) {
      params_[id] = std::make_unique<DoubleParam>(id, value);
    }
    void set(const std::string& id, const std::string& value) {
      params_[id] = std::make_unique<StringParam>(id, value);
    }
    int getIntValue(const std::string& id) const {
     return dynamic_cast<IntParam*>(params_.at(id).get())->value_;
    }
    double getDoubleValue(const std::string& id) const {
      return dynamic_cast<DoubleParam* > (params_.at(id).get())->value_;
    }
    const std::string& getStringValue(const std::string& id) const {
      return dynamic_cast<StringParam*> (params_.at(id).get())->value_;
    }
    const Param* get(const std::string& id) const {
      return params_.at(id).get();
    }
    bool hasParams() { return params_.size() > 0; }

    const std::map<std::string, std::unique_ptr <Param>>& params() const {
      return params_;
    };


  
  private:
    std::map<std::string, std::unique_ptr <Param>> params_;
  };
  OrtoolsOptionsManager optionsManager_;
};


// This macro is useful to automatically throw an error if a function in the 
// solver API does not return a valid errorcode. In this mock driver, we define it 
// ourselves, normally this constant would be defined in the solver's API.
#define ORTOOLS_RETCODE_OK 0
#define ORTOOLS_CCALL( call ) do { if (int e = (call) != ORTOOLS_RETCODE_OK) \
  throw std::runtime_error( \
    fmt::format("  Call failed: '{}' with code {}", #call, e )); } while (0)

} // namespace mp

#endif // ORTOOLSCOMMON_H
