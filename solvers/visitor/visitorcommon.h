#ifndef VISITORCOMMON_H
#define VISITORCOMMON_H

#include <string>

#include "mp/backend-to-model-api.h"

extern "C" {
// TODO Typically import here the solver's C API headers
//  #include "visitor.h"
}

/// Instead, faking a typical solver namespace and defs:
namespace Solver {

  enum ATTRIBS {
    NVARS_INT,
    NVARS_CONT,

    NCONS,
    NCONS_TYPE,

    NOBJS,
    ISQOBJ

  };
  enum TYPE {
    CONS_LIN,
    CONS_QUAD,
    CONS_QUAD_CONE,
    CONS_QUAD_CONE_ROTATED,
    CONS_INDIC,
    CONS_SOS,

    CONS_MAX,
    CONS_MIN,
    CONS_ABS,
    CONS_AND,
    CONS_OR,

    CONS_EXP,
    CONS_EXPA,
    CONS_LOG,
    CONS_LOGA,

    CONS_POW,
    CONS_SIN,
    CONS_COS,
    CONS_TAN,

    CONS_PL
  };
  
  const int NTYPES=20;
  class SolverModel {
    
    int nEntities_[NTYPES];

    std::vector<bool> vars_;
    int nobj_ = 0;
    bool quadObj_ = false;
  public:

    int getAttribute(ATTRIBS a, TYPE t) {
      if (a == NCONS)
      {
        int n = 0;
        for (int i = 0; i < NTYPES; i++)
          n += nEntities_[i];
        return n;
      }
      if (a==NCONS_TYPE)
        return nEntities_[t];
      if (a == NVARS_CONT)
        return getNumVars(false);
      if (a == NVARS_INT)
        return getNumVars(true);
      if (a == NOBJS)
        return nobj_;
      if (a == ISQOBJ)
        return quadObj_;

    }
    void allocateVars(int nvars) {
      vars_.resize(nvars);
    }
    void setVariable(int index, bool integer) {
      vars_[index]=integer;
    }
    std::size_t getNumVars(bool integer) {
      std::size_t count = 0;
      for (auto v : vars_)
        if (v == integer) count++;
      return count;
    }
    void addEntity(TYPE type) {
      nEntities_[type]++;
    }

    void addObjective() {
      nobj_++;
    }
    void addQuadTerms() {
      quadObj_ = true;
    }
  };

  SolverModel* CreateSolverModel(); 
}


/// The below would go into actual ...common.h:

#include "mp/format.h"

namespace mp {

/// Information shared by both
/// `VisitorBackend` and `VisitorModelAPI`
struct VisitorCommonInfo {
  // TODO provide accessors to the solver's in-memory model/environment
  //visitor_env* env() const { return env_; }
  Solver::SolverModel* lp() const { return lp_; }

  // TODO provide accessors to the solver's in-memory model/environment
  //void set_env(visitor_env* e) { env_ = e; }
  void set_lp(Solver::SolverModel* lp) { lp_ = lp; }


private:
  // TODO provide accessors to the solver's in-memory model/environment
  //visitor_env*      env_ = NULL;
  Solver::SolverModel*      lp_ = NULL;

};


/// Common API for Visitor classes
class VisitorCommon :
    public Backend2ModelAPIConnector<VisitorCommonInfo> {
public:
  /// These methods access Visitor options. Used by AddSolverOption()
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
  int getIntAttr(Solver::ATTRIBS name, Solver::TYPE subtype=Solver::CONS_ABS) const;
  double getDblAttr(const char* name) const;

  int NumLinCons() const;
  int NumVars() const;
  int NumObjs() const;
  int NumQPCons() const;
  int NumSOSCons() const;
  int NumIndicatorCons() const;

protected:
  // TODO if desirable, provide function to create the solver's environment
  // with own license
  // int (*createEnv) (solver_env**) = nullptr;
  
};


/// Convenience macro
// TODO This macro is useful to automatically throw an error if a function in the 
// solver API does not return a valid errorcode. In this mock driver, we define it 
// ourselves, normally this constant would be defined in the solver's API.
#define VISITOR_RETCODE_OK 0
#define VISITOR_CCALL( call ) do { if (int e = (call) != VISITOR_RETCODE_OK) \
  throw std::runtime_error( \
    fmt::format("  Call failed: '{}' with code {}", #call, e )); } while (0)

} // namespace mp

#endif // VISITORCOMMON_H
