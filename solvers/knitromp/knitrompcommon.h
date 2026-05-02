#ifndef KNITROMPCOMMON_H
#define KNITROMPCOMMON_H

#include <string>
#include <functional>

#include "mp/backend-to-model-api.h"

//extern "C" {
// TODO Typically import here the solver's C API headers
// Here we use a cplus plus stub instead
   #include "knitro.h"
//}

#include "mp/format.h"

namespace mp {

class NonlinearConstraintData;
/// Information shared by both
/// `KnitrompBackend` and `KnitrompModelAPI`
struct KnitrompCommonInfo {

  KN_context_ptr lp() const { return lp_; };
  KN_context_ptr& lp_ref() { return lp_; }
  void set_lp(KN_context_ptr lp) { lp_ = lp; }

  void set_format_model(std::function<void(fmt::MemoryWriter &)> fn) { 
      formatModelFn = std::move(fn);
  }
  void call_format_model(fmt::MemoryWriter &w) {
      if (formatModelFn) {
          formatModelFn(w);
      }
  }
  bool useJacobian = false; // If to compute gradient/hessian or let knitro do its
  bool useHessian = false; // If to compute hessian or let knitro do it
  int printProblem = 0;

 
private:
    KN_context_ptr      lp_ = NULL;
    // Function to print constraints (set by ModelAPI)
    std::function<void(fmt::MemoryWriter&)> formatModelFn = nullptr;

};


/// Common API for Knitromp classes
class KnitrompCommon :
    public Backend2ModelAPIConnector<KnitrompCommonInfo> {
public:
  /// These methods access Knitromp options. Used by AddSolverOption()
  void GetSolverOption(int key, int& value) const;
  void SetSolverOption(int key, int value);
  void GetSolverOption(int key, double& value) const;
  void SetSolverOption(int key, double value);
  void GetSolverOption(int key, std::string& value) const;
  void SetSolverOption(int key, const std::string& value);
  static constexpr double Infinity() { return KN_INFINITY;  }
  static constexpr double MinusInfinity() { return -KN_INFINITY; }

protected:
  //int getIntAttr(Solver::ATTRIBS name, Solver::ConsType subtype=Solver::ConsType::CONS_LIN) const;
  double getDblAttr(const char* name) const;

  int NumVars() const;
  int NumCons() const;
  int NumObjs() const;
};


#define KNITROMP_RETCODE_OK 0
#define KNITROMP_CCALL( call ) do { if (int e = (call) != KNITROMP_RETCODE_OK) \
  throw std::runtime_error( \
    fmt::format("  Call failed: '{}' with code {}", #call, e )); } while (0)

} // namespace mp

#endif // KNITROMPCOMMON_H
