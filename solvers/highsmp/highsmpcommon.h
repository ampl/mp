#ifndef HIGHSCOMMON_H
#define HIGHSCOMMON_H

#include <string>
#include <memory> // For std::unique_ptr
#include <vector>

extern "C" {
  #include "interfaces/highs_c_api.h"
}

#include "mp/backend-to-model-api.h"
#include "mp/format.h"
#include "mp/arrayref.h"
#include "highsmploader.h"

namespace mp{ 


  /// Class to store the objectives as they are added from the Model API.
  /// Need a shared reference because priorities and other properties are known
  /// only in the backend, and all info must be set at the same time
  /// (currently in HighsBackend::InputExtras)
  class AccObjectives {
    int numVars_;
    std::vector<double> coeffs;
    std::vector<HighsInt> senses;
    bool hadNativeMultiObj_ = false;
    bool clearedOnce_, hadEmulatedMultiObj_ = false;
    std::vector<double> weight, offset, reltol, abstol;
    std::vector<int> priority;
    HighsLoader* parent_;
  public:
    void setLoader(HighsLoader * l) {
      parent_ = l;
    }
    AccObjectives() : hadNativeMultiObj_(false), 
                      clearedOnce_(false),
                      hadEmulatedMultiObj_(false) {    
    }
    void setNumVars(int numVars) {
      numVars_ = numVars;
    }
    void add(const ::std::vector<int>& indices,
      const ::std::vector<double>& c, bool max) {
      if (clearedOnce_) hadEmulatedMultiObj_ = true;
      if (senses.size() > 0) hadNativeMultiObj_ = true;
      coeffs.resize(coeffs.size() + numVars_);
      for (auto i = 0; i < c.size(); i++)
        coeffs[senses.size() * numVars_ + indices[i]] = c[i];
      senses.push_back(max ? kHighsObjSenseMaximize : kHighsObjSenseMinimize);
    }
    void setInHighs(void* highs) const;
    void setAllInHighs(void* highs) const;
    void setWeights(ArrayRef<double> w);
    void setOffsets(ArrayRef<double> o);
    void setRelTols(ArrayRef<double> r);
    void setAbsTols(ArrayRef<double> r);
    void setPriorities(ArrayRef<int> p);
    int numObjs() const { return senses.size(); }
    bool hadNativeMultiObj() { return hadNativeMultiObj_; }
    bool hadEmulatedMultiObj() { return hadEmulatedMultiObj_; }
    void clear() {
      coeffs.clear();
      senses.clear();
      weight.clear();
      offset.clear();
      reltol.clear();
      abstol.clear();
      priority.clear();
      clearedOnce_ = true;
      
    }
  };


/// Information shared by both
/// `HighsBackend` and `HighsModelAPI`
  struct HighsCommonInfo {
    void* lp() const { return lp_; }
    void set_lp(void* lp) { lp_ = lp; }
    HighsCommonInfo() { }

    AccObjectives& accObjectives() {
      return *accobjs_;
    }
    HighsLoader& loader() const { return *loader_; }
    void setLoader(std::shared_ptr<HighsLoader> l) {
      loader_ = std::move(l);
      
    }
    void setObjectiveAccumulator(std::shared_ptr<AccObjectives> acc) {
      accobjs_ = std::move(acc);
      accobjs_->setLoader(loader_.get());
    }
  private:
    std::shared_ptr<AccObjectives> accobjs_;
    void* lp_ = nullptr;
    mutable std::shared_ptr <HighsLoader> loader_;
  };


/// Common API for Highs classes
  class HighsCommon :
    public Backend2ModelAPIConnector<HighsCommonInfo> {
  public:
    /// These methods �ess Highs options. Used by AddSolverOption()
    void GetSolverOption(const char* key, int& value) const;
    void SetSolverOption(const char* key, int value);
    void GetSolverOption(const char* key, double& value) const;
    void SetSolverOption(const char* key, double value);
    void GetSolverOption(const char* key, std::string& value) const;
    void SetSolverOption(const char* key, const std::string& value);

    double myinf = 0;
    double Infinity() {
      if (!myinf) myinf = Highs_getInfinity(lp());
      return myinf;
    }
    double MinusInfinity() { return -Infinity(); }

  protected:
    void OpenSolver();
    void CloseSolver();

    int64_t getInt64Attr(const char* name)  const;
    int getIntAttr(const char* name) const;
    double getDblAttr(const char* name) const;

    int NumLinCons() const;
    int NumVars() const;
    int NumObjs();

  };

#define HIGHS_CCALL( call ) do { \
  int e = (call); \
  if (e != kHighsStatusOk && e != kHighsStatusWarning) \
    throw std::runtime_error( \
      fmt::format("  Error {} for call {}", e, #call ).c_str()); \
  if (e == kHighsStatusWarning) \
    std::printf("%s\n", \
      fmt::format("  Warning code {} for call {}", e, #call ).c_str()); \
} while (0)

} // namespace mp

#endif // HIGHSCOMMON_H
