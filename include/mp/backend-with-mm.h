#ifndef BACKEND_WITH_MM_H
#define BACKEND_WITH_MM_H

/// Backend with a Model Manager

#include <algorithm>
#include <memory>

#include "mp/backend-base.h"
#include "mp/model-mgr-base.h"

namespace mp {

/// Backends using a separate Model Manager
/// could derive from this
class BackendWithModelManager
    : public BasicBackend {
public:
  /// Initialize backend, incl. solver options
  /// @param argv: the command-line arguments, NULL-terminated
  void Init(char** argv) override {
    BasicBackend::Init(argv);
    InitMetaInfoAndOptions();
    GetMM().InitOptions();
  }

  void ReportError(int solve_result, fmt::CStringRef msg) override {
    HandleSolution(solve_result, msg, 0, 0, 0.0);
  }


protected:
  /// Chance for the Backend to note base IO filename
  virtual void SetBasename(const std::string& filename_base) override {
    GetMM().SetBasename(filename_base);
  }

  /// Deriving backends can use this
  virtual void InitMetaInfoAndOptions() { }

  const BasicModelManager& GetMM() const {
    assert(p_model_mgr_);
    return *p_model_mgr_;
  }
  BasicModelManager& GetMM() {
    assert(p_model_mgr_);
    return *p_model_mgr_;
  }
  using PMM = std::unique_ptr<BasicModelManager>;
  void SetMM(PMM pmm) {
    p_model_mgr_ = std::move(pmm);
  }


protected:
  /// Access to ModelManager's interface
  virtual void HandleSolution(int status, fmt::CStringRef msg,
      const double *x, const double *y, double obj) {
    GetMM().SetSolutionFileName(GetOverridenSolutionFile());
    GetMM().HandleSolution(status, msg, x, y, obj);
  }

  virtual void HandleFeasibleSolution(fmt::CStringRef msg,
      const double *x, const double *y, double obj) {
    GetMM().HandleFeasibleSolution(msg, x, y, obj);
  }

  /// Variables' initial values
  virtual ArrayRef<double> InitialValues() {
    return GetMM().InitialValues();
  }
  /// Variables' initial values: sparsity
  virtual ArrayRef<int> InitialValuesSparsity() {
    return GetMM().InitialValuesSparsity();
  }

  /// Initial dual values
  virtual ArrayRef<double> InitialDualValues() {
    return GetMM().InitialDualValues();
  }
  /// Initial dual values
  virtual ArrayRef<int> InitialDualValuesSparsity() {
    return GetMM().InitialDualValuesSparsity();
  }



  /// Read unpresolved suffix
  template <class N>
  ArrayRef<N> ReadSuffix(const SuffixDef<N>& suf) {
    return GetMM().ReadSuffix(suf);
  }

  virtual ArrayRef<int> ReadIntSuffix(const SuffixDef<int>& suf) {
    return GetMM().ReadSuffix(suf);
  }

  virtual ArrayRef<double> ReadDblSuffix(const SuffixDef<double>& suf) {
    return GetMM().ReadSuffix(suf);
  }

  virtual size_t GetSuffixSize(int kind) {
    return GetMM().GetSuffixSize(kind);
  }

  /// Record suffix values which are written into .sol
  /// by HandleSolution()
  /// Does nothing if vector empty
  virtual void ReportSuffix(const SuffixDef<int>& suf,
                    ArrayRef<int> values) {
    GetMM().ReportSuffix(suf, values);
  }
  virtual void ReportSuffix(const SuffixDef<double>& suf,
                    ArrayRef<double> values) {
    GetMM().ReportSuffix(suf, values);
  }
  virtual void ReportIntSuffix(const SuffixDef<int>& suf,
                       ArrayRef<int> values) {
    GetMM().ReportSuffix(suf, values);
  }
  virtual void ReportDblSuffix(const SuffixDef<double>& suf,
                       ArrayRef<double> values) {
    GetMM().ReportSuffix(suf, values);
  }

  /// Report single value for all elements of the suffix
  template <class N>
  void ReportSingleSuffix(const SuffixDef<N>& suf,
                          N value) {
    std::vector<N> values(
          GetMM().GetSuffixSize(suf.kind()), value);
    GetMM().ReportSuffix(suf, values);
  }

  /// Access underlying model instance: integrality flags
  const std::vector<bool>& IsVarInt() const {
    return GetMM().IsVarInt();
  }

  /// Underlying instance info: has integer variables?
  bool HasIntVars() const {
    const auto& ivi=IsVarInt();
    return ivi.end()!=std::find(ivi.begin(), ivi.end(), true);
  }


private:
  PMM p_model_mgr_;
};

} // namespace mp

#endif // BACKEND_WITH_MM_H
