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
    HandleSolution(solve_result,   // prints if no -AMPL flag
                   GetWarnings() + msg.c_str(),
                   0, 0, 0.0);
    std::fflush(stdout);
    std::fflush(stderr);
    if (!ampl_flag())
      std::exit(EXIT_FAILURE);
  }


protected:
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

  virtual void HandleFeasibleSolution(
      int solve_code, fmt::CStringRef msg,
      const double *x, const double *y, double obj) {
    GetMM().HandleFeasibleSolution(solve_code, msg, x, y, obj);
  }


public:
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


  

  SuffixSet& GetSuffixes(suf::Kind kind) {
    return GetMM().GetSuffixes(kind);
  }
  /// Get suffix names
  std::set<std::string> GetSuffixNames() {
    return GetMM().GetSuffixNames();
  }

  /// Read unpresolved suffix
  template <class N>
  ArrayRef<N> ReadSuffix(const SuffixDef<N>& suf) {
    return GetMM().ReadSuffix(suf);
  }

  virtual ArrayRef<int> ReadIntSuffix(const SuffixDef<int>& suf) {
    return GetMM().ReadSuffix(suf);
  }

  /// @param fint: if not NULL,
  ///   is set to 1 iff the suffix was integer.
  virtual ArrayRef<double> ReadDblSuffix(
      const SuffixDef<double>& suf, int *fint=nullptr) {
    return GetMM().ReadSuffix(suf, fint);
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
  /// Report dbl suffix
  virtual void ReportSuffix(const SuffixDef<double>& suf,
                    ArrayRef<double> values) {
    GetMM().ReportSuffix(suf, values);
  }

  /// Report int or dbl suffix from int or dbl data.
  /// kind & suf::FLOAT determines the resulting kind.
  template <class Array>
  void ReportSuffix(const std::string& name,
                    int kind,
                    const Array& vals,
                    const std::string& suf_table = {}) {
    if (kind & suf::FLOAT) {
      ReportDblSuffix( {name, kind, suf_table}, vals);
    } else {
      ReportIntSuffix( {name, kind, suf_table},
                      std::vector<int>{ vals.begin(), vals.end() });
    }
  }

  /// Report int suffix
  virtual void ReportIntSuffix(const SuffixDef<int>& suf,
                       ArrayRef<int> values) {
    GetMM().ReportSuffix(suf, values);
  }
  /// Report dbl suffix
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

  /// Get obj option setter
  BasicObjOptionSetter* GetObjOptionSetter() {
    return GetMM().GetObjOptionSetter();
  }

  /// Access original (NL) model instance:
  /// integrality flags. Used in solution rounding.
  const std::vector<bool>& IsVarInt() const {
    return GetMM().IsVarInt();
  }

  /// Solver-facing instance info:
  /// has unfixed integer variables?
  bool HasUnfixedIntVars() const {
    return GetMM().HasUnfixedIntVars();
  }


private:
  PMM p_model_mgr_;
};

} // namespace mp

#endif // BACKEND_WITH_MM_H
