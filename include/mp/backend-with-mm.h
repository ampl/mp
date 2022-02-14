#ifndef BACKEND_WITH_MM_H
#define BACKEND_WITH_MM_H

/// Backend with a Model Manager
#include <memory>

#include "mp/backend-base.h"
#include "mp/model-mgr-base.h"

namespace mp {

/// Backends using a separate Model Manager
/// could derive from this
class BackendWithModelManager
    : public BasicBackend {
public:
  /// Constructs a BasicBackend object.
  /// date:  The solver date in YYYYMMDD format.
  /// flags: Bitwise OR of zero or more of the following values
  ///          MULTIPLE_SOL
  ///          MULTIPLE_OBJ
  BackendWithModelManager(fmt::CStringRef name, fmt::CStringRef long_name,
               long date, int flags) :
    BasicBackend(name, long_name, date , flags) { }

  /// Initialize backend, incl. solver options
  /// @param argv: the command-line arguments, NULL-terminated
  void Init(char** argv) override {
    BasicBackend::Init(argv);
    GetMM().InitOptions();
    InitMetaInfoAndOptions();
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

private:
  PMM p_model_mgr_;
};

} // namespace mp

#endif // BACKEND_WITH_MM_H
