#ifndef OBJ_OPTION_BASE_H
#define OBJ_OPTION_BASE_H

#include <vector>

namespace mp {

/// Abstract obj option setter
class BasicObjOptionSetter {
public:
  /// Destroy
  virtual ~BasicObjOptionSetter() { }
  /// Get vector of multi-obj passes with options
  virtual std::vector<int> GetPassesWithOptions() = 0;
  /// Set solver options for given pass
  virtual void SetOptionsForMultiobjPass(int iPass) = 0;

  /// Read-only access to a specific pass's suffix-provided option value,
  /// without applying it anywhere. Unlike SetOptionsForMultiobjPass(),
  /// which writes a whole pass's options into the live environment
  /// (correct for real native solver parameters, since it can be scoped
  /// to a per-objective sub-environment), this is for options that have
  /// no native equivalent (e.g. mip:plateau* options) and so must be
  /// captured by the caller itself, per pass, rather than relied upon to
  /// persist anywhere between passes.
  /// @return true if pass iPass has a suffix-provided value for optname
  ///   (written to val), false otherwise (caller should keep its own
  ///   default).
  virtual bool GetPassOptionValueDbl(
      int iPass, const char* optname, double& val) const = 0;
  /// Same as GetPassOptionValueDbl(), for int-valued options.
  virtual bool GetPassOptionValueInt(
      int iPass, const char* optname, int& val) const = 0;
};

}  // namespace mp

#endif // OBJ_OPTION_BASE_H
