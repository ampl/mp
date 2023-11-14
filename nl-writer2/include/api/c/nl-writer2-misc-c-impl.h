/**
 * "Implementation" of NL Writer utils:
 * C++ classes "wrapping" them
 * in order to interface for NLWriter2, NLReader2
 */
#ifndef NLWRITER2MISCCIMPL_H
#define NLWRITER2MISCCIMPL_H

#include "api/c/nl-writer2-misc-c.h"

#include "mp/nl-writer2-misc.h"

namespace mp {

/// Wrap NLUtils_C into a C++ class,
/// in order to interface it for NLWriter2, NLReader2
class NLUtils_C_Impl
    : public NLUtils {
public:
  /// Construct
  NLUtils_C_Impl(NLUtils_C* pu)
    : nlu_c_(*pu) { }

  /// log message
  void log_message(const char* format, ...) override {
    assert(NLU().log_message);
    va_list args;
    va_start (args, format);
    NLU().log_message (NLU().p_user_data_, format, args);
    va_end (args);
  }
  /// log warning
  void log_warning(const char* format, ...) override {
    assert(NLU().log_warning);
    va_list args;
    va_start (args, format);
    NLU().log_warning (NLU().p_user_data_, format, args);
    va_end (args);
  }
  /// Override this to your error handler.
  virtual void myexit(const std::string& msg) override {
    assert(NLU().myexit);
    NLU().myexit(NLU().p_user_data_, msg.c_str());
  }

protected:
  const NLUtils_C& NLU() const { return nlu_c_; }

private:
  const NLUtils_C nlu_c_;
};

}  // namespace mp

#endif // NLWRITER2MISCCIMPL_H
