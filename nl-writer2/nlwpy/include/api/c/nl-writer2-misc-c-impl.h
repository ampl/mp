/*
 C API: "Implementation" of NL Writer utils:
 C++ classes "wrapping" them
 in order to interface for NLWriter2, SOLReader2.

 Copyright (C) 2024 AMPL Optimization Inc.

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization Inc disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Gleb Belov
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
  NLUtils_C_Impl(NLW2_NLUtils_C* pu)
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
  const NLW2_NLUtils_C& NLU() const { return nlu_c_; }

private:
  const NLW2_NLUtils_C nlu_c_;
};

}  // namespace mp

#endif // NLWRITER2MISCCIMPL_H
