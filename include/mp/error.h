/*
 Exception classes and assertions

 Copyright (C) 2013 AMPL Optimization Inc

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

 Author: Victor Zverovich
 */

#ifndef MP_ERROR_H_
#define MP_ERROR_H_

#include <cstdlib>

#include "mp/format.h"

namespace mp {

#ifndef MP_ASSERT
  /// Debug assert
  #define MP_ASSERT(condition, message) \
    assert((condition) && (message))
#endif

/// Assert with throw, only for Debug
#ifdef NDEBUG
  #define MP_ASSERT__RAISE(condition, message) \
    do { } while(0)
#else
  #define MP_ASSERT__RAISE(condition, message) \
    do { if (!(condition)) MP_RAISE(message); } while(0)
#endif

/// Assert even for Release
#define MP_ASSERT_ALWAYS(condition, message) \
  do { if (!(condition)) MP_RAISE(message); } while(0)

/// Use this to raise UnsupportedError
#define MP_UNSUPPORTED(name) \
  throw MakeUnsupportedError( name )

/// Raise error
#define MP_RAISE(msg) throw mp::Error(msg)

/// Raise with a solve result code (see sol::Status),
/// which the driver reports in the .sol file
#define MP_RAISE_WITH_CODE(exit_code, msg) \
  throw mp::Error(msg, mp::SolveCode(exit_code))

/// Raise infeasibility
#define MP_INFEAS(msg) \
  MP_RAISE_WITH_CODE(200, std::string("Model infeasible: ") + msg)

/// Silence unused parameter warnings
#define MP_UNUSED(x) (void)(x)

/// Solve result code carried by an Error (see sol::Status).
/// A distinct type, so that an int meant as a message argument
/// can never be taken for a code.
struct SolveCode {
  explicit SolveCode(int v) : value(v) { }
  int value;
};

/// A general error.
///
/// Its exit code is the solve result code that a driver reports
/// for it (see sol::Status), or -1 (the default) for an error
/// without one, which drivers report as a failure (500).
/// Only Error(msg, SolveCode) sets a code.
class Error : public fmt::internal::RuntimeError {
  int exit_code_ = -1;
 protected:
  Error() {}

  void SetMessage(const std::string &message) {
    std::runtime_error &base = *this;
    base = std::runtime_error(message);
  }

  void init(fmt::CStringRef format_str, fmt::ArgList args) {
    SetMessage(fmt::format(format_str, args));
  }

 public:
  /// Construct from a message; no solve result code (-1)
  Error(fmt::CStringRef msg) { SetMessage(msg.c_str()); }

  /// Construct from a message and a solve result code
  Error(fmt::CStringRef msg, SolveCode code) : exit_code_(code.value)
  { SetMessage(msg.c_str()); }

  /// Deprecated: an int argument here was easily meant as a
  /// message argument. Use Error(msg, SolveCode(code)).
  [[deprecated("use mp::Error(msg, mp::SolveCode(code))")]]
  Error(fmt::CStringRef msg, int code) : exit_code_(code)
  { SetMessage(msg.c_str()); }

  /// Deprecated: construct from a format string and its arguments.
  /// Use Error(fmt::format(format_str, args...)).
  template <typename Arg1, typename... Args>
  [[deprecated("use mp::Error(fmt::format(...))")]]
  Error(fmt::CStringRef format_str, const Arg1& arg1, const Args&... args)
  { SetMessage(fmt::format(format_str, arg1, args...)); }

  /// The solve result code, or -1 if none
  int exit_code() const { return exit_code_; }
};

/// The operation is not supported by the object.
class UnsupportedError : public Error {
 public:
  FMT_VARIADIC_(char, , UnsupportedError, init, fmt::CStringRef)
};

/// Makes UnsupportedError with prefix "unsupported: ".
/// Use via macro MP_UNSUPPORTED
inline UnsupportedError MakeUnsupportedError(
    fmt::CStringRef format_str, fmt::ArgList args) {
  return UnsupportedError("unsupported: {}", fmt::format(format_str, args));
}
FMT_VARIADIC(UnsupportedError, MakeUnsupportedError, fmt::CStringRef)

/// An option error
class OptionError : public Error {
public:
  explicit OptionError(fmt::CStringRef message) : Error(message) {}
};

}  // namespace mp

#endif  // MP_ERROR_H_
