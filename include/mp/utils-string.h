/*
 String utils.

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

#ifndef MP_UTILS_STRING_H_
#define MP_UTILS_STRING_H_

#include <vector>
#include <string>

namespace mp {

/// Split string
std::vector<std::string> split_string(const char* );

/// Split string
template <class Str>
inline
std::vector<std::string> split_string(const Str& str)
{ return split_string(str.c_str()); }

/// Class logging strings to file or screen or memory etc.
class BasicLogger {
public:
  /// Destruct
  virtual ~BasicLogger() { }

  /// Is the log active and ok?
  virtual bool IsOpen() const = 0;

  /// append string.
  ///
  /// @return whether all ok.
  template<class Str>
  bool Append (const Str& s) { return Append(s.c_str()); }

  /// Append string.
  ///
  /// @return whether all ok.
  virtual bool Append(const char* ) = 0;
};

}  // namespace mp

#endif  // MP_UTILS_STRING_H_
