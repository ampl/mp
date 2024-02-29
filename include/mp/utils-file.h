/*
 File utils

 Copyright (C) 2021 AMPL Optimization Inc

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization Inc disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use	, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Gleb Belov
 */

#ifndef MP_UTILS_FILE_H_
#define MP_UTILS_FILE_H_

#include <string>
#include <memory>

#include "mp/utils-string.h"

namespace mp {

/// Class appending strings to file with given name
class BasicFileAppender
    : public BasicLogger {
public:
  /// Is the log active and ok?
  bool IsOpen() const override = 0;

  /// Open file
  virtual bool Open(const std::string& fln, bool fErase) = 0;

  /// Close file
  virtual void Close() = 0;

  /// append string.
  ///
  /// @return whether all ok.
  template<class Str>
  bool Append (const Str& s) { return Append(s.c_str()); }

  /// Append string.
  ///
  /// @return whether all ok.
  bool Append(const char* ) override = 0;
};

/// FileAppender maker
std::unique_ptr<BasicFileAppender> MakeFileAppender();

}  // namespace mp

#endif  // MP_UTILS_FILE_H_
