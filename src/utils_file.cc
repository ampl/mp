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
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Authors: Gleb Belov
 */

#include <string>
#include <algorithm>
#include <cctype>
#include <fstream>

#include "mp/utils-file.h"

namespace mp {

/// Class appending strings to file with given name
class FileAppender__fstream : public BasicFileAppender {
public:
  /// Is the log active and ok?
  bool IsOpen() const override {
    return fs_.is_open() ? fs_.good() : false;
  }

  /// Open file
  bool Open(const std::string& fln, bool fErase=false) override {
    if (fErase) {
      std::ofstream ofs(fln, std::ios::trunc | std::ios::out);
    }
    fs_.open(fln, std::ios::app | std::ios::out);
    return fs_.good();
  }

  /// Close file
  void Close() override { fs_.close(); }

  /// Append string
  bool Append(const char* s) override { fs_ << s; return fs_.good(); }


private:
  std::ofstream fs_;
};

/// FileAppender maker
std::unique_ptr<BasicFileAppender> MakeFileAppender() {
  return std::unique_ptr<BasicFileAppender>
  { new FileAppender__fstream() };
}

}  // namespace mp
