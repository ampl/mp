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

#include "mp/utils_file.h"

namespace mp {

void ProcessLines_AvoidComments(std::istream& stream,
                  std::function<void(const char*)> processor) {
  std::string line;
  while (stream.good() && !stream.eof()) {
    std::getline(stream, line);
    if (line.size()) {
      auto itfirstns = std::find_if(line.begin(), line.end(),
                                  [](char c){ return !std::isspace(c); });
      if (line.end()!=itfirstns &&
              '#'!=*itfirstns) {                  // Skip commented line
        processor(line.c_str()+(itfirstns-line.begin()));
      }
    }
  }
}

}  // namespace mp
