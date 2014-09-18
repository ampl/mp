/*
 .sol format support.

 Copyright (C) 2014 AMPL Optimization Inc

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

#include "mp/sol.h"

#include <cstdio>
#include "mp/posix.h"

template <typename T>
static void WriteArray(fmt::BufferedFile &file, int size, const T *values) {
  for (int i = 0; i < size; ++i)
    file.print("{}\n", values[i]);
}

void mp::WriteSol(fmt::StringRef filename, const SolutionRef &sol) {
  fmt::BufferedFile file(filename, "w");
  const char *line_start = sol.message.c_str();
  for (;;) {
    const char *line_end = line_start;
    while (*line_end && *line_end != '\n')
      ++line_end;
    // Replace an empty line with a line containing a single space
    // because an empty line indicates the end of message.
    if (line_end == line_start + 1)
      std::fputc(' ', file.get());
    else
      std::fwrite(line_start, 1, line_end - line_start, file.get());
    std::fputc('\n', file.get());
  }
  std::fputc('\n', file.get());
  // Write options.
  if (sol.num_options) {
    file.print("{}\n", sol.num_options);
    WriteArray(file, sol.num_options, sol.options);
  }
  // TODO: write dual values
  WriteArray(file, sol.num_vars, sol.values);
}
