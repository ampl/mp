/*
 Mock file

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

#ifndef MP_MOCK_FILE_H_
#define MP_MOCK_FILE_H_

#include "gmock/gmock.h"
#include "mp/format.h"

struct MockFile {
  MockFile() {}
  MockFile(fmt::StringRef, int) {}
  MockFile(const MockFile &) {}
  MockFile &operator=(const MockFile &) { return *this; }

  MOCK_CONST_METHOD0(descriptor, int ());
  MOCK_CONST_METHOD0(size, fmt::LongLong ());
  MOCK_CONST_METHOD2(read, std::size_t (void *buffer, std::size_t count));
};

#endif  // MP_MOCK_FILE_H_
