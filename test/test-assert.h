/*
 Test version of MP_ASSERT

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

#ifndef MP_TEST_ASSERT_H_
#define MP_TEST_ASSERT_H_

#include <stdexcept>
#include "gtest-extra.h"

class AssertionFailure : public std::logic_error {
 public:
  explicit AssertionFailure(const char *message) : std::logic_error(message) {}
};

#define MP_ASSERT(condition, message) \
  if (!(condition)) throw AssertionFailure(message);

// Expects an assertion failure.
#define EXPECT_ASSERT(stmt, message) \
  EXPECT_THROW_MSG(stmt, AssertionFailure, message)

#endif  // MP_TEST_ASSERT_H_
