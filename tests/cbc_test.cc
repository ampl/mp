/*
 CBC solver tests.

 Copyright (C) 2012 AMPL Optimization LLC

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization LLC disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Victor Zverovich
 */

#include <cstdio>
#include <algorithm>

#include "gtest/gtest.h"
#include "tests/util.h"

TEST(CBCTest, RunSolver) {
  std::remove("data/assign0.sol");
  EXPECT_EQ(0, std::system("../solvers/cbc/bin/cbc data/assign0 -AMPL"));
  std::string s = ReadFile("data/assign0.sol");
  std::replace(s.begin(), s.end(), '\r', '\n');
  EXPECT_TRUE(s.find("objective 6\n") != std::string::npos);
}
