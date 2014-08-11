/*
 Solver C API tests.

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

#include <gtest/gtest.h>
#include "solvers/util/solver-c.h"

#include <stdlib.h>

#ifdef _WIN32
# define putenv _putenv
#endif

namespace {

TEST(SolverCTest, CreateSolver) {
  ASL_Solver *s = ASL_CreateSolver(0, 0);
  EXPECT_TRUE(s != 0);
  ASL_DestroySolver(s);
}

TEST(SolverCTest, DestroyNullSolver) {
  ASL_DestroySolver(0);
}

TEST(SolverCTest, CreateSolverError) {
  ASL_Error *error = 0;
  ASL_Solver *s = ASL_CreateSolver("", &error);
  EXPECT_TRUE(!s);
  ASSERT_TRUE(error != 0);
  EXPECT_STREQ("epic fail", ASL_GetErrorMessage(error));
  ASL_DestroyError(error);
}

TEST(SolverCTest, DestroyNullError) {
  ASL_DestroyError(0);
}

TEST(SolverCTest, GetLastError) {
  ASL_Solver *s = ASL_CreateSolver(0, 0);
  EXPECT_TRUE(!ASL_GetLastError(s));
  char arg0[] = "test";
  char arg1[] = "../data/test";
  char arg2[] = "opt1=die";
  char *argv[] = {arg0, arg1, arg2, 0};
  EXPECT_EQ(-1, ASL_RunSolver(s, 2, argv));
  ASL_Error *error = ASL_GetLastError(s);
  EXPECT_TRUE(error != 0);
  EXPECT_STREQ("epic fail", ASL_GetErrorMessage(error));
  ASL_DestroyError(error);
  ASL_DestroySolver(s);
}

TEST(SolverCTest, GetOptionHeader) {
  ASL_Solver *s = ASL_CreateSolver(0, 0);
  EXPECT_STREQ("Options rock!", ASL_GetOptionHeader(s));
  ASL_DestroySolver(s);
}

TEST(SolverCTest, GetSolverOptions) {
  ASL_Solver *s = ASL_CreateSolver(0, 0);
  int num_options = ASL_GetSolverOptions(s, 0, 0);
  EXPECT_EQ(5, num_options);
  std::vector<ASL_SolverOptionInfo> options(num_options);
  EXPECT_EQ(num_options, ASL_GetSolverOptions(s, &options[0], num_options));
  EXPECT_STREQ("opt1", options[0].name);
  EXPECT_STREQ("desc1", options[0].description);
  EXPECT_EQ(0, options[0].flags);
  EXPECT_STREQ("opt2", options[1].name);
  EXPECT_STREQ("desc2", options[1].description);
  EXPECT_EQ(ASL_OPT_HAS_VALUES, options[1].flags);
  EXPECT_STREQ("timing", options[2].name);
  EXPECT_STREQ("version", options[3].name);
  EXPECT_STREQ("wantsol", options[4].name);
  ASL_DestroySolver(s);
}

TEST(SolverCTest, GetPartOfSolverOptions) {
  ASL_Solver *s = ASL_CreateSolver(0, 0);
  int num_options = ASL_GetSolverOptions(s, 0, 0);
  EXPECT_EQ(5, num_options);
  std::vector<ASL_SolverOptionInfo> options(3);
  EXPECT_EQ(num_options, ASL_GetSolverOptions(s, &options[0], 2));
  EXPECT_STREQ("opt1", options[0].name);
  EXPECT_STREQ("opt2", options[1].name);
  EXPECT_TRUE(!options[2].name);
  EXPECT_TRUE(!options[2].description);
  EXPECT_TRUE(!options[2].flags);
  EXPECT_TRUE(!options[2].option);
  ASL_DestroySolver(s);
}

TEST(SolverCTest, GetOptionValues) {
  ASL_Solver *s = ASL_CreateSolver(0, 0);
  ASL_SolverOptionInfo info[2];
  ASL_GetSolverOptions(s, info, 2);
  int num_values = ASL_GetOptionValues(s, info[1].option, 0, 0);
  EXPECT_EQ(3, num_values);
  std::vector<ASL_OptionValueInfo> values(num_values);
  EXPECT_EQ(num_values,
      ASL_GetOptionValues(s, info[1].option, &values[0], num_values));
  EXPECT_STREQ("val1", values[0].value);
  EXPECT_STREQ("valdesc1", values[0].description);
  EXPECT_STREQ("val2", values[1].value);
  EXPECT_STREQ("valdesc2", values[1].description);
  EXPECT_STREQ("val3", values[2].value);
  EXPECT_STREQ("valdesc3", values[2].description);
  ASL_DestroySolver(s);
}

TEST(SolverCTest, GetPartOfOptionValues) {
  ASL_Solver *s = ASL_CreateSolver(0, 0);
  ASL_SolverOptionInfo info[2];
  ASL_GetSolverOptions(s, info, 2);
  int num_values = ASL_GetOptionValues(s, info[1].option, 0, 0);
  EXPECT_EQ(3, num_values);
  std::vector<ASL_OptionValueInfo> values(num_values);
  EXPECT_EQ(num_values,
      ASL_GetOptionValues(s, info[1].option, &values[0], 2));
  EXPECT_STREQ("val1", values[0].value);
  EXPECT_STREQ("valdesc1", values[0].description);
  EXPECT_STREQ("val2", values[1].value);
  EXPECT_STREQ("valdesc2", values[1].description);
  EXPECT_TRUE(!values[2].value);
  EXPECT_TRUE(!values[2].description);
  ASL_DestroySolver(s);
}

TEST(SolverCTest, RunSolver) {
  ASL_Solver *s = ASL_CreateSolver(0, 0);
  char arg0[] = "test";
  char arg1[] = "../data/test";
  char *argv[] = {arg0, arg1, 0};
  EXPECT_EQ(0, ASL_RunSolver(s, 2, argv));
  ASL_DestroySolver(s);
}
}
