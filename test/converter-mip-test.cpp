/*
 Tests of a MIP interface

 Copyright (C) 2020 AMPL Optimization Inc

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

 Author: Gleb Belov <Gleb.Belov@monash.edu>
 */

#include "gtest/gtest.h"

#include "converter-mip-test.h"

// TODO adapt for interface-to-end test of any concrete backend


namespace {

using namespace mip_converter_test;

TEST(RedefsMIPTest, PureMILP__01__seemsPassedToBackendCorrectly) {
  MIPInstance milp {
    { { minimize_,
        { { 1, 1, 2.5 },
          { 0, 1, 2 } } } },
    { {  0.0,  0.0, -1.2 },
      {  1.0,  3.0,  5.2 },
      {  I_,   F_,   F_ } },
    {
      { { { 2, 0.0, 4 },
          { 0, 1,   2 } }, -infty_, 56.4 }
    }
  };
  MIPConverterTester tester;
  feedInstance(tester, milp);
  tester.ConvertModel();
  ASSERT_TRUE(tester.ObjsEqual(milp));
  ASSERT_TRUE(tester.VarBoundsEqual(milp));
  ASSERT_TRUE(tester.NConstrEqual(milp));
}

}
