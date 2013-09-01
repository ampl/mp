/*
 SMPS writer tests.

 Copyright (C) 2013 AMPL Optimization Inc

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

#include <sstream>

#include "gtest/gtest.h"

#include "solvers/smpswriter/smpswriter.h"
#include "tests/args.h"
#include "tests/util.h"

namespace {

std::vector<std::string> Split(const std::string &s) {
  std::vector<std::string> lines;
  std::stringstream ss(s);
  std::string line;
  while (std::getline(ss, line))
    lines.push_back(line);
  return lines;
}

::testing::AssertionResult AssertFilesEqual(const char *, const char *,
    const std::string &expected_file, const std::string &actual_file) {
  auto expected = Split(ReadFile(expected_file));
  auto actual = Split(ReadFile(actual_file));
  if (expected.size() > actual.size())
    actual.resize(expected.size());
  else if (actual.size() > expected.size())
    expected.resize(actual.size());
  for (std::size_t i = 0,
      n = std::min(expected.size(), actual.size()); i != n; ++i) {
    if (expected[i] != actual[i]) {
      return ::testing::AssertionFailure()
        << "Files differ in line " << i + 1 << "\n"
        << "Expected file: " << expected_file << "\n"
        << "Actual file: " << actual_file << "\n"
        << "Expected line: " << expected[i] << "\n"
        << "Actual line: " << actual[i] << "\n";
    }
  }
  return ::testing::AssertionSuccess();
}

#define EXPECT_FILES_EQ(expected_file, actual_file) \
    EXPECT_PRED_FORMAT2(AssertFilesEqual, expected_file, actual_file)

TEST(SMPSWriterTest, SMPSOutput) {
  static const char *const EXTS[] = {".cor", ".sto", ".tim"};
  static const char *const PROBLEMS[] = {
      "int-var", "random-bound", "random-con-matrix", "random-con-matrix2",
      "random-rhs", "single-scenario", "single-stage",
      "vars-not-in-stage-order", "zero-core-coefs", "zero-core-con"
  };
  int count = 0;
  for (size_t i = 0, n = sizeof(PROBLEMS) / sizeof(*PROBLEMS); i != n; ++i) {
    ampl::SMPSWriter w;
    std::string path("../data/smps/");
    path += PROBLEMS[i];
    WriteFile("test.nl", ReadFile(path + ".nl"));
    WriteFile("test.col", ReadFile(path + ".col"));
    WriteFile("test.row", ReadFile(path + ".row"));
    EXPECT_EQ(0, w.Run(Args("", "test.nl")));
    for (size_t j = 0, n = sizeof(EXTS) / sizeof(*EXTS); j != n; ++j, ++count) {
      EXPECT_FILES_EQ(
          std::string(path) + EXTS[j], std::string("test") + EXTS[j]);
    }
  }
  EXPECT_EQ(10 * 3, count);
}

TEST(SMPSWriterTest, NonlinearNotSupported) {
  ampl::SMPSWriter w;
  WriteFile("test.nl", ReadFile("../data/smps/nonlinear.nl"));
  EXPECT_THROW(w.Run(Args("", "test.nl")), ampl::Error);
}

TEST(SMPSWriterTest, MoreThan2StagesNotSupported) {
  ampl::SMPSWriter w;
  WriteFile("test.nl", ReadFile("../data/smps/three-stage.nl"));
  EXPECT_THROW(w.Run(Args("", "test.nl")), ampl::Error);
}

TEST(SMPSWriterTest, RangesNotSupported) {
  ampl::SMPSWriter w;
  WriteFile("test.nl", ReadFile("../data/smps/range-con.nl"));
  EXPECT_THROW(w.Run(Args("", "test.nl")), ampl::Error);
}

TEST(SMPSWriterTest, InconsistentProbabilities) {
  ampl::SMPSWriter w;
  WriteFile("test.nl", ReadFile("../data/smps/inconsistent-probabilities.nl"));
  WriteFile("test.row",
      ReadFile("../data/smps/inconsistent-probabilities.row"));
  WriteFile("test.col",
      ReadFile("../data/smps/inconsistent-probabilities.col"));
  EXPECT_THROW(w.Run(Args("", "test.nl")), ampl::Error);
}
}
