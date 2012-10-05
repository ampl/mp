#include "gtest/gtest.h"
#include "tests/function.h"
#include "solvers/asl.h"

namespace {

class ODBCTest : public ::testing::Test {};

TEST_F(ODBCTest, Load) {
  fun::Library lib("../tables/ampltabl.dll");
  lib.Load();
  // TODO: check if table handlers are registered
}
}
