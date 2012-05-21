#include "gtest/gtest.h"

#include <ilconcert/ilomodel.h>
#include "asl.h"
#include "nlp.h"
#include "opcode.hd"

IloExpr build_expr(expr *e);

namespace {

TEST(CPLEXCPTest, Fails) {
  expr_n n = {reinterpret_cast<efunc_n*>(OPNUM), 0.42};
  IloExpr result = build_expr(reinterpret_cast<expr*>(&n));
  EXPECT_EQ(0.42, result.getConstant());
}

}

