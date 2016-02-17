/*
 Tests of common definitions

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

#include "test-assert.h"
#include "mp/common.h"

#include <algorithm>
#include <gtest/gtest.h>

namespace ex = mp::expr;

struct ExprInfo {
  ex::Kind kind;
  const char *str;
  int opcode;
  ex::Kind first_kind;
};

// The opcodes are defined in the technical report Writing .nl Files
// (http://ampl.github.io/nlwrite.pdf).
const ExprInfo INFO[] = {
  {ex::UNKNOWN,        "unknown",                 -1, ex::UNKNOWN},

  {ex::NUMBER,         "number",                  80, ex::NUMBER},
  {ex::VARIABLE,       "variable",                82, ex::FIRST_REFERENCE},
  {ex::COMMON_EXPR,    "common expression",       -1, ex::FIRST_REFERENCE},

  // Unary expressions.
  {ex::MINUS,          "unary -",                 16, ex::FIRST_UNARY},
  {ex::ABS,            "abs",                     15, ex::FIRST_UNARY},
  {ex::FLOOR,          "floor",                   13, ex::FIRST_UNARY},
  {ex::CEIL,           "ceil",                    14, ex::FIRST_UNARY},
  {ex::SQRT,           "sqrt",                    39, ex::FIRST_UNARY},
  {ex::POW2,           "^2",                      77, ex::FIRST_UNARY},
  {ex::EXP,            "exp",                     44, ex::FIRST_UNARY},
  {ex::LOG,            "log",                     43, ex::FIRST_UNARY},
  {ex::LOG10,          "log10",                   42, ex::FIRST_UNARY},
  {ex::SIN,            "sin",                     41, ex::FIRST_UNARY},
  {ex::SINH,           "sinh",                    40, ex::FIRST_UNARY},
  {ex::COS,            "cos",                     46, ex::FIRST_UNARY},
  {ex::COSH,           "cosh",                    45, ex::FIRST_UNARY},
  {ex::TAN,            "tan",                     38, ex::FIRST_UNARY},
  {ex::TANH,           "tanh",                    37, ex::FIRST_UNARY},
  {ex::ASIN,           "asin",                    51, ex::FIRST_UNARY},
  {ex::ASINH,          "asinh",                   50, ex::FIRST_UNARY},
  {ex::ACOS,           "acos",                    53, ex::FIRST_UNARY},
  {ex::ACOSH,          "acosh",                   52, ex::FIRST_UNARY},
  {ex::ATAN,           "atan",                    49, ex::FIRST_UNARY},
  {ex::ATANH,          "atanh",                   47, ex::FIRST_UNARY},

  // Binary expressions.
  {ex::ADD,            "+",                        0, ex::FIRST_BINARY},
  {ex::SUB,            "-",                        1, ex::FIRST_BINARY},
  {ex::LESS,           "less",                     6, ex::FIRST_BINARY},
  {ex::MUL,            "*",                        2, ex::FIRST_BINARY},
  {ex::DIV,            "/",                        3, ex::FIRST_BINARY},
  {ex::TRUNC_DIV,      "div",                     55, ex::FIRST_BINARY},
  {ex::MOD,            "mod",                      4, ex::FIRST_BINARY},
  {ex::POW,            "^",                        5, ex::FIRST_BINARY},
  {ex::POW_CONST_BASE, "^",                       78, ex::FIRST_BINARY},
  {ex::POW_CONST_EXP,  "^",                       76, ex::FIRST_BINARY},
  {ex::ATAN2,          "atan2",                   48, ex::FIRST_BINARY},
  {ex::PRECISION,      "precision",               56, ex::FIRST_BINARY},
  {ex::ROUND,          "round",                   57, ex::FIRST_BINARY},
  {ex::TRUNC,          "trunc",                   58, ex::FIRST_BINARY},

  {ex::IF,             "if",                      35, ex::IF},
  {ex::PLTERM,         "piecewise-linear term",   64, ex::PLTERM},
  {ex::CALL,           "function call",           79, ex::CALL},

  {ex::MIN,            "min",                     11, ex::FIRST_VARARG},
  {ex::MAX,            "max",                     12, ex::FIRST_VARARG},
  {ex::SUM,            "sum",                     54, ex::SUM},
  {ex::NUMBEROF,       "numberof",                60, ex::NUMBEROF},
  {ex::NUMBEROF_SYM,   "symbolic numberof",       61, ex::NUMBEROF_SYM},
  {ex::COUNT,          "count",                   59, ex::COUNT},

  {ex::BOOL,           "bool",                    80, ex::NUMBER},
  {ex::NOT,            "!",                       34, ex::NOT},

  {ex::OR,             "||",                      20, ex::FIRST_BINARY_LOGICAL},
  {ex::AND,            "&&",                      21, ex::FIRST_BINARY_LOGICAL},
  {ex::IFF,            "<==>",                    73, ex::FIRST_BINARY_LOGICAL},

  {ex::LT,             "<",                       22, ex::FIRST_RELATIONAL},
  {ex::LE,             "<=",                      23, ex::FIRST_RELATIONAL},
  {ex::EQ,             "=",                       24, ex::FIRST_RELATIONAL},
  {ex::GE,             ">=",                      28, ex::FIRST_RELATIONAL},
  {ex::GT,             ">",                       29, ex::FIRST_RELATIONAL},
  {ex::NE,             "!=",                      30, ex::FIRST_RELATIONAL},

  {ex::ATLEAST,        "atleast",                 62, ex::FIRST_LOGICAL_COUNT},
  {ex::ATMOST,         "atmost",                  63, ex::FIRST_LOGICAL_COUNT},
  {ex::EXACTLY,        "exactly",                 66, ex::FIRST_LOGICAL_COUNT},
  {ex::NOT_ATLEAST,    "!atleast",                67, ex::FIRST_LOGICAL_COUNT},
  {ex::NOT_ATMOST,     "!atmost",                 68, ex::FIRST_LOGICAL_COUNT},
  {ex::NOT_EXACTLY,    "!exactly",                69, ex::FIRST_LOGICAL_COUNT},

  {ex::IMPLICATION,    "==>",                     72, ex::IMPLICATION},

  {ex::EXISTS,         "exists",                  71,
   ex::FIRST_ITERATED_LOGICAL},
  {ex::FORALL,         "forall",                  70,
   ex::FIRST_ITERATED_LOGICAL},

  {ex::ALLDIFF,        "alldiff",                 74, ex::FIRST_PAIRWISE},
  {ex::NOT_ALLDIFF,    "!alldiff",                75, ex::FIRST_PAIRWISE},

  {ex::STRING,         "string",                  81, ex::STRING},
  {ex::IFSYM,          "symbolic if",             65, ex::IFSYM}
};

const std::size_t NUM_KINDS = sizeof(INFO) / sizeof(*INFO);

TEST(CommonTest, IsValid) {
  EXPECT_TRUE(mp::internal::IsValid(ex::UNKNOWN));
  EXPECT_TRUE(mp::internal::IsValid(ex::LAST_EXPR));
  EXPECT_FALSE(mp::internal::IsValid(static_cast<ex::Kind>(-1)));
  EXPECT_FALSE(mp::internal::IsValid(static_cast<ex::Kind>(ex::LAST_EXPR + 1)));
}

TEST(CommonTest, Str) {
  EXPECT_EQ(ex::LAST_EXPR + 1u, NUM_KINDS);
  for (std::size_t i = 0; i < NUM_KINDS; ++i)
    EXPECT_STREQ(INFO[i].str, str(INFO[i].kind));
  EXPECT_ASSERT(ex::str(static_cast<ex::Kind>(-1)), "invalid expression kind");
  EXPECT_ASSERT(ex::str(static_cast<ex::Kind>(ex::LAST_EXPR + 1)),
                "invalid expression kind");
}

TEST(CommonTest, NLOpCode) {
  EXPECT_EQ(ex::LAST_EXPR + 1u, NUM_KINDS);
  for (std::size_t i = 0; i < NUM_KINDS; ++i)
    EXPECT_EQ(INFO[i].opcode, nl_opcode(INFO[i].kind));
  EXPECT_ASSERT(ex::nl_opcode(static_cast<ex::Kind>(-1)),
                "invalid expression kind");
  EXPECT_ASSERT(ex::nl_opcode(static_cast<ex::Kind>(ex::LAST_EXPR + 1)),
                "invalid expression kind");
}

TEST(CommonTest, MaxOpCode) {
  int max_opcode = INFO[0].opcode;
  for (std::size_t i = 0; i < NUM_KINDS; ++i)
    max_opcode = std::max(max_opcode, INFO[i].opcode);
  EXPECT_EQ(max_opcode, mp::internal::MAX_OPCODE);
}

struct TestExpr {
  enum { FIRST_KIND = ex::FLOOR, LAST_KIND = ex::SQRT };
};

TEST(CommonTest, Is) {
  using mp::internal::Is;
  EXPECT_TRUE(Is<TestExpr>(ex::FLOOR));
  EXPECT_TRUE(Is<TestExpr>(ex::CEIL));
  EXPECT_TRUE(Is<TestExpr>(ex::SQRT));
  EXPECT_FALSE(Is<TestExpr>(static_cast<ex::Kind>(ex::FLOOR - 1)));
  EXPECT_FALSE(Is<TestExpr>(static_cast<ex::Kind>(ex::SQRT + 1)));
}

TEST(CommonTest, GetOpCodeInfo) {
  for (std::size_t i = 0; i < NUM_KINDS; ++i) {
    int opcode = INFO[i].opcode;
    if (opcode == -1) continue;
    const mp::internal::OpCodeInfo &info = mp::internal::GetOpCodeInfo(opcode);
    ex::Kind kind = INFO[i].kind != ex::BOOL ? INFO[i].kind : ex::NUMBER;
    EXPECT_EQ(kind, info.kind) << str(info.kind);
    EXPECT_EQ(INFO[i].first_kind, info.first_kind);
  }
}
