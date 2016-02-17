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
};

// The opcodes are defined in the technical report Writing .nl Files
// (http://ampl.github.io/nlwrite.pdf).
const ExprInfo INFO[] = {
  {ex::UNKNOWN,        "unknown",                 -1},

  {ex::NUMBER,         "number",                  80},
  {ex::VARIABLE,       "variable",                82},
  {ex::COMMON_EXPR,    "common expression",       -1},

  // Unary expressions.
  {ex::MINUS,          "unary -",                 16},
  {ex::ABS,            "abs",                     15},
  {ex::FLOOR,          "floor",                   13},
  {ex::CEIL,           "ceil",                    14},
  {ex::SQRT,           "sqrt",                    39},
  {ex::POW2,           "^2",                      77},
  {ex::EXP,            "exp",                     44},
  {ex::LOG,            "log",                     43},
  {ex::LOG10,          "log10",                   42},
  {ex::SIN,            "sin",                     41},
  {ex::SINH,           "sinh",                    40},
  {ex::COS,            "cos",                     46},
  {ex::COSH,           "cosh",                    45},
  {ex::TAN,            "tan",                     38},
  {ex::TANH,           "tanh",                    37},
  {ex::ASIN,           "asin",                    51},
  {ex::ASINH,          "asinh",                   50},
  {ex::ACOS,           "acos",                    53},
  {ex::ACOSH,          "acosh",                   52},
  {ex::ATAN,           "atan",                    49},
  {ex::ATANH,          "atanh",                   47},

  // Binary expressions.
  {ex::ADD,            "+",                        0},
  {ex::SUB,            "-",                        1},
  {ex::LESS,           "less",                     6},
  {ex::MUL,            "*",                        2},
  {ex::DIV,            "/",                        3},
  {ex::TRUNC_DIV,      "div",                     55},
  {ex::MOD,            "mod",                      4},
  {ex::POW,            "^",                        5},
  {ex::POW_CONST_BASE, "^",                       78},
  {ex::POW_CONST_EXP,  "^",                       76},
  {ex::ATAN2,          "atan2",                   48},
  {ex::PRECISION,      "precision",               56},
  {ex::ROUND,          "round",                   57},
  {ex::TRUNC,          "trunc",                   58},

  {ex::IF,             "if",                      35},
  {ex::PLTERM,         "piecewise-linear term",   64},
  {ex::CALL,           "function call",           79},

  {ex::MIN,            "min",                     11},
  {ex::MAX,            "max",                     12},
  {ex::SUM,            "sum",                     54},
  {ex::NUMBEROF,       "numberof",                60},
  {ex::NUMBEROF_SYM,   "symbolic numberof",       61},
  {ex::COUNT,          "count",                   59},

  {ex::BOOL,           "bool",                    80},
  {ex::NOT,            "!",                       34},

  {ex::OR,             "||",                      20},
  {ex::AND,            "&&",                      21},
  {ex::IFF,            "<==>",                    73},

  {ex::LT,             "<",                       22},
  {ex::LE,             "<=",                      23},
  {ex::EQ,             "=",                       24},
  {ex::GE,             ">=",                      28},
  {ex::GT,             ">",                       29},
  {ex::NE,             "!=",                      30},

  {ex::ATLEAST,        "atleast",                 62},
  {ex::ATMOST,         "atmost",                  63},
  {ex::EXACTLY,        "exactly",                 66},
  {ex::NOT_ATLEAST,    "!atleast",                67},
  {ex::NOT_ATMOST,     "!atmost",                 68},
  {ex::NOT_EXACTLY,    "!exactly",                69},

  {ex::IMPLICATION,    "==>",                     72},

  {ex::EXISTS,         "exists",                  71},
  {ex::FORALL,         "forall",                  70},

  {ex::ALLDIFF,        "alldiff",                 74},
  {ex::NOT_ALLDIFF,    "!alldiff",                75},

  {ex::STRING,         "string",                  81},
  {ex::IFSYM,          "symbolic if",             65}
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
