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

#include <gtest/gtest.h>

#include "mp/common.h"

namespace ex = mp::expr;

struct ExprInfo {
  ex::Kind kind;
  const char *str;
};

const ExprInfo INFO[] = {
  {ex::UNKNOWN,     "unknown"},

  {ex::VARIABLE,    "variable"},
  {ex::COMMON_EXPR, "common expression"},

  // Unary expressions.
  {ex::MINUS,          "unary -"},
  {ex::ABS,            "abs"},
  {ex::FLOOR,          "floor"},
  {ex::CEIL,           "ceil"},
  {ex::SQRT,           "sqrt"},
  {ex::POW2,           "^2"},
  {ex::EXP,            "exp"},
  {ex::LOG,            "log"},
  {ex::LOG10,          "log10"},
  {ex::SIN,            "sin"},
  {ex::SINH,           "sinh"},
  {ex::COS,            "cos"},
  {ex::COSH,           "cosh"},
  {ex::TAN,            "tan"},
  {ex::TANH,           "tanh"},
  {ex::ASIN,           "asin"},
  {ex::ASINH,          "asinh"},
  {ex::ACOS,           "acos"},
  {ex::ACOSH,          "acosh"},
  {ex::ATAN,           "atan"},
  {ex::ATANH,          "atanh"},

  // Binary expressions.
  {ex::ADD,            "+"},
  {ex::SUB,            "-"},
  {ex::LESS,           "less"},
  {ex::MUL,            "*"},
  {ex::DIV,            "/"},
  {ex::INT_DIV,        "div"},
  {ex::MOD,            "mod"},
  {ex::POW,            "^"},
  {ex::POW_CONST_BASE, "^"},
  {ex::POW_CONST_EXP,  "^"},
  {ex::ATAN2,          "atan2"},
  {ex::PRECISION,      "precision"},
  {ex::ROUND,          "round"},
  {ex::TRUNC,          "trunc"},

  {ex::IF,              "if"},
  {ex::PLTERM,          "piecewise-linear term"},
  {ex::CALL,            "function call"},

  {ex::MIN,             "min"},
  {ex::MAX,             "max"},
  {ex::SUM,             "sum"},
  {ex::NUMBEROF,        "numberof"},
  {ex::NUMBEROF_SYM,    "symbolic numberof"},
  {ex::COUNT,           "count"},

  {ex::CONSTANT,        "constant"},

  {ex::NOT,             "!"},

  {ex::OR,              "||"},
  {ex::AND,             "&&"},
  {ex::IFF,             "<==>"},

  {ex::LT,              "<"},
  {ex::LE,              "<="},
  {ex::EQ,              "="},
  {ex::GE,              ">="},
  {ex::GT,              ">"},
  {ex::NE,              "!="},

  {ex::ATLEAST,         "atleast"},
  {ex::ATMOST,          "atmost"},
  {ex::EXACTLY,         "exactly"},
  {ex::NOT_ATLEAST,     "!atleast"},
  {ex::NOT_ATMOST,      "!atmost"},
  {ex::NOT_EXACTLY,     "!exactly"},

  {ex::IMPLICATION,     "==>"},

  {ex::EXISTS,          "exists"},
  {ex::FORALL,          "forall"},

  {ex::ALLDIFF,         "alldiff"},
  {ex::NOT_ALLDIFF,     "!alldiff"},

  {ex::STRING,          "string"},
  {ex::IFSYM,           "symbolic if"}
};

TEST(CommonTest, Str) {
  std::size_t num_kinds = sizeof(INFO) / sizeof(*INFO);
  EXPECT_EQ(ex::LAST_EXPR + 1u, num_kinds);
  for (std::size_t i = 0; i < num_kinds; ++i)
    EXPECT_STREQ(INFO[i].str, str(INFO[i].kind));
}
