/*
 A program for generating tables with information
 about expressions and NL opcodes.

 Copyright (C) 2023 AMPL Optimization Inc

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

 Author: Victor Zverovich, Gleb Belov
 */

#include <algorithm>
#include <vector>
#include <cstdlib>
#include <cassert>

#include "mp/common.h"
#include "mp/posix.h"
#include "expr-writer.h"

namespace expr = mp::expr;
namespace prec = mp::prec;


/// MP expression information.
/// Temporary storage to generate information files.
struct ExprInfo {
  const char *name;
  mp::expr::Kind kind;
  const char *start;
  int opcode;
  mp::prec::Precedence prec;
  const char *prec_str;
  const char *str;
};

#define EXPR_(name, first, opcode, precedence, str) \
  {#name, expr::name, #first, opcode, prec::precedence, #precedence, str}

#define EXPR(name, opcode, precedence, str) \
  EXPR_(name, name, opcode, precedence, str)

const ExprInfo info[] = {
  /// Not an NL opcode - MP internal
  EXPR(NUMBER,      80, PRIMARY, "number"),
  /// Not an NL opcode - MP internal
  EXPR(VARIABLE,    82, PRIMARY, "variable"),
  /// Not an NL opcode - MP internal
  EXPR(COMMON_EXPR, -1, PRIMARY, "common expression"),

  /// Unary expressions.
#define UNARY(name, opcode, precedence, str) \
  EXPR_(name, FIRST_UNARY, opcode, precedence, str)
  UNARY(MINUS, 16, UNARY,          "unary -"),
  UNARY(POW2,  77, EXPONENTIATION, "^2"),
  UNARY(FLOOR, 13, CALL,           "floor"),
  UNARY(CEIL,  14, CALL,           "ceil"),
  UNARY(ABS,   15, CALL,           "abs"),
  UNARY(TANH,  37, CALL,           "tanh"),
  UNARY(TAN,   38, CALL,           "tan"),
  UNARY(SQRT,  39, CALL,           "sqrt"),
  UNARY(SINH,  40, CALL,           "sinh"),
  UNARY(SIN,   41, CALL,           "sin"),
  UNARY(LOG10, 42, CALL,           "log10"),
  UNARY(LOG,   43, CALL,           "log"),
  UNARY(EXP,   44, CALL,           "exp"),
  UNARY(COSH,  45, CALL,           "cosh"),
  UNARY(COS,   46, CALL,           "cos"),
  UNARY(ATANH, 47, CALL,           "atanh"),
  UNARY(ATAN,  49, CALL,           "atan"),
  UNARY(ASINH, 50, CALL,           "asinh"),
  UNARY(ASIN,  51, CALL,           "asin"),
  UNARY(ACOSH, 52, CALL,           "acosh"),
  UNARY(ACOS,  53, CALL,           "acos"),

  /// Binary expressions.
#define BINARY(name, opcode, precedence, str) \
  EXPR_(name, FIRST_BINARY, opcode, precedence, str)
  BINARY(ADD,             0, ADDITIVE,       "+"),
  BINARY(SUB,             1, ADDITIVE,       "-"),
  BINARY(LESS,            6, ADDITIVE,       "less"),
  BINARY(MUL,             2, MULTIPLICATIVE, "*"),
  BINARY(DIV,             3, MULTIPLICATIVE, "/"),
  BINARY(TRUNC_DIV,      55, MULTIPLICATIVE, "div"),
  BINARY(MOD,             4, MULTIPLICATIVE, "mod"),
  BINARY(POW,             5, EXPONENTIATION, "^"),
  BINARY(POW_CONST_BASE, 78, EXPONENTIATION, "^"),
  BINARY(POW_CONST_EXP,  76, EXPONENTIATION, "^"),
  BINARY(ATAN2,          48, CALL,           "atan2"),
  BINARY(PRECISION,      56, CALL,           "precision"),
  BINARY(ROUND,          57, CALL,           "round"),
  BINARY(TRUNC,          58, CALL,           "trunc"),

  EXPR(IF,     35, CONDITIONAL,  "if"),
  EXPR(PLTERM, 64, CALL,         "piecewise-linear term"),

  /// Not an NL opcode - MP internal
  EXPR(CALL,   79, CALL,         "function call"),

  EXPR_(MIN, FIRST_VARARG, 11, CALL, "min"),
  EXPR_(MAX, FIRST_VARARG, 12, CALL, "max"),
  EXPR(SUM,          54, ITERATIVE, "sum"),
  EXPR(NUMBEROF,     60, CALL,      "numberof"),
  EXPR(NUMBEROF_SYM, 61, CALL,      "symbolic numberof"),
  EXPR(COUNT,        59, CALL,      "count"),

  /// Not an NL opcode - MP internal
  EXPR(BOOL, 80, PRIMARY, "bool"),

  EXPR(NOT, 34, NOT, "!"),

  /// Binary logical expressions.
  EXPR_(OR,  FIRST_BINARY_LOGICAL, 20, LOGICAL_OR,  "||"),
  EXPR_(AND, FIRST_BINARY_LOGICAL, 21, LOGICAL_AND, "&&"),
  EXPR_(IFF, FIRST_BINARY_LOGICAL, 73, IFF,         "<==>"),

  /// Relational expressions.
#define RELATIONAL(name, opcode, str) \
  EXPR_(name, FIRST_RELATIONAL, opcode, RELATIONAL, str)
  RELATIONAL(LT, 22, "<"),
  RELATIONAL(LE, 23, "<="),
  RELATIONAL(EQ, 24, "="),
  RELATIONAL(GE, 28, ">="),
  RELATIONAL(GT, 29, ">"),
  RELATIONAL(NE, 30, "!="),

  /// Logical count expressions.
#define LOGICAL_COUNT(name, opcode, str) \
  EXPR_(name, FIRST_LOGICAL_COUNT, opcode, CALL, str)
  LOGICAL_COUNT(ATLEAST,     62, "atleast"),
  LOGICAL_COUNT(ATMOST,      63, "atmost"),
  LOGICAL_COUNT(EXACTLY,     66, "exactly"),
  LOGICAL_COUNT(NOT_ATLEAST, 67, "!atleast"),
  LOGICAL_COUNT(NOT_ATMOST,  68, "!atmost"),
  LOGICAL_COUNT(NOT_EXACTLY, 69, "!exactly"),

  EXPR(IMPLICATION, 72, IMPLICATION, "==>"),

  EXPR_(FORALL, FIRST_ITERATED_LOGICAL, 70, CALL, "forall"),
  EXPR_(EXISTS, FIRST_ITERATED_LOGICAL, 71, CALL, "exists"),

  EXPR_(ALLDIFF,     FIRST_PAIRWISE, 74, CALL, "alldiff"),
  EXPR_(NOT_ALLDIFF, FIRST_PAIRWISE, 75, CALL, "!alldiff"),

  /// Not an NL opcode - MP internal
  EXPR(STRING, 81, PRIMARY,     "string"),

  EXPR(IFSYM,  65, CONDITIONAL, "symbolic if")
};


/// Expression info comparator
struct OpCodeLess {
  bool operator()(const ExprInfo &lhs, const ExprInfo &rhs) const {
    return lhs.opcode < rhs.opcode;
  }
};


/// Functor to generate expression / opcode files
class ExprInfoFileGenerator {
public:
  /// Construct
  ExprInfoFileGenerator(
      const ExprInfo* info, std::size_t n);

  /// Generate expr-info.cc
  void GenExprInfoFile(const char* fln);

  /// Generate nl-opcodes.h
  void GenNLOpcodesFile(const char* fln);

private:
  const ExprInfo* info_;
  std::size_t num_exprs_, num_opcodes_;
  std::vector<const ExprInfo*> opcode_info_;
};

int main(int argc, char **argv) {
  if (argc != 3)
    return 1;

  ExprInfoFileGenerator eg(info,
                           sizeof(info) / sizeof(*info));
  eg.GenExprInfoFile(argv[1]);
  eg.GenNLOpcodesFile(argv[2]);

  return 0;
}

ExprInfoFileGenerator::ExprInfoFileGenerator(
    const ExprInfo* info, std::size_t n)
  : info_(info), num_exprs_(n) {
  // Arrange information by opcodes.
  num_opcodes_ = std::max_element(
        info, info + num_exprs_, OpCodeLess())->opcode + 1;
  opcode_info_.resize(num_opcodes_);
  for (std::size_t i = 0; i < num_exprs_; ++i) {
    int opcode = info[i].opcode;
    if (opcode != -1 && !opcode_info_[opcode])
      opcode_info_[opcode] = info + i;
  }
}

void ExprInfoFileGenerator::GenExprInfoFile(const char *fln) {
  // Print a table that maps opcodes to expression kinds.
  fmt::BufferedFile f(fln, "w");
  f.print(
        "// This file is automatically generated. Do not edit!\n"
        "\n"
        "#include \"mp/common.h\"\n"
        "#include \"expr-writer.h\"\n"
        "\n"
        "const mp::internal::OpCodeInfo mp::internal::OpCodeInfo::INFO[] = {{\n"
        );
  for (std::size_t i = 0; i < num_opcodes_; ++i) {
    if (i != 0)
      f.print(",\n");
    if (const ExprInfo *ei = opcode_info_[i])
      f.print("  {{expr::{}, expr::{}}}", ei->name, ei->start);
    else
      f.print("  {{expr::UNKNOWN, expr::UNKNOWN}}");
  }
  f.print("\n}};\n\n");

  // Arrange information by kinds.
  std::vector<const ExprInfo*> expr_info(expr::LAST_EXPR + 1);
  for (std::size_t i = 0; i < num_exprs_; ++i)
    expr_info[info_[i].kind] = info + i;

  // Write expression info.
  f.print("const mp::internal::ExprInfo mp::internal::ExprInfo::INFO[] = {{\n");
  f.print("  {{-1, \"unknown\"}}");
  for (std::size_t i = 1; i <= expr::LAST_EXPR; ++i) {
    f.print(",\n");
    const ExprInfo *ei = expr_info[i];
    if (!ei) {
      fmt::print(stderr, "unknown expression kind {}", i);
      std::exit(1);
    }
    f.print("  {{{}, \"{}\"}}", ei->opcode, ei->str);
  }
  f.print("\n}};\n\n");

  // Write precedence info.
  f.print("const mp::prec::Precedence mp::expr::PrecInfo::INFO[] = {{\n");
  f.print("  prec::UNKNOWN");
  for (std::size_t i = 1; i <= expr::LAST_EXPR; ++i) {
    f.print(",\n");
    const ExprInfo *ei = expr_info[i];
    if (!ei) {
      fmt::print(stderr, "unknown expression kind {}", i);
      std::exit(1);
    }
    f.print("  prec::{}", ei->prec_str);
  }
  f.print("\n}};\n");
}

void ExprInfoFileGenerator::GenNLOpcodesFile(const char *fln) {
  // Print a header defining IDs for opcodes + names
  fmt::BufferedFile f(fln, "w");
  f.print(
        "// This file is automatically generated. Do not edit!\n"
        "\n"
        "namespace mp {{\n"
        "\n"
        "/// Namespace nl.\n"
        "/// Defines NL opcodes to be used, e.g., by the NL writer.\n"
        "namespace nl {{\n"
        "\n"
        "/// Struct Opcode\n"
        "struct Opcode {{\n"
        "  int code;           // the opcode\n"
        "  const char* name;   // name\n"
        "}};\n\n"
        );
  for (std::size_t i = 0; i < num_opcodes_; ++i) {
    if (const ExprInfo *ei = opcode_info_[i]) {
      assert((int)i==ei->opcode);
      if (i<79)
        f.print(
              "const Opcode {:14} = {{ {}, \"{}\" }};\n",
              ei->name, ei->opcode, ei->str);
    }
    else
      f.print("       // Unused opcode: {}\n", i);
  }
  f.print("\n}}  // namespace nl\n");
  f.print("}}  // namespace mp\n\n");
}
