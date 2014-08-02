/*
 An .nl reader tests.

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

#include <climits>

#include "gtest/gtest.h"
#include "solvers/util/aslbuilder.h"
#include "solvers/util/nl.h"
#include "solvers/util/problem.h"
#include "tests/util.h"

using ampl::NLHeader;
using ampl::ReadNLString;

namespace {

class TestNLHandler {
 private:
  // Writes a comma-separated list.
  static void WriteList(fmt::Writer &w, int size, const double *values) {
    for (int i = 0; i < size; ++i) {
      if (i != 0) w << ", ";
      w << values[i];
    }
  }

 public:
  NLHeader header;
  fmt::Writer log;  // Call log.
  std::vector<std::string> obj_exprs;

  typedef std::string Expr;
  typedef std::string NumericExpr;
  typedef std::string LogicalExpr;
  typedef std::string CountExpr;
  typedef std::string Variable;

  void BeginBuild(const char *, const NLHeader &h, int) {
    header = h;
    obj_exprs.resize(h.num_objs);
    log.clear();
  }

  void SetObj(int index, ampl::obj::Type type, std::string expr) {
    log << (type == ampl::obj::MAX ? "maximize" : "minimize")
        << " o" << (index + 1) << ": " << expr;
    obj_exprs[index] = expr;
    log << ";\n";
  }

  void SetCon(int index, std::string expr) {
    log << "c" << index << ": " << expr << ";\n";
  }

  void SetFunction(
      int index, const char *name, int num_args, ampl::func::Type type) {
    // TODO
  }

  std::string MakeNumericConstant(double value) {
    return fmt::format("{}", value);
  }

  std::string MakeVariable(int index) {
    return fmt::format("x{}", index);
  }

  std::string MakeUnary(int opcode, std::string arg) {
    return fmt::format("o{}({})", opcode, arg);
  }

  std::string MakeBinary(int opcode, std::string lhs, std::string rhs) {
    return fmt::format("o{}({}, {})", opcode, lhs, rhs);
  }

  std::string MakeIf(std::string condition,
                     std::string true_expr, std::string false_expr) {
    return fmt::format("if {} then {} else {}",
                       condition, true_expr, false_expr);
  }

  std::string MakePiecewiseLinear(
      int num_breakpoints, const double *breakpoints,
      const double *slopes, std::string var) {
    fmt::Writer w;
    w << "<<";
    WriteList(w, num_breakpoints, breakpoints);
    w << "; ";
    WriteList(w, num_breakpoints + 1, slopes);
    w << ">> " << var;
    return w.str();
  }

  std::string MakeVarArg(int opcode, ampl::ArrayRef<std::string> args) {
    // TODO
    return "";
  }

  std::string MakeSum(ampl::ArrayRef<std::string> args) {
    // TODO
    return "";
  }

  std::string MakeCount(ampl::ArrayRef<std::string> args) {
    // TODO
    return "";
  }

  std::string MakeNumberOf(ampl::ArrayRef<std::string> args) {
    // TODO
    return "";
  }

  std::string MakeLogicalConstant(bool value) {
    return fmt::format("l{}", value);
  }

  std::string MakeNot(std::string arg) { return fmt::format("not {}", arg); }

  std::string MakeBinaryLogical(int opcode, std::string lhs, std::string rhs) {
    return fmt::format("o{}({}, {})", opcode, lhs, rhs);
  }

  std::string MakeRelational(int opcode, std::string lhs, std::string rhs) {
    return fmt::format("o{}({}, {})", opcode, lhs, rhs);
  }

  std::string MakeLogicalCount(int opcode, std::string lhs, std::string rhs) {
    return fmt::format("o{}({}, {})", opcode, lhs, rhs);
  }

  std::string MakeImplication(std::string condition,
                              std::string true_expr, std::string false_expr) {
    return fmt::format("{} ==> {} else {}",
                       condition, true_expr, false_expr);
  }

  std::string MakeIteratedLogical(
      int opcode, ampl::ArrayRef<std::string> args) {
    // TODO
    return "";
  }

  std::string MakeAllDiff(ampl::ArrayRef<std::string> args) {
    // TODO
    return "";
  }
};

TEST(NLTest, WriteTextHeader) {
  NLHeader header = {
    NLHeader::TEXT,
    9, {2, 3, 5, 7, 11, 13, 17, 19, 23}, 1.23,
    29, 47, 37, 41, 43, 31,
    53, 59, 67, 61, 71, 73,
    79, 83,
    89, 97, 101,
    103, 107, 109,
    113, 127, 131, 137, 139,
    149, 151,
    157, 163,
    167, 173, 179, 181, 191
  };
  fmt::Writer w;
  w << header;
  EXPECT_EQ(
      "g9 2 3 5 7 11 13 17 19 23 1.23\n"
      " 29 47 37 41 43 31\n"
      " 53 59 6 61 71 73\n"
      " 79 83\n"
      " 89 97 101\n"
      " 103 107 0 109\n"
      " 113 127 131 137 139\n"
      " 149 151\n"
      " 157 163\n"
      " 167 173 179 181 191\n",
      w.str());
}

TEST(NLTest, WriteBinaryHeader) {
  NLHeader header = {NLHeader::BINARY, 3, {11, 22, 33}};
  fmt::Writer w;
  w << header;
  EXPECT_EQ(
      "b3 11 22 33\n"
      " 0 0 0 0 0 0\n"
      " 0 0 0 0 0 0\n"
      " 0 0\n"
      " 0 0 0\n"
      " 0 0 0 0\n"
      " 0 0 0 0 0\n"
      " 0 0\n"
      " 0 0\n"
      " 0 0 0 0 0\n",
      w.str());
}

TEST(NLTest, WriteBinarySwappedHeader) {
  NLHeader header = {NLHeader::BINARY_SWAPPED};
  fmt::Writer w;
  w << header;
  EXPECT_EQ(
      "b0\n"
      " 0 0 0 0 0 0\n"
      " 0 0 0 0 0 0\n"
      " 0 0\n"
      " 0 0 0\n"
      " 0 0 2 0\n"
      " 0 0 0 0 0\n"
      " 0 0\n"
      " 0 0\n"
      " 0 0 0 0 0\n",
      w.str());
}

// Formats header as a string.
std::string FormatHeader(const NLHeader &h) {
  fmt::Writer w;
  w << h;
  return w.str();
}

// Reads a header from the specified string.
NLHeader ReadHeader(const std::string &s) {
  TestNLHandler handler;
  ReadNLString(s, handler, "(input)", true);
  return handler.header;
}

// Reads a zero header with one modified line.
NLHeader ReadHeader(int line_index, fmt::StringRef line) {
  return ReadHeader(ReplaceLine(
      FormatHeader(NLHeader()), line_index, line.c_str()));
}

TEST(NLTest, NoNewlineAtEOF) {
  TestNLHandler handler;
  ReadNLString("g\n"
    " 1 1 0\n"
    " 0 0\n"
    " 0 0\n"
    " 0 0 0\n"
    " 0 0 0 1\n"
    " 0 0 0 0 0\n"
    " 0 0\n"
    " 0 0\n"
    " 0 0 0 0 0\n"
    "k0\0h", handler);
}

TEST(NLTest, InvalidFormat) {
  EXPECT_THROW_MSG(ReadHeader(0, "x"),
      ampl::ParseError, "(input):1:1: expected format specifier");
}

TEST(NLTest, InvalidNumOptions) {
  EXPECT_EQ(0, ReadHeader(0, "ga").num_options);
  EXPECT_EQ(0, ReadHeader(0, "g-1").num_options);
  EXPECT_THROW_MSG(ReadHeader(0, "g10"),
      ampl::ParseError, "(input):1:2: too many options");
  EXPECT_THROW_MSG(ReadHeader(0,
      fmt::format("g{}", static_cast<unsigned>(INT_MAX) + 1)),
      ampl::ParseError, "(input):1:2: number is too big");
}

void CheckReadOptions(int num_options,
    int num_options_to_write, const int *options) {
  fmt::Writer w;
  w << 'g' << num_options;
  for (int i = 0; i < num_options_to_write; ++i)
    w << ' ' << options[i];
  NLHeader header = ReadHeader(0, w.c_str());
  ASSERT_EQ(num_options, header.num_options);
  int min_num_options = std::min(num_options, num_options_to_write);
  for (int i = 0; i < min_num_options; ++i)
    EXPECT_EQ(options[i], header.options[i]);
  for (int i = min_num_options; i < num_options_to_write; ++i)
    EXPECT_EQ(0, header.options[i]);
}

TEST(NLTest, ReadOptions) {
  const int options[ampl::MAX_NL_OPTIONS + 1] = {
      3, 5, 7, 11, 13, 17, 19, 23, 29, 31
  };
  for (int i = 0; i < ampl::MAX_NL_OPTIONS; ++i) {
    for (int j = 0; j < ampl::MAX_NL_OPTIONS + 1; ++j)
      CheckReadOptions(i, j, options);
  }
  EXPECT_EQ(0, ReadHeader(0, "g").num_options);
}

TEST(NLTest, ReadAMPLVBTol) {
  EXPECT_EQ(4.2, ReadHeader(0, "g2 0 3 4.2").ampl_vbtol);
  EXPECT_EQ(0, ReadHeader(0, "g2 0 0 4.2").ampl_vbtol);
  EXPECT_EQ(0, ReadHeader(0, "g2 0 3").ampl_vbtol);
}

void CheckHeader(const NLHeader &h) {
  std::string nl = FormatHeader(h);
  NLHeader actual_header = ReadHeader(nl);

  EXPECT_EQ(h.format, actual_header.format);

  EXPECT_EQ(h.num_options, actual_header.num_options);
  for (int i = 0; i < ampl::MAX_NL_OPTIONS; ++i)
    EXPECT_EQ(h.options[i], actual_header.options[i]);
  EXPECT_EQ(h.ampl_vbtol, actual_header.ampl_vbtol);

  EXPECT_EQ(h.num_vars, actual_header.num_vars);
  EXPECT_EQ(h.num_algebraic_cons, actual_header.num_algebraic_cons);
  EXPECT_EQ(h.num_objs, actual_header.num_objs);
  EXPECT_EQ(h.num_ranges, actual_header.num_ranges);
  EXPECT_EQ(h.num_eqns, actual_header.num_eqns);
  EXPECT_EQ(h.num_logical_cons, actual_header.num_logical_cons);

  EXPECT_EQ(h.num_nl_cons, actual_header.num_nl_cons);
  EXPECT_EQ(h.num_nl_objs, actual_header.num_nl_objs);
  EXPECT_EQ(h.num_compl_conds, actual_header.num_compl_conds);
  EXPECT_EQ(h.num_nl_compl_conds, actual_header.num_nl_compl_conds);
  EXPECT_EQ(h.num_compl_dbl_ineqs, actual_header.num_compl_dbl_ineqs);
  EXPECT_EQ(h.num_compl_vars_with_nz_lb,
      actual_header.num_compl_vars_with_nz_lb);

  EXPECT_EQ(h.num_nl_net_cons, actual_header.num_nl_net_cons);
  EXPECT_EQ(h.num_linear_net_cons, actual_header.num_linear_net_cons);

  EXPECT_EQ(h.num_nl_vars_in_cons, actual_header.num_nl_vars_in_cons);
  EXPECT_EQ(h.num_nl_vars_in_objs, actual_header.num_nl_vars_in_objs);
  EXPECT_EQ(h.num_nl_vars_in_both, actual_header.num_nl_vars_in_both);

  EXPECT_EQ(h.num_linear_net_vars, actual_header.num_linear_net_vars);
  EXPECT_EQ(h.num_funcs, actual_header.num_funcs);
  EXPECT_EQ(h.flags, actual_header.flags);

  EXPECT_EQ(h.num_linear_binary_vars, actual_header.num_linear_binary_vars);
  EXPECT_EQ(h.num_linear_integer_vars, actual_header.num_linear_integer_vars);
  EXPECT_EQ(h.num_nl_integer_vars_in_both,
      actual_header.num_nl_integer_vars_in_both);
  EXPECT_EQ(h.num_nl_integer_vars_in_cons,
      actual_header.num_nl_integer_vars_in_cons);
  EXPECT_EQ(h.num_nl_integer_vars_in_objs,
      actual_header.num_nl_integer_vars_in_objs);

  EXPECT_EQ(h.num_con_nonzeros, actual_header.num_con_nonzeros);
  EXPECT_EQ(h.num_obj_nonzeros, actual_header.num_obj_nonzeros);

  EXPECT_EQ(h.max_con_name_len, actual_header.max_con_name_len);
  EXPECT_EQ(h.max_var_name_len, actual_header.max_var_name_len);

  EXPECT_EQ(h.num_common_exprs_in_both, actual_header.num_common_exprs_in_both);
  EXPECT_EQ(h.num_common_exprs_in_cons, actual_header.num_common_exprs_in_cons);
  EXPECT_EQ(h.num_common_exprs_in_objs, actual_header.num_common_exprs_in_objs);
  EXPECT_EQ(h.num_common_exprs_in_cons1,
      actual_header.num_common_exprs_in_cons1);
  EXPECT_EQ(h.num_common_exprs_in_objs1,
      actual_header.num_common_exprs_in_objs1);

  if (h.num_vars == 0)
    return;  // jac0dim fails if there are no vars

  WriteFile("test.nl", nl);
  char stub[] = "test.nl";
  ASL *asl = ASL_alloc(ASL_read_fg);
  jac0dim_ASL(asl, stub, static_cast<int>(strlen(stub)));
  std::remove(stub);

  EXPECT_EQ(asl->i.ampl_options_[0], actual_header.num_options);
  for (int i = 0; i < ampl::MAX_NL_OPTIONS; ++i)
    EXPECT_EQ(asl->i.ampl_options_[i + 1], actual_header.options[i]);
  EXPECT_EQ(asl->i.ampl_vbtol_, actual_header.ampl_vbtol);

  EXPECT_EQ(asl->i.n_var_, actual_header.num_vars);
  EXPECT_EQ(asl->i.n_con_, actual_header.num_algebraic_cons);
  EXPECT_EQ(asl->i.n_obj_, actual_header.num_objs);
  EXPECT_EQ(asl->i.nranges_, actual_header.num_ranges);
  EXPECT_EQ(asl->i.n_eqn_, actual_header.num_eqns);
  EXPECT_EQ(asl->i.n_lcon_, actual_header.num_logical_cons);

  EXPECT_EQ(asl->i.nlc_, actual_header.num_nl_cons);
  EXPECT_EQ(asl->i.nlo_, actual_header.num_nl_objs);
  EXPECT_EQ(asl->i.n_cc_, actual_header.num_compl_conds);
  EXPECT_EQ(asl->i.nlcc_, actual_header.num_nl_compl_conds);
  EXPECT_EQ(asl->i.ndcc_, actual_header.num_compl_dbl_ineqs);
  EXPECT_EQ(asl->i.nzlb_, actual_header.num_compl_vars_with_nz_lb);

  EXPECT_EQ(asl->i.nlnc_, actual_header.num_nl_net_cons);
  EXPECT_EQ(asl->i.lnc_, actual_header.num_linear_net_cons);

  EXPECT_EQ(asl->i.nlvc_, actual_header.num_nl_vars_in_cons);
  EXPECT_EQ(asl->i.nlvo_, actual_header.num_nl_vars_in_objs);
  EXPECT_EQ(asl->i.nlvb_, actual_header.num_nl_vars_in_both);

  EXPECT_EQ(asl->i.nwv_, actual_header.num_linear_net_vars);
  EXPECT_EQ(asl->i.nfunc_, actual_header.num_funcs);
  EXPECT_EQ(asl->i.flags, actual_header.flags);

  EXPECT_EQ(asl->i.nbv_, actual_header.num_linear_binary_vars);
  EXPECT_EQ(asl->i.niv_, actual_header.num_linear_integer_vars);
  EXPECT_EQ(asl->i.nlvbi_, actual_header.num_nl_integer_vars_in_both);
  EXPECT_EQ(asl->i.nlvci_, actual_header.num_nl_integer_vars_in_cons);
  EXPECT_EQ(asl->i.nlvoi_, actual_header.num_nl_integer_vars_in_objs);

  EXPECT_EQ(asl->i.nzc_, actual_header.num_con_nonzeros);
  EXPECT_EQ(asl->i.nzo_, actual_header.num_obj_nonzeros);

  EXPECT_EQ(asl->i.maxrownamelen_, actual_header.max_con_name_len);
  EXPECT_EQ(asl->i.maxcolnamelen_, actual_header.max_var_name_len);

  EXPECT_EQ(asl->i.comb_, actual_header.num_common_exprs_in_both);
  EXPECT_EQ(asl->i.comc_, actual_header.num_common_exprs_in_cons);
  EXPECT_EQ(asl->i.como_, actual_header.num_common_exprs_in_objs);
  EXPECT_EQ(asl->i.comc1_, actual_header.num_common_exprs_in_cons1);
  EXPECT_EQ(asl->i.como1_, actual_header.num_common_exprs_in_objs1);

  ASL_free(&asl);
}

TEST(NLTest, ReadFullHeader) {
  NLHeader header = {
    NLHeader::BINARY,
    9, {2, 3, 5, 7, 11, 13, 17, 19, 23}, 1.23,
    29, 47, 37, 41, 43, 31,
    53, 59, 67, 61, 71, 73,
    79, 83,
    89, 97, 101,
    103, 107, 109,
    113, 127, 131, 137, 139,
    149, 151,
    157, 163,
    167, 173, 179, 181, 191
  };
  CheckHeader(header);
  NLHeader zero_header = {};
  CheckHeader(zero_header);
}

TEST(NLTest, NumComplDblIneq) {
  EXPECT_EQ(42, ReadHeader(2, " 0 0 0 0 42").num_compl_dbl_ineqs);
  EXPECT_EQ(-1, ReadHeader(2, " 0 0 70 0 42").num_compl_dbl_ineqs);
}

TEST(NLTest, Arith) {
  EXPECT_EQ(NLHeader::TEXT, ReadHeader(5, " 0 0").format);
  EXPECT_EQ(NLHeader::TEXT, ReadHeader(5, " 0 0 0").format);
  EXPECT_EQ(NLHeader::TEXT,
      ReadHeader(5, fmt::format(" 0 0 {}", Arith_Kind_ASL)).format);
  EXPECT_EQ(NLHeader::BINARY_SWAPPED,
      ReadHeader(5, fmt::format(" 0 0 {}", 3 - Arith_Kind_ASL)).format);
  EXPECT_THROW_MSG(
      ReadHeader(5, fmt::format(" 0 0 {}", 3 - Arith_Kind_ASL + 1)),
      ampl::ParseError, "(input):6:6: unrecognized binary format");
  // TODO: check if the bytes are actually swapped
}

TEST(NLTest, IncompleteHeader) {
  ReadHeader(0, "g");
  EXPECT_THROW_MSG(
      ReadHeader(0, "\n"),
      ampl::ParseError, "(input):1:1: expected format specifier");
  ReadHeader(1, " 1 0 0");
  EXPECT_THROW_MSG(
      ReadHeader(1, " 1 0"),
      ampl::ParseError, "(input):2:5: expected nonnegative integer");
  for (int i = 2; i <= 8; ++i) {
    if (i == 6)
      continue;
    ReadHeader(i, " 0 0");
    EXPECT_THROW_MSG(
        ReadHeader(i, " 0"), ampl::ParseError,
        fmt::format("(input):{}:3: expected nonnegative integer", i + 1));
  }
  for (int i = 6; i <= 9; i += 3) {
    ReadHeader(1, " 0 0 0 0 0");
    EXPECT_THROW_MSG(
        ReadHeader(i, " 0 0 0 0"), ampl::ParseError,
        fmt::format("(input):{}:9: expected nonnegative integer", i + 1));
  }
  std::string input = ReplaceLine(FormatHeader(NLHeader()), 4, " 0 0");
  ReadHeader(ReplaceLine(input, 6, " 0 0"));
  EXPECT_THROW_MSG(
      ReadHeader(ReplaceLine(input, 6, " 0")),
      ampl::ParseError, "(input):7:3: expected nonnegative integer");
}

void ReadNL(const NLHeader &header, const std::string &body) {
  TestNLHandler handler;
  ReadNLString(FormatHeader(header) + body, handler);
}

TEST(NLTest, ObjIndex) {
  NLHeader header = {};
  header.num_vars = 1;
  header.num_objs = 10;
  EXPECT_THROW_MSG(
    ReadNL(header, "O-1 0\nn0"),
    ampl::ParseError, "(input):11:2: expected nonnegative integer");
  ReadNL(header, "O0 9\nn0");
  EXPECT_THROW_MSG(
    ReadNL(header, "O10 0\nn0"),
    ampl::ParseError, "(input):11:2: objective index 10 out of bounds");
}

TEST(NLTest, ObjType) {
  NLHeader header = {};
  header.num_vars = header.num_objs = 1;
  ReadNL(header, "O0 0\nn0");
  ReadNL(header, "O0 1\nn0");
  ReadNL(header, "O0 10\nn0");
  EXPECT_THROW_MSG(
    ReadNL(header, "O0 -1\nn0"),
    ampl::ParseError, "(input):11:4: expected nonnegative integer");
}

NLHeader MakeHeader() {
  NLHeader header = {};
  header.num_vars = header.num_objs = header.num_algebraic_cons = 10;
  return header;
}

TEST(NLTest, ObjExpr) {
  TestNLHandler handler;
  NLHeader header = MakeHeader();
  header.num_objs = 2;
  ReadNLString(FormatHeader(header) + "O1 0\nn0", handler);
  EXPECT_TRUE(handler.obj_exprs[0].empty());
  EXPECT_EQ("minimize o2: 0;\n", handler.log.str());
  ReadNLString(FormatHeader(header) + "O0 1\nv0", handler);
  EXPECT_EQ("maximize o1: x0;\n", handler.log.str());
}

#define EXPECT_PARSE(expected_output, nl_body) {\
  TestNLHandler handler; \
  ReadNLString(FormatHeader(MakeHeader()) + nl_body, handler); \
  EXPECT_EQ(expected_output, handler.log.str()); \
}

template <typename Int>
void CheckParseInt(char code) {
  EXPECT_PARSE("c0: 4;\n", fmt::format("C0\n{}4.2", code));
  Int min = std::numeric_limits<Int>::min();
  EXPECT_PARSE(fmt::format("c0: {};\n", min + 0.),
               fmt::format("C0\n{}{}", code, min));
  typename safeint::MakeUnsigned<Int>::Type max =
      std::numeric_limits<Int>::max();
  EXPECT_PARSE(fmt::format("c0: {};\n", max + 0.),
               fmt::format("C0\n{}{}", code, max));
  EXPECT_THROW_MSG(
    ReadNL(MakeHeader(), fmt::format("C0\n{}{}", code, max + 1)),
    ampl::ParseError, "(input):12:2: number is too big");
}

TEST(NLTest, ParseNumericConstant) {
  EXPECT_PARSE("c0: 4.2;\n", "C0\nn4.2");
  EXPECT_PARSE("c0: -100;\n", "C0\nn-1e+2");
  CheckParseInt<short>('s');
  CheckParseInt<long>('l');
}

TEST(NLTest, ParseVariable) {
  EXPECT_PARSE("c0: x5;\n", "C0\nv5");
}

TEST(NLTest, ParseUnaryExpr) {
  EXPECT_PARSE("c0: o13(x3);\n", "C0\no13\nv3");
}

TEST(NLTest, ParseBinaryExpr) {
  EXPECT_PARSE("c0: o0(x1, 42);\n", "C0\no0\nv1\nn42");
}

TEST(NLTest, ParseIfExpr) {
  EXPECT_PARSE("c0: if l1 then x1 else x2;\n", "C0\no35\nn1\nv1\nv2");
}

TEST(NLTest, ParsePiecewiseLinearExpr) {
  EXPECT_PARSE("c0: <<0; -1, 1>> x1;\n", "C0\no64\n2\nn-1\nn0\nn1\nv1");
  EXPECT_THROW_MSG(
    ReadNL(MakeHeader(), "C0\no64\n1\nn0\nv1"),
    ampl::ParseError, "(input):13:1: too few slopes in piecewise-linear term");
}

TEST(NLTest, VarArgExpr) {
  // TODO
}

// TODO: test parsing expressions, expression hierarchy in handler

// TODO: more tests
}
