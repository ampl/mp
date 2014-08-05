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
using ampl::ParseError;

namespace {

class TestNLHandler {
 private:
  // Writes a comma-separated list.
  template <typename T>
  static void WriteList(fmt::Writer &w, int size, const T *values) {
    for (int i = 0; i < size; ++i) {
      if (i != 0) w << ", ";
      w << values[i];
    }
  }

  std::string MakeVarArg(std::string op, ampl::ArrayRef<std::string> args) {
    fmt::Writer w;
    w << op << '(';
    WriteList(w, args.size(), args.data());
    w << ')';
    return w.str();
  }

  fmt::Writer &WriteSep() {
    if (log.size() != 0)
      log << ' ';
    return log;
  }

  void WriteBounds(char type, int index, double lb, double ub) {
    WriteSep();
    if (lb != -AMPL_INFINITY && lb != ub)
      log << lb << " <= ";
    log << type << index;
    if (lb == ub)
      log << " = " << ub;
    else if (ub != AMPL_INFINITY)
      log << " <= " << ub;
    log << ';';
  }

 public:
  NLHeader header;
  fmt::Writer log;  // Call log.

  typedef std::string Expr;
  typedef std::string NumericExpr;
  typedef std::string LogicalExpr;
  typedef std::string CountExpr;
  typedef std::string Variable;

  void BeginBuild(const char *, const NLHeader &h, int) {
    header = h;
    log.clear();
  }

  void SetVarBounds(int index, double lb, double ub) {
    WriteBounds('v', index, lb, ub);
  }

  void SetConBounds(int index, double lb, double ub) {
    WriteBounds('c', index, lb, ub);
  }

  void SetComplVar(int con_index, int var_index) {
    // TODO
  }

  class LinearObjHandler {
   public:
    void AddTerm(int var_index, double coef) {
      // TODO
    }
  };

  LinearObjHandler GetLinearObjHandler(int index, int num_terms) {
    return LinearObjHandler();
  }

  void SetObj(int index, ampl::obj::Type type, std::string expr) {
    WriteSep() << (type == ampl::obj::MAX ? "maximize" : "minimize")
        << " o" << index << ": " << expr << ";";
  }

  void SetCon(int index, std::string expr) {
    WriteSep() << "c" << index << ": " << expr << ";";
  }

  void SetLogicalCon(int index, std::string expr) {
    WriteSep() << "l" << index << ": " << expr << ";";
  }

  void SetFunction(
      int index, const char *name, int num_args, ampl::func::Type type) {
    // TODO
  }

  std::string MakeNumericConstant(double value) {
    return fmt::format("{}", value);
  }

  std::string MakeVariable(int index) {
    return fmt::format("v{}", index);
  }

  std::string MakeUnary(int opcode, std::string arg) {
    return fmt::format("u{}({})", opcode, arg);
  }

  std::string MakeBinary(int opcode, std::string lhs, std::string rhs) {
    return fmt::format("b{}({}, {})", opcode, lhs, rhs);
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

  std::string MakeCall(int func_index, ampl::ArrayRef<std::string> args) {
    fmt::Writer w;
    w << 'f' << func_index << '(';
    WriteList(w, args.size(), args.data());
    w << ')';
    return w.str();
  }

  std::string MakeVarArg(int opcode, ampl::ArrayRef<std::string> args) {
    return MakeVarArg(fmt::format("v{}", opcode), args);
  }

  std::string MakeSum(ampl::ArrayRef<std::string> args) {
    return MakeVarArg("sum", args);
  }

  std::string MakeCount(ampl::ArrayRef<std::string> args) {
    return MakeVarArg("count", args);
  }

  std::string MakeNumberOf(ampl::ArrayRef<std::string> args) {
    fmt::Writer w;
    w << "numberof " << args[0] << " in (";
    WriteList(w, args.size() - 1, args.data() + 1);
    w << ')';
    return w.str();
  }

  std::string MakeLogicalConstant(bool value) {
    return fmt::format("l{}", value);
  }

  std::string MakeNot(std::string arg) { return fmt::format("not {}", arg); }

  std::string MakeBinaryLogical(int opcode, std::string lhs, std::string rhs) {
    return fmt::format("bl{}({}, {})", opcode, lhs, rhs);
  }

  std::string MakeRelational(int opcode, std::string lhs, std::string rhs) {
    return fmt::format("r{}({}, {})", opcode, lhs, rhs);
  }

  std::string MakeLogicalCount(int opcode, std::string lhs, std::string rhs) {
    return fmt::format("lc{}({}, {})", opcode, lhs, rhs);
  }

  std::string MakeImplication(std::string condition,
                              std::string true_expr, std::string false_expr) {
    return fmt::format("{} ==> {} else {}",
                       condition, true_expr, false_expr);
  }

  std::string MakeIteratedLogical(
      int opcode, ampl::ArrayRef<std::string> args) {
    return MakeVarArg(fmt::format("il{}", opcode), args);
  }

  std::string MakeAllDiff(ampl::ArrayRef<std::string> args) {
    return MakeVarArg("alldiff", args);
  }

  std::string MakeString(fmt::StringRef value) {
    return fmt::format("'{}'", std::string(value.c_str(), value.size()));
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
  ReadNLString(s, handler, "(input)");
  return handler.header;
}

// Reads a zero header with one modified line.
NLHeader ReadHeader(int line_index, fmt::StringRef line) {
  return ReadHeader(ReplaceLine(
      FormatHeader(NLHeader()), line_index, line.c_str()));
}

TEST(NLTest, NoNewlineAtEOF) {
  TestNLHandler handler;
  EXPECT_THROW_MSG(ReadNLString("g\n"
    " 1 1 0\n"
    " 0 0\n"
    " 0 0\n"
    " 0 0 0\n"
    " 0 0 0 1\n"
    " 0 0 0 0 0\n"
    " 0 0\n"
    " 0 0\n"
    " 0 0 0 0 0\n"
    "k0\0deadbeef", handler),
      ParseError, "(input):11:3: expected newline");
}

TEST(NLTest, InvalidFormat) {
  EXPECT_THROW_MSG(ReadHeader(0, "x"),
      ParseError, "(input):1:1: expected format specifier");
}

TEST(NLTest, InvalidNumOptions) {
  EXPECT_EQ(0, ReadHeader(0, "ga").num_options);
  EXPECT_EQ(0, ReadHeader(0, "g-1").num_options);
  EXPECT_THROW_MSG(ReadHeader(0, "g10"),
      ParseError, "(input):1:2: too many options");
  EXPECT_THROW_MSG(ReadHeader(0,
      fmt::format("g{}", static_cast<unsigned>(INT_MAX) + 1)),
      ParseError, "(input):1:2: number is too big");
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
      ParseError, "(input):6:6: unrecognized binary format");
  // TODO: check if the bytes are actually swapped
}

TEST(NLTest, IncompleteHeader) {
  ReadHeader(0, "g");
  EXPECT_THROW_MSG(
      ReadHeader(0, "\n"),
      ParseError, "(input):1:1: expected format specifier");
  ReadHeader(1, " 1 0 0");
  EXPECT_THROW_MSG(
      ReadHeader(1, " 1 0"),
      ParseError, "(input):2:5: expected nonnegative integer");
  for (int i = 2; i <= 8; ++i) {
    if (i == 6)
      continue;
    ReadHeader(i, " 0 0");
    EXPECT_THROW_MSG(
        ReadHeader(i, " 0"), ParseError,
        fmt::format("(input):{}:3: expected nonnegative integer", i + 1));
  }
  for (int i = 6; i <= 9; i += 3) {
    ReadHeader(1, " 0 0 0 0 0");
    EXPECT_THROW_MSG(
        ReadHeader(i, " 0 0 0 0"), ParseError,
        fmt::format("(input):{}:9: expected nonnegative integer", i + 1));
  }
  std::string input = ReplaceLine(FormatHeader(NLHeader()), 4, " 0 0");
  ReadHeader(ReplaceLine(input, 6, " 0 0"));
  EXPECT_THROW_MSG(
      ReadHeader(ReplaceLine(input, 6, " 0")),
      ParseError, "(input):7:3: expected nonnegative integer");
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
    ReadNL(header, "O-1 0\nn0\n"),
    ParseError, "(input):11:2: expected nonnegative integer");
  ReadNL(header, "O0 9\nn0\n");
  EXPECT_THROW_MSG(
    ReadNL(header, "O10 0\nn0\n"),
    ParseError, "(input):11:2: objective index 10 out of bounds");
}

TEST(NLTest, ObjType) {
  NLHeader header = {};
  header.num_vars = header.num_objs = 1;
  ReadNL(header, "O0 0\nn0\n");
  ReadNL(header, "O0 1\nn0\n");
  ReadNL(header, "O0 10\nn0\n");
  EXPECT_THROW_MSG(
    ReadNL(header, "O0 -1\nn0\n"),
    ParseError, "(input):11:4: expected nonnegative integer");
}

NLHeader MakeHeader() {
  NLHeader header = {};
  header.num_vars = 5;
  header.num_objs = 6;
  header.num_algebraic_cons = 7;
  header.num_logical_cons = 8;
  header.num_funcs = 9;
  return header;
}

TEST(NLTest, ReadObjExpr) {
  TestNLHandler handler;
  NLHeader header = MakeHeader();
  header.num_objs = 2;
  ReadNLString(FormatHeader(header) + "O1 0\nn0\n", handler);
  EXPECT_EQ("minimize o1: 0;", handler.log.str());
  ReadNLString(FormatHeader(header) + "O0 1\nv0\n", handler);
  EXPECT_EQ("maximize o0: v0;", handler.log.str());
}

#define EXPECT_READ(expected_output, nl_body) {\
  TestNLHandler handler; \
  ReadNLString(FormatHeader(MakeHeader()) + nl_body, handler); \
  EXPECT_EQ(expected_output, handler.log.str()); \
}

#define EXPECT_READ_ERROR(nl_body, error) \
  EXPECT_THROW_MSG(ReadNL(MakeHeader(), nl_body), ParseError, error);

template <typename Int>
void CheckReadInt(char code) {
  EXPECT_READ("c0: 4;", fmt::format("C0\n{}4.2\n", code));
  Int min = std::numeric_limits<Int>::min();
  EXPECT_READ(fmt::format("c0: {};", min + 0.),
               fmt::format("C0\n{}{}\n", code, min));
  typename safeint::MakeUnsigned<Int>::Type max =
      std::numeric_limits<Int>::max();
  EXPECT_READ(fmt::format("c0: {};", max + 0.),
               fmt::format("C0\n{}{}\n", code, max));
  EXPECT_READ_ERROR(fmt::format("C0\n{}{}\n", code, max + 1),
    "(input):12:2: number is too big");
}

TEST(NLTest, ReadNumericConstant) {
  EXPECT_READ("c0: 4.2;", "C0\nn4.2\n");
  EXPECT_READ("c0: -100;", "C0\nn-1e+2\n");
  CheckReadInt<short>('s');
  CheckReadInt<long>('l');
}

TEST(NLTest, ReadVariable) {
  // TODO: test variable index out of bounds
  EXPECT_READ("c0: v4;", "C0\nv4\n");
  EXPECT_READ_ERROR("C0\nv-1\n", "(input):12:2: expected nonnegative integer");
}

TEST(NLTest, ReadUnaryExpr) {
  EXPECT_READ("c0: u13(v3);", "C0\no13\nv3\n");
}

TEST(NLTest, ReadBinaryExpr) {
  EXPECT_READ("c0: b0(v1, 42);", "C0\no0\nv1\nn42\n");
}

TEST(NLTest, ReadIfExpr) {
  EXPECT_READ("c0: if l1 then v1 else v2;", "C0\no35\nn1\nv1\nv2\n");
}

TEST(NLTest, ReadPiecewiseLinearExpr) {
  EXPECT_READ("c0: <<0; -1, 1>> v1;", "C0\no64\n2\nn-1\ns0\nl1\nv1\n");
  EXPECT_READ_ERROR("C0\no64\n-1\nn0\nv1\n",
    "(input):13:1: expected nonnegative integer");
  EXPECT_READ_ERROR("C0\no64\n1\nn0\nv1\n",
    "(input):13:1: too few slopes in piecewise-linear term");
  EXPECT_READ_ERROR("C0\no64\n2\nv1\nn0\nn1\nv1\n",
    "(input):14:1: expected constant");
  EXPECT_READ_ERROR("C0\no64\n2\nn-1\nv0\nn1\nv1\n",
    "(input):15:1: expected constant");
  EXPECT_READ_ERROR("C0\no64\n2\nn-1\nn0\nn1\nn1\n",
    "(input):17:1: expected variable");
}

TEST(NLTest, ReadCallExpr) {
  EXPECT_READ("c0: f1(v1, 0);", "C0\nf1 2\nv1\nn0\n");
  EXPECT_READ_ERROR("C0\nf-1 1\nn0\n",
    "(input):12:2: expected nonnegative integer");
  EXPECT_READ_ERROR("C0\nf10 1\nn0\n",
    "(input):12:2: function index 10 out of bounds");
  EXPECT_READ_ERROR("C0\nf1 1\nx\n", "(input):13:1: expected expression");
}

TEST(NLTest, ReadVarArgExpr) {
  EXPECT_READ("c0: v11(v4, 5, v1);", "C0\no11\n3\nv4\nn5\nv1\n");
  EXPECT_READ("c0: v12(v4);", "C0\no12\n1\nv4\n");
  EXPECT_READ_ERROR("C0\no12\n0\n" , "(input):13:1: too few arguments");
}

TEST(NLTest, ReadSumExpr) {
  EXPECT_READ("c0: sum(v4, 5, v1);", "C0\no54\n3\nv4\nn5\nv1\n");
  EXPECT_READ_ERROR("C0\no54\n2\nv4\nn5\n", "(input):13:1: too few arguments");
}

TEST(NLTest, ReadCountExpr) {
  EXPECT_READ("c0: count(l1, r24(v1, 42), l0);",
               "C0\no59\n3\nn1\no24\nv1\nn42\nn0\n");
  EXPECT_READ("c0: count(l1);", "C0\no59\n1\nn1\n");
  EXPECT_READ_ERROR("C0\no59\n0\n", "(input):13:1: too few arguments");
}

TEST(NLTest, ReadNumberOfExpr) {
  EXPECT_READ("c0: numberof v4 in (5, v1);", "C0\no60\n3\nv4\nn5\nv1\n");
  EXPECT_READ("c0: numberof v4 in ();", "C0\no60\n1\nv4\n");
  EXPECT_READ_ERROR( "C0\no60\n0\n", "(input):13:1: too few arguments");
}

TEST(NLTest, ReadLogicalConstant) {
  EXPECT_READ("l0: l0;", "L0\nn0\n");
  EXPECT_READ("l0: l1;", "L0\nn1\n");
  EXPECT_READ("l0: l1;", "L0\nn4.2\n");
  EXPECT_READ("l0: l1;", "L0\ns1\n");
  EXPECT_READ("l0: l1;", "L0\nl1\n");
}

TEST(NLTest, ReadNotExpr) {
  EXPECT_READ("l0: not l0;", "L0\no34\nn0\n");
}

TEST(NLTest, ReadBinaryLogicalExpr) {
  EXPECT_READ("l0: bl20(l1, l0);", "L0\no20\nn1\nn0\n");
}

TEST(NLTest, ReadRelationalExpr) {
  EXPECT_READ("l0: r23(v1, 0);", "L0\no23\nv1\nn0\n");
}

TEST(NLTest, ReadLogicalCountExpr) {
  EXPECT_READ("l0: lc63(v1, count(l1));", "L0\no63\nv1\no59\n1\nn1\n");
  EXPECT_READ_ERROR("L0\no63\nv1\nn0\n",
    "(input):14:1: expected count expression");
  EXPECT_READ_ERROR("L0\no63\nv1\no16\nn0\n",
    "(input):14:2: expected count expression opcode");
}

TEST(NLTest, ReadImplicationExpr) {
  EXPECT_READ("l0: l1 ==> l0 else l1;", "L0\no72\nn1\nn0\nn1\n");
}

TEST(NLTest, ReadIteratedLogicalExpr) {
  EXPECT_READ("l0: il71(l1, l0, l1);", "L0\no71\n3\nn1\nn0\nn1\n");
  EXPECT_READ_ERROR("L0\no71\n2\nn1\nn0\n", "(input):13:1: too few arguments");
}

TEST(NLTest, ReadAllDiffExpr) {
  EXPECT_READ("l0: alldiff(v4, 5, v1);", "L0\no74\n3\nv4\nn5\nv1\n");
  EXPECT_READ_ERROR("L0\no74\n2\nv4\nn5\n", "(input):13:1: too few arguments");
}

TEST(NLTest, ReadStringLiteral) {
  EXPECT_READ("c0: f1('');", "C0\nf1 1\nh0:\n");
  EXPECT_READ("c0: f1('abc');", "C0\nf1 1\nh3:abc\n");
  EXPECT_READ("c0: f1('ab\nc');", "C0\nf1 1\nh4:ab\nc\n");
  EXPECT_READ_ERROR("C0\nf1 1\nh3:ab",
        "(input):13:6: unexpected end of file in string");
  EXPECT_READ_ERROR("C0\nf1 1\nh3:a\n",
        "(input):14:1: unexpected end of file in string");
  EXPECT_READ_ERROR("C0\nf1 1\nh3:abc", "(input):13:7: expected newline");
  EXPECT_READ_ERROR( "C0\nf1 1\nh3:ab\n", "(input):14:1: expected newline");
}

TEST(NLTest, ReadInvalidOpCode) {
  EXPECT_READ_ERROR("C0\no-1\n", "(input):12:2: expected nonnegative integer");
  EXPECT_READ_ERROR("C0\no82\n", "(input):12:2: invalid opcode 82");
}

TEST(NLTest, ReadInvalidNumericExpr) {
  EXPECT_READ_ERROR("C0\nx\n", "(input):12:1: expected expression");
  EXPECT_READ_ERROR("C0\no22\nv1\nn0\n",
        "(input):12:2: expected numeric expression opcode");
}

TEST(NLTest, ReadInvalidLogicalExpr) {
  EXPECT_READ_ERROR("L0\nx\n", "(input):12:1: expected logical expression");
  EXPECT_READ_ERROR("L0\no0\nv1\nn0\n",
        "(input):12:2: expected logical expression opcode");
}

struct TestNLHandler2 {
  typedef struct TestExpr {} Expr;
  typedef struct TestNumericExpr : TestExpr {} NumericExpr;
  typedef struct TestVariable : TestNumericExpr {} Variable;
  typedef struct TestCountExpr : TestNumericExpr {} CountExpr;
  typedef struct TestLogicalExpr : TestExpr {} LogicalExpr;

  void BeginBuild(const char *, const NLHeader &, int) {}

  void SetVarBounds(int, double, double) {}
  void SetConBounds(int, double, double) {}
  void SetComplVar(int, int) {}

  struct LinearObjHandler {
    void AddTerm(int, double) {}
  };
  LinearObjHandler GetLinearObjHandler(int, int) { return LinearObjHandler(); }

  void SetObj(int, ampl::obj::Type, TestNumericExpr) {}
  void SetCon(int, TestNumericExpr) {}
  void SetLogicalCon(int, TestLogicalExpr) {}

  void SetFunction(int, const char *, int, ampl::func::Type) {}

  TestNumericExpr MakeNumericConstant(double) { return TestNumericExpr(); }
  TestVariable MakeVariable(int) { return TestVariable(); }
  TestNumericExpr MakeUnary(int, TestNumericExpr) { return TestNumericExpr(); }

  TestNumericExpr MakeBinary(int, TestNumericExpr, TestNumericExpr) {
    return TestNumericExpr();
  }

  TestNumericExpr MakeIf(TestLogicalExpr, TestNumericExpr, TestNumericExpr) {
    return TestNumericExpr();
  }

  TestNumericExpr MakePiecewiseLinear(
      int, const double *, const double *, TestVariable) {
    return TestNumericExpr();
  }

  TestNumericExpr MakeCall(int, ampl::ArrayRef<TestExpr>) {
    return TestNumericExpr();
  }

  TestNumericExpr MakeVarArg(int, ampl::ArrayRef<TestNumericExpr>) {
    return TestNumericExpr();
  }

  TestNumericExpr MakeSum(ampl::ArrayRef<TestNumericExpr>) {
    return TestNumericExpr();
  }

  TestCountExpr MakeCount(ampl::ArrayRef<TestLogicalExpr>) {
    return TestCountExpr();
  }

  TestNumericExpr MakeNumberOf(ampl::ArrayRef<TestNumericExpr>) {
    return TestNumericExpr();
  }

  TestLogicalExpr MakeLogicalConstant(bool) { return TestLogicalExpr(); }
  TestLogicalExpr MakeNot(TestLogicalExpr) { return TestLogicalExpr(); }

  TestLogicalExpr MakeBinaryLogical(int, TestLogicalExpr, TestLogicalExpr) {
    return TestLogicalExpr();
  }

  TestLogicalExpr MakeRelational(int, TestNumericExpr, TestNumericExpr) {
    return TestLogicalExpr();
  }

  TestLogicalExpr MakeLogicalCount(int, TestNumericExpr, TestCountExpr) {
    return TestLogicalExpr();
  }

  TestLogicalExpr MakeImplication(
      TestLogicalExpr, TestLogicalExpr, TestLogicalExpr) {
    return TestLogicalExpr();
  }

  TestLogicalExpr MakeIteratedLogical(int, ampl::ArrayRef<TestLogicalExpr>) {
    return TestLogicalExpr();
  }

  TestLogicalExpr MakeAllDiff(ampl::ArrayRef<TestNumericExpr>) {
    return TestLogicalExpr();
  }

  TestExpr MakeString(fmt::StringRef) { return TestExpr(); }
};

// Test that the .nl reader accepts expression class hierarchy rather than
// a single expression type.
TEST(NLTest, ExprHierarchy) {
  TestNLHandler2 handler;
  ReadNLString(FormatHeader(MakeHeader()), handler);
}

TEST(NLTest, ReadVarBounds) {
  EXPECT_READ("11 <= v0; v1 <= 22; v2 = 33; v3; 44 <= v4 <= 55;",
              "b\n2 11\n1 22\n4 33\n3\n0 44 55\n");
  EXPECT_READ_ERROR("b\n-1\n", "(input):12:1: expected nonnegative integer");
  EXPECT_READ_ERROR("b\n5 1\n", "(input):12:1: invalid bound type");
  EXPECT_READ_ERROR("b\n2 11\n1 22\n4 33\n3\n",
                    "(input):16:1: expected nonnegative integer");
}

TEST(NLTest, ReadConBounds) {
  EXPECT_READ("11 <= c0; c1 <= 22; c2 = 33; c3; 44 <= c4 <= 55; c5 = 0; c6;",
              "r\n2 11\n1 22\n4 33\n3\n0 44 55\n5 0 4\n5 3 1\n");
  EXPECT_READ_ERROR("r\n-1\n", "(input):12:1: expected nonnegative integer");
  EXPECT_READ_ERROR("r\n6 1\n", "(input):12:1: invalid bound type");
  EXPECT_READ_ERROR("r\n2 11\n1 22\n4 33\n3\n",
                    "(input):16:1: expected nonnegative integer");
  // TODO: test complementarity
}

// TODO: test reading linear obj. expression

// TODO: more tests
}
