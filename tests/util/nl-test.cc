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
#include "solvers/util/nl.h"
#include "tests/util.h"

using ampl::NLHeader;
using ampl::ReadError;
using ampl::ReadNLString;

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
    double infinity = std::numeric_limits<double>::infinity();
    if (lb != -infinity && lb != ub)
      log << lb << " <= ";
    log << type << index;
    if (lb == ub)
      log << " = " << ub;
    else if (ub != infinity)
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

  void SetComplement(int con_index, int var_index, int flags) {
    WriteSep().write("c{} complements v{} {};", con_index, var_index, flags);
  }

  class LinearExprHandler {
   private:
    std::string str_;
    fmt::Writer &log_;

   public:
    explicit LinearExprHandler(fmt::Writer &log) : log_(log) {}
    ~LinearExprHandler() { log_ << str_ << ';'; }
    void AddTerm(int var_index, double coef) {
      if (!str_.empty())
        str_ += " + ";
      str_ += fmt::format("{} * v{}", coef, var_index);
    }
  };

  LinearExprHandler GetLinearVarHandler(int index, int num_terms) {
    WriteSep().write("v{} {}: ", index, num_terms);
    return LinearExprHandler(log);
  }
  LinearExprHandler GetLinearObjHandler(int index, int num_terms) {
    WriteSep().write("o{} {}: ", index, num_terms);
    return LinearExprHandler(log);
  }
  LinearExprHandler GetLinearConHandler(int index, int num_terms) {
    WriteSep().write("c{} {}: ", index, num_terms);
    return LinearExprHandler(log);
  }

  class ColumnSizeHandler {
   private:
    fmt::Writer &log_;

   public:
    explicit ColumnSizeHandler(fmt::Writer &log) : log_(log) {}
    ~ColumnSizeHandler() { log_ << ';'; }
    void Add(int offset) { log_ << ' ' << offset; };
  };
  ColumnSizeHandler GetColumnSizeHandler() {
    log << "sizes:";
    return ColumnSizeHandler(log);
  }

  void SetVar(int index, std::string expr, int position) {
    WriteSep().write("v{}/{} = {};", index, position, expr);
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

  void SetInitialValue(int var_index, double value) {
    WriteSep().write("v{} := {};", var_index, value);
  }

  void SetInitialDualValue(int con_index, double value) {
    WriteSep().write("c{} := {};", con_index, value);
  }

  void SetFunction(
      int index, fmt::StringRef name, int num_args, ampl::func::Type type) {
    WriteSep().write("f{}: {} {} {};", index,
                     std::string(name.c_str(), name.size()), num_args, type);
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
      ReadError, "(input):11:3: expected newline");
}

TEST(NLTest, InvalidFormat) {
  EXPECT_THROW_MSG(ReadHeader(0, "x"),
      ReadError, "(input):1:1: expected format specifier");
}

TEST(NLTest, InvalidNumOptions) {
  EXPECT_EQ(0, ReadHeader(0, "ga").num_options);
  EXPECT_EQ(0, ReadHeader(0, "g-1").num_options);
  EXPECT_THROW_MSG(ReadHeader(0, "g10"),
      ReadError, "(input):1:2: too many options");
  EXPECT_THROW_MSG(ReadHeader(0,
      fmt::format("g{}", static_cast<unsigned>(INT_MAX) + 1)),
      ReadError, "(input):1:2: number is too big");
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
      ReadError, "(input):6:6: unrecognized binary format");
  // TODO: check if the bytes are actually swapped
}

TEST(NLTest, IncompleteHeader) {
  ReadHeader(0, "g");
  EXPECT_THROW_MSG(
      ReadHeader(0, "\n"),
      ReadError, "(input):1:1: expected format specifier");
  ReadHeader(1, " 1 0 0");
  EXPECT_THROW_MSG(
      ReadHeader(1, " 1 0"),
      ReadError, "(input):2:5: expected nonnegative integer");
  for (int i = 2; i <= 8; ++i) {
    if (i == 6)
      continue;
    ReadHeader(i, " 0 0");
    EXPECT_THROW_MSG(
        ReadHeader(i, " 0"), ReadError,
        fmt::format("(input):{}:3: expected nonnegative integer", i + 1));
  }
  for (int i = 6; i <= 9; i += 3) {
    ReadHeader(1, " 0 0 0 0 0");
    EXPECT_THROW_MSG(
        ReadHeader(i, " 0 0 0 0"), ReadError,
        fmt::format("(input):{}:9: expected nonnegative integer", i + 1));
  }
  std::string input = ReplaceLine(FormatHeader(NLHeader()), 4, " 0 0");
  ReadHeader(ReplaceLine(input, 6, " 0 0"));
  EXPECT_THROW_MSG(
      ReadHeader(ReplaceLine(input, 6, " 0")),
      ReadError, "(input):7:3: expected nonnegative integer");
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
    ReadError, "(input):11:2: expected nonnegative integer");
  ReadNL(header, "O0 9\nn0\n");
  EXPECT_THROW_MSG(
    ReadNL(header, "O10 0\nn0\n"),
    ReadError, "(input):11:2: objective index 10 out of bounds");
}

TEST(NLTest, ObjType) {
  NLHeader header = {};
  header.num_vars = header.num_objs = 1;
  ReadNL(header, "O0 0\nn0\n");
  ReadNL(header, "O0 1\nn0\n");
  ReadNL(header, "O0 10\nn0\n");
  EXPECT_THROW_MSG(
    ReadNL(header, "O0 -1\nn0\n"),
    ReadError, "(input):11:4: expected nonnegative integer");
}

NLHeader MakeHeader() {
  NLHeader header = {};
  header.num_vars = 5;
  header.num_objs = 6;
  header.num_algebraic_cons = 7;
  header.num_logical_cons = 8;
  header.num_funcs = 9;
  header.num_common_exprs_in_objs = 1;
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
  EXPECT_THROW_MSG(ReadNL(MakeHeader(), nl_body), ReadError, error);

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
  EXPECT_READ("c0: v4;", "C0\nv4\n");
  EXPECT_READ("c0: v5;", "C0\nv5\n");
  EXPECT_READ_ERROR("C0\nv-1\n", "(input):12:2: expected nonnegative integer");
  EXPECT_READ_ERROR("C0\nv6\n", "(input):12:2: variable index 6 out of bounds");
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
  EXPECT_READ("c0: <<0; -1, 1>> v1;", "C0\no64\n2\nn-1.0\ns0\nl1\nv1\n");
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
  const char input[] = "C0\nf1 1\nh1:\0\n";
  const char output[] = "c0: f1('\0');";
  EXPECT_READ(std::string(output, sizeof(output) - 1),
              std::string(input, sizeof(input) - 1));
  EXPECT_READ_ERROR("C0\nf1 1\nh3:ab",
        "(input):13:6: unexpected end of file in string");
  EXPECT_READ_ERROR("C0\nf1 1\nh3:a\n",
        "(input):14:1: unexpected end of file in string");
  EXPECT_READ_ERROR("C0\nf1 1\nh3:abc", "(input):13:7: expected newline");
  EXPECT_READ_ERROR("C0\nf1 1\nh3:ab\n", "(input):14:1: expected newline");
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
  void SetComplement(int, int, int) {}

  struct LinearExprHandler {
    void AddTerm(int, double) {}
  };
  LinearExprHandler GetLinearVarHandler(int, int) {
    return LinearExprHandler();
  }
  LinearExprHandler GetLinearObjHandler(int, int) {
    return LinearExprHandler();
  }
  LinearExprHandler GetLinearConHandler(int, int) {
    return LinearExprHandler();
  }

  struct ColumnSizeHandler {
    void Add(int) {}
  };
  ColumnSizeHandler GetColumnSizeHandler() { return ColumnSizeHandler(); }

  void SetVar(int, TestNumericExpr, int) {}
  void SetObj(int, ampl::obj::Type, TestNumericExpr) {}
  void SetCon(int, TestNumericExpr) {}
  void SetLogicalCon(int, TestLogicalExpr) {}

  void SetInitialValue(int, double) {}
  void SetInitialDualValue(int, double) {}

  void SetFunction(int, fmt::StringRef, int, ampl::func::Type) {}

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
  EXPECT_READ("1.1 <= v0; v1 <= 22; v2 = 33; v3; 44 <= v4 <= 55;",
              "b\n2 1.1\n1 22\n4 33\n3\n0 44 55\n");
  EXPECT_READ_ERROR("b\n-1\n", "(input):12:1: expected nonnegative integer");
  EXPECT_READ_ERROR("b\n5 1\n", "(input):12:1: invalid bound type");
  EXPECT_READ_ERROR("b\n2 11\n1 22\n4 33\n3\n",
                    "(input):16:1: expected nonnegative integer");
}

TEST(NLTest, ReadConBounds) {
  EXPECT_READ("1.1 <= c0; c1 <= 22; c2 = 33; c3; 44 <= c4 <= 55; "
              "c5 complements v1 3; c6 complements v4 2;",
              "r\n2 1.1\n1 22\n4 33\n3\n0 44 55\n5 7 2\n5 2 5\n");
  EXPECT_READ_ERROR("r\n-1\n", "(input):12:1: expected nonnegative integer");
  EXPECT_READ_ERROR("r\n6 1\n", "(input):12:1: invalid bound type");
  EXPECT_READ_ERROR("r\n2 11\n1 22\n4 33\n3\n",
                    "(input):16:1: expected nonnegative integer");
  EXPECT_READ_ERROR("r\n5 1 0\n",
                    "(input):12:5: variable index -1 out of bounds");
  EXPECT_READ_ERROR("r\n5 1 6\n",
                    "(input):12:5: variable index 5 out of bounds");
}

TEST(NLTest, ReadLinearObjExpr) {
  EXPECT_READ("o0 2: 1.3 * v1 + 5 * v3;", "G0 2\n1 1.3\n3 5\n");
  EXPECT_READ("o5 5: 1 * v1 + 1 * v2 + 1 * v3 + 1 * v4 + 1 * v5;",
              "G5 5\n1 1\n2 1\n3 1\n4 1\n5 1\n");
  EXPECT_READ_ERROR("G-1", "(input):11:2: expected nonnegative integer");
  EXPECT_READ_ERROR("G6", "(input):11:2: objective index 6 out of bounds");
  EXPECT_READ_ERROR("G0 0",
    "(input):11:4: number of linear terms 0 out of bounds");
  EXPECT_READ_ERROR("G0 6",
    "(input):11:4: number of linear terms 6 out of bounds");
  EXPECT_READ_ERROR("G0 1\n-1 0\n",
    "(input):12:1: expected nonnegative integer");
  EXPECT_READ_ERROR("G0 1\n6 0\n",
    "(input):12:1: variable index 6 out of bounds");
}

TEST(NLTest, ReadLinearConExpr) {
  EXPECT_READ("c0 2: 1.3 * v1 + 5 * v3;", "J0 2\n1 1.3\n3 5\n");
  EXPECT_READ("c5 5: 1 * v1 + 1 * v2 + 1 * v3 + 1 * v4 + 1 * v5;",
              "J5 5\n1 1\n2 1\n3 1\n4 1\n5 1\n");
  EXPECT_READ_ERROR("J-1", "(input):11:2: expected nonnegative integer");
  EXPECT_READ_ERROR("J8", "(input):11:2: constraint index 8 out of bounds");
  EXPECT_READ_ERROR("J0 0",
    "(input):11:4: number of linear terms 0 out of bounds");
  EXPECT_READ_ERROR("J0 6",
    "(input):11:4: number of linear terms 6 out of bounds");
  EXPECT_READ_ERROR("J0 1\n-1 0\n",
    "(input):12:1: expected nonnegative integer");
  EXPECT_READ_ERROR("J0 1\n6 0\n",
    "(input):12:1: variable index 6 out of bounds");
}

TEST(NLTest, ReadJacobianColumns) {
  EXPECT_READ("sizes: 1 2 2 4;", "k4\n1\n3\n5\n9\n");
  EXPECT_READ("sizes: 1 2 2 4;", "K4\n1\n2\n2\n4\n");
  EXPECT_READ_ERROR("k3\n", "(input):11:2: expected 4");
  EXPECT_READ_ERROR("k4\n-1\n", "(input):12:1: expected nonnegative integer");
}

TEST(NLTest, ReadInitialValues) {
  EXPECT_READ("v4 := 1.1; v3 := 0; v2 := 1; v1 := 2; v0 := 3;",
              "x5\n4 1.1\n3 0\n2 1\n1 2\n0 3\n");
  EXPECT_READ_ERROR("x6\n", "(input):11:2: too many initial values");
  EXPECT_READ_ERROR("x1\n-1 0\n", "(input):12:1: expected nonnegative integer");
  EXPECT_READ_ERROR("x1\n5 0\n",
                    "(input):12:1: variable index 5 out of bounds");
  EXPECT_READ_ERROR("x2\n4 1.1\n\n",
                    "(input):13:1: expected nonnegative integer");
}

TEST(NLTest, ReadInitialDualValues) {
  EXPECT_READ("c4 := 1.1; c3 := 0; c2 := 1; c1 := 2; "
              "c0 := 3; c5 := 1; c6 := 2;",
              "d7\n4 1.1\n3 0\n2 1\n1 2\n0 3\n5 1\n6 2\n");
  EXPECT_READ_ERROR("d8\n", "(input):11:2: too many initial values");
  EXPECT_READ_ERROR("d1\n-1 0\n", "(input):12:1: expected nonnegative integer");
  EXPECT_READ_ERROR("d1\n7 0\n",
                    "(input):12:1: constraint index 7 out of bounds");
  EXPECT_READ_ERROR("d2\n4 1.1\n\n",
                    "(input):13:1: expected nonnegative integer");
}

TEST(NLTest, ReadFunction) {
  EXPECT_READ("f0: foo 2 1;", "F0 1 2 foo\n");
  EXPECT_READ("f0:  2 1;", "F0 1 2 \n");
  EXPECT_READ("f0: foo -1 0;", "F0 0 -1 foo\n");
  EXPECT_READ_ERROR("F-1 0 0 f\n",
                    "(input):11:2: expected nonnegative integer");
  EXPECT_READ_ERROR("F9 0 0 f\n",
                    "(input):11:2: function index 9 out of bounds");
  EXPECT_READ_ERROR("F0 -1 0 f\n",
                    "(input):11:4: expected nonnegative integer");
  EXPECT_READ_ERROR("F0 2 0 f\n", "(input):11:4: invalid function type");
}

TEST(NLTest, ReadDefinedVars) {
  EXPECT_READ("v5/1 = b2(v0, 42);", "V5 0 1\no2\nv0\nn42\n");
  EXPECT_READ("v5 2: 2 * v1 + 3 * v0; v5/1 = 0;", "V5 2 1\n1 2.0\n0 3\nn0\n");
  EXPECT_READ_ERROR("V4 0 1\nv0\n",
                    "(input):11:2: defined variable index 4 out of bounds");
  EXPECT_READ_ERROR("V6 0 1\nv0\n",
                    "(input):11:2: defined variable index 6 out of bounds");
}

// TODO: test suffixes

}  // namespace
