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
#include <cstring>

#include "gtest-extra.h"
#include "mp/nl.h"
#include "util.h"

using mp::NLHeader;
using mp::ReadError;
using mp::ReadNLString;
using mp::internal::TextReader;
using mp::internal::BinaryReader;
namespace expr = mp::expr;

namespace {

TEST(ReaderBaseTest, ReadChar) {
  mp::internal::ReaderBase rb(" b", "test");
  EXPECT_EQ(' ', rb.ReadChar());
  EXPECT_EQ('b', rb.ReadChar());
  EXPECT_EQ(0, rb.ReadChar());
  TextReader reader("abc\nd", "test");
  reader.ReadTillEndOfLine();
  // ReadTillEndOfLine doesn't change token location
  EXPECT_THROW_MSG(reader.ReportError("oops"), mp::ReadError,
                   "test:1:1: oops");
  // while ReadChar does.
  reader.ReadChar();
  EXPECT_THROW_MSG(reader.ReportError("oops"), mp::ReadError,
                   "test:2:1: oops");
}

TEST(TextReaderTest, ReportError) {
  TextReader reader("x\n", "test");
  EXPECT_THROW_MSG(reader.ReportError("a{}c", 'b'), mp::ReadError,
                   "test:1:1: abc");
}

TEST(TextReaderTest, ReadTillEndOfLine) {
  TextReader reader("ab cde\nfg h", "test");
  reader.ReadTillEndOfLine();
  EXPECT_THROW_MSG(reader.ReadTillEndOfLine(), mp::ReadError,
                   "test:2:5: expected newline");
}

TEST(TextReaderTest, ReadInt) {
  TextReader reader("  11   -22 ", "test");
  EXPECT_EQ(11, reader.ReadInt<int>());
  EXPECT_EQ(-22, reader.ReadInt<long>());
  EXPECT_THROW_MSG(reader.ReadInt<int>(), mp::ReadError,
                   "test:1:12: expected integer");
}

TEST(TextReaderTest, ReadUInt) {
  TextReader reader("  11   -22 ", "test");
  EXPECT_EQ(11, reader.ReadUInt());
  EXPECT_THROW_MSG(reader.ReadUInt(), mp::ReadError,
                   "test:1:8: expected unsigned integer");
}

TEST(TextReaderTest, ReadDouble) {
  TextReader reader("  11   -2.2 ", "test");
  EXPECT_EQ(11, reader.ReadDouble());
  EXPECT_EQ(-2.2, reader.ReadDouble());
  EXPECT_THROW_MSG(reader.ReadDouble(), mp::ReadError,
                   "test:1:13: expected double");
}

TEST(TextReaderTest, ReadName) {
  TextReader reader("  abc  \n  ", "test");
  fmt::StringRef name = reader.ReadName();
  EXPECT_EQ("abc", std::string(name.c_str(), name.size()));
  EXPECT_THROW_MSG(reader.ReadName(), mp::ReadError, "test:1:8: expected name");
  reader.ReadTillEndOfLine();
  EXPECT_THROW_MSG(reader.ReadName(), mp::ReadError, "test:2:3: expected name");
}

TEST(TextReaderTest, ReadString) {
  TextReader reader("  3:a\nb\nc\n2x\n2:de \n2:fg", "test");
  fmt::StringRef name = reader.ReadString();
  EXPECT_EQ("a\nb", std::string(name.c_str(), name.size()));
  EXPECT_THROW_MSG(reader.ReadString(), mp::ReadError,
                   "test:3:1: expected unsigned integer");
  reader.ReadTillEndOfLine();
  EXPECT_THROW_MSG(reader.ReadString(), mp::ReadError,
                   "test:4:2: expected ':'");
  reader.ReadTillEndOfLine();
  EXPECT_THROW_MSG(reader.ReadString(), mp::ReadError,
                   "test:5:5: expected newline");
  reader.ReadTillEndOfLine();
  EXPECT_THROW_MSG(reader.ReadString(), mp::ReadError,
                   "test:6:5: expected newline");
}

// Formats header as a string.
std::string FormatHeader(const NLHeader &h) {
  fmt::Writer w;
  w << h;
  return w.str();
}

// Reads a header from the specified string.
NLHeader ReadHeader(const std::string &s) {
  TextReader reader(s, "(input)");
  NLHeader header = NLHeader();
  reader.ReadHeader(header);
  return header;
}

// Reads a zero header with one modified line.
NLHeader ReadHeader(int line_index, fmt::StringRef line) {
  return ReadHeader(ReplaceLine(
      FormatHeader(NLHeader()), line_index, line.c_str()));
}

TEST(TextReaderTest, InvalidFormat) {
  EXPECT_THROW_MSG(ReadHeader(0, "x"),
      ReadError, "(input):1:1: expected format specifier");
}

TEST(TextReaderTest, InvalidNumOptions) {
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

TEST(TextReaderTest, ReadOptions) {
  const int options[mp::MAX_NL_OPTIONS + 1] = {
      3, 5, 7, 11, 13, 17, 19, 23, 29, 31
  };
  for (int i = 0; i < mp::MAX_NL_OPTIONS; ++i) {
    for (int j = 0; j < mp::MAX_NL_OPTIONS + 1; ++j)
      CheckReadOptions(i, j, options);
  }
  EXPECT_EQ(0, ReadHeader(0, "g").num_options);
}

TEST(TextReaderTest, ReadAMPLVBTol) {
  EXPECT_EQ(4.2, ReadHeader(0, "g2 0 3 4.2").ampl_vbtol);
  EXPECT_EQ(0, ReadHeader(0, "g2 0 0 4.2").ampl_vbtol);
  EXPECT_EQ(0, ReadHeader(0, "g2 0 3").ampl_vbtol);
}

TEST(TextReaderTest, NumComplDblIneq) {
  EXPECT_EQ(42, ReadHeader(2, " 0 0 0 0 42").num_compl_dbl_ineqs);
  EXPECT_EQ(-1, ReadHeader(2, " 0 0 70 0 42").num_compl_dbl_ineqs);
}

TEST(TextReaderTest, ReadArithKind) {
  EXPECT_EQ(NLHeader::TEXT, ReadHeader(5, " 0 0").format);
  EXPECT_EQ(NLHeader::TEXT, ReadHeader(5, " 0 0 0").format);
  EXPECT_EQ(NLHeader::TEXT,
      ReadHeader(5, fmt::format(" 0 0 {}", mp::arith::LAST)).format);
  EXPECT_THROW_MSG(
      ReadHeader(5, fmt::format(" 0 0 {}", mp::arith::LAST + 1)),
      ReadError, "(input):6:6: unknown floating-point arithmetic kind");
}

TEST(TextReaderTest, IncompleteHeader) {
  ReadHeader(0, "g");
  EXPECT_THROW_MSG(
      ReadHeader(0, "\n"),
      ReadError, "(input):1:1: expected format specifier");
  ReadHeader(1, " 1 0 0");
  EXPECT_THROW_MSG(
      ReadHeader(1, " 1 0"),
      ReadError, "(input):2:5: expected unsigned integer");
  for (int i = 2; i <= 8; ++i) {
    if (i == 6)
      continue;
    ReadHeader(i, " 0 0");
    EXPECT_THROW_MSG(
        ReadHeader(i, " 0"), ReadError,
        fmt::format("(input):{}:3: expected unsigned integer", i + 1));
  }
  for (int i = 6; i <= 9; i += 3) {
    ReadHeader(1, " 0 0 0 0 0");
    EXPECT_THROW_MSG(
        ReadHeader(i, " 0 0 0 0"), ReadError,
        fmt::format("(input):{}:9: expected unsigned integer", i + 1));
  }
  std::string input = ReplaceLine(FormatHeader(NLHeader()), 4, " 0 0");
  ReadHeader(ReplaceLine(input, 6, " 0 0"));
  EXPECT_THROW_MSG(
      ReadHeader(ReplaceLine(input, 6, " 0")),
      ReadError, "(input):7:3: expected unsigned integer");
}

#define CHECK_INT_OVERFLOW(field, col) { \
  NLHeader h = NLHeader(); \
  h.format = NLHeader::TEXT; \
  h.num_vars = INT_MAX; \
  h.field = 1; \
  fmt::Writer w; \
  w << h; \
  NLHeader actual = NLHeader(); \
  EXPECT_THROW_MSG(TextReader(w.str(), "in").ReadHeader(actual), \
                   ReadError, fmt::format("in:10:{}: integer overflow", col)); \
}

TEST(TextReaderTest, ReadHeaderIntegerOverflow) {
  CHECK_INT_OVERFLOW(num_common_exprs_in_both, 2);
  CHECK_INT_OVERFLOW(num_common_exprs_in_cons, 4);
  CHECK_INT_OVERFLOW(num_common_exprs_in_objs, 6);
  CHECK_INT_OVERFLOW(num_common_exprs_in_single_cons, 8);
  CHECK_INT_OVERFLOW(num_common_exprs_in_single_objs, 10);
}

class TestReaderBase : public mp::internal::ReaderBase {
 public:
   TestReaderBase(fmt::StringRef data, fmt::StringRef name)
     : ReaderBase(data, name) {}
};

class TestBinaryReader : private TestReaderBase, public BinaryReader {
 public:
  TestBinaryReader(fmt::StringRef data, fmt::StringRef name = "test")
    : TestReaderBase(data, name),
      BinaryReader(static_cast<TestReaderBase&>(*this)) {}
};

TEST(BinaryReaderTest, ReportError) {
  TestBinaryReader reader("x");
  EXPECT_THROW_MSG(reader.ReportError("a{}c", 'b'), mp::BinaryReadError,
                   "test:offset 0: abc");
}

TEST(BinaryReaderTest, ReadInt) {
  int data[] = {11, -22};
  std::size_t size = sizeof(data);
  TestBinaryReader reader(fmt::StringRef(reinterpret_cast<char*>(data), size));
  EXPECT_EQ(11, reader.ReadInt<int>());
  EXPECT_EQ(-22, reader.ReadInt<int>());
  EXPECT_THROW_MSG(reader.ReadInt<int>(), mp::BinaryReadError,
                   fmt::format("test:offset {}: unexpected end of file", size));
}

TEST(BinaryReaderTest, ReadLong) {
  long data[] = {42};
  std::size_t size = sizeof(data);
  TestBinaryReader reader(fmt::StringRef(reinterpret_cast<char*>(data), size));
  EXPECT_EQ(42, reader.ReadInt<long>());
}

TEST(BinaryReaderTest, ReadUInt) {
  int data[] = {11, -22};
  std::size_t size = sizeof(data);
  TestBinaryReader reader(fmt::StringRef(reinterpret_cast<char*>(data), size));
  EXPECT_EQ(11, reader.ReadUInt());
  std::string message =
      fmt::format("test:offset {}: expected unsigned integer", sizeof(int));
  EXPECT_THROW_MSG(reader.ReadUInt(), mp::BinaryReadError, message);
}

TEST(BinaryReaderTest, ReadDouble) {
  double data[] = {11, -2.2};
  std::size_t size = sizeof(data);
  TestBinaryReader reader(fmt::StringRef(reinterpret_cast<char*>(data), size));
  EXPECT_EQ(11, reader.ReadDouble());
  EXPECT_EQ(-2.2, reader.ReadDouble());
  EXPECT_THROW_MSG(reader.ReadDouble(), mp::BinaryReadError,
                   fmt::format("test:offset {}: unexpected end of file", size));
}

void TestReadString(fmt::StringRef (BinaryReader::*read)()) {
  int data[] = {3, 0, 0};
  std::strcpy(reinterpret_cast<char*>(data + 1), "abc");
  {
    TestBinaryReader reader(
          fmt::StringRef(reinterpret_cast<char*>(data), sizeof(int) + 3));
    fmt::StringRef name = (reader.*read)();
    EXPECT_EQ("abc", std::string(name.c_str(), name.size()));
  }
  {
    std::size_t size = sizeof(int) + 2;
    TestBinaryReader reader(
          fmt::StringRef(reinterpret_cast<char*>(data), size));
    std::string message =
        fmt::format("test:offset {}: unexpected end of file", size);
    EXPECT_THROW_MSG((reader.*read)(), mp::BinaryReadError, message);
  }
  data[0] = -1;
  {
    TestBinaryReader reader(
          fmt::StringRef(reinterpret_cast<char*>(data), sizeof(data)));
    EXPECT_THROW_MSG((reader.*read)(), mp::BinaryReadError,
                     "test:offset 0: expected unsigned integer");
  }
}

TEST(BinaryReaderTest, ReadName) {
  TestReadString(&BinaryReader::ReadName);
}

TEST(BinaryReaderTest, ReadString) {
  TestReadString(&BinaryReader::ReadString);
}

TEST(NLTest, ArithKind) {
  namespace arith = mp::arith;
  EXPECT_GE(arith::GetKind(), arith::UNKNOWN);
  EXPECT_LE(arith::GetKind(), arith::LAST);
  EXPECT_TRUE(arith::IsIEEE(arith::IEEE_BIG_ENDIAN));
  EXPECT_TRUE(arith::IsIEEE(arith::IEEE_BIG_ENDIAN));
  EXPECT_FALSE(arith::IsIEEE(arith::UNKNOWN));
  EXPECT_FALSE(arith::IsIEEE(arith::IBM));
  EXPECT_FALSE(arith::IsIEEE(arith::VAX));
  EXPECT_FALSE(arith::IsIEEE(arith::CRAY));
}

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

  std::string MakeVarArg(std::string op, const std::vector<std::string> &args) {
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
  fmt::Writer log;  // Call log.

  typedef std::string Expr;
  typedef std::string NumericExpr;
  typedef std::string LogicalExpr;
  typedef std::string CountExpr;
  typedef std::string Variable;

  void SetHeader(const NLHeader &) { log.clear(); }

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
  typedef LinearExprHandler LinearObjHandler;
  typedef LinearExprHandler LinearConHandler;

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
    void Add(int offset) { log_ << ' ' << offset; }
  };
  ColumnSizeHandler GetColumnSizeHandler() {
    log << "sizes:";
    return ColumnSizeHandler(log);
  }

  void SetVar(int index, std::string expr, int position) {
    WriteSep().write("v{}/{} = {};", index, position, expr);
  }

  void SetObj(int index, mp::obj::Type type, std::string expr) {
    WriteSep() << (type == mp::obj::MAX ? "maximize" : "minimize")
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
      int index, fmt::StringRef name, int num_args, mp::func::Type type) {
    WriteSep().write("f{}: {} {} {};", index,
                     std::string(name.c_str(), name.size()), num_args, type);
  }

  class SuffixHandler {
   private:
    fmt::Writer &log_;
    bool first_;

   public:
    explicit SuffixHandler(fmt::Writer &log) : log_(log), first_(true) {}
    ~SuffixHandler() { log_ << ';'; }

    void SetValue(int index, int value) {
      if (!first_)
        log_ << ',';
      first_ = false;
      log_.write(" i{} = {}", index, value);
    }
    void SetValue(int index, double value) {
      if (!first_)
        log_ << ',';
      first_ = false;
      log_.write(" d{} = {}", index, value);
    }
  };

  SuffixHandler AddSuffix(int kind, int num_values, fmt::StringRef name) {
    WriteSep().write("suffix {}:{}:{}:", name, kind, num_values);
    return SuffixHandler(log);
  }

  std::string MakeNumericConstant(double value) {
    return fmt::format("{}", value);
  }

  std::string MakeVariable(int index) {
    return fmt::format("v{}", index);
  }

  std::string MakeUnary(expr::Kind kind, std::string arg) {
    return fmt::format("u{}({})", opcode(kind), arg);
  }

  std::string MakeBinary(expr::Kind kind, std::string lhs, std::string rhs) {
    return fmt::format("b{}({}, {})", opcode(kind), lhs, rhs);
  }

  std::string MakeIf(std::string condition,
                     std::string true_expr, std::string false_expr) {
    return fmt::format("if {} then {} else {}",
                       condition, true_expr, false_expr);
  }

  class PLTermHandler {
   private:
    std::vector<double> slopes_;
    std::vector<double> breakpoints_;

    friend class TestNLHandler;

   public:
    void AddSlope(double slope) { slopes_.push_back(slope); }
    void AddBreakpoint(double breakpoint) {
      breakpoints_.push_back(breakpoint);
    }
  };

  PLTermHandler BeginPLTerm(int) {
    return PLTermHandler();
  }
  std::string EndPLTerm(PLTermHandler h, std::string var) {
    fmt::Writer w;
    w << "<<";
    WriteList(w, h.breakpoints_.size(), h.breakpoints_.data());
    w << "; ";
    WriteList(w, h.slopes_.size(), h.slopes_.data());
    w << ">> " << var;
    return w.str();
  }

  class CallArgHandler {
   private:
    int func_index_;
    std::vector<std::string> args_;
    friend class TestNLHandler;

   public:
    explicit CallArgHandler(int func_index) : func_index_(func_index) {}
    void AddArg(std::string arg) { args_.push_back(arg); }
  };

  CallArgHandler BeginCall(int func_index, int) {
    return CallArgHandler(func_index);
  }

  std::string EndCall(CallArgHandler h) {
    fmt::Writer w;
    w << 'f' << h.func_index_ << '(';
    WriteList(w, h.args_.size(), h.args_.data());
    w << ')';
    return w.str();
  }

  class ArgHandler {
   private:
    std::string name_;
    std::vector<std::string> args_;
    friend class TestNLHandler;

   public:
    explicit ArgHandler(std::string name) : name_(name) {}
    void AddArg(std::string arg) { args_.push_back(arg); }
  };

  typedef ArgHandler NumericArgHandler;
  typedef ArgHandler LogicalArgHandler;

  ArgHandler BeginVarArg(expr::Kind kind, int) {
    return ArgHandler(fmt::format("v{}", opcode(kind)));
  }
  std::string EndVarArg(ArgHandler h) { return MakeVarArg(h.name_, h.args_); }

  ArgHandler BeginSum(int) { return ArgHandler("sum"); }
  std::string EndSum(ArgHandler h) { return MakeVarArg(h.name_, h.args_); }

  ArgHandler BeginCount(int) { return ArgHandler("count"); }
  std::string EndCount(ArgHandler h) { return MakeVarArg(h.name_, h.args_); }

  ArgHandler BeginNumberOf(int) { return ArgHandler("numberof"); }
  std::string EndNumberOf(ArgHandler h) {
    fmt::Writer w;
    w << "numberof " << h.args_[0] << " in (";
    WriteList(w, h.args_.size() - 1, h.args_.data() + 1);
    w << ')';
    return w.str();
  }

  std::string MakeLogicalConstant(bool value) {
    return fmt::format("l{}", value);
  }

  std::string MakeNot(std::string arg) { return fmt::format("not {}", arg); }

  std::string MakeBinaryLogical(
      expr::Kind kind, std::string lhs, std::string rhs) {
    return fmt::format("bl{}({}, {})", opcode(kind), lhs, rhs);
  }

  std::string MakeRelational(
      expr::Kind kind, std::string lhs, std::string rhs) {
    return fmt::format("r{}({}, {})", opcode(kind), lhs, rhs);
  }

  std::string MakeLogicalCount(
      expr::Kind kind, std::string lhs, std::string rhs) {
    return fmt::format("lc{}({}, {})", opcode(kind), lhs, rhs);
  }

  std::string MakeImplication(std::string condition,
                              std::string true_expr, std::string false_expr) {
    return fmt::format("{} ==> {} else {}",
                       condition, true_expr, false_expr);
  }

  ArgHandler BeginIteratedLogical(expr::Kind kind, int) {
    return ArgHandler(fmt::format("il{}", opcode(kind)));
  }
  std::string EndIteratedLogical(ArgHandler h) {
    return MakeVarArg(h.name_, h.args_);
  }

  typedef ArgHandler AllDiffArgHandler;

  ArgHandler BeginAllDiff(int) { return ArgHandler("alldiff"); }
  std::string EndAllDiff(ArgHandler h) { return MakeVarArg(h.name_, h.args_); }

  std::string MakeStringLiteral(fmt::StringRef value) {
    return fmt::format("'{}'", std::string(value.c_str(), value.size()));
  }
};

TEST(NLTest, WriteTextHeader) {
  NLHeader h = MakeTestHeader();
  fmt::Writer w;
  w << h;
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
  NLHeader header = NLHeader();
  header.format = NLHeader::BINARY;
  header.num_options = 3;
  for (int i = 0; i < header.num_options; ++i)
    header.options[i] = 11 * (i + 1);
  header.arith_kind = mp::arith::CRAY;
  fmt::Writer w;
  w << header;
  EXPECT_EQ(
      "b3 11 22 33\n"
      " 0 0 0 0 0 0\n"
      " 0 0 0 0 0 0\n"
      " 0 0\n"
      " 0 0 0\n"
      " 0 0 5 0\n"
      " 0 0 0 0 0\n"
      " 0 0\n"
      " 0 0\n"
      " 0 0 0 0 0\n",
      w.str());
}

struct TestNLHandler2 {
  typedef struct TestExpr {} Expr;
  typedef struct TestNumericExpr : TestExpr {} NumericExpr;
  typedef struct TestVariable : TestNumericExpr {} Variable;
  typedef struct TestCountExpr : TestNumericExpr {} CountExpr;
  typedef struct TestLogicalExpr : TestExpr {} LogicalExpr;

  void SetHeader(const NLHeader &) {}

  void SetVarBounds(int, double, double) {}
  void SetConBounds(int, double, double) {}
  void SetComplement(int, int, int) {}

  struct LinearObjHandler {
    void AddTerm(int, double) {}
  };
  struct LinearConHandler {
    void AddTerm(int, double) {}
  };
  LinearConHandler GetLinearVarHandler(int, int) {
    return LinearConHandler();
  }
  LinearObjHandler GetLinearObjHandler(int, int) {
    return LinearObjHandler();
  }
  LinearConHandler GetLinearConHandler(int, int) {
    return LinearConHandler();
  }

  struct ColumnSizeHandler {
    void Add(int) {}
  };
  ColumnSizeHandler GetColumnSizeHandler() { return ColumnSizeHandler(); }

  void SetVar(int, TestNumericExpr, int) {}
  void SetObj(int, mp::obj::Type, TestNumericExpr) {}
  void SetCon(int, TestNumericExpr) {}
  void SetLogicalCon(int, TestLogicalExpr) {}

  void SetInitialValue(int, double) {}
  void SetInitialDualValue(int, double) {}

  void SetFunction(int, fmt::StringRef, int, mp::func::Type) {}

  struct SuffixHandler {
    void SetValue(int, double) {}
  };
  SuffixHandler AddSuffix(int, int, fmt::StringRef) { return SuffixHandler(); }

  TestNumericExpr MakeNumericConstant(double) { return TestNumericExpr(); }
  TestVariable MakeVariable(int) { return TestVariable(); }
  TestNumericExpr MakeUnary(expr::Kind, TestNumericExpr) {
    return TestNumericExpr();
  }

  TestNumericExpr MakeBinary(expr::Kind, TestNumericExpr, TestNumericExpr) {
    return TestNumericExpr();
  }

  TestNumericExpr MakeIf(TestLogicalExpr, TestNumericExpr, TestNumericExpr) {
    return TestNumericExpr();
  }

  struct PLTermHandler {
    void AddSlope(double) {}
    void AddBreakpoint(double) {}
  };

  PLTermHandler BeginPLTerm(int) { return PLTermHandler(); }
  TestNumericExpr EndPLTerm(PLTermHandler, TestVariable) {
    return TestNumericExpr();
  }

  struct CallArgHandler {
    void AddArg(TestNumericExpr) {}
    void AddArg(TestExpr) {}
  };

  CallArgHandler BeginCall(int, int) { return CallArgHandler(); }
  TestNumericExpr EndCall(CallArgHandler) { return TestNumericExpr(); }

  struct NumericArgHandler {
    void AddArg(NumericExpr) {}
  };

  struct LogicalArgHandler {
    void AddArg(LogicalExpr) {}
  };

  NumericArgHandler BeginVarArg(expr::Kind, int) {
    return NumericArgHandler();
  }
  TestNumericExpr EndVarArg(NumericArgHandler) {
    return TestNumericExpr();
  }

  NumericArgHandler BeginSum(int) { return NumericArgHandler(); }
  TestNumericExpr EndSum(NumericArgHandler) { return TestNumericExpr(); }

  LogicalArgHandler BeginCount(int) { return LogicalArgHandler(); }
  TestCountExpr EndCount(LogicalArgHandler) { return TestCountExpr(); }

  NumericArgHandler BeginNumberOf(int) { return NumericArgHandler(); }
  TestNumericExpr EndNumberOf(NumericArgHandler) { return TestNumericExpr(); }

  TestLogicalExpr MakeLogicalConstant(bool) { return TestLogicalExpr(); }
  TestLogicalExpr MakeNot(TestLogicalExpr) { return TestLogicalExpr(); }

  TestLogicalExpr MakeBinaryLogical(
      expr::Kind, TestLogicalExpr, TestLogicalExpr) {
    return TestLogicalExpr();
  }

  TestLogicalExpr MakeRelational(expr::Kind, TestNumericExpr, TestNumericExpr) {
    return TestLogicalExpr();
  }

  TestLogicalExpr MakeLogicalCount(expr::Kind, TestNumericExpr, TestCountExpr) {
    return TestLogicalExpr();
  }

  TestLogicalExpr MakeImplication(
      TestLogicalExpr, TestLogicalExpr, TestLogicalExpr) {
    return TestLogicalExpr();
  }

  LogicalArgHandler BeginIteratedLogical(expr::Kind, int) {
    return LogicalArgHandler();
  }
  TestLogicalExpr EndIteratedLogical(LogicalArgHandler) {
    return TestLogicalExpr();
  }

  struct AllDiffArgHandler {
    void AddArg(NumericExpr) {}
  };

  AllDiffArgHandler BeginAllDiff(int) { return AllDiffArgHandler(); }
  TestLogicalExpr EndAllDiff(AllDiffArgHandler) { return TestLogicalExpr(); }

  TestExpr MakeStringLiteral(fmt::StringRef) { return TestExpr(); }
};

NLHeader MakeHeader() {
  NLHeader header = NLHeader();
  header.num_vars = 5;
  header.num_objs = 6;
  header.num_algebraic_cons = 7;
  header.num_logical_cons = 8;
  header.num_funcs = 9;
  header.num_common_exprs_in_objs = 1;
  return header;
}

// Test that the .nl reader accepts expression class hierarchy rather than
// a single expression type.
TEST(NLTest, ExprHierarchy) {
  TestNLHandler2 handler;
  ReadNLString(FormatHeader(MakeHeader()), handler);
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

std::string ReadNL(std::string body) {
  TestNLHandler handler;
  ReadNLString(FormatHeader(MakeHeader()) + body, handler);
  return handler.log.str();
}

#define EXPECT_READ(expected_output, nl_body) \
  EXPECT_EQ(expected_output, ReadNL(nl_body))

#define EXPECT_READ_ERROR(nl_body, error) \
  EXPECT_THROW_MSG(ReadNL(nl_body), ReadError, error)

TEST(NLTest, ReadObj) {
  EXPECT_READ("minimize o1: 0;", "O1 0\nn0\n");
  EXPECT_READ("maximize o0: v0;", "O0 1\nv0\n");
  EXPECT_READ("maximize o5: v0;", "O5 10\nv0\n");
  EXPECT_READ_ERROR("O0 -1\nn0\n", "(input):11:4: expected unsigned integer");
  EXPECT_READ_ERROR("O-1 0\nn0\n", "(input):11:2: expected unsigned integer");
  EXPECT_READ_ERROR("O6 0\nn0\n", "(input):11:2: integer 6 out of bounds");
}

template <typename Int>
void CheckReadInt(char code) {
  EXPECT_READ("c0: 4;", fmt::format("C0\n{}4.2\n", code));
  Int min = std::numeric_limits<Int>::min();
  EXPECT_READ(fmt::format("c0: {};", min + 0.),
               fmt::format("C0\n{}{}\n", code, min));
  typename mp::MakeUnsigned<Int>::Type max = std::numeric_limits<Int>::max();
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
  EXPECT_READ_ERROR("C0\nv-1\n", "(input):12:2: expected unsigned integer");
  EXPECT_READ_ERROR("C0\nv6\n", "(input):12:2: integer 6 out of bounds");
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
    "(input):13:1: expected unsigned integer");
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
    "(input):12:2: expected unsigned integer");
  EXPECT_READ_ERROR("C0\nf10 1\nn0\n",
                    "(input):12:2: integer 10 out of bounds");
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
    "(input):14:2: expected count expression");
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
  EXPECT_READ("l0: alldiff(v4);", "L0\no74\n1\nv4\n");
  EXPECT_READ_ERROR("L0\no74\n0\n", "(input):13:1: too few arguments");
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
  EXPECT_READ_ERROR("C0\no-1\n", "(input):12:2: expected unsigned integer");
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

TEST(NLTest, ReadVarBounds) {
  EXPECT_READ("1.1 <= v0; v1 <= 22; v2 = 33; v3; 44 <= v4 <= 55;",
              "b\n21.1\n1 22\n4 33\n3\n0 44 55\n");
  EXPECT_READ_ERROR("b\n-1\n", "(input):12:1: expected bound");
  EXPECT_READ_ERROR("b\n5 1\n", "(input):12:1: expected bound");
  EXPECT_READ_ERROR("b\n2 11\n1 22\n4 33\n3\n", "(input):16:1: expected bound");
}

TEST(NLTest, ReadConBounds) {
  EXPECT_READ("1.1 <= c0; c1 <= 22; c2 = 33; c3; 44 <= c4 <= 55; "
              "c5 complements v1 3; c6 complements v4 2;",
              "r\n21.1\n1 22\n4 33\n3\n0 44 55\n5 7 2\n5 2 5\n");
  EXPECT_READ_ERROR("r\n-1\n", "(input):12:1: expected bound");
  EXPECT_READ_ERROR("r\n6 1\n", "(input):12:1: expected bound");
  EXPECT_READ_ERROR("r\n2 11\n1 22\n4 33\n3\n", "(input):16:1: expected bound");
  EXPECT_READ_ERROR("r\n5 1 0\n", "(input):12:5: integer 0 out of bounds");
  EXPECT_READ_ERROR("r\n5 1 6\n", "(input):12:5: integer 6 out of bounds");
  // Check that there is no overflow for largest possible var index.
  TestNLHandler handler;
  NLHeader header = NLHeader();
  header.num_vars = INT_MAX;
  header.num_algebraic_cons = 1;
  ReadNLString(FormatHeader(header) + fmt::format("r\n5 1 {}\n", INT_MAX),
               handler);
  EXPECT_EQ(fmt::format("c0 complements v{} 1;", INT_MAX - 1),
            handler.log.str());
}

TEST(NLTest, ReadLinearObjExpr) {
  EXPECT_READ("o0 2: 1.3 * v1 + 5 * v3;", "G0 2\n1 1.3\n3 5\n");
  EXPECT_READ("o5 5: 1 * v1 + 1 * v2 + 1 * v3 + 1 * v4 + 1 * v5;",
              "G5 5\n1 1\n2 1\n3 1\n4 1\n5 1\n");
  EXPECT_READ_ERROR("G-1", "(input):11:2: expected unsigned integer");
  EXPECT_READ_ERROR("G6", "(input):11:2: integer 6 out of bounds");
  EXPECT_READ_ERROR("G0 0", "(input):11:4: integer 0 out of bounds");
  EXPECT_READ_ERROR("G0 7", "(input):11:4: integer 7 out of bounds");
  EXPECT_READ_ERROR("G0 1\n-1 0\n", "(input):12:1: expected unsigned integer");
  EXPECT_READ_ERROR("G0 1\n6 0\n", "(input):12:1: integer 6 out of bounds");
}

TEST(NLTest, ReadLinearConExpr) {
  EXPECT_READ("c0 2: 1.3 * v1 + 5 * v3;", "J0 2\n1 1.3\n3 5\n");
  EXPECT_READ("c5 5: 1 * v1 + 1 * v2 + 1 * v3 + 1 * v4 + 1 * v5;",
              "J5 5\n1 1\n2 1\n3 1\n4 1\n5 1\n");
  EXPECT_READ_ERROR("J-1", "(input):11:2: expected unsigned integer");
  EXPECT_READ_ERROR("J8", "(input):11:2: integer 8 out of bounds");
  EXPECT_READ_ERROR("J0 0", "(input):11:4: integer 0 out of bounds");
  EXPECT_READ_ERROR("J0 7", "(input):11:4: integer 7 out of bounds");
  EXPECT_READ_ERROR("J0 1\n-1 0\n", "(input):12:1: expected unsigned integer");
  EXPECT_READ_ERROR("J0 1\n6 0\n", "(input):12:1: integer 6 out of bounds");
}

TEST(NLTest, ReadColumnSizes) {
  EXPECT_READ("sizes: 1 2 2 4;", "k4\n1\n3\n5\n9\n");
  EXPECT_READ("sizes: 1 2 2 4;", "K4\n1\n2\n2\n4\n");
  EXPECT_READ_ERROR("k3\n", "(input):11:2: expected 4");
  EXPECT_READ_ERROR("k4\n-1\n", "(input):12:1: expected unsigned integer");
  EXPECT_READ_ERROR("k4\n2\n1\n", "(input):13:1: invalid column offset");
}

TEST(NLTest, ReadInitialValues) {
  EXPECT_READ("v4 := 1.1; v3 := 0; v2 := 1; v1 := 2; v0 := 3;",
              "x5\n4 1.1\n3 0\n2 1\n1 2\n0 3\n");
  EXPECT_READ_ERROR("x6\n", "(input):11:2: too many initial values");
  EXPECT_READ_ERROR("x1\n-1 0\n", "(input):12:1: expected unsigned integer");
  EXPECT_READ_ERROR("x1\n5 0\n", "(input):12:1: integer 5 out of bounds");
  EXPECT_READ_ERROR("x2\n4 1.1\n\n", "(input):13:1: expected unsigned integer");
}

TEST(NLTest, ReadInitialDualValues) {
  EXPECT_READ("c4 := 1.1; c3 := 0; c2 := 1; c1 := 2; "
              "c0 := 3; c5 := 1; c6 := 2;",
              "d7\n4 1.1\n3 0\n2 1\n1 2\n0 3\n5 1\n6 2\n");
  EXPECT_READ_ERROR("d8\n", "(input):11:2: too many initial values");
  EXPECT_READ_ERROR("d1\n-1 0\n", "(input):12:1: expected unsigned integer");
  EXPECT_READ_ERROR("d1\n7 0\n", "(input):12:1: integer 7 out of bounds");
  EXPECT_READ_ERROR("d2\n4 1.1\n\n", "(input):13:1: expected unsigned integer");
}

TEST(NLTest, ReadFunction) {
  EXPECT_READ("f0: foo 2 1;", "F0 1 2 foo\n");
  EXPECT_READ("f0: foo -1 0;", "F0 0 -1 foo\n");
  EXPECT_READ_ERROR("F0 1 2 \n", "(input):11:8: expected name");
  EXPECT_READ_ERROR("F-1 0 0 f\n", "(input):11:2: expected unsigned integer");
  EXPECT_READ_ERROR("F9 0 0 f\n", "(input):11:2: integer 9 out of bounds");
  EXPECT_READ_ERROR("F0 -1 0 f\n", "(input):11:4: expected unsigned integer");
  EXPECT_READ_ERROR("F0 2 0 f\n", "(input):11:4: invalid function type");
}

TEST(NLTest, ReadDefinedVars) {
  EXPECT_READ("v5/1 = b2(v0, 42);", "V5 0 1\no2\nv0\nn42\n");
  EXPECT_READ("v5 2: 2 * v1 + 3 * v0; v5/1 = 0;", "V5 2 1\n1 2.0\n0 3\nn0\n");
  EXPECT_READ_ERROR("V4 0 1\nv0\n", "(input):11:2: integer 4 out of bounds");
  EXPECT_READ_ERROR("V6 0 1\nv0\n", "(input):11:2: integer 6 out of bounds");
}

TEST(NLTest, ReadSuffix) {
  EXPECT_READ("suffix foo:0:5: i0 = 3, i1 = 2, i2 = 1, i3 = 2, i4 = 3;",
              "S0 5 foo\n0 3\n1 2\n2 1\n3 2\n4 3\n");
  EXPECT_READ_ERROR("S-1 1 foo\n", "(input):11:2: expected unsigned integer");
  EXPECT_READ_ERROR("S8 1 foo\n", "(input):11:2: invalid suffix kind");
  EXPECT_READ_ERROR("S0 0 foo\n", "(input):11:4: integer 0 out of bounds");
  EXPECT_READ_ERROR("S0 6 foo\n", "(input):11:4: integer 6 out of bounds");
}

TEST(NLTest, InvalidSegmentType) {
  EXPECT_READ_ERROR("?", "(input):11:1: invalid segment type");
  EXPECT_READ_ERROR(std::string("C0\nn4.2\n") + '\0',
                    "(input):13:1: invalid segment type");
}
}  // namespace
