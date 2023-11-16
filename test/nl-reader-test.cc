/*
 nl reader tests

 Copyright (C) 2013 - 2016 AMPL Optimization Inc

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

// Include the source file to test the implementation.
#include "../src/nl-reader.cc"

#include <cmath>
#include <climits>
#include <cstring>

#include "gmock/gmock.h"
#include "gtest-extra.h"
#include "mock-file.h"
#include "mock-problem-builder.h"
#include "util.h"

using mp::NLHeader;
using mp::ReadError;
using mp::ReadNLString;
using mp::internal::BinaryReader;
using mp::internal::EndiannessConverter;
using mp::internal::NLProblemBuilder;
namespace expr = mp::expr;

typedef mp::internal::TextReader<> TextReader;

using testing::_;
using testing::Field;
using testing::Return;
using testing::ReturnRef;
using testing::StrictMock;
using testing::Throw;

// Suppress a spurious warning
// "'<anonymous>' may be used uninitialized in this function" in GCC 4.4.
#if FMT_GCC_VERSION == 404
# pragma GCC diagnostic ignored "-Wuninitialized"
#endif

class MockNLHandler : public mp::NLHandler<MockNLHandler, TestExpr> {
 public:
  MOCK_METHOD1(OnUnhandled, void (const char *name));
  MOCK_METHOD1(OnHeader, void (const NLHeader &h));
  MOCK_METHOD3(OnVarBounds, void (int index, double lb, double ub));
  MOCK_METHOD3(OnObj, void (int index, mp::obj::Type type, NumericExpr expr));
  MOCK_CONST_METHOD1(NeedObj, bool (int obj_index));
  MOCK_METHOD2(OnLinearObjExpr,
               LinearObjHandler (int obj_index, int num_linear_terms));
  MOCK_METHOD0(EndInput, void ());
};

TEST(NLStringRefTest, ConstructFromCString) {
  const char *data = "foo";
  const mp::NLStringRef s(data);
  EXPECT_EQ(data, s.c_str());
  EXPECT_EQ(3u, s.size());
}

TEST(NLStringRefTest, ConstructFromCStringAndSize) {
  const char *data = "foo\0bar";
  const mp::NLStringRef s(data, 7);
  EXPECT_EQ(data, s.c_str());
  EXPECT_EQ(7u, s.size());
  EXPECT_ASSERT(mp::NLStringRef("foo", 2), "string not null-terminated");
}

TEST(NLStringRefTest, ConstructFromStdString) {
  std::string str = "foo";
  const mp::NLStringRef s(str);
  EXPECT_EQ(0, s.c_str()[str.size()]);
  EXPECT_EQ(str.c_str(), s.c_str());
  EXPECT_EQ(str.size(), s.size());
}

TEST(ReaderBaseTest, ReadChar) {
  struct TestReaderBase : mp::internal::ReaderBase {
    TestReaderBase(mp::NLStringRef data, fmt::CStringRef name)
         : ReaderBase(data, name) {}
  };
  TestReaderBase rb(" b", "test");
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

class MockLocale {
 public:
  MOCK_METHOD1(strtod, double (const char *&str));
};

TEST(TextReaderTest, ReadDoubleUsesCLocale) {
  const char *str = "4.2";
  mp::internal::TextReader<MockLocale> reader(str, "test.nl");
  EXPECT_CALL(reader.locale(), strtod(str))
    .WillOnce(testing::DoAll(testing::SetArgReferee<0>(str + 1), Return(1.23)));
  EXPECT_EQ(1.23, reader.ReadDouble());
#ifdef FMT_LOCALE
  TextReader default_reader(str, "test.nl");
  fmt::Locale &locale = default_reader.locale();
  mp::internal::Unused(&locale);
#endif
}

TEST(TextReaderTest, ReadName) {
  TextReader reader("  abc  \n  ", "test");
  fmt::StringRef name = reader.ReadName();
  EXPECT_EQ("abc", name.to_string());
  EXPECT_THROW_MSG(reader.ReadName(), mp::ReadError, "test:1:8: expected name");
  reader.ReadTillEndOfLine();
  EXPECT_THROW_MSG(reader.ReadName(), mp::ReadError, "test:2:3: expected name");
}

TEST(TextReaderTest, ReadString) {
  TextReader reader("  3:a\nb\nc\n2x\n2:de \n2:fg", "test");
  fmt::StringRef name = reader.ReadString();
  EXPECT_EQ("a\nb", name.to_string());
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

// Formats header and variable bounds as a string.
std::string FormatHeader(const NLHeader &h, bool var_bounds = true) {
  fmt::MemoryWriter w;
  w << h;
  if (var_bounds) {
    w << "b\n";
    for (int i = 0; i < h.num_vars; ++i)
      w << "1 0\n";
  }
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
  return ReadHeader(ReplaceLine(FormatHeader(NLHeader()), line_index, line));
}

TEST(TextReaderTest, InvalidFormat) {
  EXPECT_THROW_MSG(ReadHeader(0, "x"),
      ReadError, "(input):1:1: expected format specifier");
}

TEST(TextReaderTest, InvalidNumOptions) {
  EXPECT_EQ(0, ReadHeader(0, "ga").num_ampl_options);
  EXPECT_EQ(0, ReadHeader(0, "g-1").num_ampl_options);
  EXPECT_THROW_MSG(ReadHeader(0, "g10"),
      ReadError, "(input):1:2: too many options");
  EXPECT_THROW_MSG(ReadHeader(0,
      fmt::format("g{}", static_cast<unsigned>(INT_MAX) + 1)),
      ReadError, "(input):1:2: number is too big");
}

void CheckReadOptions(int num_options,
    int num_options_to_write, const int *options) {
  fmt::MemoryWriter w;
  w << 'g' << num_options;
  for (int i = 0; i < num_options_to_write; ++i)
    w << ' ' << options[i];
  NLHeader header = ReadHeader(0, w.c_str());
  ASSERT_EQ(num_options, header.num_ampl_options);
  int min_num_options = std::min(num_options, num_options_to_write);
  for (int i = 0; i < min_num_options; ++i)
    EXPECT_EQ(options[i], header.ampl_options[i]);
  for (int i = min_num_options; i < num_options_to_write; ++i)
    EXPECT_EQ(0, header.ampl_options[i]);
}

TEST(TextReaderTest, ReadOptions) {
  const int options[mp::MAX_AMPL_OPTIONS + 1] = {
      3, 5, 7, 11, 13, 17, 19, 23, 29, 31
  };
  for (int i = 0; i < mp::MAX_AMPL_OPTIONS; ++i) {
    for (int j = 0; j < mp::MAX_AMPL_OPTIONS + 1; ++j)
      CheckReadOptions(i, j, options);
  }
  EXPECT_EQ(0, ReadHeader(0, "g").num_ampl_options);
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
  fmt::MemoryWriter w; \
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
   TestReaderBase(mp::NLStringRef data, fmt::CStringRef name)
     : ReaderBase(data, name) {}

   TestReaderBase &get() { return *this; }
};

template <typename Converter = mp::internal::IdentityConverter>
class TestBinaryReader :
    private TestReaderBase, public BinaryReader<Converter> {
 public:
  TestBinaryReader(mp::NLStringRef data, fmt::CStringRef name = "test")
    : TestReaderBase(data, name), BinaryReader<Converter>(get()) {}
};

TEST(BinaryReaderTest, ReportError) {
  TestBinaryReader<> reader("x");
  EXPECT_THROW_MSG(reader.ReportError("a{}c", 'b'), mp::BinaryReadError,
                   "test:offset 0: abc");
}

template <typename T>
void ChangeEndianness(T &value) {
  char *data = reinterpret_cast<char*>(&value);
  std::reverse(data, data + sizeof(T));
}

TEST(BinaryReaderTest, ReadInt) {
  int data[] = {11, -22, 0};
  std::size_t size = sizeof(data) - sizeof(int);
  mp::NLStringRef str(reinterpret_cast<char*>(data), size);
  TestBinaryReader<> reader(str);
  EXPECT_EQ(11, reader.ReadInt<int>());
  EXPECT_EQ(-22, reader.ReadInt<int>());
  EXPECT_THROW_MSG(reader.ReadInt<int>(), mp::BinaryReadError,
                   fmt::format("test:offset {}: unexpected end of file", size));
  ChangeEndianness(data[0]);
  EXPECT_EQ(11, TestBinaryReader<EndiannessConverter>(str).ReadInt<int>());
}

TEST(BinaryReaderTest, ReadLong) {
  long data[] = {42, 0};
  mp::NLStringRef str(reinterpret_cast<char*>(data), sizeof(data) - 1);
  TestBinaryReader<> reader(str);
  EXPECT_EQ(42, reader.ReadInt<long>());
  ChangeEndianness(data[0]);
  EXPECT_EQ(42, TestBinaryReader<EndiannessConverter>(str).ReadInt<long>());
}

TEST(BinaryReaderTest, ReadUInt) {
  int data[] = {11, -22, 0};
  mp::NLStringRef str(reinterpret_cast<char*>(data), sizeof(data) - 1);
  TestBinaryReader<> reader(str);
  EXPECT_EQ(11, reader.ReadUInt());
  std::string message =
      fmt::format("test:offset {}: expected unsigned integer", sizeof(int));
  EXPECT_THROW_MSG(reader.ReadUInt(), mp::BinaryReadError, message);
  ChangeEndianness(data[0]);
  EXPECT_EQ(11, TestBinaryReader<EndiannessConverter>(str).ReadUInt());
}

TEST(BinaryReaderTest, ReadDouble) {
  double data[] = {11, -2.2, 0};
  std::size_t size = sizeof(data) - 1;
  mp::NLStringRef str(reinterpret_cast<char*>(data), size);
  TestBinaryReader<> reader(str);
  EXPECT_EQ(11, reader.ReadDouble());
  EXPECT_EQ(-2.2, reader.ReadDouble());
  EXPECT_THROW_MSG(reader.ReadDouble(), mp::BinaryReadError,
                   fmt::format("test:offset {}: unexpected end of file", size));
  ChangeEndianness(data[0]);
  EXPECT_EQ(11, TestBinaryReader<EndiannessConverter>(str).ReadDouble());
}

template <typename Converter>
void TestReadString(fmt::StringRef (BinaryReader<Converter>::*read)(),
                    bool change_endianness = false) {
  int data[] = {3, 0, 0};
  if (change_endianness)
    ChangeEndianness(data[0]);
  char *str = reinterpret_cast<char*>(data + 1);
  std::strcpy(str, "abc");
  {
    TestBinaryReader<Converter> reader(
          mp::NLStringRef(reinterpret_cast<char*>(data), sizeof(int) + 3));
    fmt::StringRef name = (reader.*read)();
    EXPECT_EQ("abc", name.to_string());
  }
  {
    str[2] = 0;
    std::size_t size = sizeof(int) + 2;
    TestBinaryReader<Converter> reader(
          mp::NLStringRef(reinterpret_cast<char*>(data), size));
    std::string message =
        fmt::format("test:offset {}: unexpected end of file", size);
    EXPECT_THROW_MSG((reader.*read)(), mp::BinaryReadError, message);
  }
  data[0] = -1;
  if (change_endianness)
    ChangeEndianness(data[0]);
  {
    TestBinaryReader<Converter> reader(
          mp::NLStringRef(reinterpret_cast<char*>(data), sizeof(int) + 3));
    EXPECT_THROW_MSG((reader.*read)(), mp::BinaryReadError,
                     "test:offset 0: expected unsigned integer");
  }
}

TEST(BinaryReaderTest, ReadName) {
  TestReadString(&BinaryReader<>::ReadName);
  TestReadString(&BinaryReader<EndiannessConverter>::ReadName, true);
}

TEST(BinaryReaderTest, ReadString) {
  TestReadString(&BinaryReader<>::ReadString);
  TestReadString(&BinaryReader<EndiannessConverter>::ReadString, true);
}

TEST(NLReaderTest, ArithKind) {
  namespace arith = mp::arith;
  EXPECT_GE(arith::GetKind(), arith::UNKNOWN);
  EXPECT_LE(arith::GetKind(), arith::LAST);
  EXPECT_TRUE(arith::IsIEEE(arith::IEEE_LITTLE_ENDIAN));
  EXPECT_TRUE(arith::IsIEEE(arith::IEEE_BIG_ENDIAN));
  EXPECT_FALSE(arith::IsIEEE(arith::UNKNOWN));
  EXPECT_FALSE(arith::IsIEEE(arith::IBM));
  EXPECT_FALSE(arith::IsIEEE(arith::VAX));
  EXPECT_FALSE(arith::IsIEEE(arith::CRAY));
}

class TestNLHandler {
 private:
  // Writes a comma-separated list.
  template <typename Array>
  static void WriteList(fmt::Writer &w, const Array &values) {
    for (std::size_t i = 0, n = values.size(); i < n; ++i) {
      if (i != 0)
        w << ", ";
      w << values[i];
    }
  }

  std::string MakeVarArg(std::string op, const std::vector<std::string> &args) {
    fmt::MemoryWriter w;
    w << op << '(';
    WriteList(w, args);
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
    double infinity = INFINITY;
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
  fmt::MemoryWriter log;  // Call log.

  typedef std::string Expr;
  typedef std::string NumericExpr;
  typedef std::string LogicalExpr;
  typedef std::string CountExpr;
  typedef std::string Reference;

  void OnHeader(const NLHeader &) { log.clear(); }

  bool NeedObj(int) const { return true; }
  int resulting_obj_index(int obj_index) const { return obj_index; }

  void OnVarBounds(int index, double lb, double ub) {
    WriteBounds('v', index, lb, ub);
  }

  void OnConBounds(int index, double lb, double ub) {
    WriteBounds('c', index, lb, ub);
  }

  void OnComplementarity(int con_index, int var_index, mp::ComplInfo info) {
    WriteSep().write("{} <= c{} <= {} complements v{};",
                     info.con_lb(), con_index, info.con_ub(), var_index);
  }

  class LinearExprHandler {
   private:
    std::string str_;
    fmt::Writer &log_;
    bool terminate_;

   public:
    explicit LinearExprHandler(fmt::Writer &log, bool terminate = true)
      : log_(log), terminate_(terminate) {}
    ~LinearExprHandler() {
      log_ << str_;
      if (terminate_)
        log_ << ';';
    }

    void AddTerm(int var_index, double coef) {
      if (!str_.empty())
        str_ += " + ";
      str_ += fmt::format("{} * v{}", coef, var_index);
    }
  };
  typedef LinearExprHandler LinearObjHandler;
  typedef LinearExprHandler LinearConHandler;

  LinearExprHandler OnLinearObjExpr(int index, int num_terms) {
    WriteSep().write("o{} {}: ", index, num_terms);
    return LinearExprHandler(log);
  }
  LinearExprHandler OnLinearConExpr(int index, int num_terms) {
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
  ColumnSizeHandler OnColumnSizes() {
    WriteSep() << "sizes:";
    return ColumnSizeHandler(log);
  }

  void OnObj(int index, mp::obj::Type type, std::string expr) {
    WriteSep() << (type == mp::obj::MAX ? "maximize" : "minimize")
        << " o" << index << ": " << expr << ";";
  }

  void OnAlgebraicCon(int index, std::string expr) {
    WriteSep() << "c" << index << ": " << expr << ";";
  }

  void OnLogicalCon(int index, std::string expr) {
    WriteSep() << "l" << index << ": " << expr << ";";
  }

  LinearExprHandler BeginCommonExpr(int index, int num_terms) {
    WriteSep().write("e{} {}: ", index, num_terms);
    return LinearExprHandler(log, false);
  }
  void EndCommonExpr(int, std::string expr, int position) {
    log.write(" + {} {};", expr, position);
  }

  void OnInitialValue(int var_index, double value) {
    WriteSep().write("v{} := {};", var_index, value);
  }

  void OnInitialDualValue(int con_index, double value) {
    WriteSep().write("c{} := {};", con_index, value);
  }

  void OnFunction(
      int index, fmt::StringRef name, int num_args, mp::func::Type type) {
    WriteSep().write("f{}: {} {} {};", index, name, num_args, type);
  }

  template <typename T>
  class SuffixHandler {
   private:
    fmt::Writer &log_;
    bool first_;

   public:
    explicit SuffixHandler(fmt::Writer &log) : log_(log), first_(true) {}
    ~SuffixHandler() { log_ << ';'; }

    void SetValue(int index, T value) {
      if (!first_)
        log_ << ',';
      first_ = false;
      log_.write(" i{} = {}", index, value);
    }
  };

  typedef SuffixHandler<int> IntSuffixHandler;

  IntSuffixHandler OnIntSuffix(fmt::StringRef name, int kind, int num_values) {
    WriteSep().write("suffix {}:{}:{}:", name, kind, num_values);
    return IntSuffixHandler(log);
  }

  typedef SuffixHandler<double> DblSuffixHandler;

  DblSuffixHandler OnDblSuffix(fmt::StringRef name, int kind, int num_values) {
    WriteSep().write("suffix {}:{}:{}:", name, kind, num_values);
    return DblSuffixHandler(log);
  }

  std::string OnNumber(double value) {
    return fmt::format("{}", value);
  }

  std::string OnVariableRef(int index) {
    return fmt::format("v{}", index);
  }

  std::string OnCommonExprRef(int index) {
    return fmt::format("e{}", index);
  }

  std::string OnUnary(expr::Kind kind, std::string arg) {
    return fmt::format("u{}({})", nl_opcode(kind), arg);
  }

  std::string OnBinary(expr::Kind kind, std::string lhs, std::string rhs) {
    return fmt::format("b{}({}, {})", nl_opcode(kind), lhs, rhs);
  }

  std::string OnIf(std::string condition,
                   std::string then_expr, std::string else_expr) {
    return fmt::format("if {} then {} else {}",
                       condition, then_expr, else_expr);
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
    fmt::MemoryWriter w;
    w << "<<";
    WriteList(w, h.breakpoints_);
    w << "; ";
    WriteList(w, h.slopes_);
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
    fmt::MemoryWriter w;
    w << 'f' << h.func_index_ << '(';
    WriteList(w, h.args_);
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

  typedef ArgHandler VarArgHandler;

  VarArgHandler BeginVarArg(expr::Kind kind, int) {
    return VarArgHandler(fmt::format("v{}", nl_opcode(kind)));
  }
  std::string EndVarArg(VarArgHandler h) {
    return MakeVarArg(h.name_, h.args_);
  }

  ArgHandler BeginSum(int) { return ArgHandler("sum"); }
  std::string EndSum(ArgHandler h) { return MakeVarArg(h.name_, h.args_); }

  typedef ArgHandler NumberOfArgHandler;

  ArgHandler BeginNumberOf(int, std::string arg0) {
    return ArgHandler("numberof " + arg0 + " in ");
  }
  std::string EndNumberOf(ArgHandler h) { return MakeVarArg(h.name_, h.args_); }

  typedef ArgHandler SymbolicArgHandler;

  ArgHandler BeginSymbolicNumberOf(int, std::string arg0) {
    return ArgHandler("numberof " + arg0 + " in ");
  }
  std::string EndSymbolicNumberOf(ArgHandler h) {
    return MakeVarArg(h.name_, h.args_);
  }

  typedef ArgHandler CountArgHandler;

  ArgHandler BeginCount(int) { return ArgHandler("count"); }
  std::string EndCount(ArgHandler h) { return MakeVarArg(h.name_, h.args_); }

  std::string OnBool(bool value) { return fmt::format("l{:d}", value); }
  std::string OnNot(std::string arg) { return fmt::format("not {}", arg); }

  std::string OnBinaryLogical(
      expr::Kind kind, std::string lhs, std::string rhs) {
    return fmt::format("bl{}({}, {})", nl_opcode(kind), lhs, rhs);
  }

  std::string OnRelational(expr::Kind kind, std::string lhs, std::string rhs) {
    return fmt::format("r{}({}, {})", nl_opcode(kind), lhs, rhs);
  }

  std::string OnLogicalCount(
      expr::Kind kind, std::string lhs, std::string rhs) {
    return fmt::format("lc{}({}, {})", nl_opcode(kind), lhs, rhs);
  }

  std::string OnImplication(std::string condition,
                            std::string then_expr, std::string else_expr) {
    return fmt::format("{} ==> {} else {}", condition, then_expr, else_expr);
  }

  typedef ArgHandler LogicalArgHandler;

  ArgHandler BeginIteratedLogical(expr::Kind kind, int) {
    return ArgHandler(fmt::format("il{}", nl_opcode(kind)));
  }
  std::string EndIteratedLogical(ArgHandler h) {
    return MakeVarArg(h.name_, h.args_);
  }

  typedef ArgHandler PairwiseArgHandler;

  ArgHandler BeginPairwise(expr::Kind kind, int) {
    return ArgHandler(kind == expr::ALLDIFF ? "alldiff" : "!alldiff");
  }
  std::string EndPairwise(ArgHandler h) { return MakeVarArg(h.name_, h.args_); }

  std::string OnString(fmt::StringRef value) {
    return fmt::format("'{}'", value.data() ? value : "");
  }

  std::string OnSymbolicIf(std::string condition,
                           std::string then_expr, std::string else_expr) {
    return OnIf(condition, then_expr, else_expr);
  }

  void EndInput() {}
};

TEST(NLReaderTest, WriteTextHeader) {
  NLHeader h = MakeTestHeader();
  fmt::MemoryWriter w;
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

TEST(NLReaderTest, WriteBinaryHeader) {
  NLHeader header = NLHeader();
  header.format = NLHeader::BINARY;
  header.num_ampl_options = 3;
  for (int i = 0; i < header.num_ampl_options; ++i)
    header.ampl_options[i] = 11 * (i + 1);
  header.arith_kind = (NLArithKind)mp::arith::CRAY;
  fmt::MemoryWriter w;
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
  typedef struct TestReference : TestNumericExpr {} Reference;
  typedef struct TestCountExpr : TestNumericExpr {} CountExpr;
  typedef struct TestLogicalExpr : TestExpr {} LogicalExpr;

  void OnHeader(const NLHeader &) {}

  bool NeedObj(int) const { return true; }
  int resulting_obj_index(int obj_index) const { return obj_index; }

  void OnVarBounds(int, double, double) {}
  void OnConBounds(int, double, double) {}
  void OnComplementarity(int, int, mp::ComplInfo) {}

  struct LinearObjHandler {
    void AddTerm(int, double) {}
  };

  LinearObjHandler OnLinearObjExpr(int, int) {
    return LinearObjHandler();
  }

  struct LinearConHandler {
    void AddTerm(int, double) {}
  };

  LinearConHandler OnLinearConExpr(int, int) {
    return LinearConHandler();
  }

  struct LinearExprHandler {
    void AddTerm(int, double) {}
  };

  LinearExprHandler BeginCommonExpr(int, int) {
    return LinearExprHandler();
  }

  void EndCommonExpr(int, TestNumericExpr, int) {}

  struct ColumnSizeHandler {
    void Add(int) {}
  };
  ColumnSizeHandler OnColumnSizes() { return ColumnSizeHandler(); }

  void OnObj(int, mp::obj::Type, TestNumericExpr) {}
  void OnAlgebraicCon(int, TestNumericExpr) {}
  void OnLogicalCon(int, TestLogicalExpr) {}

  void OnInitialValue(int, double) {}
  void OnInitialDualValue(int, double) {}

  void OnFunction(int, fmt::StringRef, int, mp::func::Type) {}

  template <typename T>
  struct SuffixHandler {
    void SetValue(int, T) {}
  };

  typedef SuffixHandler<int> IntSuffixHandler;

  IntSuffixHandler OnIntSuffix(fmt::StringRef, int, int) {
    return IntSuffixHandler();
  }

  typedef SuffixHandler<double> DblSuffixHandler;

  DblSuffixHandler OnDblSuffix(fmt::StringRef, int, int) {
    return DblSuffixHandler();
  }

  TestNumericExpr OnNumber(double) { return TestNumericExpr(); }
  TestReference OnVariableRef(int) { return TestReference(); }
  TestReference OnCommonExprRef(int) { return TestReference(); }
  TestNumericExpr OnUnary(expr::Kind, TestNumericExpr) {
    return TestNumericExpr();
  }

  TestNumericExpr OnBinary(expr::Kind, TestNumericExpr, TestNumericExpr) {
    return TestNumericExpr();
  }

  TestNumericExpr OnIf(TestLogicalExpr, TestNumericExpr, TestNumericExpr) {
    return TestNumericExpr();
  }

  struct PLTermHandler {
    void AddSlope(double) {}
    void AddBreakpoint(double) {}
  };

  PLTermHandler BeginPLTerm(int) { return PLTermHandler(); }
  TestNumericExpr EndPLTerm(PLTermHandler, TestNumericExpr) {
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

  struct VarArgHandler {
    void AddArg(NumericExpr) {}
  };

  VarArgHandler BeginVarArg(expr::Kind, int) { return VarArgHandler(); }
  TestNumericExpr EndVarArg(VarArgHandler) { return TestNumericExpr(); }

  NumericArgHandler BeginSum(int) { return NumericArgHandler(); }
  TestNumericExpr EndSum(NumericArgHandler) { return TestNumericExpr(); }

  struct NumberOfArgHandler {
    void AddArg(NumericExpr) {}
  };

  NumberOfArgHandler BeginNumberOf(int, TestNumericExpr) {
    return NumberOfArgHandler();
  }
  TestNumericExpr EndNumberOf(NumberOfArgHandler) { return TestNumericExpr(); }

  struct SymbolicArgHandler {
    void AddArg(TestExpr) {}
  };

  SymbolicArgHandler BeginSymbolicNumberOf(int, TestExpr) {
    return SymbolicArgHandler();
  }
  TestNumericExpr EndSymbolicNumberOf(SymbolicArgHandler) {
    return TestNumericExpr();
  }

  struct CountArgHandler {
    void AddArg(LogicalExpr) {}
  };

  CountArgHandler BeginCount(int) { return CountArgHandler(); }
  TestCountExpr EndCount(CountArgHandler) { return TestCountExpr(); }

  TestLogicalExpr OnBool(bool) { return TestLogicalExpr(); }
  TestLogicalExpr OnNot(TestLogicalExpr) { return TestLogicalExpr(); }

  TestLogicalExpr OnBinaryLogical(
      expr::Kind, TestLogicalExpr, TestLogicalExpr) {
    return TestLogicalExpr();
  }

  TestLogicalExpr OnRelational(expr::Kind, TestNumericExpr, TestNumericExpr) {
    return TestLogicalExpr();
  }

  TestLogicalExpr OnLogicalCount(expr::Kind, TestNumericExpr, TestCountExpr) {
    return TestLogicalExpr();
  }

  TestLogicalExpr OnImplication(
      TestLogicalExpr, TestLogicalExpr, TestLogicalExpr) {
    return TestLogicalExpr();
  }

  struct LogicalArgHandler {
    void AddArg(LogicalExpr) {}
  };

  LogicalArgHandler BeginIteratedLogical(expr::Kind, int) {
    return LogicalArgHandler();
  }
  TestLogicalExpr EndIteratedLogical(LogicalArgHandler) {
    return TestLogicalExpr();
  }

  struct PairwiseArgHandler {
    void AddArg(NumericExpr) {}
  };

  PairwiseArgHandler BeginPairwise(expr::Kind, int) {
    return PairwiseArgHandler();
  }
  TestLogicalExpr EndPairwise(PairwiseArgHandler) { return TestLogicalExpr(); }

  TestExpr OnString(fmt::StringRef) { return TestExpr(); }

  TestExpr OnSymbolicIf(TestLogicalExpr, TestExpr, TestExpr) {
    return TestExpr();
  }

  void EndInput() {}
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
TEST(NLReaderTest, ExprHierarchy) {
  TestNLHandler2 handler;
  ReadNLString(FormatHeader(MakeHeader()), handler);
}

TEST(NLReaderTest, NoNewlineAtEOF) {
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
    "b\0deadbeef", handler),
      ReadError, "(input):11:2: expected newline");
}

std::string ReadNL(std::string body, bool var_bounds = true) {
  TestNLHandler handler;
  ReadNLString(FormatHeader(MakeHeader(), var_bounds) + body, handler);
  return handler.log.str();
}

#define EXPECT_READ(expected_output, nl_body) \
  EXPECT_EQ(std::string("v0 <= 0; v1 <= 0; v2 <= 0; v3 <= 0; v4 <= 0; ") + \
            expected_output, ReadNL(nl_body))

#define EXPECT_READ_ERROR(nl_body, error) \
  EXPECT_THROW_MSG(ReadNL(nl_body), ReadError, error)

TEST(NLReaderTest, ReadObj) {
  EXPECT_READ("minimize o1: ;", "O1 0\nn0\n");
  EXPECT_READ("maximize o0: v0;", "O0 1\nv0\n");
  EXPECT_READ("maximize o5: v0;", "O5 10\nv0\n");
  EXPECT_READ_ERROR("O0 -1\nn0\n", "(input):17:4: expected unsigned integer");
  EXPECT_READ_ERROR("O-1 0\nn0\n", "(input):17:2: expected unsigned integer");
  EXPECT_READ_ERROR("O6 0\nn0\n", "(input):17:2: integer 6 out of bounds");
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
    "(input):18:2: number is too big");
}

TEST(NLReaderTest, ReadCon) {
  EXPECT_READ("c0: ;", "C0\nn0\n");
  EXPECT_READ("c0: u13(0);", "C0\no13\nn0\n");
}

TEST(NLReaderTest, ReadNumericConstant) {
  EXPECT_READ("c0: 4.2;", "C0\nn4.2\n");
  EXPECT_READ("c0: -100;", "C0\nn-1e+2\n");
  CheckReadInt<short>('s');
  if (sizeof(double) == 2 * sizeof(int))
    CheckReadInt<int>('l');
  else
    CheckReadInt<long>('l');
}

TEST(NLReaderTest, ReadVariable) {
  EXPECT_READ("c0: v4;", "C0\nv4\n");
  EXPECT_READ_ERROR("C0\nv-1\n", "(input):18:2: expected unsigned integer");
  EXPECT_READ_ERROR("C0\nv6\n", "(input):18:2: integer 6 out of bounds");
}

TEST(NLReaderTest, ReadCommonExprRef) {
  EXPECT_READ("c0: e0;", "C0\nv5\n");
}

TEST(NLReaderTest, ReadUnaryExpr) {
  EXPECT_READ("c0: u13(v3);", "C0\no13\nv3\n");
}

TEST(NLReaderTest, ReadBinaryExpr) {
  EXPECT_READ("c0: b0(v1, 42);", "C0\no0\nv1\nn42\n");
}

TEST(NLReaderTest, ReadIfExpr) {
  EXPECT_READ("c0: if l1 then v1 else v2;", "C0\no35\nn1\nv1\nv2\n");
}

TEST(NLReaderTest, ReadPiecewiseLinearExpr) {
  EXPECT_READ("c0: <<0; -1, 1>> v1;", "C0\no64\n2\nn-1.0\ns0\nl1\nv1\n");
  EXPECT_READ_ERROR("C0\no64\n-1\nn0\nv1\n",
    "(input):19:1: expected unsigned integer");
  EXPECT_READ_ERROR("C0\no64\n1\nn0\nv1\n",
    "(input):19:1: too few slopes in piecewise-linear term");
  EXPECT_READ_ERROR("C0\no64\n2\nv1\nn0\nn1\nv1\n",
    "(input):20:1: expected constant");
  EXPECT_READ_ERROR("C0\no64\n2\nn-1\nv0\nn1\nv1\n",
    "(input):21:1: expected constant");
  EXPECT_READ_ERROR("C0\no64\n2\nn-1\nn0\nn1\nn1\n",
    "(input):23:1: expected reference");
}

TEST(NLReaderTest, ReadCallExpr) {
  EXPECT_READ("c0: f1(v1, 0);", "C0\nf1 2\nv1\nn0\n");
  EXPECT_READ_ERROR("C0\nf-1 1\nn0\n",
    "(input):18:2: expected unsigned integer");
  EXPECT_READ_ERROR("C0\nf10 1\nn0\n",
                    "(input):18:2: integer 10 out of bounds");
  EXPECT_READ_ERROR("C0\nf1 1\nx\n", "(input):19:1: expected expression");
}

TEST(NLReaderTest, ReadVarArgExpr) {
  EXPECT_READ("c0: v11(v4, 5, v1);", "C0\no11\n3\nv4\nn5\nv1\n");
  EXPECT_READ("c0: v12(v4);", "C0\no12\n1\nv4\n");
  EXPECT_READ_ERROR("C0\no12\n0\n" , "(input):19:1: too few arguments");
}

TEST(NLReaderTest, ReadSumExpr) {
  EXPECT_READ("c0: sum(v4, 5, v1);", "C0\no54\n3\nv4\nn5\nv1\n");
  EXPECT_READ_ERROR("C0\no54\n2\nv4\nn5\n", "(input):19:1: too few arguments");
}

TEST(NLReaderTest, ReadCountExpr) {
  EXPECT_READ("c0: count(l1, r24(v1, 42), l0);",
               "C0\no59\n3\nn1\no24\nv1\nn42\nn0\n");
  EXPECT_READ("c0: count(l1);", "C0\no59\n1\nn1\n");
  EXPECT_READ_ERROR("C0\no59\n0\n", "(input):19:1: too few arguments");
}

TEST(NLReaderTest, ReadNumberOfExpr) {
  EXPECT_READ("c0: numberof v4 in (5, v1);", "C0\no60\n3\nv4\nn5\nv1\n");
  EXPECT_READ("c0: numberof v4 in ();", "C0\no60\n1\nv4\n");
  EXPECT_READ_ERROR("C0\no60\n0\n", "(input):19:1: too few arguments");
}

TEST(NLReaderTest, ReadSymbolicNumberOfExpr) {
  EXPECT_READ("c0: numberof 'a' in ('b', 42);",
              "C0\no61\n3\nh1:a\nh1:b\nn42\n");
  EXPECT_READ("c0: numberof 42 in ('b', 'c');",
              "C0\no61\n3\nn42\nh1:b\nh1:c\n");
  EXPECT_READ_ERROR("C0\no61\n0\n", "(input):19:1: too few arguments");
  EXPECT_READ_ERROR("C0\no61\n1\nx", "(input):20:1: expected expression");
}

TEST(NLReaderTest, ReadLogicalConstant) {
  EXPECT_READ("l0: l0;", "L0\nn0\n");
  EXPECT_READ("l0: l1;", "L0\nn1\n");
  EXPECT_READ("l0: l1;", "L0\nn4.2\n");
  EXPECT_READ("l0: l1;", "L0\ns1\n");
  EXPECT_READ("l0: l1;", "L0\nl1\n");
}

TEST(NLReaderTest, ReadNotExpr) {
  EXPECT_READ("l0: not l0;", "L0\no34\nn0\n");
}

TEST(NLReaderTest, ReadBinaryLogicalExpr) {
  EXPECT_READ("l0: bl20(l1, l0);", "L0\no20\nn1\nn0\n");
}

TEST(NLReaderTest, ReadRelationalExpr) {
  EXPECT_READ("l0: r23(v1, 0);", "L0\no23\nv1\nn0\n");
}

TEST(NLReaderTest, ReadLogicalCountExpr) {
  EXPECT_READ("l0: lc63(v1, count(l1));", "L0\no63\nv1\no59\n1\nn1\n");
  EXPECT_READ_ERROR("L0\no63\nv1\nn0\n",
    "(input):20:1: expected count expression");
  EXPECT_READ_ERROR("L0\no63\nv1\no16\nn0\n",
    "(input):20:2: expected count expression");
}

TEST(NLReaderTest, ReadImplicationExpr) {
  EXPECT_READ("l0: l1 ==> l0 else l1;", "L0\no72\nn1\nn0\nn1\n");
}

TEST(NLReaderTest, ReadIteratedLogicalExpr) {
  EXPECT_READ("l0: il71(l1, l0, l1);", "L0\no71\n3\nn1\nn0\nn1\n");
  EXPECT_READ_ERROR("L0\no71\n2\nn1\nn0\n", "(input):19:1: too few arguments");
}

TEST(NLReaderTest, ReadAllDiffExpr) {
  EXPECT_READ("l0: alldiff(v4, 5, v1);", "L0\no74\n3\nv4\nn5\nv1\n");
  EXPECT_READ("l0: !alldiff(v4);", "L0\no75\n1\nv4\n");
  EXPECT_READ_ERROR("L0\no74\n0\n", "(input):19:1: too few arguments");
}

TEST(NLReaderTest, ReadStringLiteral) {
  EXPECT_READ("c0: f1('');", "C0\nf1 1\nh0:\n");
  EXPECT_READ("c0: f1('abc');", "C0\nf1 1\nh3:abc\n");
  EXPECT_READ("c0: f1('ab\nc');", "C0\nf1 1\nh4:ab\nc\n");
  const char input[] = "C0\nf1 1\nh1:\0\n";
  const char output[] = "c0: f1('\0');";
  EXPECT_READ(std::string(output, sizeof(output) - 1),
              std::string(input, sizeof(input) - 1));
  EXPECT_READ_ERROR("C0\nf1 1\nh3:ab",
        "(input):19:6: unexpected end of file in string");
  EXPECT_READ_ERROR("C0\nf1 1\nh3:a\n",
        "(input):20:1: unexpected end of file in string");
  EXPECT_READ_ERROR("C0\nf1 1\nh3:abc", "(input):19:7: expected newline");
  EXPECT_READ_ERROR("C0\nf1 1\nh3:ab\n", "(input):20:1: expected newline");
}

TEST(NLReaderTest, ReadSymbolicIfExpr) {
  EXPECT_READ("c0: f1(if l1 then v1 else 'abc');",
              "C0\nf1 1\no65\nn1\nv1\nh3:abc\n");
  EXPECT_READ("c0: f1(if l1 then 'abc' else 42);",
              "C0\nf1 1\no65\nn1\nh3:abc\nn42\n");
  EXPECT_READ_ERROR("C0\nf1 1\no65\nx",
                    "(input):20:1: expected logical expression");
  EXPECT_READ_ERROR("C0\nf1 1\no65\nn1\nx",
                    "(input):21:1: expected expression");
}

TEST(NLReaderTest, ReadInvalidOpCode) {
  EXPECT_READ_ERROR("C0\no-1\n", "(input):18:2: expected unsigned integer");
  EXPECT_READ_ERROR("C0\no83\n", "(input):18:2: invalid opcode 83");
}

TEST(NLReaderTest, ReadInvalidNumericExpr) {
  EXPECT_READ_ERROR("C0\nx\n", "(input):18:1: expected expression");
  EXPECT_READ_ERROR("C0\no22\nv1\nn0\n",
        "(input):18:2: expected numeric expression opcode");
}

TEST(NLReaderTest, ReadInvalidLogicalExpr) {
  EXPECT_READ_ERROR("L0\nx\n", "(input):18:1: expected logical expression");
  EXPECT_READ_ERROR("L0\no0\nv1\nn0\n",
        "(input):18:2: expected logical expression opcode");
}

TEST(NLReaderTest, ReadVarBounds) {
  EXPECT_THROW_MSG(ReadNL("", false), ReadError,
                   "(input):11:1: segment 'b' missing");
  EXPECT_THROW_MSG(ReadNL("b\n"), ReadError,
                   "(input):17:1: duplicate 'b' segment");
  EXPECT_EQ("1.1 <= v0; v1 <= 22; v2 = 33; v3; 44 <= v4 <= 55;",
            ReadNL("b\n21.1\n1 22\n4 33\n3\n0 44 55\n", false));
  EXPECT_THROW_MSG(ReadNL("b\n-1\n", false), ReadError,
                   "(input):12:1: expected bound");
  EXPECT_THROW_MSG(ReadNL("b\n5 1\n", false), ReadError,
                   "(input):12:1: COMPL bound type is invalid for variables");
  EXPECT_THROW_MSG(ReadNL("b\n2 11\n1 22\n4 33\n3\n", false), ReadError,
                   "(input):16:1: expected bound");
}

TEST(NLReaderTest, ReadConBounds) {
  EXPECT_READ("1.1 <= c0; c1 <= 22; c2 = 33; c3; 44 <= c4 <= 55; "
              "-inf <= c5 <= inf complements v1; "
              "-inf <= c6 <= 0 complements v4;",
              "r\n21.1\n1 22\n4 33\n3\n0 44 55\n5 7 2\n5 2 5\n");
  EXPECT_READ_ERROR("r\n-1\n", "(input):18:1: expected bound");
  EXPECT_READ_ERROR("r\n6 1\n", "(input):18:1: expected bound");
  EXPECT_READ_ERROR("r\n2 11\n1 22\n4 33\n3\n", "(input):22:1: expected bound");
  EXPECT_READ_ERROR("r\n5 1 0\n", "(input):18:5: integer 0 out of bounds");
  EXPECT_READ_ERROR("r\n5 1 6\n", "(input):18:5: integer 6 out of bounds");
  // Check that there is no overflow for largest possible var index.
  TestNLHandler handler;
  NLHeader header = NLHeader();
  header.num_vars = INT_MAX;
  header.num_algebraic_cons = 1;
  std::string input = fmt::format("r\n5 1 {}\n", INT_MAX);
  TextReader text_reader(input, "(input");
  typedef mp::internal::NLReader<TextReader, TestNLHandler> NLReader;
  NLReader reader(text_reader, header, handler, 0);
  reader.ReadBounds<NLReader::AlgebraicConHandler>();
  EXPECT_EQ(fmt::format("0 <= c0 <= inf complements v{};", INT_MAX - 1),
            handler.log.str());
}

TEST(NLReaderTest, ReadLinearObjExpr) {
  EXPECT_READ("o0 2: 1.3 * v1 + 5 * v3;", "G0 2\n1 1.3\n3 5\n");
  EXPECT_READ("o5 4: 1 * v1 + 1 * v2 + 1 * v3 + 1 * v4;",
              "G5 4\n1 1\n2 1\n3 1\n4 1\n");
  EXPECT_READ_ERROR("G-1", "(input):17:2: expected unsigned integer");
  EXPECT_READ_ERROR("G6", "(input):17:2: integer 6 out of bounds");
  EXPECT_READ_ERROR("G0 0", "(input):17:4: integer 0 out of bounds");
  EXPECT_READ_ERROR("G0 6", "(input):17:4: integer 6 out of bounds");
  EXPECT_READ_ERROR("G0 1\n-1 0\n", "(input):18:1: expected unsigned integer");
  EXPECT_READ_ERROR("G0 1\n5 0\n", "(input):18:1: integer 5 out of bounds");
}

// Test that handler's OnLinearObjExpr is not called if NeedObj returns false.
TEST(NLReaderTest, SkipObj) {
  testing::NiceMock<MockNLHandler> handler;
  EXPECT_CALL(handler, NeedObj(0)).WillOnce(Return(false));
  EXPECT_CALL(handler, OnLinearObjExpr(_, _)).Times(0);
  auto header = NLHeader();
  header.num_vars = header.num_objs = 1;
  ReadNLString(FormatHeader(header) + "G0 1\n0 1\n", handler);
}

// Test that handler's OnLinearObjExpr is called if NeedObj returns true.
TEST(NLReaderTest, PassObj) {
  testing::NiceMock<MockNLHandler> handler;
  EXPECT_CALL(handler, NeedObj(0)).WillOnce(Return(true));
  EXPECT_CALL(handler, OnLinearObjExpr(0, 1)).
      WillOnce(Return(MockNLHandler::LinearObjHandler()));
  auto header = NLHeader();
  header.num_vars = header.num_objs = 1;
  ReadNLString(FormatHeader(header) + "G0 1\n0 1\n", handler);
}

TEST(NLReaderTest, ReadLinearConExpr) {
  EXPECT_READ("c0 2: 1.3 * v1 + 5 * v3;", "J0 2\n1 1.3\n3 5\n");
  EXPECT_READ("c5 4: 1 * v1 + 1 * v2 + 1 * v3 + 1 * v4;",
              "J5 4\n1 1\n2 1\n3 1\n4 1\n");
  EXPECT_READ_ERROR("J-1", "(input):17:2: expected unsigned integer");
  EXPECT_READ_ERROR("J8", "(input):17:2: integer 8 out of bounds");
  EXPECT_READ_ERROR("J0 0", "(input):17:4: integer 0 out of bounds");
  EXPECT_READ_ERROR("J0 6", "(input):17:4: integer 6 out of bounds");
  EXPECT_READ_ERROR("J0 1\n-1 0\n", "(input):18:1: expected unsigned integer");
  EXPECT_READ_ERROR("J0 1\n5 0\n", "(input):18:1: integer 5 out of bounds");
}

TEST(NLReaderTest, ReadColumnSizes) {
  EXPECT_READ("sizes: 1 2 2 4;", "k4\n1\n3\n5\n9\n");
  EXPECT_READ("sizes: 1 2 2 4;", "K4\n1\n2\n2\n4\n");
  EXPECT_READ_ERROR("k3\n", "(input):17:2: expected 4");
  EXPECT_READ_ERROR("k4\n-1\n", "(input):18:1: expected unsigned integer");
  EXPECT_READ_ERROR("k4\n2\n1\n", "(input):19:1: invalid column offset");
}

TEST(NLReaderTest, ReadInitialValues) {
  EXPECT_READ("v4 := 1.1; v3 := 0; v2 := 1; v1 := 2; v0 := 3;",
              "x5\n4 1.1\n3 0\n2 1\n1 2\n0 3\n");
  EXPECT_READ_ERROR("x6\n", "(input):17:2: too many initial values");
  EXPECT_READ_ERROR("x1\n-1 0\n", "(input):18:1: expected unsigned integer");
  EXPECT_READ_ERROR("x1\n5 0\n", "(input):18:1: integer 5 out of bounds");
  EXPECT_READ_ERROR("x2\n4 1.1\n\n", "(input):19:1: expected unsigned integer");
}

TEST(NLReaderTest, ReadInitialDualValues) {
  EXPECT_READ("c4 := 1.1; c3 := 0; c2 := 1; c1 := 2; "
              "c0 := 3; c5 := 1; c6 := 2;",
              "d7\n4 1.1\n3 0\n2 1\n1 2\n0 3\n5 1\n6 2\n");
  EXPECT_READ_ERROR("d8\n", "(input):17:2: too many initial values");
  EXPECT_READ_ERROR("d1\n-1 0\n", "(input):18:1: expected unsigned integer");
  EXPECT_READ_ERROR("d1\n7 0\n", "(input):18:1: integer 7 out of bounds");
  EXPECT_READ_ERROR("d2\n4 1.1\n\n", "(input):19:1: expected unsigned integer");
}

TEST(NLReaderTest, ReadFunction) {
  EXPECT_READ("f0: foo 2 1;", "F0 1 2 foo\n");
  EXPECT_READ("f0: foo -1 0;", "F0 0 -1 foo\n");
  EXPECT_READ_ERROR("F0 1 2 \n", "(input):17:8: expected name");
  EXPECT_READ_ERROR("F-1 0 0 f\n", "(input):17:2: expected unsigned integer");
  EXPECT_READ_ERROR("F9 0 0 f\n", "(input):17:2: integer 9 out of bounds");
  EXPECT_READ_ERROR("F0 -1 0 f\n", "(input):17:4: expected unsigned integer");
  EXPECT_READ_ERROR("F0 2 0 f\n", "(input):17:4: invalid function type");
}

TEST(NLReaderTest, ReadCommonExpr) {
  EXPECT_READ("e0 0:  + b2(v0, 42) 1;", "V5 0 1\no2\nv0\nn42\n");
  EXPECT_READ("e0 2: 2 * v1 + 3 * v0 + 0 1;", "V5 2 1\n1 2.0\n0 3\nn0\n");
  EXPECT_READ_ERROR("V4 0 1\nv0\n", "(input):17:2: integer 4 out of bounds");
  EXPECT_READ_ERROR("V6 0 1\nv0\n", "(input):17:2: integer 6 out of bounds");
}

TEST(NLReaderTest, ReadSuffix) {
  EXPECT_READ("suffix foo:0:5: i0 = 3, i1 = 2, i2 = 1, i3 = 2, i4 = 3;",
              "S0 5 foo\n0 3\n1 2\n2 1\n3 2\n4 3\n");
  EXPECT_READ_ERROR("S-1 1 foo\n", "(input):17:2: expected unsigned integer");
  EXPECT_READ_ERROR("S8 1 foo\n", "(input):17:2: invalid suffix kind");
  EXPECT_READ_ERROR("S0 0 foo\n", "(input):17:4: integer 0 out of bounds");
  EXPECT_READ_ERROR("S0 6 foo\n", "(input):17:4: integer 6 out of bounds");
}

TEST(NLReaderTest, InvalidSegmentType) {
  EXPECT_READ_ERROR("?", "(input):17:1: invalid segment type");
  EXPECT_READ_ERROR(std::string("C0\nn4.2\n") + '\0',
                    "(input):19:1: invalid segment type");
}

TEST(NLReaderTest, ErrorOnNonexistentNLFile) {
  TestNLHandler handler;
  EXPECT_THROW(mp::ReadNLFile("nonexistent", handler), fmt::SystemError);
}

void CheckReadFile(std::string nl) {
  const char *filename = "test.nl";
  WriteFile(filename, nl);
  TestNLHandler handler;
  mp::ReadNLFile(filename, handler);
  EXPECT_EQ("v0 <= 0; v1 <= 0; v2 <= 0; v3 <= 0; v4 <= 0; c0: 4.2;",
            handler.log.str());
  mp::internal::NLFileReader<> reader;
  reader.Read(filename, handler, 0);
}

TEST(NLReaderTest, ReadNLFile) {
  std::string header = FormatHeader(MakeHeader());
  std::string nl = header + "C0\nn4.2";
  std::size_t page_size = fmt::getpagesize();
  EXPECT_LT(nl.size() + 1, page_size);
  CheckReadFile(nl + "\n");
  for (std::size_t i = nl.size(); i < page_size - 2; ++i)
    nl.push_back(' ');
  EXPECT_EQ(page_size, nl.size() + 2);
  CheckReadFile(nl + "\n");
}

TEST(NLReaderTest, ReadNLFileMultipleOfPageSize) {
  std::string header = FormatHeader(MakeHeader());
  std::string nl = header + "C0\nn4.2";
  std::size_t page_size = fmt::getpagesize();
  for (std::size_t i = nl.size(); i < page_size - 1; ++i)
    nl.push_back(' ');
  EXPECT_EQ(page_size, nl.size() + 1);
  CheckReadFile(nl + "\n");
}

struct Cancel {};

TEST(NLReaderTest, FileTooBig) {
  fmt::ULongLong max_size = std::numeric_limits<std::size_t>::max();
  fmt::ULongLong max_long_long = std::numeric_limits<fmt::LongLong>::max();
  mp::internal::NLFileReader<MockFile> reader;
  TestNLHandler handler;
  if (max_size < max_long_long) {
    EXPECT_CALL(reader.file(), size()).WillOnce(Return(max_size));
    EXPECT_CALL(reader.file(), descriptor()).WillOnce(Throw(Cancel()));
    EXPECT_THROW(reader.Read("test", handler, 0), Cancel);
    EXPECT_CALL(reader.file(), size()).WillOnce(Return(max_size + 1));
    EXPECT_THROW_MSG(reader.Read("test", handler, 0),
                     mp::Error, "file test is too big");
  } else {
    EXPECT_CALL(reader.file(), size()).WillOnce(Return(max_long_long));
    EXPECT_CALL(reader.file(), descriptor()).WillOnce(Throw(Cancel()));
    EXPECT_THROW(reader.Read("test", handler, 0), Cancel);
  }
}

struct TestNLHandler3 : mp::NLHandler<TestNLHandler3, int> {
  void OnUnhandled(const char *) {}
};

TEST(NLReaderTest, NLHandler) {
  TestNLHandler3 handler;
  ReadNLString(FormatHeader(MakeHeader()) + "C0\nn4.2\n", handler);
}

TEST(NLReaderTest, ReadBoundsFirst) {
  testing::NiceMock<MockNLHandler> handler;
  testing::InSequence dummy;
  EXPECT_CALL(handler, OnHeader(_));
  EXPECT_CALL(handler, OnVarBounds(0, 1, 2));
  EXPECT_CALL(handler, OnObj(_, _, _));
  auto header = NLHeader();
  header.num_vars = header.num_objs = 1;
  ReadNLString(FormatHeader(header, false) + "O0 0\no2\nv0\nv0\nb\n0 1 2\n",
               handler, "", mp::READ_BOUNDS_FIRST);
}

// Count the number of variable references in all nonlinear expressions.
struct VarCounter : mp::NullNLHandler<int> {
  int num_vars;
  VarCounter() : num_vars(0) {}
  Reference OnVariableRef(int) {
    ++num_vars;
    return Reference();
  }
};

TEST(NLReaderTest, Example) {
  WriteFile("test.nl", FormatHeader(MakeHeader()) + "O0 0\no2\nv0\nv0\n");
  VarCounter counter;
  mp::ReadNLFile("test.nl", counter);
  EXPECT_WRITE(stdout, fmt::print("The number of variable references is {}.",
                                  counter.num_vars),
               "The number of variable references is 2.");
}

TEST(NLReaderTest, EndInput) {
  testing::NiceMock<MockNLHandler> handler;
  testing::InSequence dummy;
  EXPECT_CALL(handler, OnVarBounds(0, 1, 2));
  EXPECT_CALL(handler, EndInput());
  auto header = NLHeader();
  header.num_vars = header.num_objs = 1;
  ReadNLString(FormatHeader(header, false) + "b\n0 1 2\n", handler, "");
}

struct MockNameHandler {
  MOCK_METHOD1(OnName, void (fmt::StringRef name));
};

MATCHER_P2(StringRefEq, data, size, "") {
  return arg.data() == data && arg.size() == static_cast<std::size_t>(size);
}

TEST(NLReaderTest, ReadNames) {
  MockNameHandler handler;
  const char data[] = "ab\nc d\n";
  testing::InSequence sequence;
  EXPECT_CALL(handler, OnName(StringRefEq(data, 2)));
  EXPECT_CALL(handler, OnName(StringRefEq(data + 3, 3)));
  mp::internal::ReadNames("test", data, handler);
}

TEST(NLReaderTest, ReadNamesError) {
  MockNameHandler handler;
  const char data[] = "abc\ndef";
  EXPECT_CALL(handler, OnName(StringRefEq(data, 3)));
  EXPECT_THROW_MSG(mp::internal::ReadNames("test", data, handler),
                   mp::Error, "test:2:4: missing newline");
}

MATCHER_P(GetStringRef, names, "") {
  names->push_back(arg);
  return true;
}

TEST(NLReaderTest, NameReader) {
  TempFile file("test.col");
  WriteFile(file.name(), "abc\ndef\n");
  MockNameHandler handler;
  std::vector<fmt::StringRef> names;
  EXPECT_CALL(handler, OnName(GetStringRef(&names))).Times(2);
  mp::internal::NameReader reader;
  reader.Read(file.name(), handler);
  // Names should be valid after the call to Read.
  EXPECT_EQ("abc", names[0].to_string());
  EXPECT_EQ("def", names[1].to_string());
}

struct DummyNLHandler : public mp::NLHandler<DummyNLHandler, int> {};

TEST(NLReaderTest, OnUnhandled) {
  DummyNLHandler handler;
  EXPECT_THROW_MSG(handler.OnUnhandled("foo"), mp::UnsupportedError,
                   "unsupported: foo");
  EXPECT_THROW_MSG(handler.OnHeader(mp::NLHeader()), mp::UnsupportedError,
                   "unsupported: NL header");
}

#define EXPECT_UNHANDLED(kind, call) { \
    MockNLHandler handler; \
    mp::NLHandler<MockNLHandler, TestExpr> &base = handler; \
    EXPECT_CALL(handler, OnUnhandled(testing::StrEq(kind))); \
    base.call; \
  }

TEST(NLReaderTest, OnUnhandledDispatch) {
  EXPECT_UNHANDLED("NL header", OnHeader(mp::NLHeader()));
  EXPECT_UNHANDLED("objective", OnObj(0, mp::obj::MAX, TestExpr()));
  EXPECT_UNHANDLED("nonlinear constraint", OnAlgebraicCon(0, TestExpr()));
  EXPECT_UNHANDLED("logical constraint", OnLogicalCon(0, TestExpr()));
  EXPECT_UNHANDLED("common expression", BeginCommonExpr(0, 0));
  EXPECT_UNHANDLED("common expression", EndCommonExpr(0, TestExpr(), 0));
  EXPECT_UNHANDLED("complementarity constraint",
                   OnComplementarity(0, 0, mp::ComplInfo(0)));
  EXPECT_UNHANDLED("linear objective", OnLinearObjExpr(0, 0));
  EXPECT_UNHANDLED("linear constraint", OnLinearConExpr(0, 0));
  EXPECT_UNHANDLED("linear common expression", OnLinearCommonExpr(0, 0));
  EXPECT_UNHANDLED("variable bounds", OnVarBounds(0, 0, 0));
  EXPECT_UNHANDLED("constraint bounds", OnConBounds(0, 0, 0));
  EXPECT_UNHANDLED("initial value", OnInitialValue(0, 0));
  EXPECT_UNHANDLED("initial dual value", OnInitialDualValue(0, 0));
  EXPECT_UNHANDLED("column sizes", OnColumnSizes());
  EXPECT_UNHANDLED("function", OnFunction(0, "", 0, mp::func::NUMERIC));
  EXPECT_UNHANDLED("integer suffix", OnIntSuffix("",  mp::suf::OBJ, 0));
  EXPECT_UNHANDLED("double suffix", OnDblSuffix("",  mp::suf::OBJ, 0));
  EXPECT_UNHANDLED("number", OnNumber(0));
  EXPECT_UNHANDLED("variable reference", OnVariableRef(0));
  EXPECT_UNHANDLED("common expression reference", OnCommonExprRef(0));
  EXPECT_UNHANDLED("unary expression", OnUnary(expr::MINUS, TestExpr()));
  EXPECT_UNHANDLED("binary expression",
                   OnBinary(expr::ADD, TestExpr(), TestExpr()));
  EXPECT_UNHANDLED("if expression", OnIf(TestExpr(), TestExpr(), TestExpr()));
  EXPECT_UNHANDLED("piecewise-linear term", BeginPLTerm(0));
  EXPECT_UNHANDLED("piecewise-linear term",
                   EndPLTerm(MockNLHandler::PLTermHandler(), TestExpr()));
  EXPECT_UNHANDLED("call expression", BeginCall(0, 0));
  EXPECT_UNHANDLED("call expression", EndCall(MockNLHandler::CallArgHandler()));
  EXPECT_UNHANDLED("vararg expression", BeginVarArg(mp::expr::MIN, 0));
  EXPECT_UNHANDLED("vararg expression",
                   EndVarArg(MockNLHandler::VarArgHandler()));
  EXPECT_UNHANDLED("summation", BeginSum(0));
  EXPECT_UNHANDLED("summation", EndSum(MockNLHandler::NumericArgHandler()));
  EXPECT_UNHANDLED("count expression", BeginCount(0));
  EXPECT_UNHANDLED("count expression",
                   EndCount(MockNLHandler::CountArgHandler()));
  EXPECT_UNHANDLED("numberof expression", BeginNumberOf(0, TestExpr()));
  EXPECT_UNHANDLED("numberof expression",
                   EndNumberOf(MockNLHandler::NumberOfArgHandler()));
  EXPECT_UNHANDLED("symbolic numberof expression",
                   BeginSymbolicNumberOf(0, TestExpr()));
  EXPECT_UNHANDLED("symbolic numberof expression",
                   EndSymbolicNumberOf(MockNLHandler::SymbolicArgHandler()));
  EXPECT_UNHANDLED("bool", OnBool(true));
  EXPECT_UNHANDLED("logical not", OnNot(TestExpr()));
  EXPECT_UNHANDLED("binary logical expression",
                   OnBinaryLogical(mp::expr::OR, TestExpr(), TestExpr()));
  EXPECT_UNHANDLED("relational expression",
                   OnRelational(mp::expr::LT, TestExpr(), TestExpr()));
  EXPECT_UNHANDLED("logical count expression",
                   OnLogicalCount(mp::expr::ATMOST, TestExpr(), TestExpr()));
  EXPECT_UNHANDLED("implication expression",
                   OnImplication(TestExpr(), TestExpr(), TestExpr()));
  EXPECT_UNHANDLED("iterated logical expression",
                   BeginIteratedLogical(mp::expr::EXISTS, 0));
  EXPECT_UNHANDLED("iterated logical expression",
                   EndIteratedLogical(MockNLHandler::LogicalArgHandler()));
  EXPECT_UNHANDLED("pairwise expression",
                   BeginPairwise(mp::expr::ALLDIFF, 0));
  EXPECT_UNHANDLED("pairwise expression",
                   EndPairwise(MockNLHandler::PairwiseArgHandler()));
  EXPECT_UNHANDLED("string", OnString(""));
  EXPECT_UNHANDLED("symbolic if expression",
                   OnSymbolicIf(TestExpr(), TestExpr(), TestExpr()));
}

TEST(NLReaderTest, NullNLHandler) {
  mp::NullNLHandler<int> handler;
  handler.OnHeader(mp::NLHeader());
}

// Checks if NLProblemBuilder forwards arguments passed to
// adapter_func to ProblemBuilder's builder_func.
#define EXPECT_FORWARD(adapter_func, builder_func, args) { \
  StrictMock<MockProblemBuilder> builder; \
  EXPECT_CALL(builder, builder_func args); \
  NLProblemBuilder<MockProblemBuilder> adapter(builder); \
  adapter.adapter_func args; \
}

// Version of EXPECT_FORWARD for methods with a return value.
#define EXPECT_FORWARD_RET(adapter_func, builder_func, args, result) { \
  StrictMock<MockProblemBuilder> builder; \
  EXPECT_CALL(builder, builder_func args).WillOnce(Return(result)); \
  NLProblemBuilder<MockProblemBuilder> adapter(builder); \
  adapter.adapter_func args; \
}

class NLProblemBuilderTest : public ::testing::Test {
 protected:
  StrictMock<MockProblemBuilder> builder;
  NLProblemBuilder<MockProblemBuilder> adapter;

  NLProblemBuilderTest() : adapter(builder) {}
};

TEST_F(NLProblemBuilderTest, Forward) {
  // Use the same StringRef object in arguments, because StringRef objects
  // are compared as pointers and string literals they point to may not
  // be merged.
  fmt::StringRef str("foo");
  EXPECT_FORWARD_RET(OnIntSuffix, AddIntSuffix, (str, mp::suf::OBJ, 11),
                     TestSuffixHandler<0>(ID));
  EXPECT_FORWARD_RET(OnDblSuffix, AddDblSuffix, (str, mp::suf::OBJ, 11),
                     TestSuffixHandler<1>(ID));

  EXPECT_FORWARD_RET(OnNumber, MakeNumericConstant,
                     (2.2), TestNumericExpr(ID));
  EXPECT_FORWARD_RET(OnVariableRef, MakeVariable, (33), TestReference(ID));
  EXPECT_FORWARD_RET(OnCommonExprRef, MakeCommonExpr, (33), TestReference(ID));
  EXPECT_FORWARD_RET(OnUnary, MakeUnary, (expr::ABS, TestNumericExpr(ID)),
                     TestNumericExpr(ID2));
  EXPECT_FORWARD_RET(OnBinary, MakeBinary,
                     (expr::ADD, TestNumericExpr(ID), TestNumericExpr(ID2)),
                     TestNumericExpr(ID3));
  EXPECT_FORWARD_RET(OnIf, MakeIf,
                     (TestLogicalExpr(ID), TestNumericExpr(ID2),
                      TestNumericExpr(ID3)), TestNumericExpr(ID4));

  EXPECT_FORWARD_RET(BeginPLTerm, BeginPLTerm, (44), TestPLTermBuilder(ID));
  EXPECT_FORWARD_RET(EndPLTerm, EndPLTerm,
                     (TestPLTermBuilder(ID), TestReference(ID2)),
                     TestNumericExpr(ID3));

  EXPECT_FORWARD_RET(BeginVarArg, BeginIterated, (expr::MAX, 77),
                     TestIteratedExprBuilder(ID));
  EXPECT_FORWARD_RET(EndVarArg, EndIterated, (TestIteratedExprBuilder(ID)),
                     TestNumericExpr(ID2));

  EXPECT_FORWARD_RET(BeginSum, BeginSum, (88),
                     TestIteratedExprBuilder(ID));
  EXPECT_FORWARD_RET(EndSum, EndSum, (TestIteratedExprBuilder(ID)),
                     TestNumericExpr(ID2));

  EXPECT_FORWARD_RET(BeginNumberOf, BeginNumberOf, (11, TestNumericExpr(ID)),
                     TestNumberOfExprBuilder(ID2));
  EXPECT_FORWARD_RET(EndNumberOf, EndNumberOf, (TestNumberOfExprBuilder(ID)),
                     TestNumericExpr(ID2));

  EXPECT_FORWARD_RET(BeginSymbolicNumberOf, BeginSymbolicNumberOf,
                     (11, TestExpr(ID)), TestSymbolicNumberOfExprBuilder(ID2));
  EXPECT_FORWARD_RET(EndSymbolicNumberOf, EndSymbolicNumberOf,
                     (TestSymbolicNumberOfExprBuilder(ID)),
                     TestNumericExpr(ID2));

  EXPECT_FORWARD_RET(BeginCount, BeginCount, (99), TestCountExprBuilder(ID));
  EXPECT_FORWARD_RET(EndCount, EndCount, (TestCountExprBuilder(ID)),
                     TestCountExpr(ID2));

  EXPECT_FORWARD_RET(OnBool, MakeLogicalConstant, (true), TestLogicalExpr(ID));
  EXPECT_FORWARD_RET(OnNot, MakeNot, (TestLogicalExpr(ID)),
                     TestLogicalExpr(ID2));
  EXPECT_FORWARD_RET(OnBinaryLogical, MakeBinaryLogical,
                     (expr::OR, TestLogicalExpr(ID), TestLogicalExpr(ID2)),
                     TestLogicalExpr(ID3));
  EXPECT_FORWARD_RET(OnRelational, MakeRelational,
                     (expr::LT, TestNumericExpr(ID), TestNumericExpr(ID2)),
                     TestLogicalExpr(ID3));
  EXPECT_FORWARD_RET(OnLogicalCount, MakeLogicalCount,
                     (expr::ATLEAST, TestNumericExpr(ID), TestCountExpr(ID2)),
                     TestLogicalExpr(ID3));
  EXPECT_FORWARD_RET(OnImplication, MakeImplication,
                     (TestLogicalExpr(ID), TestLogicalExpr(ID2),
                      TestLogicalExpr(ID3)), TestLogicalExpr(ID4));

  EXPECT_FORWARD_RET(BeginIteratedLogical, BeginIteratedLogical,
                     (expr::EXISTS, 22), TestIteratedLogicalExprBuilder(ID));
  EXPECT_FORWARD_RET(EndIteratedLogical, EndIteratedLogical,
                     (TestIteratedLogicalExprBuilder(ID)),
                     TestLogicalExpr(ID2));

  EXPECT_FORWARD_RET(BeginPairwise, BeginPairwise, (expr::ALLDIFF, 33),
                     TestPairwiseExprBuilder(ID));
  EXPECT_FORWARD_RET(EndPairwise, EndPairwise,
                     (TestPairwiseExprBuilder(ID)), TestLogicalExpr(ID2));

  EXPECT_FORWARD_RET(OnString, MakeStringLiteral, (str), TestExpr(ID));

  EXPECT_FORWARD_RET(OnSymbolicIf, MakeSymbolicIf,
                     (TestLogicalExpr(ID), TestExpr(ID2), TestExpr(ID3)),
                     TestExpr(ID4));
}

TEST_F(NLProblemBuilderTest, OnHeader) {
  NLHeader h;
  EXPECT_CALL(builder, SetInfo(testing::Ref(h)));
  adapter.OnHeader(h);
}

TEST_F(NLProblemBuilderTest, AddVars) {
  auto header = mp::NLHeader();
  header.num_vars = 42;
  EXPECT_CALL(builder, SetInfo(testing::Ref(header)));
  EXPECT_CALL(builder, AddVars(header.num_vars, mp::var::CONTINUOUS));
  adapter.OnHeader(header);
  header.num_linear_integer_vars = 42;
  EXPECT_CALL(builder, SetInfo(testing::Ref(header)));
  EXPECT_CALL(builder, AddVars(header.num_vars, mp::var::INTEGER));
  adapter.OnHeader(header);
  header.num_linear_integer_vars = 10;
  EXPECT_CALL(builder, SetInfo(testing::Ref(header)));
  EXPECT_CALL(builder, AddVars(header.num_linear_integer_vars,
                               mp::var::INTEGER));
  EXPECT_CALL(builder, AddVars(header.num_vars - header.num_linear_integer_vars,
                               mp::var::CONTINUOUS));
  adapter.OnHeader(header);
}

TEST_F(NLProblemBuilderTest, OnVarBounds) {
  MockProblemBuilder::MutVariable var;
  EXPECT_CALL(builder, var(66)).WillOnce(ReturnRef(var));
  EXPECT_CALL(var, set_lb(7.7));
  EXPECT_CALL(var, set_ub(8.8));
  adapter.OnVarBounds(66, 7.7, 8.8);
}

TEST_F(NLProblemBuilderTest, AddObjs_allByDefault) {
  auto header = mp::NLHeader();
  header.num_objs = 42;
  EXPECT_CALL(builder, SetInfo(testing::Ref(header)));
  EXPECT_CALL(builder, AddObjs(42));
  adapter.OnHeader(header);
}

TEST_F(NLProblemBuilderTest, OnObj) {
  auto expr = TestNumericExpr(ID);
  StrictMock<MockProblemBuilder::Objective> obj;
  EXPECT_CALL(builder, obj(42)).WillOnce(ReturnRef(obj));
  EXPECT_CALL(obj, set_type(mp::obj::MAX));
  EXPECT_CALL(obj, set_nonlinear_expr(expr));
  adapter.OnObj(42, mp::obj::MAX, expr);
}

TEST_F(NLProblemBuilderTest, OnLinearObjExpr) {
  MockProblemBuilder::Objective obj;
  EXPECT_CALL(builder, obj(42)).WillOnce(ReturnRef(obj));
  auto obj_builder = TestLinearObjBuilder(ID);
  EXPECT_CALL(obj, set_linear_expr(11)).WillOnce(Return(obj_builder));
  EXPECT_EQ(obj_builder, adapter.OnLinearObjExpr(42, 11));
}

TEST_F(NLProblemBuilderTest, AddAlgebraicCons) {
  auto header = mp::NLHeader();
  header.num_algebraic_cons = 42;
  EXPECT_CALL(builder, SetInfo(testing::Ref(header)));
  EXPECT_CALL(builder, AddAlgebraicCons(header.num_algebraic_cons));
  adapter.OnHeader(header);
}

TEST_F(NLProblemBuilderTest, OnAlgebraicCon) {
  MockProblemBuilder::AlgebraicCon con;
  EXPECT_CALL(builder, algebraic_con(42)).WillOnce(ReturnRef(con));
  auto expr = TestNumericExpr(ID);
  EXPECT_CALL(con, set_nonlinear_expr(expr));
  adapter.OnAlgebraicCon(42, expr);
}

TEST_F(NLProblemBuilderTest, OnLinearConExpr) {
  MockProblemBuilder::AlgebraicCon con;
  EXPECT_CALL(builder, algebraic_con(42)).WillOnce(ReturnRef(con));
  auto con_builder = TestLinearConBuilder(ID);
  EXPECT_CALL(con, set_linear_expr(11)).WillOnce(Return(con_builder));
  EXPECT_EQ(con_builder, adapter.OnLinearConExpr(42, 11));
}

TEST_F(NLProblemBuilderTest, OnConBounds) {
  MockProblemBuilder::AlgebraicCon con;
  EXPECT_CALL(builder, algebraic_con(42)).WillOnce(ReturnRef(con));
  EXPECT_CALL(con, set_lb(11));
  EXPECT_CALL(con, set_ub(22));
  adapter.OnConBounds(42, 11, 22);
}

TEST_F(NLProblemBuilderTest, OnInitialValue) {
  MockProblemBuilder::MutVariable var;
  EXPECT_CALL(builder, var(33)).WillOnce(ReturnRef(var));
  EXPECT_CALL(var, set_value(4.4));
  adapter.OnInitialValue(33, 4.4);
}

TEST_F(NLProblemBuilderTest, OnInitialDualValue) {
  MockProblemBuilder::AlgebraicCon con;
  EXPECT_CALL(builder, algebraic_con(33)).WillOnce(ReturnRef(con));
  EXPECT_CALL(con, set_dual(4.4));
  adapter.OnInitialDualValue(33, 4.4);
}

TEST_F(NLProblemBuilderTest, AddLogicalCons) {
  auto header = mp::NLHeader();
  header.num_logical_cons = 42;
  EXPECT_CALL(builder, SetInfo(testing::Ref(header)));
  EXPECT_CALL(builder, AddLogicalCons(header.num_logical_cons));
  adapter.OnHeader(header);
}

TEST_F(NLProblemBuilderTest, OnLogicalCon) {
  MockProblemBuilder::LogicalCon con;
  EXPECT_CALL(builder, logical_con(42)).WillOnce(ReturnRef(con));
  auto expr = TestLogicalExpr(ID);
  EXPECT_CALL(con, set_expr(expr));
  adapter.OnLogicalCon(42, expr);
}

MATCHER_P2(MatchComplInfo, lb, ub, "") {
  return arg.con_lb() == lb && arg.con_ub() == ub;
}

TEST_F(NLProblemBuilderTest, OnComplementarity) {
  double inf = INFINITY;
  EXPECT_CALL(builder, SetComplementarity(66, 77, MatchComplInfo(0, inf)));
  adapter.OnComplementarity(66, 77, mp::ComplInfo(1));
}

TEST_F(NLProblemBuilderTest, AddCommonExprs) {
  auto header = mp::NLHeader();
  header.num_common_exprs_in_both = 42;
  EXPECT_CALL(builder, SetInfo(testing::Ref(header)));
  EXPECT_CALL(builder, AddCommonExprs(header.num_common_exprs_in_both));
  adapter.OnHeader(header);
}

TEST_F(NLProblemBuilderTest, BeginCommonExpr) {
  MockProblemBuilder::MutCommonExpr common_expr;
  EXPECT_CALL(builder, common_expr(42)).WillOnce(ReturnRef(common_expr));
  auto expr_builder = TestLinearExprBuilder(ID);
  EXPECT_CALL(common_expr, set_linear_expr(11)).WillOnce(Return(expr_builder));
  EXPECT_EQ(expr_builder, adapter.BeginCommonExpr(42, 11));
}

TEST_F(NLProblemBuilderTest, EndCommonExpr) {
  MockProblemBuilder::MutCommonExpr common_expr;
  auto expr = TestNumericExpr(ID);
  EXPECT_CALL(builder, common_expr(42)).WillOnce(ReturnRef(common_expr));
  EXPECT_CALL(common_expr, set_nonlinear_expr(expr));
  EXPECT_CALL(common_expr, set_position(11));
  adapter.EndCommonExpr(42, expr, 11);
}

TEST_F(NLProblemBuilderTest, OnColumnSizes) {
  // OnColumnSizes is ignored by default.
  adapter.OnColumnSizes();
}

TEST_F(NLProblemBuilderTest, AddFunctions) {
  auto header = mp::NLHeader();
  header.num_funcs = 42;
  EXPECT_CALL(builder, SetInfo(testing::Ref(header)));
  EXPECT_CALL(builder, AddFunctions(header.num_funcs));
  adapter.OnHeader(header);
}

TEST_F(NLProblemBuilderTest, OnFunction) {
  auto func = TestFunction(ID);
  fmt::StringRef name("f");
  EXPECT_CALL(builder, DefineFunction(42, name, 11, mp::func::SYMBOLIC)).
      WillOnce(Return(func));
  adapter.OnFunction(42, name, 11, mp::func::SYMBOLIC);
}

TEST_F(NLProblemBuilderTest, OnCall) {
  auto func = TestFunction(ID);
  auto call_builder = TestCallExprBuilder(ID);
  EXPECT_CALL(builder, function(42)).WillOnce(Return(func));
  EXPECT_CALL(builder, BeginCall(func, 11)).WillOnce(Return(call_builder));
  adapter.BeginCall(42, 11);
  EXPECT_CALL(builder, EndCall(call_builder)).
      WillOnce(Return(TestNumericExpr(ID)));
  adapter.EndCall(call_builder);
}
