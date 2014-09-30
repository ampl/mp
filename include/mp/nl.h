/*
 .nl format support.

 .nl is a format for representing optimization problems such as linear,
 quadratic, nonlinear, complementarity and constraint programming problems
 in discrete or continuous variables. It is described in the technical report
 "Writing .nl Files" (http://www.cs.sandia.gov/~dmgay/nlwrite.pdf).

 This is a complete reusable C++ implementation of an .nl reader.

 Usage:
   // Read an .nl file:
   ReadNLFile(filename, handler);

   // Read an .nl string:
   ReadNLString(nl_string, handler);

 where handler is an object that receives notifications of problem
 components. See NLHandler and ProblemBuilderToNLAdapter for examples of
 handler classes.

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

#ifndef MP_NL_H_
#define MP_NL_H_

#include "mp/error.h"
#include "mp/os.h"
#include "mp/problem-base.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <limits>
#include <string>

namespace mp {

using fmt::internal::MakeUnsigned;

// A read error with location information.
class ReadError : public Error {
 private:
  std::string filename_;
  int line_;
  int column_;

 public:
  ReadError(fmt::StringRef filename,
      int line, int column, fmt::StringRef message)
  : Error(message), filename_(filename), line_(line), column_(column) {}
  ~ReadError() throw() {}

  const std::string &filename() const { return filename_; }
  int line() const { return line_; }
  int column() const { return column_; }
};

// A read error with information about offset in a binary input.
class BinaryReadError : public Error {
 private:
  std::string filename_;
  std::size_t offset_;

 public:
  BinaryReadError(
      fmt::StringRef filename, std::size_t offset, fmt::StringRef message)
  : Error(message), filename_(filename), offset_(offset) {}
  ~BinaryReadError() throw() {}

  const std::string &filename() const { return filename_; }
  std::size_t offset() const { return offset_; }
};

enum {
  // Maximum number of options in .nl and .sol formats.
  MAX_NL_OPTIONS = 9,
  VBTOL_OPTION   = 1,
  READ_VBTOL     = 3
};

namespace arith {
// Floating-point arithmetic kind.
enum Kind {
  UNKNOWN            = 0,
  // Standard IEEE-754 floating-point:
  IEEE_LITTLE_ENDIAN = 1,
  IEEE_BIG_ENDIAN    = 2,
  // Historical floating-point formats:
  IBM                = 3,
  VAX                = 4,
  CRAY               = 5,
  LAST               = CRAY
};

// Returns floating-point arithmetic kind used on the current system.
Kind GetKind();

inline bool IsIEEE(arith::Kind k) {
  return k == IEEE_LITTLE_ENDIAN || k == IEEE_BIG_ENDIAN;
}
}

// .nl file header.
struct NLHeader : ProblemInfo {
  // .nl file format.
  enum Format { TEXT = 0, BINARY = 1 };
  Format format;

  int num_options;
  int options[MAX_NL_OPTIONS];

  // Extra info for writing solution.
  double ampl_vbtol;

  // Floating-point arithmetic kind used with binary format to check
  // if an .nl file is written using a compatible representation of
  // floating-point numbers. It is not used with text format and normally
  // set to arith::UNKNOWN there.
  arith::Kind arith_kind;

  // Flags: 1 = want output suffixes.
  int flags;

  NLHeader()
    : ProblemInfo(), format(TEXT), num_options(0), ampl_vbtol(0),
      arith_kind(arith::UNKNOWN), flags(0) {
    std::fill(options, options + MAX_NL_OPTIONS, 0);
  }
};

// Writes NLHeader in the .nl file format.
fmt::Writer &operator<<(fmt::Writer &w, const NLHeader &h);

// An .nl handler that ignores all input.
template <typename ExprT>
class NLHandler {
 protected:
  ~NLHandler() {}

 public:
  typedef ExprT Expr;
  typedef Expr NumericExpr;
  typedef Expr LogicalExpr;
  typedef Expr CountExpr;
  typedef Expr Variable;

  // Receives notification of an .nl header.
  void OnHeader(const NLHeader &h) { MP_UNUSED(h); }

  // Receives notification of an objective type and the nonlinear part of
  // an objective expression.
  void OnObj(int index, obj::Type type, NumericExpr expr) {
    MP_UNUSED3(index, type, expr);
  }

  // Receives notification of the nonlinear part of an algebraic constraint
  // expression.
  void OnAlgebraicCon(int index, NumericExpr expr) {
    MP_UNUSED2(index, expr);
  }

  // Receives notification of a logical constraint expression.
  void OnLogicalCon(int index, LogicalExpr expr) {
    MP_UNUSED2(index, expr);
  }

  // Receives notification of the nonlinear part of a defined variable
  // expression.
  void OnDefinedVar(int index, NumericExpr expr, int position) {
    MP_UNUSED3(index, expr, position);
  }

  // Receives notification of a complementarity relation.
  void OnComplement(int con_index, int var_index, int flags) {
    MP_UNUSED3(con_index, var_index, flags);
  }

  struct LinearExprHandler {
    void AddTerm(int var_index, double coef) { MP_UNUSED2(var_index, coef); }
  };

  typedef LinearExprHandler LinearObjHandler;

  // Receives notification of the linear part of an objective expression.
  LinearObjHandler OnLinearObjExpr(int obj_index, int num_linear_terms) {
    MP_UNUSED2(obj_index, num_linear_terms);
    return LinearObjHandler();
  }

  typedef LinearExprHandler LinearConHandler;

  // Receives notification of the linear part of a constraint expression.
  LinearConHandler OnLinearConExpr(int con_index, int num_linear_terms) {
    MP_UNUSED2(con_index, num_linear_terms);
    return LinearConHandler();
  }

  typedef LinearExprHandler LinearVarHandler;

  // Receives notification of the linear part of a defined variable expression.
  LinearVarHandler OnLinearVarExpr(int var_index, int num_linear_terms) {
    MP_UNUSED2(var_index, num_linear_terms);
    return LinearVarHandler();
  }

  // Receives notification of variable bounds.
  void OnVarBounds(int index, double lb, double ub) {
    MP_UNUSED3(index, lb, ub);
  }

  // Receives notification of constraint bounds (ranges).
  void OnConBounds(int index, double lb, double ub) {
    MP_UNUSED3(index, lb, ub);
  }

  // Receives notification of the initial value for a variable.
  void OnInitialValue(int var_index, double value) {
    MP_UNUSED2(var_index, value);
  }

  // Receives notification of the initial value for a dual variable.
  void OnInitialDualValue(int con_index, double value) {
    MP_UNUSED2(con_index, value);
  }

  struct ColumnSizeHandler {
    void Add(int size) { MP_UNUSED(size); }
  };

  // Receives notification of Jacobian column sizes.
  ColumnSizeHandler OnColumnSizes() { return ColumnSizeHandler(); }

  // Receives notification of a function.
  void OnFunction(int index, fmt::StringRef name,
                  int num_args, func::Type type) {
    MP_UNUSED3(index, name, num_args); MP_UNUSED(type);
  }

  struct SuffixHandler {
    void SetValue(int index, int value) { MP_UNUSED2(index, value); }
    void SetValue(int index, double value) { MP_UNUSED2(index, value); }
  };

  // Receives notification of a suffix.
  SuffixHandler OnSuffix(int kind, int num_values, fmt::StringRef name) {
    MP_UNUSED3(kind, num_values, name);
    return SuffixHandler();
  }

  struct ArgHandler {
    void AddArg(Expr arg) { MP_UNUSED(arg); }
  };

  typedef ArgHandler NumericArgHandler;
  typedef ArgHandler LogicalArgHandler;
  typedef ArgHandler CallArgHandler;

  // Receives notification of a numeric constant in a nonlinear expression.
  NumericExpr OnNumericConstant(double value) {
    MP_UNUSED(value);
    return NumericExpr();
  }

  // Receives notification of a variable in a nonlinear expression.
  Variable OnVariable(int var_index) {
    MP_UNUSED(var_index);
    return Variable();
  }

  // Receives notification of a unary expression.
  NumericExpr OnUnary(expr::Kind kind, NumericExpr arg) {
    MP_UNUSED2(kind, arg);
    return NumericExpr();
  }

  // Receives notification of a binary expression.
  NumericExpr OnBinary(expr::Kind kind, NumericExpr lhs, NumericExpr rhs) {
    MP_UNUSED3(kind, lhs, rhs);
    return NumericExpr();
  }

  // Receives notification of an if expression.
  NumericExpr OnIf(LogicalExpr condition,
      NumericExpr true_expr, NumericExpr false_expr) {
    MP_UNUSED3(condition, true_expr, false_expr);
    return NumericExpr();
  }

  struct PLTermHandler {
    void AddSlope(double slope) { MP_UNUSED(slope); }
    void AddBreakpoint(double breakpoint) { MP_UNUSED(breakpoint); }
  };

  // Receives notification of the beginning of a piecewise-linear term.
  PLTermHandler BeginPLTerm(int num_breakpoints) {
    MP_UNUSED(num_breakpoints);
    return PLTermHandler();
  }
  // Receives notification of the end of a piecewise-linear term.
  NumericExpr EndPLTerm(PLTermHandler handler, Variable var) {
    MP_UNUSED2(handler, var);
    return NumericExpr();
  }

  // Receives notification of the beginning of a call expression.
  CallArgHandler BeginCall(int func_index, int num_args) {
    MP_UNUSED2(func_index, num_args);
    return CallArgHandler();
  }
  // Receives notification of the end of a call expression.
  NumericExpr EndCall(CallArgHandler handler) {
    MP_UNUSED(handler);
    return NumericExpr();
  }

  // Receives notification of the beginning of a vararg expression (min or max).
  NumericArgHandler BeginVarArg(expr::Kind kind, int num_args) {
    MP_UNUSED2(kind, num_args);
    return NumericArgHandler();
  }
  // Receives notification of the end of a vararg expression (min or max).
  NumericExpr EndVarArg(NumericArgHandler handler) {
    MP_UNUSED(handler);
    return NumericExpr();
  }

  // Receives notification of the beginning of a sum expression.
  NumericArgHandler BeginSum(int num_args) {
    MP_UNUSED(num_args);
    return NumericArgHandler();
  }
  // Receives notification of the end of a sum expression.
  NumericExpr EndSum(NumericArgHandler handler) {
    MP_UNUSED(handler);
    return NumericExpr();
  }

  // Receives notification of the beginning of a count expression.
  LogicalArgHandler BeginCount(int num_args) {
    MP_UNUSED(num_args);
    return LogicalArgHandler();
  }
  // Receives notification of the end of a count expression.
  CountExpr EndCount(LogicalArgHandler handler) {
    MP_UNUSED(handler);
    return NumericExpr();
  }

  // Receives notification of the beginning of a numberof expression.
  NumericArgHandler BeginNumberOf(int num_args) {
    MP_UNUSED(num_args);
    return NumericArgHandler();
  }
  // Receives notification of the end of a numberof expression.
  NumericExpr EndNumberOf(NumericArgHandler handler) {
    MP_UNUSED(handler);
    return NumericExpr();
  }

  // Receives notification of a logical constant.
  LogicalExpr OnLogicalConstant(bool value) {
    MP_UNUSED(value);
    return LogicalExpr();
  }

  // Receives notification of a logical not expression.
  LogicalExpr OnNot(LogicalExpr arg) {
    MP_UNUSED(arg);
    return LogicalExpr();
  }

  // Receives notification of a binary logical expression.
  LogicalExpr OnBinaryLogical(
      expr::Kind kind, LogicalExpr lhs, LogicalExpr rhs) {
    MP_UNUSED3(kind, lhs, rhs);
    return LogicalExpr();
  }

  // Receives notification of a relational expression.
  LogicalExpr OnRelational(expr::Kind kind, NumericExpr lhs, NumericExpr rhs) {
    MP_UNUSED3(kind, lhs, rhs);
    return LogicalExpr();
  }

  // Receives notification of a logical count expression.
  LogicalExpr OnLogicalCount(expr::Kind kind, NumericExpr lhs, CountExpr rhs) {
    MP_UNUSED3(kind, lhs, rhs);
    return LogicalExpr();
  }

  // Receives notification of an implication expression.
  LogicalExpr OnImplication(
      LogicalExpr condition, LogicalExpr true_expr, LogicalExpr false_expr) {
    MP_UNUSED3(condition, true_expr, false_expr);
    return LogicalExpr();
  }

  // Receives notification of the beginning of an iterated logical expression.
  LogicalArgHandler BeginIteratedLogical(expr::Kind kind, int num_args) {
    MP_UNUSED2(kind, num_args);
    return LogicalArgHandler();
  }
  // Receives notification of the end of an iterated logical expression.
  LogicalExpr EndIteratedLogical(LogicalArgHandler handler) {
    MP_UNUSED(handler);
    return LogicalExpr();
  }

  typedef ArgHandler AllDiffArgHandler;

  // Receives notification of the beginning of an alldiff expression.
  AllDiffArgHandler BeginAllDiff(int num_args) {
    MP_UNUSED(num_args);
    return AllDiffArgHandler();
  }
  // Receives notification of the end of an alldiff expression.
  LogicalExpr EndAllDiff(AllDiffArgHandler handler) {
    MP_UNUSED(handler);
    return LogicalExpr();
  }

  // Receives notification of a string literal.
  Expr OnStringLiteral(fmt::StringRef value) {
    MP_UNUSED(value);
    return Expr();
  }
};

// Adapts ProblemBuilder for use as an .nl handler.
template <typename ProblemBuilder>
class ProblemBuilderToNLAdapter {
 private:
  ProblemBuilder &builder_;

 public:
  explicit ProblemBuilderToNLAdapter(ProblemBuilder &builder)
    : builder_(builder) {}

  ProblemBuilder &builder() { return builder_; }

  typedef typename ProblemBuilder::Expr Expr;
  typedef typename ProblemBuilder::NumericExpr NumericExpr;
  typedef typename ProblemBuilder::LogicalExpr LogicalExpr;
  typedef typename ProblemBuilder::CountExpr CountExpr;
  typedef typename ProblemBuilder::Variable Variable;

  // Receives notification of an .nl header.
  void OnHeader(const NLHeader &h) { builder_.SetInfo(h); }

  // Receives notification of an objective type and the nonlinear part of
  // an objective expression.
  void OnObj(int index, obj::Type type, NumericExpr expr) {
    builder_.SetObj(index, type, expr);
  }

  // Receives notification of the nonlinear part of an algebraic constraint
  // expression.
  void OnAlgebraicCon(int index, NumericExpr expr) {
    builder_.SetCon(index, expr);
  }

  // Receives notification of a logical constraint expression.
  void OnLogicalCon(int index, LogicalExpr expr) {
    builder_.SetLogicalCon(index, expr);
  }

  // Receives notification of the nonlinear part of a defined variable
  // expression.
  void OnDefinedVar(int index, NumericExpr expr, int position) {
    builder_.SetVar(index, expr, position);
  }

  // Receives notification of a complementarity relation.
  void OnComplement(int con_index, int var_index, int flags) {
    builder_.SetComplement(con_index, var_index, flags);
  }

  typedef typename ProblemBuilder::LinearObjBuilder LinearObjHandler;

  // Receives notification of the linear part of an objective expression.
  LinearObjHandler OnLinearObjExpr(int obj_index, int num_linear_terms) {
    return builder_.GetLinearObjBuilder(obj_index, num_linear_terms);
  }

  typedef typename ProblemBuilder::LinearConBuilder LinearConHandler;

  // Receives notification of the linear part of a constraint expression.
  LinearConHandler OnLinearConExpr(int con_index, int num_linear_terms) {
    return builder_.GetLinearConBuilder(con_index, num_linear_terms);
  }

  typedef typename ProblemBuilder::LinearVarBuilder LinearVarHandler;

  // Receives notification of the linear part of a defined variable expression.
  LinearVarHandler OnLinearVarExpr(int var_index, int num_linear_terms) {
    return builder_.GetLinearVarBuilder(var_index, num_linear_terms);
  }

  // Receives notification of variable bounds.
  void OnVarBounds(int index, double lb, double ub) {
    builder_.SetVarBounds(index, lb, ub);
  }

  // Receives notification of constraint bounds (ranges).
  void OnConBounds(int index, double lb, double ub) {
    builder_.SetConBounds(index, lb, ub);
  }

  // Receives notification of the initial value for a variable.
  void OnInitialValue(int var_index, double value) {
    builder_.SetInitialValue(var_index, value);
  }

  // Receives notification of the initial value for a dual variable.
  void OnInitialDualValue(int con_index, double value) {
    builder_.SetInitialDualValue(con_index, value);
  }

  typedef typename ProblemBuilder::ColumnSizeHandler ColumnSizeHandler;

  // Receives notification of Jacobian column sizes.
  ColumnSizeHandler OnColumnSizes() {
    return builder_.GetColumnSizeHandler();
  }

  // Receives notification of a function.
  void OnFunction(int index, fmt::StringRef name,
                  int num_args, func::Type type) {
    builder_.SetFunction(index, name, num_args, type);
  }

  typedef typename ProblemBuilder::SuffixHandler SuffixHandler;

  // Receives notification of a suffix.
  SuffixHandler OnSuffix(int kind, int num_values, fmt::StringRef name) {
    return builder_.AddSuffix(kind, num_values, name);
  }

  typedef typename ProblemBuilder::NumericArgHandler NumericArgHandler;
  typedef typename ProblemBuilder::LogicalArgHandler LogicalArgHandler;
  typedef typename ProblemBuilder::CallArgHandler CallArgHandler;

  // Receives notification of a numeric constant in a nonlinear expression.
  NumericExpr OnNumericConstant(double value) {
    return builder_.MakeNumericConstant(value);
  }

  // Receives notification of a variable in a nonlinear expression.
  Variable OnVariable(int var_index) {
    return builder_.MakeVariable(var_index);
  }

  // Receives notification of a unary expression.
  NumericExpr OnUnary(expr::Kind kind, NumericExpr arg) {
    return builder_.MakeUnary(kind, arg);
  }

  // Receives notification of a binary expression.
  NumericExpr OnBinary(expr::Kind kind, NumericExpr lhs, NumericExpr rhs) {
    return builder_.MakeBinary(kind, lhs, rhs);
  }

  // Receives notification of an if expression.
  NumericExpr OnIf(LogicalExpr condition,
      NumericExpr true_expr, NumericExpr false_expr) {
    return builder_.MakeIf(condition, true_expr, false_expr);
  }

  typedef typename ProblemBuilder::PLTermHandler PLTermHandler;

  // Receives notification of the beginning of a piecewise-linear term.
  PLTermHandler BeginPLTerm(int num_breakpoints) {
    return builder_.BeginPLTerm(num_breakpoints);
  }
  // Receives notification of the end of a piecewise-linear term.
  NumericExpr EndPLTerm(PLTermHandler handler, Variable var) {
    return builder_.EndPLTerm(handler, var);
  }

  // Receives notification of the beginning of a call expression.
  CallArgHandler BeginCall(int func_index, int num_args) {
    return builder_.BeginCall(func_index, num_args);
  }
  // Receives notification of the end of a call expression.
  NumericExpr EndCall(CallArgHandler handler) {
    return builder_.EndCall(handler);
  }

  // Receives notification of the beginning of a vararg expression (min or max).
  NumericArgHandler BeginVarArg(expr::Kind kind, int num_args) {
    return builder_.BeginVarArg(kind, num_args);
  }
  // Receives notification of the end of a vararg expression (min or max).
  NumericExpr EndVarArg(NumericArgHandler handler) {
    return builder_.EndVarArg(handler);
  }

  // Receives notification of the beginning of a sum expression.
  NumericArgHandler BeginSum(int num_args) {
    return builder_.BeginSum(num_args);
  }
  // Receives notification of the end of a sum expression.
  NumericExpr EndSum(NumericArgHandler handler) {
    return builder_.EndSum(handler);
  }

  // Receives notification of the beginning of a count expression.
  LogicalArgHandler BeginCount(int num_args) {
    return builder_.BeginCount(num_args);
  }
  // Receives notification of the end of a count expression.
  CountExpr EndCount(LogicalArgHandler handler) {
    return builder_.EndCount(handler);
  }

  // Receives notification of the beginning of a numberof expression.
  NumericArgHandler BeginNumberOf(int num_args) {
    return builder_.BeginNumberOf(num_args);
  }
  // Receives notification of the end of a numberof expression.
  NumericExpr EndNumberOf(NumericArgHandler handler) {
    return builder_.EndNumberOf(handler);
  }

  // Receives notification of a logical constant.
  LogicalExpr OnLogicalConstant(bool value) {
    return builder_.MakeLogicalConstant(value);
  }

  // Receives notification of a logical not expression.
  LogicalExpr OnNot(LogicalExpr arg) {
    return builder_.MakeNot(arg);
  }

  // Receives notification of a binary logical expression.
  LogicalExpr OnBinaryLogical(
      expr::Kind kind, LogicalExpr lhs, LogicalExpr rhs) {
    return builder_.MakeBinaryLogical(kind, lhs, rhs);
  }

  // Receives notification of a relational expression.
  LogicalExpr OnRelational(
      expr::Kind kind, NumericExpr lhs, NumericExpr rhs) {
    return builder_.MakeRelational(kind, lhs, rhs);
  }

  // Receives notification of a logical count expression.
  LogicalExpr OnLogicalCount(
      expr::Kind kind, NumericExpr lhs, CountExpr rhs) {
    return builder_.MakeLogicalCount(kind, lhs, rhs);
  }

  // Receives notification of an implication expression.
  LogicalExpr OnImplication(
      LogicalExpr condition, LogicalExpr true_expr, LogicalExpr false_expr) {
    return builder_.MakeImplication(condition, true_expr, false_expr);
  }

  // Receives notification of the beginning of an iterated logical expression.
  LogicalArgHandler BeginIteratedLogical(expr::Kind kind, int num_args) {
    return builder_.BeginIteratedLogical(kind, num_args);
  }
  // Receives notification of the end of an iterated logical expression.
  LogicalExpr EndIteratedLogical(LogicalArgHandler handler) {
    return builder_.EndIteratedLogical(handler);
  }

  typedef typename ProblemBuilder::AllDiffArgHandler AllDiffArgHandler;

  // Receives notification of the beginning of an alldiff expression.
  AllDiffArgHandler BeginAllDiff(int num_args) {
    return builder_.BeginAllDiff(num_args);
  }
  // Receives notification of the end of an alldiff expression.
  LogicalExpr EndAllDiff(AllDiffArgHandler handler) {
    return builder_.EndAllDiff(handler);
  }

  // Receives notification of a string literal.
  Expr OnStringLiteral(fmt::StringRef value) {
    return builder_.MakeStringLiteral(value);
  }
};

namespace internal {

class ReaderBase {
 protected:
  const char *ptr_, *start_, *end_;
  const char *token_;  // start of the current token
  std::string name_;

  ~ReaderBase() {}

 public:
  ReaderBase(fmt::StringRef data, fmt::StringRef name);

  char ReadChar() {
    token_ = ptr_;
    return *ptr_++;
  }

  bool IsEOF() const { return ptr_ == end_ + 1; }
};

class TextReader : public ReaderBase {
 private:
  const char *line_start_;
  int line_;

  // Reads an integer without a sign.
  // Int: signed or unsigned integer type.
  template <typename Int>
  bool ReadIntWithoutSign(Int& value) {
    char c = *ptr_;
    if (c < '0' || c > '9')
      return false;
    typedef typename MakeUnsigned<Int>::Type UInt;
    UInt result = 0;
    do {
      UInt new_result = result * 10 + (c - '0');
      if (new_result < result)
        ReportError("number is too big");
      result = new_result;
      c = *++ptr_;
    } while (c >= '0' && c <= '9');
    UInt max = std::numeric_limits<Int>::max();
    if (result > max)
      ReportError("number is too big");
    value = result;
    return true;
  }

  template <typename Int>
  bool DoReadOptionalInt(Int &value) {
    SkipSpace();
    char sign = *ptr_;
    if (sign == '+' || sign == '-')
      ++ptr_;
    typedef typename MakeUnsigned<Int>::Type UInt;
    UInt result = 0;
    if (!ReadIntWithoutSign<UInt>(result))
      return false;
    UInt max = std::numeric_limits<Int>::max();
    if (result > max && !(sign == '-' && result == max + 1))
      ReportError("number is too big");
    value = sign != '-' ? result : 0 - result;
    return true;
  }

  // Reads a nonnegative integer and checks that adding it to accumulator
  // doesn't overflow.
  int ReadUInt(int &accumulator) {
    int value = ReadUInt();
    if (accumulator > std::numeric_limits<int>::max() - value)
      ReportError("integer overflow");
    accumulator += value;
    return value;
  }

  template <typename Int>
  Int ReadUInt() {
    SkipSpace();
    Int value = 0;
    if (!ReadIntWithoutSign(value))
      ReportError("expected unsigned integer");
    return value;
  }

  bool ReadOptionalInt(int &value) { return DoReadOptionalInt(value); }

  bool ReadOptionalUInt(int &value) {
    SkipSpace();
    return ReadIntWithoutSign(value);
  }

  bool ReadOptionalDouble(double &value);

  void DoReportError(
      const char *loc, fmt::StringRef format_str,
      const fmt::ArgList &args = fmt::ArgList());

  void SkipSpace() {
    while (std::isspace(*ptr_) && *ptr_ != '\n')
      ++ptr_;
    token_ = ptr_;
  }

 public:
  TextReader(fmt::StringRef data, fmt::StringRef name);

  void ReportError(fmt::StringRef format_str, const fmt::ArgList &args) {
    DoReportError(token_, format_str, args);
  }
  FMT_VARIADIC(void, ReportError, fmt::StringRef)

  void ReadTillEndOfLine() {
    while (char c = *ptr_) {
      ++ptr_;
      if (c == '\n') {
        line_start_ = ptr_;
        ++line_;
        return;
      }
    }
    DoReportError(ptr_, "expected newline");
  }

  template <typename Int>
  Int ReadInt() {
    Int value = 0;
    if (!DoReadOptionalInt(value))
      ReportError("expected integer");
    return value;
  }

  int ReadUInt() { return ReadUInt<int>(); }

  double ReadDouble() {
    SkipSpace();
    char *end = 0;
    double value = 0;
    if (*ptr_ != '\n')
      value = std::strtod(ptr_, &end);
    if (!end || ptr_ == end)
      ReportError("expected double");
    ptr_ = end;
    return value;
  }

  fmt::StringRef ReadString();

  // Reads a function or suffix name.
  fmt::StringRef ReadName();

  // Reads an .nl file header. The header is always in text format, so this
  // function doesn't have a counterpart in BinaryReader.
  void ReadHeader(NLHeader &header);
};

class BinaryReader : public ReaderBase {
 private:
  // Reads length chars.
  const char *Read(int length) {
    if (end_ - ptr_ < length) {
      token_ = end_;
      ReportError("unexpected end of file");
    }
    const char *start = ptr_;
    ptr_ += length;
    return start;
  }

 public:
  explicit BinaryReader(const ReaderBase &base) : ReaderBase(base) {}

  void ReportError(fmt::StringRef format_str, const fmt::ArgList &args);
  FMT_VARIADIC(void, ReportError, fmt::StringRef)

  void ReadTillEndOfLine() {
    // Do nothing.
  }

  template <typename Int>
  Int ReadInt() {
    token_ = ptr_;
    return *reinterpret_cast<const Int*>(Read(sizeof(Int)));
  }

  int ReadUInt() {
    int value = ReadInt<int>();
    if (value < 0)
      ReportError("expected unsigned integer");
    return value;
  }

  double ReadDouble() {
    token_ = ptr_;
    return *reinterpret_cast<const double*>(Read(sizeof(double)));
  }

  fmt::StringRef ReadString() {
    int length = ReadUInt();
    return fmt::StringRef(length != 0 ? Read(length) : 0, length);
  }

  // Reads a function or suffix name.
  fmt::StringRef ReadName() { return ReadString(); }
};

// An .nl file reader.
// Handler: a class implementing the ProblemHandler concept that receives
//          notifications of problem components
template <typename Reader, typename Handler>
class NLReader {
 private:
  Reader &reader_;
  NLHeader &header_;
  Handler &handler_;
  int total_num_vars_;  // Total number of variables including defined ones.

  typedef typename Handler::NumericExpr NumericExpr;
  typedef typename Handler::LogicalExpr LogicalExpr;
  typedef typename Handler::Variable Variable;

  double ReadConstant(char code);
  double ReadConstant() { return ReadConstant(reader_.ReadChar()); }

  // Reads a nonnegative integer and checks that it is less than ub.
  // ub is unsigned so that it can hold value INT_MAX + 1u.
  int ReadUInt(unsigned ub) {
    int value = reader_.ReadUInt();
    unsigned unsigned_value = value;
    if (unsigned_value >= ub)
      reader_.ReportError("integer {} out of bounds", value);
    return value;
  }

  // Reads a nonnegative integer and checks that it is in the range [lb, ub).
  int ReadUInt(unsigned lb, unsigned ub) {
    int value = reader_.ReadUInt();
    unsigned unsigned_value = value;
    if (unsigned_value < lb || unsigned_value >= ub)
      reader_.ReportError("integer {} out of bounds", value);
    return value;
  }

  // Minimum number of arguments for an iterated expression that has a
  // binary counterpart. Examples: sum (+), forall (&&), exists (||).
  enum {MIN_ITER_ARGS = 3};

  int ReadNumArgs(int min_args = MIN_ITER_ARGS) {
    int num_args = reader_.ReadUInt();
    if (num_args < min_args)
      reader_.ReportError("too few arguments");
    return num_args;
  }

  Variable DoReadVariable() {
    int var_index = ReadUInt(total_num_vars_);
    reader_.ReadTillEndOfLine();
    return handler_.OnVariable(var_index);
  }

  Variable ReadVariable() {
    if (reader_.ReadChar() != 'v')
      reader_.ReportError("expected variable");
    return DoReadVariable();
  }

  template <typename ExprReader, typename ArgHandler>
  void ReadArgs(int num_args, ArgHandler &arg_handler) {
    reader_.ReadTillEndOfLine();
    ExprReader expr_reader;
    for (int i = 0; i < num_args; ++i)
      arg_handler.AddArg(expr_reader.Read(*this));
  }

  int ReadOpCode() {
    int opcode = reader_.ReadUInt();
    if (opcode > expr::MAX_OPCODE)
      reader_.ReportError("invalid opcode {}", opcode);
    reader_.ReadTillEndOfLine();
    return opcode;
  }

  typename Handler::CountExpr ReadCountExpr() {
    int num_args = ReadNumArgs(1);
    typename Handler::LogicalArgHandler args = handler_.BeginCount(num_args);
    ReadArgs<LogicalExprReader>(num_args, args);
    return handler_.EndCount(args);
  }

  // Helper structs to provide a uniform interface to Read{Numeric,Logical}Expr
  // since it is not possible to overload on expression type as NumericExpr
  // and LogicalExpr can be the same type.
  struct NumericExprReader {
    typedef NumericExpr Expr;
    Expr Read(NLReader &r) const { return r.ReadNumericExpr(); }
  };
  struct LogicalExprReader {
    typedef LogicalExpr Expr;
    Expr Read(NLReader &r) const { return r.ReadLogicalExpr(); }
  };

  // A helper struct used to make sure that the arguments to a binary
  // expression are read in the correct order and avoid errors of the form:
  //   MakeBinary(opcode, ReadNumericExpr(), ReadNumericExpr())
  // The above code is incorrect as the order of evaluation of arguments is
  // unspecified.
  template <typename ExprReader = NumericExprReader>
  struct BinaryArgReader {
    typename ExprReader::Expr lhs;
    typename ExprReader::Expr rhs;
    BinaryArgReader(NLReader &r)
      : lhs(ExprReader().Read(r)), rhs(ExprReader().Read(r)) {}
  };

  // Reads a numeric expression.
  NumericExpr ReadNumericExpr() { return ReadNumericExpr(reader_.ReadChar()); }
  NumericExpr ReadNumericExpr(char code);
  NumericExpr ReadNumericExpr(int opcode);

  // Reads a logical expression.
  LogicalExpr ReadLogicalExpr();
  LogicalExpr ReadLogicalExpr(int opcode);

  enum ItemType { VAR, OBJ, CON, PROB };

  template <ItemType T>
  class ItemHandler {
   protected:
    NLReader &reader_;

   public:
    static const ItemType TYPE = T;
    explicit ItemHandler(NLReader &r) : reader_(r) {}
  };

  struct VarHandler : ItemHandler<VAR> {
    explicit VarHandler(NLReader &r) : ItemHandler<VAR>(r) {}

    int num_items() const { return this->reader_.header_.num_vars; }

    void SetBounds(int index, double lb, double ub) {
      this->reader_.handler_.OnVarBounds(index, lb, ub);
    }
    void SetInitialValue(int index, double value) {
      this->reader_.handler_.OnInitialValue(index, value);
    }
  };

  struct ObjHandler : ItemHandler<OBJ> {
    explicit ObjHandler(NLReader &r) : ItemHandler<OBJ>(r) {}

    int num_items() const { return this->reader_.header_.num_objs; }

    typename Handler::LinearObjHandler OnLinearExpr(int index, int num_terms) {
      return this->reader_.handler_.OnLinearObjExpr(index, num_terms);
    }
  };

  struct AlgebraicConHandler : ItemHandler<CON> {
    explicit AlgebraicConHandler(NLReader &r) : ItemHandler<CON>(r) {}

    int num_items() const { return this->reader_.header_.num_algebraic_cons; }

    typename Handler::LinearConHandler OnLinearExpr(int index, int num_terms) {
      return this->reader_.handler_.OnLinearConExpr(index, num_terms);
    }

    void SetBounds(int index, double lb, double ub) {
      this->reader_.handler_.OnConBounds(index, lb, ub);
    }
    void SetInitialValue(int index, double value) {
      this->reader_.handler_.OnInitialDualValue(index, value);
    }
  };

  struct ConHandler : ItemHandler<CON> {
    explicit ConHandler(NLReader &r) : ItemHandler<CON>(r) {}

    int num_items() const {
      return this->reader_.header_.num_algebraic_cons +
          this->reader_.header_.num_logical_cons;
    }
  };

  struct ProblemHandler : ItemHandler<PROB> {
    explicit ProblemHandler(NLReader &r) : ItemHandler<PROB>(r) {}

    // An .nl file always contains one problem.
    int num_items() const { return 1; }
  };

  // Reads the linear part of an objective or constraint expression.
  template <typename LinearHandler>
  void ReadLinearExpr();

  template <typename LinearHandler>
  void ReadLinearExpr(int num_terms, LinearHandler linear_expr);

  // Reads variable or constraint bounds.
  template <typename BoundHandler>
  void ReadBounds();

  // Reads column sizes, numbers of nonzeros in the first num_var − 1
  // columns of the Jacobian sparsity matrix.
  template <bool CUMULATIVE>
  void ReadColumnSizes();

  // Reads initial values for primal or dual variables.
  template <typename ValueHandler>
  void ReadInitialValues();

  template <typename SuffixHandler>
  void ReadSuffix(int kind);

 public:
  NLReader(Reader &reader, NLHeader &header, Handler &handler)
    : reader_(reader), header_(header), handler_(handler), total_num_vars_(0) {}

  void Read();
};

template <typename Reader, typename Handler>
double NLReader<Reader, Handler>::ReadConstant(char code) {
  double value = 0;
  switch (code) {
  case 'n':
    value = reader_.ReadDouble();
    break;
  case 's':
    value = reader_.template ReadInt<short>();
    break;
  case 'l':
    value = reader_.template ReadInt<long>();
    break;
  default:
    reader_.ReportError("expected constant");
  }
  reader_.ReadTillEndOfLine();
  return value;
}

template <typename Reader, typename Handler>
typename Handler::NumericExpr
    NLReader<Reader, Handler>::ReadNumericExpr(char code) {
  switch (code) {
  case 'f': {
    int func_index = ReadUInt(header_.num_funcs);
    int num_args = reader_.ReadUInt();
    reader_.ReadTillEndOfLine();
    typename Handler::CallArgHandler args =
        handler_.BeginCall(func_index, num_args);
    for (int i = 0; i < num_args; ++i) {
      char c = reader_.ReadChar();
      switch (c) {
      case 'h':
        args.AddArg(handler_.OnStringLiteral(reader_.ReadString()));
        break;
      case 'o': {
        int opcode = ReadOpCode();
        if (opcode == expr::IFSYM) {
          // TODO: read symbolic if
        }
        args.AddArg(ReadNumericExpr(opcode));
        break;
      }
      default:
        args.AddArg(ReadNumericExpr(c));
      }
    }
    return handler_.EndCall(args);
  }
  case 'n': case 'l': case 's':
    return handler_.OnNumericConstant(ReadConstant(code));
  case 'o':
    return ReadNumericExpr(ReadOpCode());
  case 'v':
    return DoReadVariable();
  default:
    reader_.ReportError("expected expression");
  }
  return NumericExpr();
}

template <typename Reader, typename Handler>
typename Handler::NumericExpr
    NLReader<Reader, Handler>::ReadNumericExpr(int opcode) {
  const expr::OpCodeInfo &info = expr::GetOpCodeInfo(opcode);
  expr::Kind kind = info.kind;
  switch (info.first_kind) {
  case expr::FIRST_UNARY:
    return handler_.OnUnary(kind, ReadNumericExpr());
  case expr::FIRST_BINARY: {
    BinaryArgReader<> args(*this);
    return handler_.OnBinary(kind, args.lhs, args.rhs);
  }
  case expr::IF: {
    LogicalExpr condition = ReadLogicalExpr();
    NumericExpr true_expr = ReadNumericExpr();
    NumericExpr false_expr = ReadNumericExpr();
    return handler_.OnIf(condition, true_expr, false_expr);
  }
  case expr::PLTERM: {
    int num_slopes = reader_.ReadUInt();
    if (num_slopes <= 1)
      reader_.ReportError("too few slopes in piecewise-linear term");
    reader_.ReadTillEndOfLine();
    typename Handler::PLTermHandler pl_handler =
        handler_.BeginPLTerm(num_slopes - 1);
    for (int i = 0; i < num_slopes - 1; ++i) {
      pl_handler.AddSlope(ReadConstant());
      pl_handler.AddBreakpoint(ReadConstant());
    }
    pl_handler.AddSlope(ReadConstant());
    return handler_.EndPLTerm(pl_handler, ReadVariable());
  }
  case expr::FIRST_VARARG: {
    int num_args = ReadNumArgs(1);
    typename Handler::NumericArgHandler args =
        handler_.BeginVarArg(kind, num_args);
    ReadArgs<NumericExprReader>(num_args, args);
    return handler_.EndVarArg(args);
  }
  case expr::SUM: {
    int num_args = ReadNumArgs();
    typename Handler::NumericArgHandler args = handler_.BeginSum(num_args);
    ReadArgs<NumericExprReader>(num_args, args);
    return handler_.EndSum(args);
  }
  case expr::COUNT:
    return ReadCountExpr();
  case expr::NUMBEROF: {
    int num_args = ReadNumArgs(1);
    typename Handler::NumericArgHandler args = handler_.BeginNumberOf(num_args);
    ReadArgs<NumericExprReader>(num_args, args);
    return handler_.EndNumberOf(args);
  }
  case expr::NUMBEROF_SYM: {
    // TODO: read symbolic numberof
  }
  default:
    reader_.ReportError("expected numeric expression opcode");
  }
  return NumericExpr();
}

template <typename Reader, typename Handler>
typename Handler::LogicalExpr NLReader<Reader, Handler>::ReadLogicalExpr() {
  switch (char c = reader_.ReadChar()) {
  case 'n': case 'l': case 's':
    return handler_.OnLogicalConstant(ReadConstant(c) != 0);
  case 'o':
    return ReadLogicalExpr(ReadOpCode());
  }
  reader_.ReportError("expected logical expression");
  return LogicalExpr();
}

template <typename Reader, typename Handler>
typename Handler::LogicalExpr
    NLReader<Reader, Handler>::ReadLogicalExpr(int opcode) {
  const expr::OpCodeInfo &info = expr::GetOpCodeInfo(opcode);
  expr::Kind kind = info.kind;
  switch (info.first_kind) {
  case expr::NOT:
    return handler_.OnNot(ReadLogicalExpr());
  case expr::FIRST_BINARY_LOGICAL: {
    BinaryArgReader<LogicalExprReader> args(*this);
    return handler_.OnBinaryLogical(kind, args.lhs, args.rhs);
  }
  case expr::FIRST_RELATIONAL: {
    BinaryArgReader<> args(*this);
    return handler_.OnRelational(kind, args.lhs, args.rhs);
  }
  case expr::FIRST_LOGICAL_COUNT: {
    NumericExpr lhs = ReadNumericExpr();
    char c = reader_.ReadChar();
    if (c != 'o' || expr::GetOpCodeInfo(ReadOpCode()).kind != expr::COUNT)
      reader_.ReportError("expected count expression");
    return handler_.OnLogicalCount(kind, lhs, ReadCountExpr());
  }
  case expr::IMPLICATION: {
    LogicalExpr condition = ReadLogicalExpr();
    LogicalExpr true_expr = ReadLogicalExpr();
    LogicalExpr false_expr = ReadLogicalExpr();
    return handler_.OnImplication(condition, true_expr, false_expr);
  }
  case expr::FIRST_ITERATED_LOGICAL: {
    int num_args = ReadNumArgs();
    typename Handler::LogicalArgHandler args =
        handler_.BeginIteratedLogical(kind, num_args);
    ReadArgs<LogicalExprReader>(num_args, args);
    return handler_.EndIteratedLogical(args);
  }
  case expr::ALLDIFF: {
    int num_args = ReadNumArgs(1);
    typename Handler::AllDiffArgHandler args = handler_.BeginAllDiff(num_args);
    ReadArgs<NumericExprReader>(num_args, args);
    return handler_.EndAllDiff(args);
  }
  default:
    reader_.ReportError("expected logical expression opcode");
  }
  return LogicalExpr();
}

template <typename Reader, typename Handler>
template <typename LinearHandler>
void NLReader<Reader, Handler>::ReadLinearExpr() {
  LinearHandler lh(*this);
  int index = ReadUInt(lh.num_items());
  int num_terms = ReadUInt(1, total_num_vars_ + 1u);
  reader_.ReadTillEndOfLine();
  ReadLinearExpr(num_terms, lh.OnLinearExpr(index, num_terms));
}

template <typename Reader, typename Handler>
template <typename LinearHandler>
void NLReader<Reader, Handler>::ReadLinearExpr(
    int num_terms, LinearHandler linear_expr) {
  for (int i = 0; i < num_terms; ++i) {
    int var_index = ReadUInt(total_num_vars_);
    double coef = reader_.ReadDouble();
    reader_.ReadTillEndOfLine();
    linear_expr.AddTerm(var_index, coef);
  }
}

template <typename Reader, typename Handler>
template <typename BoundHandler>
void NLReader<Reader, Handler>::ReadBounds() {
  enum BoundType {
    RANGE,     // Both lower and upper bounds: l <= body <= u.
    UPPER,     // Only upper bound: body <= u.
    LOWER,     // Only lower bound: l <= body.
    FREE,      // No constraints on body (free variable or constraint).
    CONSTANT,  // Equal to constant: body = c.
    COMPL      // Body complements variable v[i - 1].
  };
  reader_.ReadTillEndOfLine();
  double lb = 0, ub = 0;
  BoundHandler bh(*this);
  int num_bounds = bh.num_items();
  double infinity = std::numeric_limits<double>::infinity();
  for (int i = 0; i < num_bounds; ++i) {
    switch (reader_.ReadChar() - '0') {
    case RANGE:
      lb = reader_.ReadDouble();
      ub = reader_.ReadDouble();
      break;
    case UPPER:
      lb = -infinity;
      ub = reader_.ReadDouble();
      break;
    case LOWER:
      lb = reader_.ReadDouble();
      ub = infinity;
      break;
    case FREE:
      lb = -infinity;
      ub =  infinity;
      break;
    case CONSTANT:
      lb = ub = reader_.ReadDouble();
      break;
    case COMPL:
      if (BoundHandler::TYPE == CON) {
        int flags = reader_.template ReadInt<int>();
        int var_index = reader_.ReadUInt();
        // Don't use NLReader::ReadUInt(int, int) as num_vars + 1 may overflow.
        if (var_index == 0 || var_index > header_.num_vars)
          reader_.ReportError("integer {} out of bounds", var_index);
        --var_index;
        int mask = comp::INF_LB | comp::INF_UB;
        handler_.OnComplement(i, var_index, flags & mask);
        reader_.ReadTillEndOfLine();
        continue;
      }
      // Fall through as COMPL bound type is invalid for variables.
    default:
      reader_.ReportError("expected bound");
    }
    reader_.ReadTillEndOfLine();
    bh.SetBounds(i, lb, ub);
  }
}

template <typename Reader, typename Handler>
template <bool CUMULATIVE>
void NLReader<Reader, Handler>::ReadColumnSizes() {
  int num_sizes = header_.num_vars - 1;
  if (reader_.ReadUInt() != num_sizes)
    reader_.ReportError("expected {}", num_sizes);
  reader_.ReadTillEndOfLine();
  typename Handler::ColumnSizeHandler size_handler = handler_.OnColumnSizes();
  int prev_size = 0;
  for (int i = 0; i < num_sizes; ++i) {
    int size = reader_.ReadUInt();
    if (CUMULATIVE) {
      if (size < prev_size)
        reader_.ReportError("invalid column offset");
      size -= prev_size;
      prev_size += size;
    }
    size_handler.Add(size);
    reader_.ReadTillEndOfLine();
  }
}

template <typename Reader, typename Handler>
template <typename ValueHandler>
void NLReader<Reader, Handler>::ReadInitialValues() {
  int num_values = reader_.ReadUInt();
  ValueHandler vh(*this);
  if (num_values > vh.num_items())
    reader_.ReportError("too many initial values");
  reader_.ReadTillEndOfLine();
  for (int i = 0; i < num_values; ++i) {
    int index = ReadUInt(vh.num_items());
    vh.SetInitialValue(index, reader_.ReadDouble());
    reader_.ReadTillEndOfLine();
  }
}

template <typename Reader, typename Handler>
template <typename SuffixHandler>
void NLReader<Reader, Handler>::ReadSuffix(int kind) {
  SuffixHandler sh(*this);
  int num_items = sh.num_items();
  int num_values = ReadUInt(1, num_items + 1);
  fmt::StringRef name = reader_.ReadName();
  reader_.ReadTillEndOfLine();
  bool is_float = (kind & suf::FLOAT) != 0;
  typename Handler::SuffixHandler
      suffix_handler = handler_.OnSuffix(kind, num_values, name);
  for (int i = 0; i < num_values; ++i) {
    int index = ReadUInt(num_items);
    if (is_float)
      suffix_handler.SetValue(index, reader_.ReadDouble());
    else
      suffix_handler.SetValue(index, reader_.template ReadInt<int>());
    reader_.ReadTillEndOfLine();
  }
}

template <typename Reader, typename Handler>
void NLReader<Reader, Handler>::Read() {
  // TextReader::ReadHeader checks that this doesn't overflow.
  total_num_vars_ = header_.num_vars +
      header_.num_common_exprs_in_both +
      header_.num_common_exprs_in_cons +
      header_.num_common_exprs_in_objs +
      header_.num_common_exprs_in_single_cons +
      header_.num_common_exprs_in_single_objs;
  for (;;) {
    char c = reader_.ReadChar();
    switch (c) {
    case 'C': {
      // Nonlinear part of an algebraic constraint body.
      int index = ReadUInt(header_.num_algebraic_cons);
      reader_.ReadTillEndOfLine();
      handler_.OnAlgebraicCon(index, ReadNumericExpr());
      break;
    }
    case 'L': {
      // Logical constraint expression.
      int index = ReadUInt(header_.num_logical_cons);
      reader_.ReadTillEndOfLine();
      handler_.OnLogicalCon(index, ReadLogicalExpr());
      break;
    }
    case 'O': {
      // Objective type and nonlinear part of an objective expression.
      int index = ReadUInt(header_.num_objs);
      int obj_type = reader_.ReadUInt();
      reader_.ReadTillEndOfLine();
      handler_.OnObj(index, obj_type != 0 ? obj::MAX : obj::MIN,
                     ReadNumericExpr());
      break;
    }
    case 'V': {
      // Defined variable definition (must precede V, C, L, O segments
      // where used).
      int var_index = ReadUInt(header_.num_vars, total_num_vars_);
      int num_linear_terms = reader_.ReadUInt();
      int position = reader_.ReadUInt();
      reader_.ReadTillEndOfLine();
      if (num_linear_terms != 0) {
        ReadLinearExpr(num_linear_terms,
            handler_.OnLinearVarExpr(var_index, num_linear_terms));
      }
      handler_.OnDefinedVar(var_index, ReadNumericExpr(), position);
      break;
    }
    case 'F': {
      // Imported function description.
      int index = ReadUInt(header_.num_funcs);
      int type = reader_.ReadUInt();
      if (type != func::NUMERIC && type != func::SYMBOLIC)
        reader_.ReportError("invalid function type");
      int num_args = reader_.template ReadInt<int>();
      fmt::StringRef name = reader_.ReadName();
      reader_.ReadTillEndOfLine();
      handler_.OnFunction(index, name, num_args, static_cast<func::Type>(type));
      break;
    }
    case 'G':
      // Linear part of an objective expression & gradient sparsity.
      ReadLinearExpr<ObjHandler>();
      break;
    case 'J':
      // Jacobian sparsity & linear terms in constraints.
      ReadLinearExpr<AlgebraicConHandler>();
      break;
    case 'S': {
      // Suffix values.
      int kind = reader_.ReadUInt();
      if (kind > (suf::MASK | suf::FLOAT))
        reader_.ReportError("invalid suffix kind");
      switch (kind & suf::MASK) {
      case suf::VAR:
        ReadSuffix<VarHandler>(kind);
        break;
      case suf::CON:
        ReadSuffix<ConHandler>(kind);
        break;
      case suf::OBJ:
        ReadSuffix<ObjHandler>(kind);
        break;
      case suf::PROBLEM:
        ReadSuffix<ProblemHandler>(kind);
        break;
      }
      break;
    }
    case 'b':
      // Bounds on variables.
      ReadBounds<VarHandler>();
      break;
    case 'r':
      // Bounds on algebraic constraint bodies ("ranges").
      ReadBounds<AlgebraicConHandler>();
      break;
    case 'K':
      // Jacobian sparsity & linear constraint term matrix column sizes
      // (must precede all J segments).
      ReadColumnSizes<false>();
      break;
    case 'k':
      // Jacobian sparsity & linear constraint term matrix cumulative column
      // sizes (must precede all J segments).
      ReadColumnSizes<true>();
      break;
    case 'x':
      // Primal initial guess.
      ReadInitialValues<VarHandler>();
      break;
    case 'd':
      // Dual initial guess.
      ReadInitialValues<AlgebraicConHandler>();
      break;
    case '\0':
      if (reader_.IsEOF())
        return;
      // Fall through.
    default:
      reader_.ReportError("invalid segment type");
    }
  }
}

// An .nl file.
class NLFile {
 private:
  fmt::File file_;
  std::size_t size_;
  std::size_t rounded_size_;  // Size rounded up to a multiple of page size.

 public:
  explicit NLFile(fmt::StringRef filename)
    : file_(filename, fmt::File::RDONLY) {
    fmt::LongLong file_size = file_.size();
    assert(file_size >= 0);
    fmt::ULongLong unsigned_file_size = file_size;
    // Check if file size fits in size_t.
    size_ = static_cast<std::size_t>(unsigned_file_size);
    if (size_ != unsigned_file_size)
      throw Error("file {} is too big", filename);
    // Round size up to a multiple of page_size. The remainded of the last
    // partial page is zero-filled both on POSIX and Windows so the resulting
    // memory buffer is zero terminated.
    std::size_t page_size = fmt::getpagesize();
    std::size_t remainder = size_ % page_size;
    rounded_size_ = remainder != 0 ? (size_ + page_size - remainder) : size_;
  }

  // Reads the file into an array.
  void Read(fmt::internal::Array<char, 1> &array);

  std::size_t size() { return size_; }
  std::size_t rounded_size() { return rounded_size_; }

  fmt::File &get() { return file_; }
};
}  // namespace internal

// Reads a string containing an optimization problem in .nl format.
// name: Name to be used when reporting errors.
template <typename Handler>
void ReadNLString(fmt::StringRef str, Handler &handler,
                  fmt::StringRef name = "(input)") {
  internal::TextReader reader(str, name);
  NLHeader header = NLHeader();
  reader.ReadHeader(header);
  handler.OnHeader(header);
  switch (header.format) {
  case NLHeader::TEXT:
    internal::NLReader<internal::TextReader, Handler>(
          reader, header, handler).Read();
    break;
  case NLHeader::BINARY: {
    arith::Kind arith_kind = arith::GetKind();
    if (arith_kind != header.arith_kind) {
      if (IsIEEE(arith_kind) && IsIEEE(header.arith_kind))
        ; // TODO: use binary reader & swap bytes
      else
        throw ReadError(name, 0, 0, "unsupported floating-point arithmetic");
    } else {
      internal::BinaryReader bin_reader(reader);
      internal::NLReader<internal::BinaryReader, Handler>(
            bin_reader, header, handler).Read();
    }
    break;
  }
  }
}

// Reads an .nl file.
template <typename Handler>
void ReadNLFile(fmt::StringRef filename, Handler &handler) {
  internal::NLFile file(filename);
  std::size_t size = file.size(), rounded_size = file.rounded_size();
  if (size == rounded_size) {
    // Don't use mmap, because the file size is a multiple of the page size
    // and therefore the mmap'ed buffer won't be zero terminated.
    fmt::internal::Array<char, 1> array;
    file.Read(array);
    return ReadNLString(fmt::StringRef(&array[0], size), handler, filename);
  }
  MemoryMappedFile mapped_file(file.get(), rounded_size);
  ReadNLString(fmt::StringRef(mapped_file.start(), size), handler, filename);
}
}  // namespace mp

#endif  // MP_NL_H_
