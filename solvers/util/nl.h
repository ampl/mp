/*
 .nl file support.

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

#ifndef SOLVERS_UTIL_NL_H_
#define SOLVERS_UTIL_NL_H_

#include "solvers/arith.h"
#include "solvers/util/error.h"
#include "solvers/util/os.h"
#include "solvers/util/problem-base.h"
#include "solvers/util/safeint.h"
#include "solvers/opcode.hd"

#include <cctype>
#include <limits>

namespace ampl {

class ParseError : public Error {
 private:
  std::string filename_;
  int line_;
  int column_;

 public:
  ParseError(fmt::StringRef filename,
      int line, int column, fmt::StringRef message)
  : Error(message), filename_(filename), line_(line), column_(column) {}
  ~ParseError() throw() {}

  const std::string &filename() const { return filename_; }
  int line() const { return line_; }
  int column() const { return column_; }
};

class TextReader;

enum {
  MAX_NL_OPTIONS = 9,
  VBTOL_OPTION   = 1,
  READ_VBTOL     = 3
};

// NL file header.
struct NLHeader {
  // NL file format.
  enum Format { TEXT = 0, BINARY = 1, BINARY_SWAPPED = 2 };
  Format format;

  int num_options;
  int options[MAX_NL_OPTIONS];

  // Extra info for writing solution.
  double ampl_vbtol;

  // Total number of variables.
  int num_vars;

  // Number of algebraic constraints including ranges and equality constraints.
  // It doesn't include logical constraints.
  int num_algebraic_cons;

  // Total number of objectives.
  int num_objs;

  // Number of ranges (constraints with -Infinity < LHS < RHS < Infinity).
  int num_ranges;

  // Number of equality constraints or -1 if unknown (AMPL prior to 19970627).
  int num_eqns;

  // Number of logical constraints.
  int num_logical_cons;

  // Nonlinear and complementarity information
  // -----------------------------------------

  // Total number of nonlinear constraints.
  int num_nl_cons;

  // Total number of nonlinear objectives.
  int num_nl_objs;

  // Total number of complementarity conditions.
  int num_compl_conds;

  // Number of nonlinear complementarity conditions.
  int num_nl_compl_conds;

  // Number of complementarities involving double inequalities
  // (for ASL_cc_simplify).
  int num_compl_dbl_ineqs;

  // Number of complemented variables with a nonzero lower bound
  // (for ASL_cc_simplify).
  int num_compl_vars_with_nz_lb;

  // Information about network constraints
  // -------------------------------------

  // Number of nonlinear network constraints.
  int num_nl_net_cons;

  // Number of linear network constraints.
  int num_linear_net_cons;

  // Information about nonlinear variables
  // -------------------------------------

  // Number of nonlinear variables in constraints including nonlinear
  // variables in both constraints and objectives.
  int num_nl_vars_in_cons;

  // Number of nonlinear variables in objectives including nonlinear
  // variables in both constraints and objectives.
  int num_nl_vars_in_objs;

  // Number of nonlinear variables in both constraints and objectives.
  int num_nl_vars_in_both;

  // Miscellaneous
  // -------------

  // Number of linear network variables (arcs).
  int num_linear_net_vars;

  // Number of functions.
  int num_funcs;

  // Flags: 1 = want output suffixes.
  int flags;

  // Information about discrete variables
  // ------------------------------------

  // Number of linear binary variables.
  int num_linear_binary_vars;

  // Number of linear non-binary integer variables.
  int num_linear_integer_vars;

  // Number of integer nonlinear variables in both constraints and objectives.
  int num_nl_integer_vars_in_both;

  // Number of integer nonlinear variables just in constraints.
  int num_nl_integer_vars_in_cons;

  // Number of integer nonlinear variables just in objectives.
  int num_nl_integer_vars_in_objs;

  // Information about nonzeros
  // --------------------------

  // Number of nonzeros in constraints' Jacobian.
  int num_con_nonzeros;


  // Number of nonzeros in all objective gradients.
  int num_obj_nonzeros;

  // Information about names
  // -----------------------

  // Length of longest constraint name (if stub.row exists).
  int max_con_name_len;

  // Length of longest variable name (if stub.col exists).
  int max_var_name_len;

  // Information about common expressions
  // ------------------------------------

  int num_common_exprs_in_both;
  int num_common_exprs_in_cons;
  int num_common_exprs_in_objs;
  int num_common_exprs_in_cons1;
  int num_common_exprs_in_objs1;
};

fmt::Writer &operator<<(fmt::Writer &w, const NLHeader &h);

// A reference to an immutable array.
template <typename T>
class ArrayRef {
 private:
  const T *data_;
  std::size_t size_;

 public:
  ArrayRef(const T *data, std::size_t size) : data_(data), size_(size) {}

  template <typename U>
  ArrayRef(ArrayRef<U> other) : data_(other.data()), size_(other.size()) {}

  template <std::size_t SIZE>
  ArrayRef(const T (&data)[SIZE]) : data_(data), size_(SIZE) {}

  const T *data() const { return data_; }
  std::size_t size() const { return size_; }

  const T &operator[](std::size_t i) const { return data_[i]; }
};

template <typename T>
ArrayRef<T> MakeArrayRef(const T *data, std::size_t size) {
  return ArrayRef<T>(data, size);
}

class TextReader {
 private:
  const char *ptr_;
  const char *line_start_;
  const char *token_;  // start of the current token
  std::string name_;
  int line_;

  template <typename Int>
  typename safeint::MakeUnsigned<Int>::Type DoReadUInt() {
    char c = *ptr_;
    typedef typename safeint::MakeUnsigned<Int>::Type UInt;
    UInt value = 0;
    do {
      UInt new_value = value * 10 + (c - '0');
      if (new_value < value)
        ReportParseError("number is too big");
      value = new_value;
      c = *++ptr_;
    } while (c >= '0' && c <= '9');
    if (value > std::numeric_limits<Int>::max())
      ReportParseError("number is too big");
    return value;
  }

 public:
  TextReader(fmt::StringRef name, const char *ptr)
  : ptr_(ptr), line_start_(ptr), token_(ptr), name_(name), line_(1) {}

  void ReportParseError(fmt::StringRef format_str, const fmt::ArgList &args) {
    int column = static_cast<int>(token_ - line_start_ + 1);
    fmt::Writer w;
    w.write(format_str, args);
    throw ampl::ParseError(name_, line_, column,
        fmt::format("{}:{}:{}: {}", name_, line_, column,
            fmt::StringRef(w.c_str(), w.size())));
  }
  FMT_VARIADIC(void, ReportParseError, fmt::StringRef)

  char ReadChar() {
    token_ = ptr_;
    return *ptr_++;
  }

  void SkipSpace() {
    while (std::isspace(*ptr_) && *ptr_ != '\n')
      ++ptr_;
    token_ = ptr_;
  }

  void ReadTillEndOfLine() {
    while (char c = *ptr_) {
      ++ptr_;
      if (c == '\n') {
        line_start_ = ptr_;
        ++line_;
        return;
      }
    }
  }

  template <typename Int>
  Int ReadUInt() {
    SkipSpace();
    char c = *ptr_;
    if (c < '0' || c > '9')
      ReportParseError("expected nonnegative integer");
    return DoReadUInt<Int>();
  }

  int ReadUInt() { return ReadUInt<int>(); }

  bool ReadOptionalUInt(int &value) {
    SkipSpace();
    char c = *ptr_;
    bool has_value = c >= '0' && c <= '9';
    if (has_value)
      value = DoReadUInt<int>();
    return has_value;
  }

  template <typename Int>
  Int ReadInt() {
    char sign = *ptr_;
    if (sign == '+' || sign == '-')
      ++ptr_;
    typedef typename safeint::MakeUnsigned<Int>::Type UInt;
    UInt result = ReadUInt<UInt>();
    UInt max = std::numeric_limits<Int>::max();
    if (result > max && !(sign == '-' && result == max + 1))
      ReportParseError("number is too big");
    return sign != '-' ? result : 0 - result;
  }

  double ReadDouble() {
    SkipSpace();
    char *end = 0;
    double value = 0;
    if (*ptr_ != '\n')
      value = strtod(ptr_, &end);
    if (!end || ptr_ == end)
      ReportParseError("expected double");
    ptr_ = end;
    return value;
  }

  bool ReadOptionalDouble(double &value) {
    SkipSpace();
    if (*ptr_ == '\n')
      return false;
    char *end = 0;
    value = strtod(ptr_, &end);
    bool has_value = ptr_ != end;
    ptr_ = end;
    return has_value;
  }
};

// An .nl file reader.
// Handler is a class that receives notifications of the content of a file.
template <typename Reader, typename Handler>
class NLReader {
 private:
  Reader &reader_;
  NLHeader &header_;
  Handler &handler_;

  // Minimum number of arguments for an iterated expression that has a
  // binary counterpart. Examples: sum (+), forall (&&), exists (||).
  enum {MIN_ITER_ARGS = 3};

  typedef typename Handler::NumericExpr NumericExpr;
  typedef typename Handler::LogicalExpr LogicalExpr;
  typedef typename Handler::Variable Variable;

  double ReadNumericConstant(char code);
  double ReadNumericConstant() {
    return ReadNumericConstant(reader_.ReadChar());
  }

  Variable DoReadVariable() {
    // TODO: variable index can be greater than num_vars
    int var_index = reader_.ReadUInt();
    if (var_index >= header_.num_vars)
      reader_.ReportParseError("variable index {} out of bounds", var_index);
    reader_.ReadTillEndOfLine();
    return handler_.MakeVariable(var_index);
  }

  Variable ReadVariable() {
    if (reader_.ReadChar() != 'v')
      reader_.ReportParseError("expected variable");
    return DoReadVariable();
  }

  typename Handler::CountExpr ReadCountExpr() {
    return handler_.MakeCount(ReadArgs<LogicalExprReader>(*this, 1));
  }

  // Helper structs to provide a uniform interface to Read{Numeric,Logical}Expr
  // since it is not possible to overload on expression type as NumericExpr
  // and LogicalExpr can be the same type.
  struct NumericExprReader {
    typedef NumericExpr Expr;
    Expr Read(NLReader &r) const { return r.ReadNumericExpr(); }
    Expr Read(NLReader &r, int opcode) const {
      return r.ReadNumericExpr(opcode);
    }
  };
  struct LogicalExprReader {
    typedef LogicalExpr Expr;
    Expr Read(NLReader &r) const { return r.ReadLogicalExpr(); }
    Expr Read(NLReader &r, int opcode) const {
      return r.ReadLogicalExpr(opcode);
    }
  };
  struct CountExprReader {
    typedef typename Handler::CountExpr Expr;
    Expr Read(NLReader &r, int opcode) const {
      if (opcode != OPCOUNT)
        r.reader_.ReportParseError("expected count expression");
      return r.ReadCountExpr();
    }
  };

  template <typename ExprReader = NumericExprReader>
  class ReadArgs {
   private:
    typedef typename ExprReader::Expr Expr;
    fmt::internal::Array<Expr, 10> args_;

   public:
    ReadArgs(NLReader &r, int min_args = MIN_ITER_ARGS) {
      int num_args = r.reader_.ReadUInt();
      if (num_args < min_args)
        r.reader_.ReportParseError("too few arguments");
      args_.resize(num_args);
      ExprReader expr_reader;
      for (int i = 0; i < num_args; ++i)
        args_[i] = expr_reader.Read(r);
    }

    operator ArrayRef<Expr>() const {
      return ArrayRef<Expr>(&args_[0], args_.size());
    }
  };

  // A helper struct used to make sure that the arguments to a binary
  // expression are read in the correct order and avoid errors of the form:
  //   MakeBinary(opcode, ReadNumericExpr(), ReadNumericExpr())
  // The above code is incorrect as the order of evaluation of arguments is
  // unspecified.
  template <typename ExprReader = NumericExprReader>
  struct BinaryArgReader {
    LogicalExpr lhs;
    LogicalExpr rhs;
    BinaryArgReader(NLReader &r)
      : lhs(ExprReader().Read(r)), rhs(ExprReader().Read(r)) {}
  };

  template <typename ExprReader>
  typename ExprReader::Expr ReadExpr() {
    int opcode = reader_.ReadUInt();
    if (opcode > N_OPS)
      reader_.ReportParseError("invalid opcode {}", opcode);
    reader_.ReadTillEndOfLine();
    return ExprReader().Read(*this, opcode);
  }

  // Reads a numeric expression.
  NumericExpr ReadNumericExpr() { return ReadNumericExpr(reader_.ReadChar()); }
  NumericExpr ReadNumericExpr(char code);
  NumericExpr ReadNumericExpr(int opcode);

  // Reads a logical expression.
  LogicalExpr ReadLogicalExpr();
  LogicalExpr ReadLogicalExpr(int opcode);

  // Reads a linear expression.
  void ReadLinearExpr(int num_terms);

  // Reads bounds.
  void ReadBounds();

  // Read the column offsets, the cumulative sums of the numbers of
  // nonzeros in the first num_var − 1 columns of the Jacobian matrix.
  void ReadColumnOffsets();

 public:
  NLReader(Reader &reader, NLHeader &header, Handler &handler)
    : reader_(reader), header_(header), handler_(handler) {}

  void Read();
};

template <typename Reader, typename Handler>
double NLReader<Reader, Handler>::ReadNumericConstant(char code) {
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
    reader_.ReportParseError("expected numeric constant");
  }
  reader_.ReadTillEndOfLine();
  return value;
}

template <typename Reader, typename Handler>
typename Handler::NumericExpr
    NLReader<Reader, Handler>::ReadNumericExpr(char code) {
  NumericExpr expr;
  switch (code) {
  case 'f': {
    int func_index = reader_.ReadUInt();
    if (func_index >= header_.num_funcs)
      reader_.ReportParseError("function index {} out of bounds", func_index);
    int num_args = reader_.ReadUInt();
    reader_.ReadTillEndOfLine();
    fmt::internal::Array<typename Handler::Expr, 10> args(num_args);
    for (int i = 0; i < num_args; ++i) {
      char c = reader_.ReadChar();
      if (c == 'h') {
        // TODO: read string
      } else {
        args[i] = ReadNumericExpr(c);
      }
    }
    // TODO
    //expr = handler_.MakeCall(func_index, args);
    break;
  }
  case 'n': case 'l': case 's':
    return handler_.MakeNumericConstant(ReadNumericConstant(code));
  case 'o':
    return ReadExpr<NumericExprReader>();
  case 'v':
    return DoReadVariable();
  default:
    reader_.ReportParseError("expected numeric expression");
  }
  reader_.ReadTillEndOfLine();
  return expr;
}

template <typename Reader, typename Handler>
typename Handler::NumericExpr
    NLReader<Reader, Handler>::ReadNumericExpr(int opcode) {
  using fmt::internal::Array;
  switch (expr::kind(opcode)) {
  case expr::UNARY:
    return handler_.MakeUnary(opcode, ReadNumericExpr());
  case expr::BINARY: {
    BinaryArgReader<> args(*this);
    return handler_.MakeBinary(opcode, args.lhs, args.rhs);
  }
  case expr::IF: {
    LogicalExpr condition = ReadLogicalExpr();
    NumericExpr true_expr = ReadNumericExpr();
    NumericExpr false_expr = ReadNumericExpr();
    return handler_.MakeIf(condition, true_expr, false_expr);
  }
  case expr::PLTERM: {
    int num_slopes = reader_.ReadUInt();
    if (num_slopes <= 1)
      reader_.ReportParseError("too few slopes in piecewise-linear term");
    reader_.ReadTillEndOfLine();
    Array<double, 10> breakpoints(num_slopes - 1);
    Array<double, 10> slopes(num_slopes);
    for (int i = 0; i < num_slopes - 1; ++i) {
      slopes[i] = ReadNumericConstant();
      breakpoints[i] = ReadNumericConstant();
    }
    slopes[num_slopes - 1] = ReadNumericConstant();
    return handler_.MakePiecewiseLinear(
          num_slopes - 1, &breakpoints[0], &slopes[0], ReadVariable());
  }
  case expr::VARARG:
    return handler_.MakeVarArg(opcode, ReadArgs<>(*this, 1));
  case expr::SUM:
    return handler_.MakeSum(ReadArgs<>(*this));
  case expr::COUNT:
    return ReadCountExpr();
  case expr::NUMBEROF:
    return handler_.MakeNumberOf(ReadArgs<>(*this, 1));
  default:
    reader_.ReportParseError("expected numeric expression");
  }
  return NumericExpr();
}

template <typename Reader, typename Handler>
typename Handler::LogicalExpr NLReader<Reader, Handler>::ReadLogicalExpr() {
  switch (char c = reader_.ReadChar()) {
  case 'n': case 'l': case 's':
    return handler_.MakeLogicalConstant(ReadNumericConstant(c) != 0);
  case 'o':
    return ReadExpr<LogicalExprReader>();
  }
  reader_.ReportParseError("expected logical expression");
  return LogicalExpr();
}

template <typename Reader, typename Handler>
typename Handler::LogicalExpr
    NLReader<Reader, Handler>::ReadLogicalExpr(int opcode) {
  switch (expr::kind(opcode)) {
  case expr::NOT:
    return handler_.MakeNot(ReadLogicalExpr());
  case expr::BINARY_LOGICAL: {
    BinaryArgReader<LogicalExprReader> args(*this);
    return handler_.MakeBinaryLogical(opcode, args.lhs, args.rhs);
  }
  case expr::RELATIONAL: {
    BinaryArgReader<> args(*this);
    return handler_.MakeRelational(opcode, args.lhs, args.rhs);
  }
  case expr::LOGICAL_COUNT: {
    NumericExpr lhs = ReadNumericExpr();
    char c = reader_.ReadChar();
    if (c != 'o')
      reader_.ReportParseError("expected count expression");
    return handler_.MakeLogicalCount(opcode, lhs, ReadExpr<CountExprReader>());
  }
  case expr::IMPLICATION: {
    LogicalExpr condition = ReadLogicalExpr();
    LogicalExpr true_expr = ReadNumericExpr();
    LogicalExpr false_expr = ReadNumericExpr();
    return handler_.MakeImplication(condition, true_expr, false_expr);
  }
  case expr::ITERATED_LOGICAL:
    return handler_.MakeIteratedLogical(
          opcode, ReadArgs<LogicalExprReader>(*this));
  case expr::ALLDIFF:
    return handler_.MakeAllDiff(ReadArgs<>(*this));
  default:
    reader_.ReportParseError("expected logical expression");
  }
  return LogicalExpr();
}

template <typename Reader, typename Handler>
void NLReader<Reader, Handler>::ReadLinearExpr(int num_terms) {
  for (int i = 0; i < num_terms; ++i) {
    reader_.ReadUInt();
    reader_.ReadUInt(); // TODO: read double
    // TODO
    reader_.ReadTillEndOfLine();
  }
}

template <typename Reader, typename Handler>
void NLReader<Reader, Handler>::ReadBounds() {
  for (int i = 0; i < header_.num_vars; ++i) {
    reader_.ReadUInt();
    reader_.ReadUInt(); // TODO: read double
    // TODO
    reader_.ReadTillEndOfLine();
  }
}

template <typename Reader, typename Handler>
void NLReader<Reader, Handler>::ReadColumnOffsets() {
  int count = reader_.ReadUInt(); // TODO
  reader_.ReadTillEndOfLine();
  for (int i = 0; i < count; ++i) {
    reader_.ReadUInt(); // TODO
    reader_.ReadTillEndOfLine();
  }
}

template <typename Reader, typename Handler>
void NLReader<Reader, Handler>::Read() {
  for (;;) {
    char c = reader_.ReadChar();
    switch (c) {
    case '\0':
      // TODO: check for end of input
      return;
    case 'C': {
      int index = reader_.ReadUInt();
      if (index >= header_.num_algebraic_cons)
        reader_.ReportParseError("constraint index {} out of bounds", index);
      reader_.ReadTillEndOfLine();
      handler_.SetCon(index, ReadNumericExpr());
      break;
    }
    case 'F': {
      int index = reader_.ReadUInt();
      if (index >= header_.num_funcs)
        reader_.ReportParseError("function index {} out of bounds", index);
      int type = reader_.ReadUInt();
      if (type != func::NUMERIC && type != func::SYMBOLIC) {
        if (type < 0 || type > 6)
          reader_.ReportParseError("invalid function type");
        // Ignore function of unsupported type.
        break;
      }
      int num_args = reader_.ReadUInt();
      const char *name = 0; // TODO: read name
      reader_.ReadTillEndOfLine();
      handler_.SetFunction(index, name, num_args,
                           static_cast<func::Type>(type));
      break;
    }
    case 'L': {
      int index = reader_.ReadUInt();
      if (index >= header_.num_logical_cons) {
        reader_.ReportParseError(
              "logical constraint index {} out of bounds", index);
      }
      reader_.ReadTillEndOfLine();
      ReadLogicalExpr();
      // TODO: send to handler
      break;
    }
    case 'V':
      // TODO
      break;
    case 'G': {
      int index = reader_.ReadUInt();
      if (index >= header_.num_objs)
        reader_.ReportParseError("objective index {} out of bounds", index);
      int num_terms = reader_.ReadUInt(); // TODO: check
      reader_.ReadTillEndOfLine();
      ReadLinearExpr(num_terms);
      // TODO: read gradient!
      break;
    }
    case 'J':
      // TODO: read Jacobian matrix
      break;
    case 'O': {
      int index = reader_.ReadUInt();
      if (index >= header_.num_objs)
        reader_.ReportParseError("objective index {} out of bounds", index);
      int obj_type = reader_.ReadUInt();
      reader_.ReadTillEndOfLine();
      handler_.SetObj(index, obj_type != 0 ? obj::MAX : obj::MIN,
                      ReadNumericExpr());
      break;
    }
    case 'S':
      // TODO: read suffix
      break;
    case 'r':
      // TODO: read RHS
      break;
    case 'b':
      reader_.ReadTillEndOfLine();
      ReadBounds();
      break;
    case 'k':
      ReadColumnOffsets();
      break;
    case 'x':
      // TODO
      break;
    case 'd':
      // TODO
      break;
    default:
      reader_.ReportParseError("invalid segment type '{}'", c);
    }
  }
}

// Reads a string containing a problem in .nl format.
// name: Name to be used when reporting errors.
// header_only: true to read the header only, false to read the whole file
template <typename Handler>
void ReadNLString(fmt::StringRef str, Handler &handler,
                  fmt::StringRef name = "(input)", bool header_only = false) {
  TextReader reader(name, str.c_str());

  // Read the format (text or binary).
  NLHeader header = NLHeader();
  switch (reader.ReadChar()) {
  case 'g':
    break;
  case 'b':
    header.format = NLHeader::BINARY;
    break;
  default:
    reader.ReportParseError("expected format specifier");
    break;
  }

  // Read options.
  reader.ReadOptionalUInt(header.num_options);
  if (header.num_options > MAX_NL_OPTIONS)
    reader.ReportParseError("too many options");
  for (int i = 0; i < header.num_options; ++i) {
    // TODO: can option values be negative?
    if (!reader.ReadOptionalUInt(header.options[i]))
      break;
  }
  if (header.options[VBTOL_OPTION] == READ_VBTOL)
    reader.ReadOptionalDouble(header.ampl_vbtol);
  reader.ReadTillEndOfLine();

  // Read problem dimensions.
  header.num_vars = reader.ReadUInt();
  header.num_algebraic_cons = reader.ReadUInt();
  header.num_objs = reader.ReadUInt();
  header.num_eqns = -1;
  if (reader.ReadOptionalUInt(header.num_ranges) &&
      reader.ReadOptionalUInt(header.num_eqns)) {
      reader.ReadOptionalUInt(header.num_logical_cons);
  }
  reader.ReadTillEndOfLine();

  // Read the nonlinear and complementarity information.
  header.num_nl_cons = reader.ReadUInt();
  header.num_nl_objs = reader.ReadUInt();
  bool all_compl =
      reader.ReadOptionalUInt(header.num_compl_conds) &&
      reader.ReadOptionalUInt(header.num_nl_compl_conds) &&
      reader.ReadOptionalUInt(header.num_compl_dbl_ineqs) &&
      reader.ReadOptionalUInt(header.num_compl_vars_with_nz_lb);
  header.num_compl_conds += header.num_nl_compl_conds;
  if (header.num_compl_conds > 0 && !all_compl)
    header.num_compl_dbl_ineqs = -1;
  reader.ReadTillEndOfLine();

  // Read the information about network constraints.
  header.num_nl_net_cons = reader.ReadUInt();
  header.num_linear_net_cons = reader.ReadUInt();
  reader.ReadTillEndOfLine();

  // Read the information about nonlinear variables.
  header.num_nl_vars_in_cons = reader.ReadUInt();
  header.num_nl_vars_in_objs = reader.ReadUInt();
  header.num_nl_vars_in_both = -1;
  reader.ReadOptionalUInt(header.num_nl_vars_in_both);
  reader.ReadTillEndOfLine();

  header.num_linear_net_vars = reader.ReadUInt();
  header.num_funcs = reader.ReadUInt();
  int arith = 0;
  if (reader.ReadOptionalUInt(arith)) {
    if (arith != Arith_Kind_ASL && arith != 0) {
      bool swap_bytes = false;
#if defined(IEEE_MC68k) || defined(IEEE_8087)
      swap_bytes = arith > 0 && arith + Arith_Kind_ASL == 3;
#endif
      if (!swap_bytes)
        reader.ReportParseError("unrecognized binary format");
      header.format = NLHeader::BINARY_SWAPPED;
      // TODO: swap bytes
    }
    reader.ReadOptionalUInt(header.flags);
  }
  reader.ReadTillEndOfLine();

  // Read the information about discrete variables.
  header.num_linear_binary_vars = reader.ReadUInt();
  header.num_linear_integer_vars = reader.ReadUInt();
  if (header.num_nl_vars_in_both >= 0) {  // ampl versions >= 19930630
    header.num_nl_integer_vars_in_both = reader.ReadUInt();
    header.num_nl_integer_vars_in_cons = reader.ReadUInt();
    header.num_nl_integer_vars_in_objs = reader.ReadUInt();
  }
  reader.ReadTillEndOfLine();

  // Read the information about nonzeros.
  header.num_con_nonzeros = reader.ReadUInt();
  header.num_obj_nonzeros = reader.ReadUInt();
  reader.ReadTillEndOfLine();

  // Read the information about names.
  header.max_con_name_len = reader.ReadUInt();
  header.max_var_name_len = reader.ReadUInt();
  reader.ReadTillEndOfLine();

  // Read the information about common expressions.
  header.num_common_exprs_in_both = reader.ReadUInt();
  header.num_common_exprs_in_cons = reader.ReadUInt();
  header.num_common_exprs_in_objs = reader.ReadUInt();
  header.num_common_exprs_in_cons1 = reader.ReadUInt();
  header.num_common_exprs_in_objs1 = reader.ReadUInt();
  reader.ReadTillEndOfLine();

  handler.BeginBuild(name.c_str(), header, 0);
  if (header_only)
    return;

  if (header.format != NLHeader::TEXT) {
    // TODO: use binary reader
  }
  NLReader<TextReader, Handler>(reader, header, handler).Read();
}

// Reads an .nl file.
template <typename Handler>
void ReadNLFile(fmt::StringRef filename, Handler &h) {
  MemoryMappedFile file(filename);
  // TODO: use a buffer instead of mmap if mmap is not available or the
  //       file length is a multiple of the page size
  std::size_t size = static_cast<std::size_t>(file.size());
  // Check if file size fits in size_t.
  if (size != file.size())
    throw Error("file {} is too big", filename);
  ReadNLString(fmt::StringRef(file.start(), size), h, filename);
}
}  // namespace ampl

#endif  // SOLVERS_UTIL_NL_H_
