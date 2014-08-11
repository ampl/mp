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

#ifndef MP_NL_H_
#define MP_NL_H_

#include "mp/error.h"
#include "mp/os.h"
#include "mp/problem-base.h"
#include "mp/safeint.h"

#include <cctype>
#include <cstdlib>
#include <limits>
#include <string>

namespace mp {

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

enum {
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
// The .nl file format is described in the technical report
// "Writing .nl Files" (http://www.cs.sandia.gov/~dmgay/nlwrite.pdf).
struct NLHeader {
  // .nl file format.
  enum Format { TEXT = 0, BINARY = 1 };
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

  // Floating-point arithmetic kind used with binary format to check
  // if an .nl file is written using a compatible representation of
  // floating-point numbers. It is not used with text format and normally
  // set to arith::UNKNOWN there.
  arith::Kind arith_kind;

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

  // Number of common expressions that only appear in a single constraint
  // and don't appear in objectives.
  int num_common_exprs_in_single_cons;

  // Number of common expressions that only appear in a single objective
  // and don't appear in constraints.
  int num_common_exprs_in_single_objs;
};

// Writes NLHeader in the .nl file format.
fmt::Writer &operator<<(fmt::Writer &w, const NLHeader &h);

class TextReader {
 private:
  const char *ptr_, *end_;
  const char *line_start_;
  const char *token_;  // start of the current token
  std::string name_;
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
    if (result > std::numeric_limits<Int>::max())
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

  char ReadChar() {
    token_ = ptr_;
    return *ptr_++;
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
    DoReportError(ptr_, "expected newline");
  }

  bool ReadOptionalInt(int &value) { return DoReadOptionalInt(value); }

  bool ReadOptionalUInt(int &value) {
    SkipSpace();
    return ReadIntWithoutSign(value);
  }

  template <typename Int>
  Int ReadInt() {
    Int value = 0;
    if (!DoReadOptionalInt(value))
      ReportError("expected integer");
    return value;
  }

  int ReadUInt() {
    SkipSpace();
    int value = 0;
    if (!ReadIntWithoutSign(value))
      ReportError("expected unsigned integer");
    return value;
  }

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

  bool ReadOptionalDouble(double &value);

  // Reads a function or suffix name.
  fmt::StringRef ReadName();

  fmt::StringRef ReadString();

  // Reads an .nl file header. The header is always in text format, so this
  // function doesn't have a counterpart in BinaryReader.
  void ReadHeader(NLHeader &header);
};

// An .nl file reader.
// Handler is a class that receives notifications of the content of a file.
template <typename Reader, typename Handler>
class NLReader {
 private:
  Reader &reader_;
  NLHeader &header_;
  Handler &handler_;
  int total_num_vars_;  // Total number of variables including defined ones.

  // Minimum number of arguments for an iterated expression that has a
  // binary counterpart. Examples: sum (+), forall (&&), exists (||).
  enum {MIN_ITER_ARGS = 3};

  typedef typename Handler::NumericExpr NumericExpr;
  typedef typename Handler::LogicalExpr LogicalExpr;
  typedef typename Handler::Variable Variable;

  double ReadConstant(char code);
  double ReadConstant() { return ReadConstant(reader_.ReadChar()); }

  // Reads a nonnegative integer and checks that it is less than ub.
  int ReadUInt(int ub) {
    int index = reader_.ReadUInt();
    if (index >= ub)
      reader_.ReportError("integer {} out of bounds", index);
    return index;
  }

  // Reads a nonnegative integer and checks that it is in the range [lb, ub).
  int ReadUInt(int lb, int ub) {
    int index = reader_.ReadUInt();
    if (index < lb || index >= ub)
      reader_.ReportError("integer {} out of bounds", index);
    return index;
  }

  Variable DoReadVariable() {
    int var_index = ReadUInt(total_num_vars_);
    reader_.ReadTillEndOfLine();
    return handler_.MakeVariable(var_index);
  }

  Variable ReadVariable() {
    if (reader_.ReadChar() != 'v')
      reader_.ReportError("expected variable");
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
        r.reader_.ReportError("expected count expression opcode");
      return r.ReadCountExpr();
    }
  };

  template <typename ExprReader = NumericExprReader>
  class ReadArgs {
   private:
    typedef typename ExprReader::Expr Expr;
    // TODO: pass arguments directly to handler
    fmt::internal::Array<Expr, 10> args_;

   public:
    ReadArgs(NLReader &r, int min_args = MIN_ITER_ARGS);

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
    typename ExprReader::Expr lhs;
    typename ExprReader::Expr rhs;
    BinaryArgReader(NLReader &r)
      : lhs(ExprReader().Read(r)), rhs(ExprReader().Read(r)) {}
  };

  template <typename ExprReader>
  typename ExprReader::Expr ReadExpr() {
    int opcode = reader_.ReadUInt();
    if (opcode >= N_OPS)
      reader_.ReportError("invalid opcode {}", opcode);
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
      this->reader_.handler_.SetVarBounds(index, lb, ub);
    }
    void SetInitialValue(int index, double value) {
      this->reader_.handler_.SetInitialValue(index, value);
    }
  };

  struct ObjHandler : ItemHandler<OBJ> {
    explicit ObjHandler(NLReader &r) : ItemHandler<OBJ>(r) {}

    int num_items() const { return this->reader_.header_.num_objs; }

    typename Handler::LinearExprHandler
        GetLinearExprHandler(int index, int num_terms) {
      return this->reader_.handler_.GetLinearObjHandler(index, num_terms);
    }
  };

  struct AlgebraicConHandler : ItemHandler<CON> {
    explicit AlgebraicConHandler(NLReader &r) : ItemHandler<CON>(r) {}

    int num_items() const { return this->reader_.header_.num_algebraic_cons; }

    typename Handler::LinearExprHandler
        GetLinearExprHandler(int index, int num_terms) {
      return this->reader_.handler_.GetLinearConHandler(index, num_terms);
    }

    void SetBounds(int index, double lb, double ub) {
      this->reader_.handler_.SetConBounds(index, lb, ub);
    }
    void SetInitialValue(int index, double value) {
      this->reader_.handler_.SetInitialDualValue(index, value);
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

  void ReadLinearExpr(
      int num_terms, typename Handler::LinearExprHandler linear_expr);

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
template <typename ExprReader>
NLReader<Reader, Handler>::ReadArgs<ExprReader>::ReadArgs(
    NLReader &r, int min_args) {
  int num_args = r.reader_.ReadUInt();
  if (num_args < min_args)
    r.reader_.ReportError("too few arguments");
  r.reader_.ReadTillEndOfLine();
  args_.resize(num_args);
  ExprReader expr_reader;
  for (int i = 0; i < num_args; ++i)
    args_[i] = expr_reader.Read(r);
}

template <typename Reader, typename Handler>
typename Handler::NumericExpr
    NLReader<Reader, Handler>::ReadNumericExpr(char code) {
  switch (code) {
  case 'f': {
    int func_index = ReadUInt(header_.num_funcs);
    int num_args = reader_.ReadUInt();
    reader_.ReadTillEndOfLine();
    typedef typename Handler::Expr Expr;
    fmt::internal::Array<Expr, 10> args(num_args); // TODO
    for (int i = 0; i < num_args; ++i) {
      char c = reader_.ReadChar();
      args[i] = c == 'h' ?
            handler_.MakeStringLiteral(reader_.ReadString()) :
            args[i] = ReadNumericExpr(c);
    }
    return handler_.MakeCall(func_index, ArrayRef<Expr>(&args[0], args.size()));
  }
  case 'n': case 'l': case 's':
    return handler_.MakeNumericConstant(ReadConstant(code));
  case 'o':
    return ReadExpr<NumericExprReader>();
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
      reader_.ReportError("too few slopes in piecewise-linear term");
    reader_.ReadTillEndOfLine();
    using fmt::internal::Array; // TODO
    Array<double, 10> breakpoints(num_slopes - 1);
    Array<double, 10> slopes(num_slopes);
    for (int i = 0; i < num_slopes - 1; ++i) {
      slopes[i] = ReadConstant();
      breakpoints[i] = ReadConstant();
    }
    slopes[num_slopes - 1] = ReadConstant();
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
    reader_.ReportError("expected numeric expression opcode");
  }
  return NumericExpr();
}

template <typename Reader, typename Handler>
typename Handler::LogicalExpr NLReader<Reader, Handler>::ReadLogicalExpr() {
  switch (char c = reader_.ReadChar()) {
  case 'n': case 'l': case 's':
    return handler_.MakeLogicalConstant(ReadConstant(c) != 0);
  case 'o':
    return ReadExpr<LogicalExprReader>();
  }
  reader_.ReportError("expected logical expression");
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
      reader_.ReportError("expected count expression");
    return handler_.MakeLogicalCount(opcode, lhs, ReadExpr<CountExprReader>());
  }
  case expr::IMPLICATION: {
    LogicalExpr condition = ReadLogicalExpr();
    LogicalExpr true_expr = ReadLogicalExpr();
    LogicalExpr false_expr = ReadLogicalExpr();
    return handler_.MakeImplication(condition, true_expr, false_expr);
  }
  case expr::ITERATED_LOGICAL:
    return handler_.MakeIteratedLogical(
          opcode, ReadArgs<LogicalExprReader>(*this));
  case expr::ALLDIFF:
    return handler_.MakeAllDiff(ReadArgs<>(*this));
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
  int num_terms = ReadUInt(1, total_num_vars_); // TODO: inclusive
  reader_.ReadTillEndOfLine();
  ReadLinearExpr(num_terms, lh.GetLinearExprHandler(index, num_terms));
}

template <typename Reader, typename Handler>
void NLReader<Reader, Handler>::ReadLinearExpr(
    int num_terms, typename Handler::LinearExprHandler linear_expr) {
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
    RANGE,  // Both lower and upper bounds: l <= body <= u.
    UPPER,  // Only upper bound: body <= u.
    LOWER,  // Only lower bound: l <= body.
    FREE,   // No constraints on body (free variable or constraint).
    CONST,  // Equal to constant: body = c.
    COMPL   // Body complements variable v[i - 1].
  };
  reader_.ReadTillEndOfLine();
  double lb = 0, ub = 0;
  BoundHandler bh(*this);
  int num_bounds = bh.num_items();
  double infinity = std::numeric_limits<double>::infinity();
  for (int i = 0; i < num_bounds; ++i) {
    switch (reader_.ReadUInt()) {
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
    case CONST:
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
        handler_.SetComplement(i, var_index, flags & mask);
        reader_.ReadTillEndOfLine();
        continue;
      }
      // Fall through as COMPL bound type is invalid for variables.
    default:
      reader_.ReportError("invalid bound type");
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
  typename Handler::ColumnSizeHandler
      size_handler = handler_.GetColumnSizeHandler();
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
      suffix_handler = handler_.AddSuffix(kind, num_values, name);
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
      handler_.SetCon(index, ReadNumericExpr());
      break;
    }
    case 'L': {
      // Logical constraint expression.
      int index = ReadUInt(header_.num_logical_cons);
      reader_.ReadTillEndOfLine();
      handler_.SetLogicalCon(index, ReadLogicalExpr());
      break;
    }
    case 'O': {
      // Objective type and nonlinear part of an objective expression.
      int index = ReadUInt(header_.num_objs);
      int obj_type = reader_.ReadUInt();
      reader_.ReadTillEndOfLine();
      handler_.SetObj(index, obj_type != 0 ? obj::MAX : obj::MIN,
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
            handler_.GetLinearVarHandler(var_index, num_linear_terms));
      }
      handler_.SetVar(var_index, ReadNumericExpr(), position);
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
      handler_.SetFunction(index, name, num_args,
                           static_cast<func::Type>(type));
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
      case suf::PROB:
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
      // TODO: check for end of input
      return;
    default:
      // TODO: handle unprintable chars
      reader_.ReportError("invalid segment type '{}'", c);
    }
  }
}

// Reads a string containing an optimization problem in .nl format.
// name: Name to be used when reporting errors.
template <typename Handler>
void ReadNLString(fmt::StringRef str, Handler &handler,
                  fmt::StringRef name = "(input)") {
  TextReader reader(str, name);
  NLHeader header = NLHeader();
  reader.ReadHeader(header);
  handler.BeginBuild(name.c_str(), header, 0);
  switch (header.format) {
  case NLHeader::TEXT:
    NLReader<TextReader, Handler>(reader, header, handler).Read();
    break;
  case NLHeader::BINARY: {
    arith::Kind arith_kind = arith::GetKind();
    if (arith_kind != header.arith_kind) {
      if (IsIEEE(arith_kind) && IsIEEE(header.arith_kind))
        ; // TODO: use binary reader & swap bytes
      else
        throw ReadError(name, 0, 0, "unsupported floating-point arithmetic");
    } else {
      // TODO: use binary reader
    }
    break;
  }
  }
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
}  // namespace mp

#endif  // MP_NL_H_
