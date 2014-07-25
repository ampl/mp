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

#include <cctype>
#include <climits>

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

class TextReader {
 private:
  const char *ptr_;
  const char *line_start_;
  const char *token_;  // start of the current token
  std::string name_;
  int line_;

  int DoReadUInt() {
    char c = *ptr_;
    unsigned value = 0;
    do {
      unsigned new_value = value * 10 + (c - '0');
      if (new_value < value)
        ReportParseError("number is too big");
      value = new_value;
      c = *++ptr_;
    } while (c >= '0' && c <= '9');
    if (value > INT_MAX)
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

  int ReadUInt() {
    SkipSpace();
    char c = *ptr_;
    if (c < '0' || c > '9')
      ReportParseError("expected nonnegative integer");
    return DoReadUInt();
  }

  bool ReadOptionalUInt(int &value) {
    SkipSpace();
    char c = *ptr_;
    bool has_value = c >= '0' && c <= '9';
    if (has_value)
      value = DoReadUInt();
    return has_value;
  }

  int ReadLong() {
    char sign = *ptr_;
    if (sign == '+' || sign == '-')
      ++ptr_;
    int result = ReadUInt();
    return sign == '-' ? -result : result;
  }

  int ReadShort() { return ReadLong(); }

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

// An nl reader.
// Handler is a class that receives notifications about
// read constructs.
template <typename Handler>
class NLReader {
 private:
  NLHeader header_;
  Handler &handler_;

  // Reads an expression.
  typename Handler::NumericExpr ReadExpr(TextReader &reader);

  // Reads a linear expression.
  void ReadLinearExpr(TextReader &reader, int num_terms);

  // Reads bounds.
  void ReadBounds(TextReader &reader);

  // Read the column offsets, the cumulative sums of the numbers of
  // nonzeros in the first num_var − 1 columns of the Jacobian matrix.
  void ReadColumnOffsets(TextReader &reader);

 public:
  explicit NLReader(Handler &h) : handler_(h) {}

  // Reads a problem from a file.
  void ReadFile(fmt::StringRef filename);

  // Reads a problem from a string.
  // name: Name to be used when reporting errors.
  // header_only: true to read the header only, false to read the whole file
  void ReadString(fmt::StringRef str,
      fmt::StringRef name = "(input)", bool header_only = false);
};

template <typename Hanlder>
typename Hanlder::NumericExpr NLReader<Hanlder>::ReadExpr(TextReader &reader) {
  typename Hanlder::NumericExpr expr;
  switch (reader.ReadChar()) {
  case 'f':
    // TODO: implement function
    break;
  case 'h':
    // TODO: implement string
    break;
  case 's':
    expr = handler_.MakeNumericConstant(reader.ReadShort());
    break;
  case 'l':
    expr = handler_.MakeNumericConstant(reader.ReadLong());
    break;
  case 'n':
    expr = handler_.MakeNumericConstant(reader.ReadDouble());
    break;
  case 'o':
    // TODO: implement expression
    break;
  case 'v': {
    // TODO: variable index can be greater than num_vars
    int var_index = reader.ReadUInt();
    if (var_index >= header_.num_vars)
      reader.ReportParseError("variable index {} is out of bounds", var_index);
    expr = handler_.MakeVariable(var_index);
    break;
  }
  default:
    reader.ReportParseError("expected expression");
  }
  reader.ReadTillEndOfLine();
  return expr;
}

template <typename Hanlder>
void NLReader<Hanlder>::ReadLinearExpr(TextReader &reader, int num_terms) {
  for (int i = 0; i < num_terms; ++i) {
    reader.ReadUInt();
    reader.ReadUInt(); // TODO: read double
    // TODO
    reader.ReadTillEndOfLine();
  }
}

template <typename Hanlder>
void NLReader<Hanlder>::ReadBounds(TextReader &reader) {
  for (int i = 0; i < header_.num_vars; ++i) {
    reader.ReadUInt();
    reader.ReadUInt(); // TODO: read double
    // TODO
    reader.ReadTillEndOfLine();
  }
}

template <typename Hanlder>
void NLReader<Hanlder>::ReadColumnOffsets(TextReader &reader) {
  int count = reader.ReadUInt(); // TODO
  reader.ReadTillEndOfLine();
  for (int i = 0; i < count; ++i) {
    reader.ReadUInt(); // TODO
    reader.ReadTillEndOfLine();
  }
}

template <typename Hanlder>
void NLReader<Hanlder>::ReadFile(fmt::StringRef filename) {
  MemoryMappedFile file(filename);
  // TODO: use a buffer instead of mmap if mmap is not available or the
  //       file length is a multiple of the page size
  std::size_t size = static_cast<std::size_t>(file.size());
  // Check if file size fits in size_t.
  if (size != file.size())
    throw Error("file {} is too big", filename);
  ReadString(fmt::StringRef(file.start(), size), filename);
}

template <typename Hanlder>
void NLReader<Hanlder>::ReadString(
    fmt::StringRef str, fmt::StringRef name, bool header_only) {
  TextReader reader(name, str.c_str());

  // Read the format (text or binary).
  header_ = NLHeader();
  switch (reader.ReadChar()) {
  case 'g':
    break;
  case 'b':
    header_.format = NLHeader::BINARY;
    break;
  default:
    reader.ReportParseError("expected format specifier");
    break;
  }

  // Read options.
  reader.ReadOptionalUInt(header_.num_options);
  if (header_.num_options > MAX_NL_OPTIONS)
    reader.ReportParseError("too many options");
  for (int i = 0; i < header_.num_options; ++i) {
    // TODO: can option values be negative?
    if (!reader.ReadOptionalUInt(header_.options[i]))
      break;
  }
  if (header_.options[VBTOL_OPTION] == READ_VBTOL)
    reader.ReadOptionalDouble(header_.ampl_vbtol);
  reader.ReadTillEndOfLine();

  // Read problem dimensions.
  header_.num_vars = reader.ReadUInt();
  header_.num_algebraic_cons = reader.ReadUInt();
  header_.num_objs = reader.ReadUInt();
  header_.num_eqns = -1;
  if (reader.ReadOptionalUInt(header_.num_ranges) &&
      reader.ReadOptionalUInt(header_.num_eqns)) {
      reader.ReadOptionalUInt(header_.num_logical_cons);
  }
  reader.ReadTillEndOfLine();

  // Read the nonlinear and complementarity information.
  header_.num_nl_cons = reader.ReadUInt();
  header_.num_nl_objs = reader.ReadUInt();
  bool all_compl =
      reader.ReadOptionalUInt(header_.num_compl_conds) &&
      reader.ReadOptionalUInt(header_.num_nl_compl_conds) &&
      reader.ReadOptionalUInt(header_.num_compl_dbl_ineqs) &&
      reader.ReadOptionalUInt(header_.num_compl_vars_with_nz_lb);
  header_.num_compl_conds += header_.num_nl_compl_conds;
  if (header_.num_compl_conds > 0 && !all_compl)
    header_.num_compl_dbl_ineqs = -1;
  reader.ReadTillEndOfLine();

  // Read the information about network constraints.
  header_.num_nl_net_cons = reader.ReadUInt();
  header_.num_linear_net_cons = reader.ReadUInt();
  reader.ReadTillEndOfLine();

  // Read the information about nonlinear variables.
  header_.num_nl_vars_in_cons = reader.ReadUInt();
  header_.num_nl_vars_in_objs = reader.ReadUInt();
  header_.num_nl_vars_in_both = -1;
  reader.ReadOptionalUInt(header_.num_nl_vars_in_both);
  reader.ReadTillEndOfLine();

  header_.num_linear_net_vars = reader.ReadUInt();
  header_.num_funcs = reader.ReadUInt();
  int arith = 0;
  if (reader.ReadOptionalUInt(arith)) {
    if (arith != Arith_Kind_ASL && arith != 0) {
      bool swap_bytes = false;
#if defined(IEEE_MC68k) || defined(IEEE_8087)
      swap_bytes = arith > 0 && arith + Arith_Kind_ASL == 3;
#endif
      if (!swap_bytes)
        reader.ReportParseError("unrecognized binary format");
      header_.format = NLHeader::BINARY_SWAPPED;
      // TODO: swap bytes
    }
    reader.ReadOptionalUInt(header_.flags);
  }
  reader.ReadTillEndOfLine();

  // Read the information about discrete variables.
  header_.num_linear_binary_vars = reader.ReadUInt();
  header_.num_linear_integer_vars = reader.ReadUInt();
  if (header_.num_nl_vars_in_both >= 0) {  // ampl versions >= 19930630
    header_.num_nl_integer_vars_in_both = reader.ReadUInt();
    header_.num_nl_integer_vars_in_cons = reader.ReadUInt();
    header_.num_nl_integer_vars_in_objs = reader.ReadUInt();
  }
  reader.ReadTillEndOfLine();

  // Read the information about nonzeros.
  header_.num_con_nonzeros = reader.ReadUInt();
  header_.num_obj_nonzeros = reader.ReadUInt();
  reader.ReadTillEndOfLine();

  // Read the information about names.
  header_.max_con_name_len = reader.ReadUInt();
  header_.max_var_name_len = reader.ReadUInt();
  reader.ReadTillEndOfLine();

  // Read the information about common expressions.
  header_.num_common_exprs_in_both = reader.ReadUInt();
  header_.num_common_exprs_in_cons = reader.ReadUInt();
  header_.num_common_exprs_in_objs = reader.ReadUInt();
  header_.num_common_exprs_in_cons1 = reader.ReadUInt();
  header_.num_common_exprs_in_objs1 = reader.ReadUInt();
  reader.ReadTillEndOfLine();

  handler_.BeginBuild(name.c_str(), header_, 0);
  if (header_only)
    return;

  if (header_.format != NLHeader::TEXT) {
    // TODO: switch to binary reader
  }

  for (;;) {
    char c = reader.ReadChar();
    switch (c) {
    case '\0':
      // TODO: check for end of input
      return;
    case 'C': {
      int con_index = reader.ReadUInt();
      if (con_index >= header_.num_algebraic_cons) {
        // TODO: error: constraint index out of bounds
      }
      reader.ReadTillEndOfLine();
      ReadExpr(reader);
      break;
    }
    case 'F':
      // TODO: read functions
      break;
    case 'L': {
      int lcon_index = reader.ReadUInt();
      if (lcon_index >= header_.num_logical_cons) {
        // TODO: error: logical constraint index out of bounds
      }
      reader.ReadTillEndOfLine();
      ReadExpr(reader);
      break;
    }
    case 'V':
      // TODO
      break;
    case 'G': {
      int obj_index = reader.ReadUInt();
      if (obj_index >= header_.num_objs) {
        // TODO: error: objective index out of bounds
      }
      int num_terms = reader.ReadUInt(); // TODO: check
      reader.ReadTillEndOfLine();
      ReadLinearExpr(reader, num_terms);
      // TODO: read gradient!
      break;
    }
    case 'J':
      // TODO: read Jacobian matrix
      break;
    case 'O': {
      int obj_index = reader.ReadUInt();
      if (obj_index >= header_.num_objs) {
        reader.ReportParseError(
            "objective index {} is out of bounds", obj_index);
      }
      int obj_type = reader.ReadUInt();
      reader.ReadTillEndOfLine();
      handler_.AddObj(obj_index, obj_type != 0, ReadExpr(reader));
      break;
    }
    case 'S':
      // TODO: read suffix
      break;
    case 'r':
      // TODO: read RHS
      break;
    case 'b':
      reader.ReadTillEndOfLine();
      ReadBounds(reader);
      break;
    case 'k':
      ReadColumnOffsets(reader);
      break;
    case 'x':
      // TODO
      break;
    case 'd':
      // TODO
      break;
    default:
      reader.ReportParseError("invalid segment type '{}'", c);
    }
  }
}
}

#endif  // SOLVERS_UTIL_NL_H_
