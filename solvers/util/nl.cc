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

#include "solvers/util/nl.h"

#include "solvers/arith.h"
#include "solvers/util/os.h"
#include "solvers/asl.h"

#include <cctype>

#undef ASL_SWAP_BYTES
#if defined(IEEE_MC68k) || defined(IEEE_8087)
# define ASL_SWAP_BYTES
#endif

namespace {

class ParseErrorReporter {
 private:
  fmt::StringRef name_;
  int line_;
  int column_;

 public:
  ParseErrorReporter(fmt::StringRef name, int line, int column)
  : name_(name), line_(line), column_(column) {}

  void operator()(const fmt::Writer &w) const {
    throw ampl::ParseError(name_, line_, column_,
        fmt::Format("{}:{}:{}: {}") << name_ << line_ << column_
          << fmt::StringRef(w.c_str(), w.size()));
  }
};
}

namespace ampl {

fmt::Writer &operator<<(fmt::Writer &w, const NLHeader &h) {
  w << (h.format == NLHeader::TEXT ? 'g' : 'b') << h.num_options;
  for (int i = 0; i < h.num_options; ++i)
    w << ' ' << h.options[i];
  if (h.options[VBTOL_OPTION] == READ_VBTOL)
    w << ' ' << h.ampl_vbtol;
  w << '\n';
  w.Format(" {} {} {} {} {} {}\n")
      << h.num_vars << h.num_algebraic_cons << h.num_objs
      << h.num_ranges << h.num_eqns << h.num_logical_cons;
  w.Format(" {} {} {} {} {} {}\n")
      << h.num_nl_cons << h.num_nl_objs
      << (h.num_compl_conds - h.num_nl_compl_conds)
      << h.num_nl_compl_conds << h.num_compl_dbl_ineqs
      << h.num_compl_vars_with_nz_lb;
  w.Format(" {} {}\n")
      << h.num_nl_net_cons << h.num_linear_net_cons;
  w.Format(" {} {} {}\n")
      << h.num_nl_vars_in_cons << h.num_nl_vars_in_objs
      << h.num_nl_vars_in_both;
  w.Format(" {} {} 0 {}\n")
      << h.num_linear_net_vars << h.num_funcs << h.flags;
  w.Format(" {} {} {} {} {}\n")
      << h.num_linear_binary_vars << h.num_linear_integer_vars
      << h.num_nl_integer_vars_in_both << h.num_nl_integer_vars_in_cons
      << h.num_nl_integer_vars_in_objs;
  w.Format(" {} {}\n")
      << h.num_con_nonzeros << h.num_obj_nonzeros;
  w.Format(" {} {}\n")
      << h.max_con_name_len << h.max_var_name_len;
  w.Format(" {} {} {} {} {}\n")
      << h.num_common_exprs_in_both << h.num_common_exprs_in_cons
      << h.num_common_exprs_in_objs << h.num_common_exprs_in_cons1
      << h.num_common_exprs_in_objs1;
  return w;
}

NLHandler::~NLHandler() {}

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

  fmt::Formatter<ParseErrorReporter> ReportParseError(fmt::StringRef message) {
    fmt::Formatter<ParseErrorReporter> f(message,
        ParseErrorReporter(name_, line_, token_ - line_start_ + 1));
    return f;
  }

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

NumericExpr NLReader::ReadExpr(TextReader &reader) {
  NumericExpr expr;
  switch (reader.ReadChar()) {
  case 'f':
    // TODO: implement function
    break;
  case 'h':
    // TODO: implement string
    break;
  case 's':
    expr = factory_->CreateNumericConstant(reader.ReadShort());
    break;
  case 'l':
    expr = factory_->CreateNumericConstant(reader.ReadLong());
    break;
  case 'n':
    expr = factory_->CreateNumericConstant(reader.ReadDouble());
    break;
  case 'o':
    // TODO: implement expression
    break;
  case 'v': {
    // TODO: variable index can be greater than num_vars
    int var_index = reader.ReadUInt();
    if (var_index >= header_.num_vars) {
      reader.ReportParseError("variable index {} is out of bounds")
          << var_index;
    }
    expr = factory_->CreateVariable(var_index);
    break;
  }
  default:
    reader.ReportParseError("expected expression");
  }
  reader.ReadTillEndOfLine();
  return expr;
}

void NLReader::ReadLinearExpr(TextReader &reader, int num_terms) {
  for (int i = 0; i < num_terms; ++i) {
    reader.ReadUInt();
    reader.ReadUInt(); // TODO: read double
    // TODO
    reader.ReadTillEndOfLine();
  }
}

void NLReader::ReadBounds(TextReader &reader) {
  for (int i = 0; i < header_.num_vars; ++i) {
    reader.ReadUInt();
    reader.ReadUInt(); // TODO: read double
    // TODO
    reader.ReadTillEndOfLine();
  }
}

void NLReader::ReadColumnOffsets(TextReader &reader) {
  int count = reader.ReadUInt(); // TODO
  reader.ReadTillEndOfLine();
  for (int i = 0; i < count; ++i) {
    reader.ReadUInt(); // TODO
    reader.ReadTillEndOfLine();
  }
}

void NLReader::ReadFile(fmt::StringRef filename) {
  MemoryMappedFile file(filename);
  // TODO: use a buffer instead of mmap if mmap is not available or the
  //       file length is a multiple of the page size
  ReadString(fmt::StringRef(file.start(), file.size()), filename);
}

void NLReader::ReadString(fmt::StringRef str, fmt::StringRef name) {
  class DefaultNLHandler : public NLHandler {
   public:
    void HandleHeader(const NLHeader &) {}
    void HandleObj(int, bool, NumericExpr) {}
  } handler;
  if (!handler_)
    handler_ = &handler;

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
#ifdef ASL_SWAP_BYTES
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

  handler_->HandleHeader(header_);

  if (header_.format != NLHeader::TEXT) {
    // TODO: switch to binary reader
  }

  ExprFactory ef(header_);
  factory_ = &ef;
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
        reader.ReportParseError("objective index {} is out of bounds")
            << obj_index;
      }
      int obj_type = reader.ReadUInt();
      reader.ReadTillEndOfLine();
      handler_->HandleObj(obj_index, obj_type != 0, ReadExpr(reader));
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
      reader.ReportParseError("invalid segment type '{}'") << c;
    }
  }
}
}
