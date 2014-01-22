/*
 An nl reader.

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

#include "solvers/util/nlreader.h"

#include "solvers/util/os.h"
#include <cctype>

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

NLHandler::~NLHandler() {}

class TextReader {
 private:
  const char *ptr_;
  const char *line_start_;
  const char *token_;  // start of the current token
  std::string name_;
  int line_;

  int DoReadInt() {
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

  void ReadEndOfLine() {
    while (*ptr_ && *ptr_ != '\n')
      ++ptr_;
    ++ptr_;
  }

  int ReadInt() {
    SkipSpace();
    char c = *ptr_;
    if (c < '0' || c > '9')
      ReportParseError("expected integer");
    return DoReadInt();
  }

  bool ReadOptionalInt(int &value) {
    SkipSpace();
    char c = *ptr_;
    bool has_value = c >= '0' && c <= '9';
    if (has_value)
      value = DoReadInt();
    return has_value;
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

void NLReader::ReadExpr(TextReader &reader) {
  char c = reader.ReadChar();
  if (c == 'n')
    reader.ReadInt(); // TODO: read double
  // TODO: other types of expressions
  reader.ReadEndOfLine();
}

void NLReader::ReadLinearExpr(TextReader &reader, int num_terms) {
  for (int i = 0; i < num_terms; ++i) {
    reader.ReadInt();
    reader.ReadInt(); // TODO: read double
    // TODO
    reader.ReadEndOfLine();
  }
}

void NLReader::ReadBounds(TextReader &reader, int num_bounds) {
  for (int i = 0; i < num_bounds; ++i) {
    reader.ReadInt();
    reader.ReadInt(); // TODO: read double
    // TODO
    reader.ReadEndOfLine();
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
  } handler;
  if (!handler_)
    handler_ = &handler;

  TextReader reader(name, str.c_str());

  // TODO: always read header as text

  // Read the format (text or binary).
  bool binary = false;
  switch (char c = reader.ReadChar()) {
  case 'g':
    break;
  case 'b':
    binary = true;
    break;
  default:
    reader.ReportParseError("invalid format '{}'") << c;
    break;
  }

  // Read options.
  NLHeader header = {};
  header.num_options = reader.ReadInt();
  if (header.num_options > MAX_NL_OPTIONS)
    reader.ReportParseError("too many options");
  for (int i = 0; i < header.num_options; ++i) {
    if (!reader.ReadOptionalInt(header.options[i]))
      break;
  }
  if (header.options[VBTOL_OPTION] == READ_VBTOL)
    reader.ReadOptionalDouble(header.ampl_vbtol);
  reader.ReadEndOfLine();
  handler_->HandleHeader(header);

  header.num_vars = reader.ReadInt();
  header.num_cons = reader.ReadInt();
  header.num_objs = reader.ReadInt();
  reader.SkipSpace();
  header.num_eqns = -1;
  if (reader.ReadOptionalInt(header.num_ranges))
    reader.ReadOptionalInt(header.num_eqns);
  reader.ReadEndOfLine();

  // Read the nonlinear and complementarity information.
  header.num_nl_cons = reader.ReadInt();
  header.num_nl_objs = reader.ReadInt();
  bool all_compl = reader.ReadOptionalInt(header.num_compl_conds) &&
      reader.ReadOptionalInt(header.num_nl_compl_conds) &&
      reader.ReadOptionalInt(header.num_compl_dbl_ineq) &&
      reader.ReadOptionalInt(header.num_compl_vars_with_nz_lb);
  header.num_compl_conds += header.num_nl_compl_conds;
  if (header.num_compl_conds > 0 && !all_compl)
    header.num_compl_dbl_ineq = -1;
  reader.ReadEndOfLine();

  // Read the information about network constraints.
  header.num_nl_net_cons = reader.ReadInt();
  header.num_linear_net_cons = reader.ReadInt();
  reader.ReadEndOfLine();

  // Read the information about nonlinear variables.
  header.num_nl_vars_in_cons = reader.ReadInt();
  header.num_nl_vars_in_objs = reader.ReadInt();
  header.num_nl_vars_in_both = -1;
  reader.ReadOptionalInt(header.num_nl_vars_in_both);
  reader.ReadEndOfLine();

  header.num_linear_net_vars = reader.ReadInt();
  header.num_funcs = reader.ReadInt();
  int arith_kind = 0, flags = 0;
  if (reader.ReadOptionalInt(arith_kind))
    reader.ReadOptionalInt(flags);
  // TODO: resolve the mystery with flags
  reader.ReadEndOfLine();

  // Read the information about discrete variables.
  header.num_linear_binary_vars = reader.ReadInt();
  header.num_linear_integer_vars = reader.ReadInt();
  if (header.num_nl_vars_in_both >= 0) {  // ampl versions >= 19930630
    header.num_nl_integer_vars_in_both = reader.ReadInt();
    header.num_nl_integer_vars_in_cons = reader.ReadInt();
    header.num_nl_integer_vars_in_objs = reader.ReadInt();
  }
  reader.ReadEndOfLine();

  // Read the information about nonzeros.
  header.num_con_nonzeros = reader.ReadInt();
  header.num_obj_nonzeros = reader.ReadInt();
  reader.ReadEndOfLine();

  // Read the information about names.
  header.max_con_name_len = reader.ReadInt();
  header.max_var_name_len = reader.ReadInt();
  reader.ReadEndOfLine();

  // Read the information about common expressions.
  header.num_common_b_exprs = reader.ReadInt(); // FIXME: what is b (both?)
  header.num_common_con_exprs = reader.ReadInt();
  header.num_common_obj_exprs = reader.ReadInt();
  header.num_common_con1_exprs = reader.ReadInt();
  header.num_common_obj1_exprs = reader.ReadInt();
  reader.ReadEndOfLine();

  for (;;) {
    char c = reader.ReadChar();
    switch (c) {
    case '\0':
      // TODO: check EOF end of input
      return;
    case 'C': {
      int con_index = reader.ReadInt();
      if (con_index >= header.num_cons) {
        // TODO: error: constraint index out of bounds
      }
      reader.ReadEndOfLine();
      ReadExpr(reader);
      break;
    }
    case 'F':
      // TODO: read functions
      break;
    case 'L':
      // TODO: read logical constraints
      break;
    case 'V':
      // TODO
      break;
    case 'G': {
      int obj_index = reader.ReadInt();
      if (obj_index < 0 || obj_index >= header.num_objs) {
        // TODO: error: objective index out of bounds
      }
      int num_terms = reader.ReadInt(); // TODO: check
      reader.ReadEndOfLine();
      ReadLinearExpr(reader, num_terms);
      // TODO: read gradient!
      break;
    }
    case 'J':
      // TODO: read Jacobian matrix
      break;
    case 'O': {
      int obj_index = reader.ReadInt();
      if (obj_index < 0 || obj_index >= header.num_objs) {
        // TODO: error: objective index out of bounds
      }
      int obj_type = reader.ReadInt();
      reader.ReadEndOfLine();
      ReadExpr(reader);
      break;
    }
    case 'S':
      // TODO: read suffix
      break;
    case 'r':
      // TODO: read RHS
      break;
    case 'b':
      reader.ReadEndOfLine();
      ReadBounds(reader, header.num_vars);
      break;
    case 'k': {
      int value = reader.ReadInt();
      // TODO: what is value?
      reader.ReadEndOfLine();
      break;
    }
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
