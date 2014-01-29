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

#include "solvers/arith.h"
#include "solvers/util/os.h"
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

  void ReadEndOfLine() {
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
    reader.ReadUInt(); // TODO: read double
  // TODO: other types of expressions
  reader.ReadEndOfLine();
}

void NLReader::ReadLinearExpr(TextReader &reader, int num_terms) {
  for (int i = 0; i < num_terms; ++i) {
    reader.ReadUInt();
    reader.ReadUInt(); // TODO: read double
    // TODO
    reader.ReadEndOfLine();
  }
}

void NLReader::ReadBounds(TextReader &reader, int num_bounds) {
  for (int i = 0; i < num_bounds; ++i) {
    reader.ReadUInt();
    reader.ReadUInt(); // TODO: read double
    // TODO
    reader.ReadEndOfLine();
  }
}

void NLReader::ReadColumnOffsets(TextReader &reader, int num_vars) {
  int count = reader.ReadUInt(); // TODO
  reader.ReadEndOfLine();
  for (int i = 0; i < count; ++i) {
    reader.ReadUInt(); // TODO
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
  reader.ReadEndOfLine();

  // Read problem dimensions.
  header.num_vars = reader.ReadUInt();
  header.num_cons = reader.ReadUInt();
  header.num_objs = reader.ReadUInt();
  header.num_eqns = -1;
  if (reader.ReadOptionalUInt(header.num_ranges) &&
      reader.ReadOptionalUInt(header.num_eqns)) {
      reader.ReadOptionalUInt(header.num_logical_cons);
      // Include the number of logical constraints in the total number of
      // constraints for consistency.
      header.num_cons += header.num_logical_cons;
  }
  reader.ReadEndOfLine();

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
  reader.ReadEndOfLine();

  // Read the information about network constraints.
  header.num_nl_net_cons = reader.ReadUInt();
  header.num_linear_net_cons = reader.ReadUInt();
  reader.ReadEndOfLine();

  // Read the information about nonlinear variables.
  header.num_nl_vars_in_cons = reader.ReadUInt();
  header.num_nl_vars_in_objs = reader.ReadUInt();
  header.num_nl_vars_in_both = -1;
  reader.ReadOptionalUInt(header.num_nl_vars_in_both);
  reader.ReadEndOfLine();

  header.num_linear_net_vars = reader.ReadUInt();
  header.num_funcs = reader.ReadUInt();
  int arith = 0;
  if (reader.ReadOptionalUInt(arith)) {
    if (arith != Arith_Kind_ASL && arith != 0) {
      bool swap_bytes = false;
#ifdef ASL_SWAP_BYTES
      swap_bytes = arith > 0 && arith + Arith_Kind_ASL == 3;
#endif
      if (swap_bytes) {
        // TODO: swap bytes
      } else {
        reader.ReportParseError("unrecognized binary format");
      }
    }
    reader.ReadOptionalUInt(header.flags);
  }
  reader.ReadEndOfLine();

  // Read the information about discrete variables.
  header.num_linear_binary_vars = reader.ReadUInt();
  header.num_linear_integer_vars = reader.ReadUInt();
  if (header.num_nl_vars_in_both >= 0) {  // ampl versions >= 19930630
    header.num_nl_integer_vars_in_both = reader.ReadUInt();
    header.num_nl_integer_vars_in_cons = reader.ReadUInt();
    header.num_nl_integer_vars_in_objs = reader.ReadUInt();
  }
  reader.ReadEndOfLine();

  // Read the information about nonzeros.
  header.num_con_nonzeros = reader.ReadUInt();
  header.num_obj_nonzeros = reader.ReadUInt();
  reader.ReadEndOfLine();

  // Read the information about names.
  header.max_con_name_len = reader.ReadUInt();
  header.max_var_name_len = reader.ReadUInt();
  reader.ReadEndOfLine();

  // Read the information about common expressions.
  header.num_common_exprs_in_both = reader.ReadUInt();
  header.num_common_exprs_in_cons = reader.ReadUInt();
  header.num_common_exprs_in_objs = reader.ReadUInt();
  header.num_common_exprs_in_cons1 = reader.ReadUInt();
  header.num_common_exprs_in_objs1 = reader.ReadUInt();
  reader.ReadEndOfLine();

  handler_->HandleHeader(header);

  for (;;) {
    char c = reader.ReadChar();
    switch (c) {
    case '\0':
      // TODO: check for end of input
      return;
    case 'C': {
      int con_index = reader.ReadUInt();
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
    case 'L': {
      int lcon_index = reader.ReadUInt();
      if (lcon_index >= header.num_logical_cons) {
        // TODO: error: logical constraint index out of bounds
      }
      reader.ReadEndOfLine();
      ReadExpr(reader);
      break;
    }
    case 'V':
      // TODO
      break;
    case 'G': {
      int obj_index = reader.ReadUInt();
      if (obj_index >= header.num_objs) {
        // TODO: error: objective index out of bounds
      }
      int num_terms = reader.ReadUInt(); // TODO: check
      reader.ReadEndOfLine();
      ReadLinearExpr(reader, num_terms);
      // TODO: read gradient!
      break;
    }
    case 'J':
      // TODO: read Jacobian matrix
      break;
    case 'O': {
      int obj_index = reader.ReadUInt();
      if (obj_index >= header.num_objs) {
        // TODO: error: objective index out of bounds
      }
      int obj_type = reader.ReadUInt();
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
    case 'k':
      ReadColumnOffsets(reader, header.num_vars);
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
