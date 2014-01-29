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

ExprFactory::ExprFactory() : asl_(0) {
  for (int i = 0; i < N_OPS; ++i)
    r_ops_[i] = reinterpret_cast<efunc*>(i);
}

ExprFactory::~ExprFactory() {
  ASL_free(&asl_);
}

void ExprFactory::Init(const NLHeader &h) {
  // TODO: move to ctor and remove this method
  ASL_free(&asl_);

  asl_ = ASL_alloc(ASL_read_fg);
  ASL_fg *asl = reinterpret_cast<ASL_fg*>(asl_);
  asl->I.r_ops_ = r_ops_;

  asl->i.ampl_options_[0] = h.num_options;
  for (int i = 0; i < ampl::MAX_NL_OPTIONS; ++i)
    asl->i.ampl_options_[i + 1] = h.options[i];
  asl->i.ampl_vbtol_ = h.ampl_vbtol;

  asl->i.n_var_ = h.num_vars;
  asl->i.n_con_ = h.num_algebraic_cons;
  asl->i.n_obj_ = h.num_objs;
  asl->i.nranges_ = h.num_ranges;
  asl->i.n_eqn_ = h.num_eqns;
  asl->i.n_lcon_ = h.num_logical_cons;

  asl->i.nlc_ = h.num_nl_cons;
  asl->i.nlo_ = h.num_nl_objs;
  asl->i.n_cc_ = h.num_compl_conds;
  asl->i.nlcc_ = h.num_nl_compl_conds;
  asl->i.ndcc_ = h.num_compl_dbl_ineqs;
  asl->i.nzlb_ = h.num_compl_vars_with_nz_lb;

  asl->i.nlnc_ = h.num_nl_net_cons;
  asl->i.lnc_ = h.num_linear_net_cons;

  asl->i.nlvc_ = h.num_nl_vars_in_cons;
  asl->i.nlvo_ = h.num_nl_vars_in_objs;
  asl->i.nlvb_ = h.num_nl_vars_in_both;

  asl->i.nwv_ = h.num_linear_net_vars;
  asl->i.nfunc_ = h.num_funcs;
  asl->i.flags = h.flags;

  asl->i.nbv_ = h.num_linear_binary_vars;
  asl->i.niv_ = h.num_linear_integer_vars;
  asl->i.nlvbi_ = h.num_nl_integer_vars_in_both;
  asl->i.nlvci_ = h.num_nl_integer_vars_in_cons;
  asl->i.nlvoi_ = h.num_nl_integer_vars_in_objs;

  asl->i.nzc_ = h.num_con_nonzeros;
  asl->i.nzo_ = h.num_obj_nonzeros;

  asl->i.maxrownamelen_ = h.max_con_name_len;
  asl->i.maxcolnamelen_ = h.max_var_name_len;

  asl->i.comb_ = h.num_common_exprs_in_both;
  asl->i.comc_ = h.num_common_exprs_in_cons;
  asl->i.como_ = h.num_common_exprs_in_objs;
  asl->i.comc1_ = h.num_common_exprs_in_cons1;
  asl->i.como1_ = h.num_common_exprs_in_objs1;

  asl->i.n_var0 = asl->i.n_var1 = asl->i.n_var_;
  asl->i.n_con0 = asl->i.n_con1 = asl->i.n_con_;
  int nlv = asl->i.nlvc_;
  if (nlv < asl->i.nlvo_)
    nlv = asl->i.nlvo_;
  if (nlv <= 0)
    nlv = 1;
  asl->i.x0len_ = nlv * sizeof(double);
  asl->i.x0kind_ = ASL_first_x;
  asl->i.n_conjac_[0] = 0;
  asl->i.n_conjac_[1] = asl->i.n_con_;
  asl->i.c_vars_ = asl->i.o_vars_ =
      asl->i.n_var_;  // confusion arises otherwise

  // TODO: allocate arrays as fg_read does
  int nv1 = asl->i.n_var_ + asl->i.nsufext[ASL_Sufkind_var];
  int ncom = 0;
  int nv = nv1 + ncom;
  int nc0 = asl->i.n_con_;
  int nc = nc0 + asl->i.nsufext[ASL_Sufkind_con];
  int no = asl->i.n_obj_;
  //int nvc = asl->i.c_vars_;
  //int nvo = asl->i.o_vars_;
  int nlcon = asl->i.n_lcon_;
  int nco = nc + no + nlcon;
  asl->i.ncom0_ = asl->i.combc_ + asl->i.como_;
  asl->i.ncom1_ = asl->i.comc1_ + asl->i.como1_;
  unsigned x =
      nco * sizeof(cde) + no * sizeof(ograd*)
    + nv * (sizeof(expr_v) + 2 * sizeof(int))
    //+ asl->i.ncom0_ * sizeof(cexp)
    + asl->i.ncom1_ * sizeof(cexp1)
    //+ nfunc * sizeof(func_info*)
    //+ nvref * sizeof(int)
    + no;
  expr_v *e = asl->I.var_e_ =
      reinterpret_cast<expr_v*>(M1zapalloc_ASL(&asl_->i, x));
  for (int i = 0; i < h.num_vars; ++i, ++e) {
    e->op = r_ops_[OPVARVAL];
    e->a = i;
  }
}

NumericConstant ExprFactory::CreateNumericConstant(double value) {
  expr_n *e = reinterpret_cast<expr_n*>(mem_ASL(asl_, asl_->i.size_expr_n_));
  e->op = reinterpret_cast<efunc_n*>(r_ops_[OPNUM]);
  e->v = value;
  return Expr::Create<NumericConstant>(reinterpret_cast<expr*>(e));
}

Variable ExprFactory::CreateVariable(int var_index) {
  assert(var_index >= 0);
  return Expr::Create<Variable>(reinterpret_cast<expr*>(
      reinterpret_cast<ASL_fg*>(asl_)->I.var_e_ + var_index));
}

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
    expr = factory_.CreateNumericConstant(reader.ReadShort());
    break;
  case 'l':
    expr = factory_.CreateNumericConstant(reader.ReadLong());
    break;
  case 'n':
    expr = factory_.CreateNumericConstant(reader.ReadDouble());
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
    expr = factory_.CreateVariable(var_index);
    break;
  }
  default:
    reader.ReportParseError("expected expression");
  }
  reader.ReadEndOfLine();
  return expr;
}

void NLReader::ReadLinearExpr(TextReader &reader, int num_terms) {
  for (int i = 0; i < num_terms; ++i) {
    reader.ReadUInt();
    reader.ReadUInt(); // TODO: read double
    // TODO
    reader.ReadEndOfLine();
  }
}

void NLReader::ReadBounds(TextReader &reader) {
  for (int i = 0; i < header_.num_vars; ++i) {
    reader.ReadUInt();
    reader.ReadUInt(); // TODO: read double
    // TODO
    reader.ReadEndOfLine();
  }
}

void NLReader::ReadColumnOffsets(TextReader &reader) {
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
  reader.ReadEndOfLine();

  // Read problem dimensions.
  header_.num_vars = reader.ReadUInt();
  header_.num_algebraic_cons = reader.ReadUInt();
  header_.num_objs = reader.ReadUInt();
  header_.num_eqns = -1;
  if (reader.ReadOptionalUInt(header_.num_ranges) &&
      reader.ReadOptionalUInt(header_.num_eqns)) {
      reader.ReadOptionalUInt(header_.num_logical_cons);
  }
  reader.ReadEndOfLine();

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
  reader.ReadEndOfLine();

  // Read the information about network constraints.
  header_.num_nl_net_cons = reader.ReadUInt();
  header_.num_linear_net_cons = reader.ReadUInt();
  reader.ReadEndOfLine();

  // Read the information about nonlinear variables.
  header_.num_nl_vars_in_cons = reader.ReadUInt();
  header_.num_nl_vars_in_objs = reader.ReadUInt();
  header_.num_nl_vars_in_both = -1;
  reader.ReadOptionalUInt(header_.num_nl_vars_in_both);
  reader.ReadEndOfLine();

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
  reader.ReadEndOfLine();

  // Read the information about discrete variables.
  header_.num_linear_binary_vars = reader.ReadUInt();
  header_.num_linear_integer_vars = reader.ReadUInt();
  if (header_.num_nl_vars_in_both >= 0) {  // ampl versions >= 19930630
    header_.num_nl_integer_vars_in_both = reader.ReadUInt();
    header_.num_nl_integer_vars_in_cons = reader.ReadUInt();
    header_.num_nl_integer_vars_in_objs = reader.ReadUInt();
  }
  reader.ReadEndOfLine();

  // Read the information about nonzeros.
  header_.num_con_nonzeros = reader.ReadUInt();
  header_.num_obj_nonzeros = reader.ReadUInt();
  reader.ReadEndOfLine();

  // Read the information about names.
  header_.max_con_name_len = reader.ReadUInt();
  header_.max_var_name_len = reader.ReadUInt();
  reader.ReadEndOfLine();

  // Read the information about common expressions.
  header_.num_common_exprs_in_both = reader.ReadUInt();
  header_.num_common_exprs_in_cons = reader.ReadUInt();
  header_.num_common_exprs_in_objs = reader.ReadUInt();
  header_.num_common_exprs_in_cons1 = reader.ReadUInt();
  header_.num_common_exprs_in_objs1 = reader.ReadUInt();
  reader.ReadEndOfLine();

  handler_->HandleHeader(header_);

  if (header_.format != NLHeader::TEXT) {
    // TODO: switch to binary reader
  }

  factory_.Init(header_);
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
      reader.ReadEndOfLine();
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
      reader.ReadEndOfLine();
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
      if (obj_index >= header_.num_objs) {
        reader.ReportParseError("objective index {} is out of bounds")
            << obj_index;
      }
      int obj_type = reader.ReadUInt();
      reader.ReadEndOfLine();
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
      reader.ReadEndOfLine();
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
