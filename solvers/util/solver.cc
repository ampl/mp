/*
 A C++ interface to an AMPL solver.

 Copyright (C) 2012 AMPL Optimization Inc

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

#include "util/solver.h"

#include <cctype>
#include <cstdarg>
#include <cstdio>
#include <cstring>

#include <algorithm>
#include <stack>

#ifndef _WIN32
# include <strings.h>
# include <unistd.h>
# define AMPL_WRITE write
#else
# include <io.h>
# define AMPL_WRITE _write
# define strcasecmp _stricmp
# undef max
#endif

#include "util/clock.h"
#include "util/format.h"
#include "util/rstparser.h"
#include "getstub.h"

namespace {

const char *SkipSpaces(const char *s) {
  while (*s && isspace(*s))
    ++s;
  return s;
}

const char *SkipNonSpaces(const char *s) {
  while (*s && !isspace(*s))
    ++s;
  return s;
}

struct Deleter {
  void operator()(ampl::SolverOption* p) const { delete p; }
};

// A reStructuredText (RST) formatter.
class RSTFormatter : public rst::ContentHandler {
 private:
  fmt::Writer &writer_;
  const ampl::EnumOptionValue *values_;
  std::stack<int> indents_;
  int indent_;
  int pos_in_line_;
  bool end_block_;  // true if there was no text after the last end of block

  enum {
    LITERAL_BLOCK_INDENT = 3,
    LIST_ITEM_INDENT     = 2
  };

  void Indent() {
    if (end_block_) {
      end_block_ = false;
      EndLine();
    }
    for (; pos_in_line_ < indent_; ++pos_in_line_)
      writer_ << ' ';
  }

  // Ends the current line.
  void EndLine() {
    writer_ << '\n';
    pos_in_line_ = 0;
  }

  // Writes a string doing a text reflow.
  void Write(fmt::StringRef s);

 public:
  RSTFormatter(fmt::Writer &w, const ampl::EnumOptionValue *values, int indent)
  : writer_(w), values_(values), indent_(indent),
    pos_in_line_(0), end_block_(false) {}

  void StartBlock(rst::BlockType type);
  void EndBlock();

  void HandleText(const char *text, std::size_t size);
  void HandleDirective(const char *type);
};

void RSTFormatter::Write(fmt::StringRef s) {
  enum { MAX_LINE_LENGTH = 78 };
  const char *p = s.c_str();
  for (;;) {
    // Skip leading spaces.
    while (*p == ' ')
      ++p;
    const char *word_start = p;
    while (*p != ' ' && *p != '\n' && *p)
      ++p;
    const char *word_end = p;
    if (pos_in_line_ + (word_end - word_start) +
        (pos_in_line_ == 0 ? 0 : 1) > MAX_LINE_LENGTH) {
      EndLine();  // The word doesn't fit in the current line, start a new one.
    }
    if (pos_in_line_ == 0) {
      Indent();
    } else {
      // Separate words.
      writer_ << ' ';
      ++pos_in_line_;
    }
    std::size_t length = word_end - word_start;
    writer_ << fmt::StringRef(word_start, length);
    pos_in_line_ += static_cast<int>(length);
    if (*p == '\n') {
      EndLine();
      ++p;
    }
    if (!*p) break;
  }
}

void RSTFormatter::StartBlock(rst::BlockType type) {
  indents_.push(indent_);
  if (type == rst::LITERAL_BLOCK) {
    indent_ += LITERAL_BLOCK_INDENT;
  } else if (type == rst::LIST_ITEM) {
    Write("*");
    indent_ += LIST_ITEM_INDENT;
  }
}

void RSTFormatter::EndBlock() {
  indent_ = indents_.top();
  indents_.pop();
  end_block_ = true;
}

void RSTFormatter::HandleText(const char *text, std::size_t size) {
  Write(std::string(text, size));
  size = writer_.size();
  if (size != 0 && writer_.data()[size - 1] != '\n')
    EndLine();
}

void RSTFormatter::HandleDirective(const char *type) {
  // TODO: test
  if (std::strcmp(type, "value-table") != 0)
    ampl::ThrowError("unknown directive {}") << type;
  if (!values_)
    ampl::ThrowError("no values to format");
  std::size_t max_len = 0;
  for (const ampl::EnumOptionValue *v = values_; v->value; ++v)
    max_len = std::max(max_len, std::strlen(v->value));

  // If the values don't have descriptions indent them as list items.
  if (!values_->description)
    indent_ += LIST_ITEM_INDENT;
  for (const ampl::EnumOptionValue *v = values_; v->value; ++v) {
    Indent();
    writer_ << fmt::pad(v->value, max_len);
    if (v->description) {
      static const char SEP[] = " -";
      writer_ << SEP;
      int saved_indent = indent_;
      indent_ += max_len + sizeof(SEP);
      pos_in_line_ = indent_;
      Write(v->description);
      indent_ = saved_indent;
    }
    EndLine();
  }

  // Unindent values if necessary and end the block.
  if (!values_->description)
    indent_ -= LIST_ITEM_INDENT;
  end_block_ = true;
}
}

namespace ampl {

namespace internal {

// Formats restructured text.
void FormatRST(fmt::Writer &w, fmt::StringRef s, int indent,
    const ampl::EnumOptionValue *values) {
  RSTFormatter formatter(w, values, indent);
  rst::Parser parser(&formatter);
  parser.Parse(s.c_str());
}

const char OptionHelper<int>::TYPE_NAME[] = "int";
const char OptionHelper<double>::TYPE_NAME[] = "double";
const char OptionHelper<std::string>::TYPE_NAME[] = "string";
int OptionHelper<int>::Parse(const char *&s) {
  char *end = 0;
  long value = std::strtol(s, &end, 10);
  s = end;
  return value;
}

void OptionHelper<double>::Write(fmt::Writer &w, Arg value) {
  char buffer[32];
  g_fmt(buffer, value);
  w << buffer;
}

double OptionHelper<double>::Parse(const char *&s) {
  char *end = 0;
  double value = strtod_ASL(s, &end);
  s = end;
  return value;
}

std::string OptionHelper<std::string>::Parse(const char *&s) {
  const char *start = s;
  s = SkipNonSpaces(s);
  return std::string(start, s - start);
}
}

std::string SignalHandler::signal_message_;
const char *SignalHandler::signal_message_ptr_;
unsigned SignalHandler::signal_message_size_;
Interruptible *SignalHandler::interruptible_;

// Set stop_ to 1 initially to avoid accessing handler_ which may not be atomic.
volatile std::sig_atomic_t SignalHandler::stop_ = 1;

SignalHandler::SignalHandler(const Solver &s, Interruptible *i) {
  signal_message_ = str(fmt::Format("\n<BREAK> ({})\n") << s.name());
  signal_message_ptr_ = signal_message_.c_str();
  signal_message_size_ = static_cast<unsigned>(signal_message_.size());
  interruptible_ = i;
  std::signal(SIGINT, HandleSigInt);
  stop_ = 0;
}

void SignalHandler::HandleSigInt(int sig) {
  unsigned count = 0;
  do {
    // Use asynchronous-safe function write instead of printf!
    int result = AMPL_WRITE(1, signal_message_ptr_ + count,
        signal_message_size_ - count);
    if (result < 0) break;
    count += result;
  } while (count < signal_message_size_);
  if (stop_) {
    // Use asynchronous-safe function _exit instead of exit!
    _exit(1);
  }
  stop_ = 1;
  if (interruptible_)
    interruptible_->Interrupt();
  // Restore the handler since it might have been reset before the handler
  // is called (this is implementation defined).
  std::signal(sig, HandleSigInt);
}

void Solver::SolutionWriter::HandleFeasibleSolution(
    Problem &p, fmt::StringRef message, const double *values,
    const double *dual_values, double) {
  ++solver_->num_solutions_;
  if (solver_->solution_stub_.empty())
    return;
  fmt::Writer w;
  w << solver_->solution_stub_ << solver_->num_solutions_ << ".sol";
  write_solf_ASL(reinterpret_cast<ASL*>(p.asl_),
      const_cast<char*>(message.c_str()), const_cast<double*>(values),
      const_cast<double*>(dual_values), solver_, w.c_str());
}

void Solver::SolutionWriter::HandleSolution(
    Problem &p, fmt::StringRef message, const double *values,
    const double *dual_values, double) {
  if (solver_->need_multiple_solutions()) {
    Suffix nsol_suffix = p.suffix("nsol", ASL_Sufkind_prob);
    if (nsol_suffix)
      nsol_suffix.set_values(&solver_->num_solutions_);
  }
  write_sol_ASL(reinterpret_cast<ASL*>(p.asl_),
      const_cast<char*>(message.c_str()), const_cast<double*>(values),
      const_cast<double*>(dual_values), solver_);
}

bool Solver::OptionNameLess::operator()(
    const SolverOption *lhs, const SolverOption *rhs) const {
  return strcasecmp(lhs->name(), rhs->name()) < 0;
}

void Solver::RegisterSuffixes(Problem &p) {
  ASL *asl = reinterpret_cast<ASL*>(p.asl_);
  if (asl->i.nsuffixes == 0 && !suffixes_.empty())
    suf_declare_ASL(asl, &suffixes_[0], suffixes_.size());
}

char *Solver::PrintOptionsAndExit(Option_Info *oi, keyword *, char *) {
  Solver *solver = static_cast<Solver*>(oi);
  fmt::Writer writer;
  internal::FormatRST(writer, solver->option_header_);
  if (!solver->option_header_.empty())
    writer << '\n';
  writer << "Options:\n";
  const int DESC_INDENT = 6;
  const OptionSet &options = solver->options_;
  for (OptionSet::const_iterator i = options.begin(); i != options.end(); ++i) {
    SolverOption *opt = *i;
    writer << '\n' << opt->name() << '\n';
    internal::FormatRST(writer, opt->description(), DESC_INDENT, opt->values());
  }
  std::fwrite(writer.data(), writer.size(), 1, stdout);
  exit(0);
  return 0;
}

Solver::Solver(
    fmt::StringRef name, fmt::StringRef long_name, long date, unsigned flags)
: name_(name), has_errors_(false), num_solutions_(0), count_solutions_(false),
  read_flags_(0), timing_(false) {
  error_handler_ = this;
  output_handler_ = this;
  sol_writer_.set_solver(this);
  sol_handler_ = &sol_writer_;

  // Workaround for GCC bug 30111 that prevents value-initialization of
  // the base POD class.
  Option_Info init = {};
  Option_Info &self = *this;
  self = init;

  sname = const_cast<char*>(name_.c_str());
  if (long_name.c_str()) {
    long_name_ = long_name;
    bsname = const_cast<char*>(long_name_.c_str());
  } else {
    bsname = sname;
  }
  options_var_name_ = name_;
  options_var_name_ += "_options";
  opname = const_cast<char*>(options_var_name_.c_str());
  Option_Info::version = bsname;
  driver_date = date;

  struct VersionOption : SolverOption {
    Solver &s;
    VersionOption(Solver &s) : SolverOption("version",
        "Single-word phrase: report version details "
        "before solving the problem.", 0, true), s(s) {}

    void Write(fmt::Writer &w) {
      w << ((s.flags() & ASL_OI_show_version) != 0);
    }
    void Parse(const char *&) {
      s.Option_Info::flags |= ASL_OI_show_version;
    }
  };
  AddOption(OptionPtr(new VersionOption(*this)));

  struct WantSolOption : TypedSolverOption<int> {
    Solver &s;
    WantSolOption(Solver &s) : TypedSolverOption<int>("wantsol",
        "In a stand-alone invocation (no ``-AMPL`` on the command line), "
        "what solution information to write.  Sum of\n"
        "\n"
        "| 1 - write ``.sol`` file\n"
        "| 2 - primal variables to stdout\n"
        "| 4 - dual variables to stdout\n"
        "| 8 - suppress solution message\n"), s(s) {}

    int GetValue() const { return s.wantsol(); }
    void SetValue(int value) {
      if ((value & ~0xf) != 0)
        throw InvalidOptionValue("wantsol", value);
      s.Option_Info::wantsol = value;
    }
  };
  AddOption(OptionPtr(new WantSolOption(*this)));

  struct BoolOption : TypedSolverOption<int> {
    bool &value_;
    BoolOption(bool &value, const char *name, const char *description)
    : TypedSolverOption<int>(name, description), value_(value) {}

    int GetValue() const { return value_; }
    void SetValue(int value) {
      if (value != 0 && value != 1)
        throw InvalidOptionValue(name(), value);
      value_ = value != 0;
    }
  };
  AddOption(OptionPtr(new BoolOption(timing_, "timing",
      "0 or 1 (default 0): Whether to display timings for the run.\n")));

  if ((flags & MULTIPLE_SOL) != 0) {
    AddSuffix("nsol", 0, ASL_Sufkind_prob | ASL_Sufkind_outonly);

    AddOption(OptionPtr(new BoolOption(count_solutions_, "countsolutions",
        "0 or 1 (default 0): Whether to count the number of solutions "
        "and return it in the ``.nsol`` problem suffix.")));

    struct StringOption : TypedSolverOption<std::string> {
      std::string &value_;
      StringOption(std::string &value,
          const char *name, const char *description)
      : TypedSolverOption<std::string>(name, description), value_(value) {}

      std::string GetValue() const { return value_; }
      void SetValue(const char *value) { value_ = value; }
    };
    AddOption(OptionPtr(new StringOption(solution_stub_, "solutionstub",
        "Stub for solution files.  If ``solutionstub`` is specified, "
        "found solutions are written to files (``solutionstub & '1' & "
        "'.sol'``) ... (``solutionstub & Current.nsol & '.sol'``), where "
        "``Current.nsol`` holds the number of returned solutions.  That is, "
        "file names are obtained by appending 1, 2, ... ``Current.nsol`` to "
        "``solutionstub``.")));
  }

  cl_option_ = keyword();
  cl_option_.name = const_cast<char*>("=");
  cl_option_.desc = const_cast<char*>("show name= possibilities");
  cl_option_.kf = Solver::PrintOptionsAndExit;
  cl_option_.info = 0;
  options = &cl_option_;
  n_options = 1;
}

Solver::~Solver() {
  std::for_each(options_.begin(), options_.end(), Deleter());
}

void Solver::AddSuffix(
    const char *name, const char *table, int kind, int nextra) {
  suffixes_.push_back(SufDecl());
  SufDecl &sd = suffixes_.back();
  sd.name = const_cast<char*>(name);
  sd.table = const_cast<char*>(table);
  sd.kind = kind;
  sd.nextra = nextra;
}

bool Solver::ProcessArgs(char **&argv, Problem &p, unsigned flags) {
  ASL *asl = reinterpret_cast<ASL*>(p.asl_);
  RegisterSuffixes(p);

  // Make "-?" show some command-line options that are relevant when
  // importing functions.
  Option_Info::flags |= ASL_OI_want_funcadd;

  char *stub = getstub_ASL(asl, &argv, this);
  if (!stub) {
    usage_noexit_ASL(this, 1);
    return false;
  }
  steady_clock::time_point start = steady_clock::now();
  p.Read(stub, read_flags_);
  double read_time = GetTimeAndReset(start);
  bool result = ParseOptions(argv, flags, &p);
  if (timing_)
    Print("Input time = {:.6f}s\n") << read_time;
  return result;
}

SolverOption *Solver::FindOption(const char *name) const {
  struct DummyOption : SolverOption {
    DummyOption(const char *name) : SolverOption(name, 0) {}
    void Write(fmt::Writer &) {}
    void Parse(const char *&) {}
  };
  DummyOption option(name);
  OptionSet::const_iterator i = options_.find(&option);
  return i != options_.end() ? *i : 0;
}

void Solver::ParseOptionString(const char *s, unsigned flags) {
  bool skip = false;
  for (;;) {
    if (!*(s = SkipSpaces(s)))
      return;

    // Parse the option name.
    const char *name_start = s;
    while (*s && !std::isspace(*s) && *s != '=')
      ++s;
    fmt::internal::Array<char, 50> name;
    std::size_t name_size = s - name_start;
    name.resize(name_size + 1);
    for (std::size_t i = 0; i < name_size; ++i)
      name[i] = name_start[i];
    name[name_size] = 0;

    // Parse the option value.
    bool equal_sign = false;
    s = SkipSpaces(s);
    if (*s == '=') {
      s = SkipSpaces(s + 1);
      equal_sign = true;
    }

    nnl = 0;
    SolverOption *opt = FindOption(&name[0]);
    if (!opt) {
      if (!skip)
        HandleUnknownOption(&name[0]);
      if (equal_sign) {
        s = SkipNonSpaces(s);
      } else {
        // Skip everything until the next known option if there is no "="
        // because it is impossible to know whether the next token is an
        // option name or a value.
        // For example, if "a" in "a b c" is an unknown option, then "b"
        // can be either a value of option "a" or another option.
        skip = true;
      }
      continue;
    }

    skip = false;
    if (*s == '?') {
      char next = s[1];
      if (!next || std::isspace(next)) {
        ++s;
        if ((flags & NO_OPTION_ECHO) == 0) {
          fmt::Writer w;
          w.Format("{}=") << &name[0];
          opt->Write(w);
          w << '\n';
          Print(fmt::StringRef(w.c_str(), w.size()));
        }
        continue;
      }
    }
    if (opt->is_keyword() && equal_sign) {
      ReportError("Option \"{}\" doesn't accept argument") << &name[0];
      s = SkipNonSpaces(s);
      continue;
    }
    try {
      opt->Parse(s);
    } catch (const OptionError &e) {
      ReportError("{}") << e.what();
    }
    if ((flags & NO_OPTION_ECHO) == 0)
      Print("{}\n") << std::string(name_start, s - name_start);
  }
}

bool Solver::ParseOptions(char **argv, unsigned flags, const Problem *) {
  has_errors_ = false;
  Option_Info::flags &= ~ASL_OI_show_version;
  if (opname) {
    if (const char *s = getenv(opname))
      ParseOptionString(s, flags);
  }
  while (const char *s = *argv++)
    ParseOptionString(s, flags);
  if (this->flags() & ASL_OI_show_version)
    show_version_ASL(this);
  std::fflush(stdout);
  return !has_errors_;
}

void Solver::Solve(Problem &p) {
  num_solutions_ = 0;
  RegisterSuffixes(p);
  DoSolve(p);
}

int Solver::Run(char **argv) {
  Problem p;
  if (!ProcessArgs(argv, p))
    return 1;
  Solve(p);
  return 0;
}
}
