/*
 Utilities for writing AMPL solvers.

 Copyright (C) 2012 AMPL Optimization LLC

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization LLC disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Victor Zverovich
 */

#include "solvers/util/solver.h"

#include <cstdarg>
#include <cstdio>

#ifndef _WIN32
# include <strings.h>
# include <unistd.h>
# define AMPL_WRITE write
#else
# include <io.h>
# define AMPL_WRITE _write
# define strcasecmp _stricmp
#endif

#include "solvers/util/clock.h"
#include "solvers/util/format.h"
#include "solvers/getstub.h"

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
}

namespace ampl {

namespace internal {
std::string IndentAndWordWrap(fmt::StringRef s, int indent) {
  std::ostringstream os;
  bool new_line = true;
  int line_offset = 0;
  int start_indent = indent;
  const int MAX_LINE_LENGTH = 78;
  const char *p = s.c_str();
  for (;;) {
    const char *start = p;
    while (*p == ' ')
      ++p;
    const char *word_start = p;
    while (*p != ' ' && *p != '\n' && *p)
      ++p;
    const char *word_end = p;
    if (new_line) {
      indent = start_indent + static_cast<int>(word_start - start);
      new_line = false;
    }
    if (line_offset + (word_end - start) > MAX_LINE_LENGTH) {
      // The word doesn't fit, start a new line.
      os << '\n';
      line_offset = 0;
    }
    if (line_offset == 0) {
      // Indent the line.
      for (; line_offset < indent; ++line_offset)
        os << ' ';
      start = word_start;
    }
    os.write(start, word_end - start);
    line_offset += static_cast<int>(word_end - start);
    if (*p == '\n') {
      os << '\n';
      line_offset = 0;
      new_line = true;
      ++p;
    }
    if (!*p) break;
  }
  if (!new_line)
    os << '\n';
  return os.str();
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

void OptionHelper<double>::Format(fmt::Formatter &f, Arg value) {
  char buffer[32];
  g_fmt(buffer, value);
  f("{}") << buffer;
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

SignalHandler::SignalHandler(const BasicSolver &s, Interruptible *i) {
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

bool BasicSolver::OptionNameLess::operator()(
    const char *lhs, const char *rhs) const {
  return strcasecmp(lhs, rhs) < 0;
}

char *BasicSolver::PrintOptionsAndExit(Option_Info *oi, keyword *, char *) {
  BasicSolver *solver = static_cast<BasicSolver*>(oi);
  std::string header = internal::IndentAndWordWrap(solver->GetOptionHeader());
  if (!header.empty())
    fmt::Print("{}\n") << header;
  fmt::Print("Directives:\n");
  const int DESC_INDENT = 6;
  const OptionMap &options = solver->options_;
  for (OptionMap::const_iterator i = options.begin(); i != options.end(); ++i) {
    fmt::Print("\n{}\n{}") << i->first
        << internal::IndentAndWordWrap(i->second->description(), DESC_INDENT);
  }
  exit(0);
  return 0;
}

BasicSolver::BasicSolver(
    fmt::StringRef name, fmt::StringRef long_name, long date)
: name_(name), has_errors_(false), read_flags_(0), timing_(false) {
  error_handler_ = this;
  output_handler_ = this;
  sol_handler_ = this;

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
    BasicSolver &s;
    VersionOption(BasicSolver &s) : SolverOption("version",
        "Single-word phrase:  report version details "
        "before solving the problem.", true), s(s) {}

    void Format(fmt::Formatter &f) {
      f("{}") << ((s.flags() & ASL_OI_show_version) != 0);
    }
    void Parse(const char *&) {
      s.Option_Info::flags |= ASL_OI_show_version;
    }
  };
  AddOption(OptionPtr(new VersionOption(*this)));

  struct WantSolOption : TypedSolverOption<int> {
    BasicSolver &s;
    WantSolOption(BasicSolver &s) : TypedSolverOption<int>("wantsol",
        "In a stand-alone invocation (no -AMPL on the command line), "
        "what solution information to write.  Sum of\n"
        "      1 = write .sol file\n"
        "      2 = primal variables to stdout\n"
        "      4 = dual variables to stdout\n"
        "      8 = suppress solution message\n"), s(s) {}

    int GetValue() const { return s.wantsol(); }
    void SetValue(int value) { s.Option_Info::wantsol = value; }
  };
  AddOption(OptionPtr(new WantSolOption(*this)));

  struct TimingOption : TypedSolverOption<int> {
    BasicSolver &s;
    TimingOption(BasicSolver &s) : TypedSolverOption<int>("timing",
        "0 or 1 (default 0):  Whether to display timings for the run.\n"),
        s(s) {}

    int GetValue() const { return s.timing_; }
    void SetValue(int value) { s.timing_ = value; }
  };
  AddOption(OptionPtr(new TimingOption(*this)));

  cl_option_ = keyword();
  cl_option_.name = const_cast<char*>("=");
  cl_option_.desc = const_cast<char*>("show name= possibilities");
  cl_option_.kf = BasicSolver::PrintOptionsAndExit;
  cl_option_.info = 0;
  options = &cl_option_;
  n_options = 1;
}

BasicSolver::~BasicSolver() {
  struct Deleter {
    void operator()(std::pair<const char *const, SolverOption*> &p) {
      delete p.second;
    }
  };
  std::for_each(options_.begin(), options_.end(), Deleter());
}

bool BasicSolver::ProcessArgs(char **&argv, unsigned flags) {
  char *stub = getstub_ASL(reinterpret_cast<ASL*>(problem_.asl_), &argv, this);
  if (!stub) {
    usage_noexit_ASL(this, 1);
    return false;
  }
  steady_clock::time_point start = steady_clock::now();
  problem_.Read(stub, read_flags_);
  double read_time = GetTimeAndReset(start);
  bool result = ParseOptions(argv, flags);
  if (timing_)
    Print("Input time = {:.6f}s\n") << read_time;
  return result;
}

SolverOption *BasicSolver::FindOption(const char *name) const {
  OptionMap::const_iterator i = options_.find(name);
  if (i == options_.end())
    throw OptionError(fmt::Format("Unknown option \"{}\"") << name);
  return i->second;
}

void BasicSolver::ParseOptionString(const char *s, unsigned flags) {
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
    OptionMap::iterator i = options_.find(&name[0]);
    if (i == options_.end()) {
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
    SolverOption *opt = i->second;
    if (*s == '?') {
      char next = s[1];
      if (!next || std::isspace(next)) {
        ++s;
        if ((flags & NO_OPTION_ECHO) == 0) {
          fmt::Formatter f;
          f("{}=") << &name[0];
          opt->Format(f);
          f << '\n';
          Print(fmt::StringRef(f.c_str(), f.size()));
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

bool BasicSolver::ParseOptions(char **argv, unsigned flags) {
  has_errors_ = false;
  Option_Info::flags &= ~ASL_OI_show_version;
  if (opname) {
    if (const char *s = getenv(opname))
      ParseOptionString(s, flags);
  }
  while (const char *s = *argv++)
    ParseOptionString(s, flags);
  problem_.asl_->i.need_nl_ = nnl;
  if (this->flags() & ASL_OI_show_version)
    show_version_ASL(this);
  std::fflush(stdout);
  return !has_errors_;
}

int BasicSolver::Run(char **argv) {
  if (!ProcessArgs(argv))
    return 1;
  Solve(problem());
  return 0;
}
}
