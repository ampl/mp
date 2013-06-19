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

#include <cctype>
#include <cstdarg>
#include <cstdio>
#include <cstring>

#ifndef _WIN32
# include <unistd.h>
# define AMPL_WRITE write
#else
# include <io.h>
# define AMPL_WRITE _write
#endif

#include "solvers/util/format.h"
#include "solvers/getstub.h"

namespace {

struct KeywordNameLess {
  bool operator()(const keyword &lhs, const keyword &rhs) const {
    return std::strcmp(lhs.name, rhs.name) < 0;
  }
};

// Skips leading spaces in string s and returns true if the string contains
// non-space characters; false otherwise.
bool SkipSpaces(const char *&s) {
  for (; *s; ++s) {
    if (!std::isspace(*s))
      return true;
  }
  return false;
}
}

namespace ampl {

namespace internal {
int OptionParser<int>::operator()(Option_Info *oi, keyword *kw, char *&s) {
  keyword thiskw(*kw);
  int value = 0;
  thiskw.info = &value;
  s = I_val(oi, &thiskw, s);
  return value;
}

double OptionParser<double>::operator()(
    Option_Info *oi, keyword *kw, char *&s) {
  keyword thiskw(*kw);
  double value = 0;
  thiskw.info = &value;
  s = D_val(oi, &thiskw, s);
  return value;
}

const char* OptionParser<const char*>::operator()(
    Option_Info *, keyword *, char *&s) {
  char *end = s;
  while (*end && !isspace(*end))
    ++end;
  value_.assign(s, end - s);
  s = end;
  return value_.c_str();
}

std::string Format(fmt::StringRef s, int indent) {
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
  stop_ = 0;
  std::signal(SIGINT, HandleSigInt);
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

void BasicSolver::SortOptions() {
  if (options_sorted_) return;
  std::sort(cl_options_.begin(), cl_options_.end(), KeywordNameLess());
  options = &cl_options_[0];
  n_options = static_cast<int>(cl_options_.size());
  options_sorted_ = true;
}

char *BasicSolver::PrintOptionsAndExit(Option_Info *oi, keyword *, char *) {
  std::string header =
      internal::Format(static_cast<BasicSolver*>(oi)->GetOptionHeader());
  if (!header.empty())
    fmt::Print("{}\n") << header;
  fmt::Print("Directives:\n");
  const int DESC_INDENT = 6;
  const OptionMap &options = static_cast<BasicSolver*>(oi)->options_;
  for (OptionMap::const_iterator i = options.begin(); i != options.end(); ++i) {
    fmt::Print("\n{}\n{}") << i->first
        << internal::Format(i->second->description(), DESC_INDENT);
  }
  exit(0);
  return 0;
}

BasicSolver::BasicSolver(
    fmt::StringRef name, fmt::StringRef long_name, long date)
: name_(name), options_sorted_(false), has_errors_(false) {
  error_handler_ = this;
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

  AddBasicOption("version",
      "Single-word phrase:  report version details "
      "before solving the problem.", true, &BasicSolver::PrintVersion);
  AddBasicOption("wantsol",
      "In a stand-alone invocation (no -AMPL on the command line), "
      "what solution information to write.  Sum of\n"
      "      1 = write .sol file\n"
      "      2 = primal variables to stdout\n"
      "      4 = dual variables to stdout\n"
      "      8 = suppress solution message\n", false, &BasicSolver::SetWantSol);

  cl_options_.push_back(keyword());
  keyword &kw = cl_options_.back();
  kw.name = const_cast<char*>("=");
  kw.desc = const_cast<char*>("show name= possibilities");
  kw.kf = BasicSolver::PrintOptionsAndExit;
  kw.info = 0;
}

BasicSolver::~BasicSolver() {
  struct Deleter {
    void operator()(std::pair<const std::string, SolverOption*> &p) {
      delete p.second;
    }
  };
  std::for_each(options_.begin(), options_.end(), Deleter());
}

bool BasicSolver::ProcessArgs(char **&argv, unsigned flags) {
  SortOptions();
  char *stub = getstub_ASL(reinterpret_cast<ASL*>(problem_.asl_), &argv, this);
  if (!stub) {
    usage_noexit_ASL(this, 1);
    return false;
  }
  problem_.Read(stub);
  return ParseOptions(argv, flags);
}

void BasicSolver::ParseOptionString(const char *s, unsigned flags) {
  bool skip = false;
  for (;;) {
    if (!SkipSpaces(s))
      return;

    // Parse the option name.
    const char *name_start = s;
    while (*s && !std::isspace(*s) && *s != '=')
      ++s;
    std::string name;
    std::size_t name_size = s - name_start;
    name.resize(name_size);
    for (std::size_t i = 0; i < name_size; ++i)
      name[i] = std::tolower(name_start[i]);

    // Parse the option value.
    bool equal_sign = false;
    SkipSpaces(s);
    if (*s == '=') {
      ++s;
      SkipSpaces(s);
      equal_sign = true;
    }

    nnl = 0;
    OptionMap::iterator i = options_.find(name);
    if (i == options_.end()) {
      if (!skip)
        ReportError("Unknown option \"{}\"") << name;
      if (equal_sign) {
        skip = false;
        while (*s && !std::isspace(*s))
          ++s;
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
    if (i->second->is_keyword() && equal_sign) {
      ReportError("Option \"{}\" doesn't accept arguments") << name;
      while (*s && !std::isspace(*s))
        ++s;
    }
    // TODO: get rid of keyword
    keyword kw = {const_cast<char*>(i->first.c_str())};
    if (!i->second->Handle(*this, &kw, const_cast<char*&>(s)))
      has_errors_ = true;
    if ((flags & NO_OPTION_ECHO) == 0)
      printf("%.*s\n", static_cast<int>(s - name_start), name_start);
  }
}

bool BasicSolver::ParseOptions(char **argv, unsigned flags) {
  has_errors_ = false;
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
  double start_time = xectim_();
  if (!ProcessArgs(argv))
    return 1;

  // Reset is used to reset read_time_ even in case of exceptions.
  // Otherwise the read time from Run may affect the time reported in
  // a subsequent Solve:
  //   solver.Run(...);
  //   solver.Solve(...); // Doesn't read anything, but reports previous
  //                      // read time.
  class Reset {
   private:
    double &value_;
   public:
    Reset(double &value) : value_(value) {}
    ~Reset() { value_ = 0; }
  };
  Reset reset(read_time_ = xectim_() - start_time);
  Solve(problem());
  return 0;
}
}
