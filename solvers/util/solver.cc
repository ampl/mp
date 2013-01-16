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
#include <cstring>

#ifndef WIN32
# include <unistd.h>
#else
# include <io.h>
#define write _write
#endif

#include "solvers/util/format.h"
#include "solvers/getstub.h"

namespace {
struct KeywordNameLess {
  bool operator()(const keyword &lhs, const keyword &rhs) const {
    return std::strcmp(lhs.name, rhs.name) < 0;
  }
};
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
}

Problem::Problem() : asl_(reinterpret_cast<ASL_fg*>(ASL_alloc(ASL_read_fg))) {}

Problem::~Problem() {
  ASL_free(reinterpret_cast<ASL**>(&asl_));
}

std::string BasicSolver::signal_message_;
const char *BasicSolver::signal_message_ptr_;
unsigned BasicSolver::signal_message_size_;
volatile std::sig_atomic_t BasicSolver::stop_;

void BasicSolver::HandleSigInt(int sig) {
  unsigned count = 0;
  do {
    // Use asynchronous-safe function write instead of printf!
    int result = write(1, signal_message_ptr_ + count,
        signal_message_size_ - count);
    if (result < 0) break;
    count += result;
  } while (count < signal_message_size_);
  if (stop_) {
    // Use asynchronous-safe function _exit instead of exit!
    _exit(1);
  }
  stop_ = 1;
  // Restore the handler since it might have been reset before the handler
  // is called (this is implementation defined).
  std::signal(sig, HandleSigInt);
}

void BasicSolver::SortOptions() {
  if (options_sorted_) return;
  std::sort(keywords_.begin(), keywords_.end(), KeywordNameLess());
  keywds = &keywords_[0];
  n_keywds = static_cast<int>(keywords_.size());
  options_sorted_ = true;
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

  version_desc_ = FormatDescription(
      "Single-word phrase:  report version details "
      "before solving the problem.");
  AddKeyword("version", version_desc_.c_str(), Ver_val, 0);
  wantsol_desc_ = FormatDescription(
      "In a stand-alone invocation (no -AMPL on the command line), "
      "what solution information to write.  Sum of\n"
      "      1 = write .sol file\n"
      "      2 = primal variables to stdout\n"
      "      4 = dual variables to stdout\n"
      "      8 = suppress solution message\n");
  AddKeyword("wantsol", wantsol_desc_.c_str(), WS_val, 0);

  signal_message_ = str(fmt::Format("\n<BREAK> ({})\n") << name.c_str());
  signal_message_ptr_ = signal_message_.c_str();
  signal_message_size_ = signal_message_.size();
  stop_ = 0;
  std::signal(SIGINT, HandleSigInt);
}

BasicSolver::~BasicSolver() {}

void BasicSolver::AddKeyword(const char *name,
    const char *description, Kwfunc func, const void *info) {
  keywords_.push_back(keyword());
  keyword &kw = keywords_.back();
  kw.name = const_cast<char*>(name);
  kw.desc = const_cast<char*>(description);
  kw.kf = func;
  kw.info = const_cast<void*>(info);
}

std::string BasicSolver::FormatDescription(const char *description) {
  std::ostringstream os;
  os << '\n';
  bool new_line = true;
  int line_offset = 0;
  int indent = 0;
  const char *s = description;
  const int MAX_LINE_LENGTH = 78;
  for (;;) {
    const char *start = s;
    while (*s == ' ')
      ++s;
    const char *word_start = s;
    while (*s != ' ' && *s != '\n' && *s)
      ++s;
    const char *word_end = s;
    if (new_line) {
      indent = 6 + static_cast<int>(word_start - start);
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
    if (*s == '\n') {
      os << '\n';
      line_offset = 0;
      new_line = true;
      ++s;
    }
    if (!*s) break;
  }
  if (!new_line)
    os << '\n';
  return os.str();
}

bool BasicSolver::ReadProblem(char **&argv) {
  SortOptions();
  ASL_fg *aslfg = problem_.asl_;
  ASL *asl = reinterpret_cast<ASL*>(aslfg);
  char *stub = getstub_ASL(asl, &argv, this);
  if (!stub) {
    usage_noexit_ASL(this, 1);
    return false;
  }
  FILE *nl = jac0dim_ASL(asl, stub, static_cast<ftnlen>(std::strlen(stub)));
  aslfg->i.Uvx_ =
      static_cast<real*>(Malloc(problem_.num_vars() * sizeof(real)));
  aslfg->i.Urhsx_ =
      static_cast<real*>(Malloc(problem_.num_cons() * sizeof(real)));
  efunc *r_ops_int[N_OPS];
  for (int i = 0; i < N_OPS; ++i)
    r_ops_int[i] = reinterpret_cast<efunc*>(i);
  aslfg->I.r_ops_ = r_ops_int;
  aslfg->p.want_derivs_ = 0;
  fg_read_ASL(asl, nl, ASL_allow_CLP);
  aslfg->I.r_ops_ = 0;
  return true;
}
}
