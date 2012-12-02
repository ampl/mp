/*
 Utilities for writing AMPL solver drivers.

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

#include "solvers/util/driver.h"

#include <cstdarg>
#include <cstdio>
#include <cstring>

#include "solvers/getstub.h"

namespace {
struct KeywordNameLess {
  bool operator()(const keyword &lhs, const keyword &rhs) const {
    return std::strcmp(lhs.name, rhs.name) < 0;
  }
};
}

namespace ampl {

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

void BaseOptionInfo::Sort() {
  if (sorted_) return;
  std::sort(keywords_.begin(), keywords_.end(), KeywordNameLess());
  keywds = &keywords_[0];
  n_keywds = keywords_.size();
  sorted_ = true;
}

BaseOptionInfo::BaseOptionInfo() : Option_Info(), sorted_(false) {
  // TODO: align text
  AddKeyword("version",
      "Single-word phrase:  report version details\n"
      "before solving the problem.\n", Ver_val, 0);
  AddKeyword("wantsol",
      "In a stand-alone invocation (no -AMPL on the\n"
      "command line), what solution information to\n"
      "write.  Sum of\n"
      "      1 = write .sol file\n"
      "      2 = primal variables to stdout\n"
      "      4 = dual variables to stdout\n"
      "      8 = suppress solution message\n", WS_val, 0);
}

void BaseOptionInfo::AddKeyword(const char *name,
    const char *description, Kwfunc func, const void *info) {
  keywords_.push_back(keyword());
  keyword &kw = keywords_.back();
  kw.name = const_cast<char*>(name);
  kw.desc = const_cast<char*>(description);
  kw.kf = func;
  kw.info = const_cast<void*>(info);
}

Problem::Problem() : asl_(reinterpret_cast<ASL_fg*>(ASL_alloc(ASL_read_fg))) {}

Problem::~Problem() {
  ASL_free(reinterpret_cast<ASL**>(&asl_));
}

bool Problem::Read(char **&argv, BaseOptionInfo &oi) {
  oi.Sort();
  ASL *asl = reinterpret_cast<ASL*>(asl_);
  char *stub = getstub_ASL(asl, &argv, &oi);
  if (!stub) {
    usage_noexit_ASL(&oi, 1);
    return false;
  }
  FILE *nl = jac0dim_ASL(asl, stub, static_cast<ftnlen>(std::strlen(stub)));
  asl_->i.Uvx_ = static_cast<real*>(Malloc(num_vars() * sizeof(real)));
  asl_->i.Urhsx_ = static_cast<real*>(Malloc(num_cons() * sizeof(real)));
  efunc *r_ops_int[N_OPS];
  for (int i = 0; i < N_OPS; ++i)
    r_ops_int[i] = reinterpret_cast<efunc*>(i);
  asl_->I.r_ops_ = r_ops_int;
  asl_->p.want_derivs_ = 0;
  fg_read_ASL(asl, nl, ASL_allow_CLP);
  asl_->I.r_ops_ = 0;
  return true;
}

void Driver::ReportError(const char *format, ...) {
  has_errors_ = true;
  std::vector<char> message(500);
  std::va_list args;
  for (;;) {
    va_start(args, format);
    int n = vsnprintf(&message[0], message.size(), format, args);
    va_end(args);
    if (n >= 0 && static_cast<unsigned>(n) < message.size())
      break;
    message.resize(n >= 0 ? n + 1 : 2 * message.size());
  }
  std::fputs(&message[0], stderr);
  std::fputc('\n', stderr);
}

bool Driver::GetOptions(char **argv, BaseOptionInfo &oi) {
  has_errors_ = false;
  oi.Sort();
  return getopts_ASL(reinterpret_cast<ASL*>(problem_.asl_), argv, &oi) == 0 &&
      !has_errors_;
}
}
