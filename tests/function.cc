/*
 AMPL function testing infrastructure.

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

#include "tests/function.h"

#include <iterator>
#include <map>
#include <sstream>
#include <cstring>

#include "solvers/asl.h"
#include "tests/config.h"
#undef VOID

static void Print(std::ostream &os, double value) {
  if (!isnan(value))
    os << value;
  else os << "NaN";
}

namespace fun {

const Type GetType<void>::VALUE = VOID;
const Type GetType<int>::VALUE = INT;
const Type GetType<unsigned>::VALUE = UINT;
const Type GetType<double>::VALUE = DOUBLE;

std::ostream &operator<<(std::ostream &os, const Tuple &t) {
  os << "(";
  if (unsigned size = t.size()) {
    Print(os, t[0]);
    for (size_t i = 1; i < size; ++i) {
      os << ", ";
      Print(os, t[i]);
    }
  }
  os << ")";
  return os;
}

BitSet::BitSet(const char *s) {
  if (!s) return;
  unsigned num_args = std::strlen(s);
  store_.resize(num_args);
  for (unsigned i = 0; i < num_args; ++i) {
    char c = s[i];
    if (c == '0')
      store_[i] = false;
    else if (c == '1')
      store_[i] = true;
    else
      throw std::invalid_argument("invalid argument to BitSet");
  }
}

FunctionInfo::~FunctionInfo() {}

FunctionInfo &FunctionInfo::SetArgNames(const char *arg_names) {
  arg_names_.clear();
  std::istringstream is(arg_names);
  copy(std::istream_iterator<std::string>(is),
      std::istream_iterator<std::string>(),
      std::back_inserter< std::vector<std::string> >(arg_names_));
  return *this;
}

FunctionInfo::Result FunctionInfo::GetDerivative(
    const Function &, unsigned, const Tuple &) const {
  return Result();
}

FunctionInfo::Result FunctionInfo::GetSecondDerivative(
    const Function &, unsigned, unsigned, const Tuple &) const {
  return Result();
}

const char *Function::name() const { return fi_->name; }

int Function::nargs() const { return fi_->nargs; }

int Function::ftype() const { return fi_->ftype; }

Function::Result Function::operator()(const Tuple &args,
    int flags, const BitSet &use_deriv, void *info) const {
  unsigned num_args = args.size();
  if (fi_->nargs != static_cast<int>(num_args))
    throw std::invalid_argument("invalid number of arguments in function call");

  // Initialize the argument list.
  std::vector<double> ra(num_args);
  for (unsigned i = 0; i < num_args; ++i)
    ra[i] = args[i];
  std::vector<char> dig(use_deriv.size());
  if (!dig.empty()) {
    if (dig.size() != num_args)
      throw std::invalid_argument("invalid size of use_deriv");
    for (unsigned i = 0; i < num_args; ++i)
      dig[i] = !use_deriv[i];
  }
  arglist al = {};
  TMInfo tmi = {};
  al.ra = !ra.empty() ? &ra[0] : nullptr;
  al.nr = al.n = num_args;
  al.TMI = &tmi;
  al.AE = ae_;
  al.dig = !dig.empty() ? &dig[0] : nullptr;
  al.funcinfo = info ? info : fi_->funcinfo;

  // Allocate storage for the derivatives if needed.
  std::vector<double> derivs, hes;
  if ((flags & DERIVS) != 0) {
    derivs.resize(al.n);
    al.derivs = &derivs[0];
  }
  if ((flags & HES) == HES) {
    hes.resize(al.n * (al.n + 1) / 2);
    al.hes = &hes[0];
  }

  // Call the function and return the result.
  double value = fi_->funcp(&al);
  return Result(value, derivs, hes, al.Errmsg);
}

DerivativeBinder::DerivativeBinder(
    Function f, unsigned deriv_arg, unsigned eval_arg, const Tuple &args)
: f_(f), deriv_arg_(deriv_arg), eval_arg_(eval_arg),
  args_(args), use_deriv_(args.size(), false) {
  unsigned num_args = args_.size();
  if (deriv_arg >= num_args || eval_arg >= num_args)
    throw std::out_of_range("argument index is out of range");
  use_deriv_[deriv_arg] = true;
}

double DerivativeBinder::operator()(double x) {
  args_[eval_arg_] = x;
  Function::Result r = f_(args_, DERIVS, use_deriv_);
  return r.error() ?
      std::numeric_limits<double>::quiet_NaN() : r.deriv(deriv_arg_);
}

class LibraryImpl : AmplExports {
 private:
  std::string name_;
  std::vector<void*> tempmem_;

  typedef std::map<std::string, func_info> FunctionMap;
  FunctionMap funcs_;

  static void AddFunc_(const char *name, rfunc f,
      int type, int nargs, void *funcinfo, AmplExports *ae);

  static void AddTableHandler(
      int (*DbRead)(AmplExports *ae, TableInfo *TI),
      int (*DbWrite)(AmplExports *ae, TableInfo *TI),
      char *handler_info, int flags, void *Vinfo) {
    // Do nothing.
  }

  static void AtExit_(AmplExports *, Exitfunc *, void *) {
    // Do nothing.
  }

  struct CustomTMInfo : TMInfo {
    LibraryImpl *impl;
  };

  static void *Tempmem_(TMInfo *tmi, size_t size) {
    LibraryImpl *impl = static_cast<CustomTMInfo*>(tmi)->impl;
    impl->tempmem_.push_back(0);
    return impl->tempmem_.back() = malloc(size);
  }

 public:
  LibraryImpl(const char *name);
  ~LibraryImpl() {
    std::for_each(tempmem_.begin(), tempmem_.end(), std::ptr_fun(free));
  }

  void Load() {
    i_option_ASL = name_.c_str();
    // Use funcadd(AmplExports*) instead of func_add(ASL*) because
    // the latter doesn't load random functions.
    funcadd(this);
  }

  const func_info *GetFunction(const char *name) const {
    FunctionMap::const_iterator i = funcs_.find(name);
    return i != funcs_.end() ? &i->second : 0;
  }
};

void LibraryImpl::AddFunc_(const char *name, rfunc f,
    int type, int nargs, void *funcinfo, AmplExports *ae) {
  func_info fi = {};
  fi.name = name;
  fi.funcp = f;
  fi.ftype = type;
  fi.nargs = nargs;
  fi.funcinfo = funcinfo;
  static_cast<LibraryImpl*>(ae)->funcs_[name] = fi;
  note_libuse_ASL();
}

LibraryImpl::LibraryImpl(const char *name) : AmplExports(), name_(name) {
  Addfunc = AddFunc_;
  Add_table_handler = AddTableHandler;
  AtExit = AtExit_;
  AtReset = AtExit_;
  Tempmem = Tempmem_;
  SnprintF = snprintf;
  VsnprintF = vsnprintf;
  Fopen = fopen;
  Fclose = fclose;
  Fread = fread;
  Fseek = fseek;
  FprintF = fprintf;
  StdErr = stderr;
}

Library::Library(const char *name) : impl_(new LibraryImpl(name)) {}

void Library::Load() {
  impl_->Load();
}

const func_info *Library::GetFunction(const char *name) const {
  return impl_->GetFunction(name);
}
}
