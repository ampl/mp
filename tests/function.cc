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

using std::string;

namespace {

void Print(std::ostream &os, double value) {
  if (!isnan(value))
    os << value;
  else
    os << "NaN";
}

fun::Library *library;
}

namespace fun {

class LibraryImpl : public AmplExports, public TMInfo {
 private:
  string name_;
  std::vector<void*> tempmem_;
  static string error_;

  typedef std::map<string, func_info> FunctionMap;
  FunctionMap funcs_;

  typedef std::map<string, Handler> HandlerMap;
  HandlerMap handlers_;

  static void ReportDuplicateFunction(const string &name) {
    error_ = "duplicate function '" + name + "'";
  }

  static void AddFunc(const char *name, rfunc f,
      int type, int nargs, void *funcinfo, AmplExports *ae);

  static void AddTableHandler(
      TableHandlerFunc read, TableHandlerFunc write,
      char *handler_info, int , void *);

  static void AtExit(AmplExports *, Exitfunc *, void *) {
    // Do nothing.
  }

  static void *Tempmem(TMInfo *tmi, size_t size) {
    LibraryImpl *impl = static_cast<LibraryImpl*>(tmi);
    impl->tempmem_.push_back(0);
    return impl->tempmem_.back() = malloc(size);
  }

 public:
  explicit LibraryImpl(const char *name);
  ~LibraryImpl() {
    std::for_each(tempmem_.begin(), tempmem_.end(), std::ptr_fun(free));
  }

  void Load() {
    error_ = string();
    i_option_ASL = name_.c_str();
    // Use funcadd(AmplExports*) instead of func_add(ASL*) because
    // the latter doesn't load random functions.
    funcadd(this);
  }

  string error() const { return error_; }

  unsigned GetNumFunctions() const { return funcs_.size(); }

  const func_info *GetFunction(const char *name) const {
    FunctionMap::const_iterator i = funcs_.find(name);
    return i != funcs_.end() ? &i->second : 0;
  }

  const Handler *GetHandler(const char *name) const {
    HandlerMap::const_iterator i = handlers_.find(name);
    return i != handlers_.end() ? &i->second : 0;
  }
};

string LibraryImpl::error_;

void LibraryImpl::AddFunc(const char *name, rfunc f,
    int type, int nargs, void *funcinfo, AmplExports *ae) {
  func_info fi = {};
  fi.name = name;
  fi.funcp = f;
  fi.ftype = type;
  fi.nargs = nargs;
  fi.funcinfo = funcinfo;
  LibraryImpl *impl = static_cast<LibraryImpl*>(ae);
  if (!impl->funcs_.insert(std::make_pair(name, fi)).second)
    ReportDuplicateFunction(name);
  note_libuse_ASL();  // Make sure the library is not unloaded.
}

void LibraryImpl::AddTableHandler(
    TableHandlerFunc read, TableHandlerFunc write,
    char *handler_info, int , void *) {
  string info(handler_info);
  string name(info.substr(0, info.find('\n')));
  Handler handler(library, read, write);
  if (!library->impl()->handlers_.insert(
      std::make_pair(name, handler)).second) {
    ReportDuplicateFunction(name);
  }
  note_libuse_ASL();  // Make sure the library is not unloaded.
}

LibraryImpl::LibraryImpl(const char *name) : AmplExports(), name_(name) {
  Addfunc = AddFunc;
  Add_table_handler = AddTableHandler;
  AmplExports::AtExit = AtExit;
  AmplExports::AtReset = AtExit;
  AmplExports::Tempmem = Tempmem;
  SprintF = sprintf;  // NOLINT(runtime/printf)
  SnprintF = snprintf;
  VsnprintF = vsnprintf;
  Fopen = fopen;
  Fclose = fclose;
  Fread = fread;
  Fseek = fseek;
  PrintF = printf;
  FprintF = fprintf;
  StdErr = stderr;
  Qsortv = qsortv;
  Getenv = getenv;
}

Library::Library(const char *name) : impl_(new LibraryImpl(name)) {
  library = this;
}

Library::~Library() {
  library = 0;
}

void Library::Load() { impl_->Load(); }

string Library::error() const { return impl_->error(); }

unsigned Library::GetNumFunctions() const {
  return impl_->GetNumFunctions();
}

const func_info *Library::GetFunction(const char *name) const {
  return impl_->GetFunction(name);
}

const Handler *Library::GetHandler(const char *name) const {
  return impl_->GetHandler(name);
}

class TableImpl : public TableInfo {
 private:
  int num_rows_;
  std::vector<char*> strings_;
  std::vector<char*> colnames_;
  std::vector<DbCol> cols_;
  std::deque<double> dvals_;
  std::deque<char*> svals_;

  void AddString(std::vector<char*> *strings, const char *str);

  int AddRows(DbCol *, long nrows) {  // NOLINT(runtime/int)
    num_rows_ += nrows;
    return 0;
  }

  static int AddRows(
      TableInfo *ti, DbCol *cols, long nrows) {  // NOLINT(runtime/int)
    return static_cast<TableImpl*>(ti)->AddRows(cols, nrows);
  }

  struct Deleter {
    void operator()(char *ptr) { delete [] ptr; }
  };

 public:
  TableImpl(const char *table_name, const char *str1,
      const char *str2, const char *str3);

  ~TableImpl();

  void AddCol(const char *name);

  int num_rows() const { return num_rows_; }

  const char *GetString(int col) const {
    return svals_[col];
  }
};

void TableImpl::AddString(std::vector<char*> *strings, const char *str) {
  if (!str) return;
  char *copy = new char[std::strlen(str) + 1];
  std::strcpy(copy, str);  // NOLINT(runtime/printf)
  strings->push_back(copy);
}

TableImpl::TableImpl(const char *table_name, const char *str1,
    const char *str2, const char *str3) : TableInfo(), num_rows_(0) {
  TableInfo::AddRows = AddRows;
  tname = const_cast<char*>(table_name);
  AddString(&strings_, str1);
  AddString(&strings_, str2);
  AddString(&strings_, str3);
  nstrings = strings_.size();
  strings = &strings_[0];
}

TableImpl::~TableImpl() {
  for_each(strings_.begin(), strings_.end(), Deleter());
  for_each(colnames_.begin(), colnames_.end(), Deleter());
}

void TableImpl::AddCol(const char *name) {
  DbCol col = {};
  svals_.push_back(0);
  dvals_.push_back(0);
  col.dval = &dvals_.back();
  col.sval = &svals_.back();
  cols_.push_back(col);
  AddString(&colnames_, name);
  ++ncols;
  colnames = &colnames_[0];
  cols = &cols_[0];
}

Table::Table(const char *table_name, const char *str1,
    const char *str2, const char *str3) :
  impl_(new TableImpl(table_name, str1, str2, str3)) {
}

int Table::num_rows() const { return impl_->num_rows(); }

const char *Table::error_message() const { return impl_->Errmsg; }

void Table::AddCol(const char *name) { return impl_->AddCol(name); }

const char *Table::GetString(int col) const { return impl_->GetString(col); }

int Handler::Read(Table *t) const {
  t->impl_->TMI = lib_->impl();
  return read_(lib_->impl(), t->impl_.get());
}

int Handler::Write(Table *t) const {
  t->impl_->TMI = lib_->impl();
  return write_(lib_->impl(), t->impl_.get());
}

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
  copy(std::istream_iterator<string>(is),
      std::istream_iterator<string>(),
      std::back_inserter< std::vector<string> >(arg_names_));
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
  al.ra = !ra.empty() ? &ra[0] : nullptr;
  al.nr = al.n = num_args;
  al.TMI = lib_->impl();
  al.AE = lib_->impl();
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
}
