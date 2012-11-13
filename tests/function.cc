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
#include <cassert>
#include <cstring>

#include "solvers/asl.h"
#include "tests/config.h"
#undef VOID

using std::string;
using fun::Table;

namespace {

class ScopedTableInfo : public TableInfo {
 private:
  std::vector<char*> strings_;
  std::vector<std::string> colnames_;
  std::vector<char*> colnameptrs_;
  std::vector<DbCol> cols_;
  std::vector<double> dvals_;
  std::vector<char*> svals_;
  Table *table_;
  static char MISSING;

  void SetString(std::vector<char*> *strings, unsigned index, const char *str);
  void AddString(std::vector<char*> *strings, const char *str);

  struct Deleter {
    void operator()(char *ptr) { delete [] ptr; }
  };

 public:
  ScopedTableInfo(const Table &t, const std::string &connection_str,
      const std::string &sql = std::string());
  ~ScopedTableInfo();

  Table *GetTable() { return table_; }
  void SetTable(Table *t) { table_ = t; }

  void SetValue(unsigned row, unsigned col, const fun::Variant &v) {
    if (v.type() == fun::STRING)
      cols_[col].sval[row] = const_cast<char*>(v.string());
    else
      cols_[col].dval[row] = v.number();
  }
};

char ScopedTableInfo::MISSING;

void ScopedTableInfo::SetString(
    std::vector<char*> *strings, unsigned index, const char *str) {
  char *&oldstr = (*strings)[index];
  if (oldstr) delete [] oldstr;
  oldstr = new char[std::strlen(str) + 1];
  std::strcpy(oldstr, str);  // NOLINT(runtime/printf)
}

void ScopedTableInfo::AddString(std::vector<char*> *strings, const char *str) {
  if (!str) return;
  strings->push_back(0);
  SetString(strings, strings->size() - 1, str);
}

ScopedTableInfo::ScopedTableInfo(const Table &t,
    const std::string &connection_str, const std::string &sql) {
  // Workaround for GCC bug 30111 that prevents value-initialization of
  // the base POD class.
  TableInfo ti = {};
  static_cast<TableInfo&>(*this) = ti;

  tname = const_cast<char*>(t.name());
  AddString(&strings_, "ODBC");
  AddString(&strings_, connection_str.c_str());
  if (!sql.empty())
    AddString(&strings_, sql.c_str());
  nstrings = strings_.size();
  strings = &strings_[0];
  Missing = &MISSING;

  unsigned num_rows = std::max(t.num_rows(), 1u);
  unsigned num_values = num_rows * t.num_cols();
  svals_.resize(num_values);
  dvals_.resize(num_values);

  nrows = maxrows = t.num_rows();
  arity = 1;
  ncols = t.num_cols() - arity;
  if (t.num_cols() != 0) {
    colnames_.resize(t.num_cols());
    colnameptrs_.resize(t.num_cols());
    colnames = &colnameptrs_[0];
    cols_.reserve(t.num_cols());
    for (unsigned i = 0; i < t.num_cols(); ++i) {
      colnames_[i] = t.GetColName(i);
      colnameptrs_[i] = const_cast<char*>(colnames_[i].c_str());
      DbCol col = {};
      col.dval = &dvals_[i * num_rows];
      col.sval = &svals_[i * num_rows];
      cols_.push_back(col);
    }
    cols = &cols_[0];
  }
}

ScopedTableInfo::~ScopedTableInfo() {
  for_each(strings_.begin(), strings_.end(), Deleter());
}

void CheckResult(int result, const ScopedTableInfo &ti) {
  if (ti.Errmsg)
    throw std::runtime_error(ti.Errmsg);
  if (result != DB_Done)
    throw std::runtime_error("DbRead failed");
}

void Print(std::ostream &os, double value) {
  if (!isnan(value))
    os << value;
  else
    os << "NaN";
}

fun::Library *library;
}

namespace fun {

bool operator==(const Variant &lhs, const Variant &rhs) {
  Type type = lhs.type();
  if (type != rhs.type())
    return false;
  switch (type) {
  case VOID:
    return true;
  case DOUBLE:
    return lhs.number() == rhs.number();
  case STRING:
    return std::strcmp(lhs.string(), rhs.string()) == 0;
  case POINTER:
    return lhs.pointer() == rhs.pointer();
  case INT:
  case UINT:
    // Not supported.
    break;
  }
  throw std::runtime_error("unknown type");
}

std::ostream &operator<<(std::ostream &os, const Variant &v) {
  switch (v.type()) {
  case POINTER:
    os << v.pointer();
    break;
  case STRING:
    os << '"' << v.string() << '"';
    break;
  default:
    os << v.number();
    break;
  }
  return os;
}

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
  Strtod = strtod;
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

int Table::AddRows(DbCol *cols, long nrows) {  // NOLINT(runtime/int)
  for (long i = 0; i < nrows; ++i) {
    for (unsigned j = 0; j < num_cols(); ++j) {
      DbCol &value = cols[j];
      if (value.sval[i])
        Add(value.sval[i]);
      else
        Add(value.dval[i]);
    }
  }
  return 0;
}

int Table::AddRows(
    TableInfo *ti, DbCol *cols, long nrows) {  // NOLINT(runtime/int)
  return static_cast<ScopedTableInfo*>(ti)->GetTable()->AddRows(cols, nrows);
}

bool operator==(const Table &lhs, const Table &rhs) {
  unsigned num_rows = lhs.num_rows();
  if (num_rows != rhs.num_rows())
    return false;
  unsigned num_cols = lhs.num_cols();
  if (num_cols != rhs.num_cols())
    return false;
  if (lhs.HasColNames() != rhs.HasColNames())
    return false;
  if (!lhs.HasColNames())
    return true;
  for (unsigned j = 0; j != num_cols; ++j) {
    if (std::strcmp(lhs.GetColName(j), rhs.GetColName(j)) != 0)
      return false;
  }
  for (unsigned i = 0; i != num_rows; ++i) {
    for (unsigned j = 0; j != num_cols; ++j) {
      if (lhs(i, j) != rhs(i, j))
        return false;
    }
  }
  return true;
}

std::ostream &operator<<(std::ostream &os, const Table &t) {
  unsigned num_rows = t.num_rows();
  unsigned num_cols = t.num_cols();
  for (unsigned j = 0; j != num_cols; ++j)
    os << t.GetColName(j) << " ";
  os << "\n";
  for (unsigned i = 0; i != num_rows; ++i) {
    for (unsigned j = 0; j != num_cols; ++j)
      os << t(i, j) << " ";
    os << "\n";
  }
  return os;
}

void Handler::Read(const std::string &connection_str,
    Table *t, const std::string &sql_statement) const {
  ScopedTableInfo table(*t, connection_str, sql_statement);
  table.TMI = lib_->impl();
  table.SetTable(t);
  table.AddRows = Table::AddRows;
  CheckResult(read_(lib_->impl(), &table), table);
}

void Handler::Write(
    const std::string &connection_str, const Table &t) const {
  ScopedTableInfo table(t, connection_str);
  for (unsigned i = 0, m = t.num_rows(); i < m; ++i) {
    for (unsigned j = 0, n = t.num_cols(); j < n; ++j)
      table.SetValue(i, j, t(i, j));
  }
  table.TMI = lib_->impl();
  CheckResult(write_(lib_->impl(), &table), table);
}

const Type GetType<void>::VALUE = VOID;
const Type GetType<int>::VALUE = INT;
const Type GetType<unsigned>::VALUE = UINT;
const Type GetType<double>::VALUE = DOUBLE;
const Type GetType<const char*>::VALUE = STRING;

std::ostream &operator<<(std::ostream &os, const Tuple &t) {
  os << "(";
  if (unsigned size = t.size()) {
    Print(os, t[0].number());
    for (size_t i = 1; i < size; ++i) {
      os << ", ";
      Print(os, t[i].number());
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
    ra[i] = args[i].number();
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
  args_[eval_arg_] = Variant::FromDouble(x);
  Function::Result r = f_(args_, DERIVS, use_deriv_);
  return r.error() ?
      std::numeric_limits<double>::quiet_NaN() : r.deriv(deriv_arg_);
}
}
