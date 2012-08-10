// Function adapters, binders, numerical differentiator and other
// function-related stuff.

#include "tests/function.h"

#include <iterator>
#include <sstream>
#include <cstring>

#include "solvers/asl.h"

namespace fun {

std::ostream &operator<<(std::ostream &os, const Tuple &t) {
  os << "(";
  if (unsigned size = t.size()) {
    os << t[0];
    for (size_t i = 1; i < size; ++i)
      os << ", " << t[i];
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

void FunctionInfo::SetArgNames(const char *arg_names) {
  std::istringstream is(arg_names);
  copy(std::istream_iterator<std::string>(is),
      std::istream_iterator<std::string>(),
      std::back_inserter<std::vector<std::string>>(arg_names_));
}

double FunctionInfo::GetDerivative(unsigned, const Tuple &) {
  return std::numeric_limits<double>::quiet_NaN();
}

double FunctionInfo::GetSecondDerivative(unsigned, unsigned, const Tuple &) {
  return std::numeric_limits<double>::quiet_NaN();
}

std::string FunctionInfo::DerivativeError(
    const Function &, unsigned, const Tuple &) {
  return "";
}

std::string FunctionInfo::Derivative2Error(const Function &, const Tuple &) {
  return "";
}

const char *Function::name() const { return fi_->name; }

Result Function::operator()(const Tuple &args,
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
  al.ra = &ra[0];
  al.nr = al.n = num_args;
  al.TMI = &tmi;
  al.AE = asl_->i.ae;
  al.dig = !dig.empty() ? &dig[0] : nullptr;
  al.funcinfo = info;

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

}
