#ifndef BASIC_CONSTR_H
#define BASIC_CONSTR_H

#include <vector>

#include "mp/convert/context.h"

namespace mp {

/// Custom constraints to derive from, so that overloaded default settings work
class BasicConstraint {
public:
  int GetResultVar() const { return -1; }
};

/// A constraint whose arguments are an array of variables
template <class Args = std::vector<int> >
class VarArrayArgConstraint : public BasicConstraint {
  Args args_;
public:
  using Arguments = Args;
  VarArrayArgConstraint() { }
  VarArrayArgConstraint(Arguments&& aa) : args_(std::move(aa)) {}
  VarArrayArgConstraint(const Arguments& aa) : args_(aa) {}
  bool operator==(const VarArrayArgConstraint& vaac) const {
    return args_==vaac.args_;
  }
  const Arguments& GetArguments() const { return args_; }
  Arguments& GetArguments() { return args_; }
};

using VarArray2ArgConstraint = VarArrayArgConstraint< std::array<int, 2> >;

/// A constraint extension which defines a variable
class DefiningConstraint {
  int result_var_=-1;                // defined var can be optional
  Context ctx;
public:
  DefiningConstraint(int v=-1) : result_var_(v) {}
  bool operator==(const DefiningConstraint& dc) {
    return result_var_==dc.result_var_;
  }
  void SetResultVar(int v) { result_var_=v; }
  int GetResultVar() const { return result_var_; }
  void SetContext(Context c) { ctx=c; }
  void AddContext(Context c) { ctx.Add(c); }
  Context GetContext() const { return ctx; }
};

/// A defining constraint with the arguments and further info as parameters
template <class Args, class Id>
class CustomDefiningConstraint :
  public DefiningConstraint, public Args, public Id {
public:
  CustomDefiningConstraint() { }
  using Arguments = typename Args::Arguments;
  using DefiningConstraint::GetResultVar;
  CustomDefiningConstraint(int varr, const Arguments& args) :
     DefiningConstraint(varr), Args(args) { }
  bool operator ==(const CustomDefiningConstraint& mc) const {
    return this->GetResultVar()==mc.GetResultVar() &&
        this->GetArguments()==mc.GetArguments();
  }
};

} // namespace mp

#endif // BASIC_CONSTR_H
