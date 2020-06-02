#ifndef BASIC_CONSTR_H
#define BASIC_CONSTR_H

#include <array>
#include <vector>

#include "mp/convert/context.h"

namespace mp {

/// Custom constraints to derive from, so that overloaded default settings work
class BasicConstraint {
public:
  int GetResultVar() const { return -1; }
};

/// A constraint extension which defines a variable
class DefiningConstraint : public BasicConstraint {
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

/// Possible argument arrays for CustomDefiningConstraint
using VarArray1 = std::array<int, 1>;
using VarArray2 = std::array<int, 2>;
template <int N>
using VarArrayN = std::array<int, N>;
using VarArray = std::vector<int>;

/// A defining constraint with the arguments and further info as parameters
template <class Args, class Id>
class CustomDefiningConstraint :
  public DefiningConstraint, public Id {
  Args args_;
public:
  CustomDefiningConstraint() { }
  using Arguments = Args;
  CustomDefiningConstraint(const Arguments& args) : args_(args) { }
  CustomDefiningConstraint(Arguments&& args) : args_(std::move(args)) { }
  CustomDefiningConstraint(int varr, const Arguments& args) :
     DefiningConstraint(varr), args_(args) { }
  CustomDefiningConstraint(int varr, Arguments&& args) :
     DefiningConstraint(varr), args_(std::move(args)) { }
  using DefiningConstraint::GetResultVar;
  const Arguments& GetArguments() const { return args_; }
  Arguments& GetArguments() { return args_; }
  bool operator ==(const CustomDefiningConstraint& mc) const {
    return this->GetResultVar()==mc.GetResultVar() &&
        this->GetArguments()==mc.GetArguments();
  }
};

////////////////////////////////////////////////////////////////////////
#define DEFINE_CUSTOM_DEFINING_CONSTRAINT(Name, Args, Descr) \
struct Name ## Id { \
  static constexpr auto description_ = Descr; \
  static constexpr auto name_        = #Name; \
}; \
using Name = \
   CustomDefiningConstraint<Args, Name ## Id>


} // namespace mp

#endif // BASIC_CONSTR_H
