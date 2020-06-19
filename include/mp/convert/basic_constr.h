#ifndef BASIC_CONSTR_H
#define BASIC_CONSTR_H

#include <array>
#include <vector>

#include "mp/convert/context.h"

namespace mp {

/// Custom constraints to derive from, so that overloaded default settings work
class BasicConstraint {
public:
  static const char* GetConstraintName() { return "BasicConstraint"; }
  int GetResultVar() const { return -1; }
};

/// A constraint extension which defines a variable
class DefiningConstraint : public BasicConstraint {
  int result_var_=-1;                // defined var can be optional
  Context ctx;
public:
  static const char* GetConstraintName() { return "DefiningConstraint"; }
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

/// Possible parameter arrays
template <class Num, size_t N>
  using ParamArrayN = std::array<Num, N>;
using ParamArray0 = ParamArrayN<int, 0>;
using DblParamArray1 = ParamArrayN<double, 1>;

/// A defining constraint with the arguments and further info as parameters
template <class Args, class Params, class Id>
class CustomDefiningConstraint :
  public DefiningConstraint, public Id {
  Args args_;
  Params params_;
public:
  static const char* GetConstraintName() { return Id::name_; }
  CustomDefiningConstraint() { }
  using Arguments = Args;
  using Parameters = Params;
  /// Construct from arguments only
  CustomDefiningConstraint(const Arguments& args) : args_(args) { }
  CustomDefiningConstraint(Arguments&& args) : args_(std::move(args)) { }
  /// Construct from arguments and parameters
  /// Might need to use explicit types when using initializer lists,
  /// in order to distinguish from the next 2 constructors
  CustomDefiningConstraint(const Arguments& args, const Parameters& prm) :
    args_(args), params_(prm) { }
  CustomDefiningConstraint(Arguments&& args, Parameters&& prm) :
    args_(std::move(args)), params_(std::move(prm)) { }
  /// From resvar and arguments
  CustomDefiningConstraint(int varr, const Arguments& args) :
     DefiningConstraint(varr), args_(args) { }
  CustomDefiningConstraint(int varr, Arguments&& args) :
     DefiningConstraint(varr), args_(std::move(args)) { }
  /////////////////////////////////////////////////////////////////////
  using DefiningConstraint::GetResultVar;
  const Arguments& GetArguments() const { return args_; }
  Arguments& GetArguments() { return args_; }
  const Parameters& GetParameters() const { return params_; }
  Parameters& GetParameters() { return params_; }
  /////////////////////////////////////////////////////////////////////
  bool operator ==(const CustomDefiningConstraint& mc) const {
    return this->GetResultVar()==mc.GetResultVar() &&
        this->GetArguments()==mc.GetArguments() &&
        this->GetParameters()==mc.GetParameters();
  }
};

////////////////////////////////////////////////////////////////////////
#define DEFINE_CUSTOM_DEFINING_CONSTRAINT_WITH_PARAMS(Name, Args, Params, Descr) \
struct Name ## Id { \
  static constexpr auto description_ = Descr; \
  static constexpr auto name_        = #Name; \
}; \
using Name = CustomDefiningConstraint<Args, Params, Name ## Id>

#define DEFINE_CUSTOM_DEFINING_CONSTRAINT(Name, Args, Descr) \
    DEFINE_CUSTOM_DEFINING_CONSTRAINT_WITH_PARAMS(Name, Args, ParamArray0, Descr)

} // namespace mp

#endif // BASIC_CONSTR_H
