#ifndef BASIC_CONSTR_H
#define BASIC_CONSTR_H

#include <array>
#include <vector>
#include <cmath>
#include <utility>

#include "mp/flat/context.h"

namespace mp {

/// Custom constraints to derive from, so that overloaded default settings work
class BasicConstraint {
public:
  /// Name for messages
  static constexpr const char* GetName()
  { return "BasicConstraint"; }
  /// Whether context is meaningful here
  static constexpr bool HasContext() { return false; }
  /// Get context, if meaningful
  Context GetContext() const { return Context::CTX_NONE; }
  /// Set context, if meaningful
  void SetContext(Context ) const { }
  /// For functional constraints, result variable index
  int GetResultVar() const { return -1; }
};

/// A special constraint 'var=...', which defines a result variable
class FunctionalConstraint : public BasicConstraint {
  int result_var_=-1;                // defined var is optional
  mutable Context ctx;               // always store context
public:
  /// Name for messages
  static constexpr const char* GetName()
  { return "FunctionalConstraint"; }
  /// Constructor
  /// @param v: result variable
  FunctionalConstraint(int v=-1) : result_var_(v) {}
  /// Basic operator==
  bool operator==(const FunctionalConstraint& dc) const {
    return result_var_==dc.result_var_;
  }
  /// Get result variable
  int GetResultVar() const { return result_var_; }
  /// Set result variable
  void SetResultVar(int v) { result_var_=v; }
  /// Whether context is meaningful
  static constexpr bool HasContext() { return true; }
  /// Get it
  Context GetContext() const { return ctx; }
  /// Set it
  void SetContext(Context c) const { ctx=c; }
  /// Add context
  void AddContext(Context c) { ctx.Add(c); }
};

/// Possible argument arrays for CustomFunctionalConstraint

/// Fixed argument array of 1 element
using VarArray1 = std::array<int, 1>;
/// Fixed argument array of 2 elements
using VarArray2 = std::array<int, 2>;
/// Fixed argument array of N elements
template <int N>
using VarArrayN = std::array<int, N>;
/// Variable-size argument array
using VarArray = std::vector<int>;

/// Possible parameter arrays

/// Fixed parameter array of N elements
template <class Num, size_t N>
  using ParamArrayN = std::array<Num, N>;
/// Empty parameter array
using ParamArray0 = ParamArrayN<int, 0>;
/// Fixed parameter array of 1 double
using DblParamArray1 = ParamArrayN<double, 1>;

/// A functional constraint with given arguments
/// and further info as parameters
/// @param Args: arguments type
/// @param Params: parameters type
/// @param NumOrLogic: base class defining a numeric or logic constraint
/// @param Id: a struct with name_
template <class Args, class Params, class NumOrLogic, class Id>
class CustomFunctionalConstraint :
  public FunctionalConstraint, public NumOrLogic, public Id {
  Args args_;
  Params params_;

public:
  /// Constraint name for messages
  static constexpr const char* GetName() { return Id::name_; }
  /// Default constructor
  CustomFunctionalConstraint() = default;
  /// Arguments typedef
  using Arguments = Args;
  /// Parameters typedef
  using Parameters = Params;
  /// Construct from arguments only
  CustomFunctionalConstraint(Arguments args) noexcept :
    args_(std::move(args)) { }
  /// Construct from arguments and parameters
  /// Might need to use explicit types when using initializer lists,
  /// in order to distinguish from the next 2 constructors
  CustomFunctionalConstraint(Arguments args, Parameters prm) noexcept :
    args_(std::move(args)), params_(std::move(prm)) { }
  /// Construct from resvar and arguments
  CustomFunctionalConstraint(int varr, Arguments args) noexcept :
     FunctionalConstraint(varr), args_(std::move(args)) { }

  /////////////////////////////////////////////////////////////////////

  /// Reuse GetResultVar()
  using FunctionalConstraint::GetResultVar;
  /// Get const Arguments&
  const Arguments& GetArguments() const { return args_; }
  /// Get Arguments&
  Arguments& GetArguments() { return args_; }
  /// Get const Parameters&
  const Parameters& GetParameters() const { return params_; }
  /// Get Parameters&
  Parameters& GetParameters() { return params_; }

  /////////////////////////////////////////////////////////////////////
  /// Specific operator==
  bool operator ==(const CustomFunctionalConstraint& mc) const {
    return this->GetResultVar()==mc.GetResultVar() &&
        this->GetArguments()==mc.GetArguments() &&
        this->GetParameters()==mc.GetParameters();
  }
};


/// A base class for numerical functional constraint.
/// It provides default properties of such a constraint
class NumericFunctionalConstraint {
public:
  /// Whether the constraint is logical
  static constexpr bool IsLogical() { return false; }
  /// Apriori bounds on the result
  static constexpr std::pair<double, double>
  GetAprioriBounds() { return {-INFINITY, INFINITY}; }
};


/// A base class for logical functional constraint.
/// It provides default properties of such a constraint
class LogicalFunctionalConstraint {
public:
  /// Whether the constraint is logical
  static constexpr bool IsLogical() { return true; }
  /// Apriori bounds on the result
  static constexpr std::pair<double, double>
  GetAprioriBounds() { return {0.0, 1.0}; }
};


////////////////////////////////////////////////////////////////////////
/// Args is the argument type, e.g., array of variables, or an expression
/// Params is the parameter type, e.g., array of numbers. Can be empty
#define DEF_CUSTOM_FUNC_CONSTR_WITH_PRM(Name, Args, Params, NumLogic, Descr) \
struct Name ## Id { \
  static constexpr auto description_ = Descr; \
  static constexpr auto name_        = #Name; \
}; \
using Name = CustomFunctionalConstraint<Args, Params, NumLogic, Name ## Id>

/// Custom numeric constraint without fixed parameters
#define DEF_NUMERIC_FUNC_CONSTR(Name, Args, Descr) \
    DEF_NUMERIC_FUNC_CONSTR_WITH_PRM(Name, Args, ParamArray0, Descr)
/// Custom logical constraint without fixed parameters
#define DEF_LOGICAL_FUNC_CONSTR(Name, Args, Descr) \
    DEF_LOGICAL_FUNC_CONSTR_WITH_PRM(Name, Args, ParamArray0, Descr)

/// Custom numeric constraint with parameter data
#define DEF_NUMERIC_FUNC_CONSTR_WITH_PRM(Name, Args, Params, Descr) \
  DEF_CUSTOM_FUNC_CONSTR_WITH_PRM(Name, Args, Params, \
    NumericFunctionalConstraint, Descr)
/// Custom logical constraint with parameter data
#define DEF_LOGICAL_FUNC_CONSTR_WITH_PRM(Name, Args, Params, Descr) \
  DEF_CUSTOM_FUNC_CONSTR_WITH_PRM(Name, Args, Params, \
    NumericFunctionalConstraint, Descr)

} // namespace mp

#endif // BASIC_CONSTR_H
