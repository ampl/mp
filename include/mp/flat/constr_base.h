#ifndef BASIC_CONSTR_H
#define BASIC_CONSTR_H

#include <array>
#include <vector>
#include <string>
#include <cmath>
#include <utility>
#include <typeinfo>

#include "mp/flat/context.h"

namespace mp {

/// Custom constraints to derive from, so that overloaded default settings work
class BasicConstraint {
public:
  /// Constraint type name for messages
  static constexpr const char* GetTypeName()
  { return "BasicConstraint"; }
  /// Constraint name
  const char* GetName() const { return name_.c_str(); }
  /// Constraint name
  const char* name() const { return GetName(); }
  /// Set constraint name
  void SetName(std::string nm) { name_ = std::move(nm); }
  /// Whether context is meaningful here
  static constexpr bool UsesContext() { return false; }
  /// Get context, if meaningful
  Context GetContext() const { return Context::CTX_NONE; }
  /// Set context, if meaningful
  void SetContext(Context ) const { }
  /// Add (merge) context, if meaningful
  void AddContext(Context ) const { }
  /// For functional constraints, result variable index
  int GetResultVar() const { return -1; }


private:
  std::string name_;
};


/// A special constraint 'var=...', which defines a result variable
class FunctionalConstraint : public BasicConstraint {
  int result_var_=-1;                // defined var is optional
  mutable Context ctx;               // always store context
public:
  /// Constraint type name for messages
  static constexpr const char* GetTypeName()
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
  static constexpr bool UsesContext() { return true; }
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
/// @param Id: a struct with GetTypeName()
template <class Args, class Params, class NumOrLogic, class Id>
class CustomFunctionalConstraint :
  public FunctionalConstraint, public NumOrLogic, public Id {
  Args args_;
  Params params_;

public:
  /// Constraint type name for messages
  static const char* GetTypeName() { return Id::GetTypeName(); }
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
};


/// A base class for numerical functional constraint.
/// It provides default properties of such a constraint
class NumericFunctionalConstraintTraits {
public:
  /// Whether the constraint is logical
  static constexpr bool IsLogical() { return false; }
  /// Apriori bounds on the result
  static std::pair<double, double>
  GetAprioriBounds() { return {-INFINITY, INFINITY}; }
};


/// A base class for logical functional constraint.
/// It provides default properties of such a constraint
class LogicalFunctionalConstraintTraits {
public:
  /// Whether the constraint is logical
  static constexpr bool IsLogical() { return true; }
  /// Apriori bounds on the result
  static std::pair<double, double>
  GetAprioriBounds() { return {0.0, 1.0}; }
};


/// Helper, conditional constraint Id
template <class Con>
struct CondConId {
  static const char* description() {
    static std::string descr =
      std::string("Conditional wrapper for constraint type ") +
      typeid(Con).name();
    return descr.c_str();
  }
  static const char* GetTypeName() {
    static std::string nm =
      std::string("Conditional< ") +
      typeid(Con).name() + " >";
    return nm.c_str();
  }
};


/// A wrapper on a static constraint \a Con making it conditional
template <class Con>
class ConditionalConstraint :
    public CustomFunctionalConstraint<
    Con,
    ParamArray0,
    LogicalFunctionalConstraintTraits,
    CondConId<Con> > {
public:
  /// ConType
  using ConType = Con;

  /// Arguments
  using Arguments = ConType;

  /// Base class
  using Base = CustomFunctionalConstraint<
    Con, ParamArray0, LogicalFunctionalConstraintTraits, CondConId<Con> >;

  /// Default constructor
  ConditionalConstraint() = default;

  /// Construct from arguments only
  ConditionalConstraint(Arguments args) noexcept :
    Base(std::move(args)) { }

  /// Construct from resvar and arguments
  ConditionalConstraint(int varr, Arguments args) noexcept :
     Base(varr, std::move(args)) { }

  /// Get the wrapped constraint, const
  const ConType& GetConstraint() const { return Base::GetArguments(); }

  /// Get the wrapped constraint
  ConType& GetConstraint() { return Base::GetArguments(); }

  /// Reuse GetResultVar()
  using Base::GetResultVar;

  /// Equality
  bool operator==(const ConditionalConstraint& cc) const {
    return GetConstraint()==cc.GetConstraint();
  }

};

////////////////////////////////////////////////////////////////////////
/// Args is the argument type, e.g., array of variables, or an expression
/// Params is the parameter type, e.g., array of numbers. Can be empty
#define DEF_CUSTOM_FUNC_CONSTR_WITH_PRM(Name, Args, Params, NumLogic, Descr) \
  struct Name ## Id { \
    static constexpr const char* description() { return Descr; } \
    static constexpr const char* GetTypeName() { return #Name; } \
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
    NumericFunctionalConstraintTraits, Descr)
/// Custom logical constraint with parameter data
#define DEF_LOGICAL_FUNC_CONSTR_WITH_PRM(Name, Args, Params, Descr) \
  DEF_CUSTOM_FUNC_CONSTR_WITH_PRM(Name, Args, Params, \
    LogicalFunctionalConstraintTraits, Descr)

/// A wrapper on a static constraint making it conditional
#define DEF_CONDITIONAL_CONSTRAINT_WRAPPER(Name, StaticConName) \
  using Name = ConditionalConstraint< StaticConName >

} // namespace mp

#endif // BASIC_CONSTR_H
