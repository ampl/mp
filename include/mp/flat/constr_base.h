#ifndef BASIC_CONSTR_H
#define BASIC_CONSTR_H

#include <array>
#include <vector>
#include <string>
#include <cmath>
#include <utility>
#include <typeinfo>

#include "mp/format.h"
#include "mp/error.h"
#include "mp/common.h"

#include "mp/flat/context.h"
#include "mp/arrayref.h"

namespace mp {

/// Constraint/obj violation
struct Violation {
  double viol_;    // abs violation: >0 if really violated
  double valX_;    // value compared to
  /// Compute whether violated + relative violation.
  /// Both absolute and relative should be violated
  /// (relative only if refvar!=0.)
  std::pair<bool, double> Check(
      double epsabs, double epsrel) const {
    double violRel {0.0};
    if (viol_ > epsabs
        && (0.0==std::fabs(valX_)
            || (violRel=std::fabs(viol_/valX_))>epsrel))
      return {true, violRel};
    return {false, 0.0};
  }
};


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
    /// Has result var (is functional)?
    bool HasResultVar() const { return false; }
  /// For functional constraints, result variable index
  int GetResultVar() const { return -1; }
  /// Compute violation
  template <class VarInfo>
  Violation ComputeViolation(const VarInfo& ) const
  { return {0.0, 0.0}; }

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
  /// Has result var (is functional)?
  /// We don't store this info in the type (yet.)
  bool HasResultVar() const { return result_var_>=0; }
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
/// Fixed argument array of 3 elements
using VarArray3 = std::array<int, 3>;
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
/// Fixed parameter array of 2 double
using DblParamArray2 = ParamArrayN<double, 2>;
/// Fixed parameter array of 3 double
using DblParamArray3 = ParamArrayN<double, 3>;
/// Variable-length parameter array
using DblParamArray = std::vector<double>;


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
  /// Construct from resvar and arguments.
  /// If Arguments = VarArray1, distinguish from the previous
  /// constructor by putting first int in just one {}
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

  /// Compute violation
  template <class VarVec>
  Violation ComputeViolation(const VarVec& x) const;
};

/// Write flat expr/obj/con variables:
/// vector/array
template <class Writer, class Vec,
          typename T = std::decay_t<
              decltype(*begin(std::declval<Vec>()))> >
inline void WriteModelItem(
    Writer& wrt, const Vec& v,
    const std::vector<std::string>& vnam) {
  static_assert (
  std::is_integral_v<typename Vec::value_type>, "Variable vector: need int's");
  wrt << '[';
  int n=-1;
  for (const auto& el: v) {
    if (++n)
      wrt << ", ";
    wrt << vnam.at(el);
  }
  wrt << ']';
}

/// Write flat expr/obj/con parameters:
/// vector/array
template <class Writer, class Vec,
          typename T = std::decay_t<
              decltype(*begin(std::declval<Vec>()))> >
void WriteModelItemParameters(
    Writer& wrt, const Vec& v) {
  wrt << '[';
  int n=-1;
  for (const auto& el: v) {
    if (++n)
      wrt << ", ";
    wrt << el;
  }
  wrt << ']';
}

/// Specialize WriteModelItem() for CustomFuncCon<>
template <class Writer, class A, class P, class N, class I>
inline void WriteModelItem(
    Writer& wrt,
    const CustomFunctionalConstraint<A,P,N,I>& cfc,
    const std::vector<std::string>& vnam) {
  if (cfc.HasResultVar())   // really functional
    wrt << vnam.at(cfc.GetResultVar()) << " == ";
  wrt << cfc.GetTypeName();
  wrt << '(';
  WriteModelItem(wrt, cfc.GetArguments(), vnam);
  wrt << ", ";
  WriteModelItemParameters(wrt, cfc.GetParameters());
  wrt << ')';
}

/// Very general template to write any flat constraint
/// with name.
template <class Writer, class Con>
inline void WriteFlatCon(Writer& wrt, const Con& c,
                  const std::vector<std::string>& vnam) {
  wrt << c.name() << ": ";
  WriteModelItem(wrt, c, vnam);
}

/// Write a CustomFunctionalConstraint<>
template <class JW, class A, class P, class N, class I>
inline void WriteJSON(JW jw,
                      const CustomFunctionalConstraint<A,P,N,I>& cfc) {
  jw["res_var"] = cfc.GetResultVar();
  WriteJSON(jw["args"], cfc.GetArguments());
  WriteJSON(jw["params"], cfc.GetParameters());
}


/// Dummy template: compute result of functional constraint.
/// Should be implemented for proper functional specializations,
/// but has no sense for conic constraints, for example.
template <class Con, class VarVec>
double ComputeValue(const Con& , const VarVec& ) {
  MP_RAISE(fmt::format("ComputeValue({}) not implemented.",
                       typeid(Con).name()));
  return 0.0;
}

/// Compute violation of functional constraint.
/// Can only be used for proper functional constraints,
/// should be redefined for cones.
template <class Args, class Params,
          class NumOrLogic, class Id, class VarVec>
Violation ComputeViolation(
    const CustomFunctionalConstraint<Args, Params, NumOrLogic, Id>& c,
    const VarVec& x) {
  auto resvar = c.GetResultVar();
  if (!x.recomp_vals()) {    // solver's var values: normal check
    auto viol = x[resvar] - ComputeValue(c, x);
    switch (c.GetContext().GetValue()) {
    case Context::CTX_MIX:
      return {std::fabs(viol), x[resvar]};
    case Context::CTX_POS:
      return {viol, x[resvar]};
    case Context::CTX_NEG:
      return {-viol, x[resvar]};
    default:
      return {INFINITY, 0.0};
    }
  }
  return                              // recomputed var minus solver's
  { std::fabs(x[resvar] - x.raw(resvar))
        + std::max(0.0, x.bounds_viol(resvar)), x[resvar]};
}


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
    Con, ParamArray0,
  LogicalFunctionalConstraintTraits, CondConId<Con> >;

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

  /// Compute violation of a conditional constraint.
  /// If the subconstr is violated but should hold,
  /// return the exact gap, despite this is a logical constraint.
  /// If the subconstr holds but should not,
  /// return the opposite gap. Similar to Indicator.
  template <class VarVec>
  Violation ComputeViolation(const VarVec& x) {
    auto viol = GetConstraint().ComputeViolation(x);
    bool ccon_valid = viol.viol_<=0.0;
    bool has_arg = x[GetResultVar()] >= 0.5;
    switch (this->GetContext().GetValue()) {
    case Context::CTX_MIX:    // Viol is non-positive if holds
      if (has_arg == ccon_valid)
        return {0.0, 0.0};
      return {std::fabs(viol.viol_), viol.valX_};
    case Context::CTX_POS:
      if (has_arg <= ccon_valid)
        return {0.0, 0.0};
      return {viol.viol_, viol.valX_};
    case Context::CTX_NEG:
      if (has_arg >= ccon_valid)
        return {0.0, 0.0};
      return {-viol.viol_, viol.valX_};
    default:
      return {INFINITY, 0.0};
    }
  }
};

/// Write ConditionalCon without name.
template <class Writer, class Con>
inline void WriteModelItem(Writer& wrt,
                    const ConditionalConstraint<Con>& condc,
                    const std::vector<std::string>& vnam) {
  wrt << vnam.at(condc.GetResultVar()) << "==1 <==> ";
  wrt << '(';
  WriteModelItem(wrt, condc.GetArguments(), vnam);
  wrt << ')';
}

/// Write a readable variable definition
template <class Writer>
void WriteVar(Writer& pr, const char* name,
              double lb, double ub, var::Type ty);

/// Write a CondCon
template <class JW, class Con>
inline void WriteJSON(JW jw,
                      const ConditionalConstraint<Con>& condc) {
  jw["res_var"] = condc.GetResultVar();
  WriteJSON(jw["con"], condc.GetConstraint());
}


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

////////////////////////////////////////////////////////////////////////
/// STATIC CONSTRAINTS
/// Workaround: defining as functional  constraint (result unused)
/// Could be solved by a mix-in parent
#define DEF_STATIC_CONSTR(Name, Args, Descr) \
    DEF_LOGICAL_FUNC_CONSTR(Name, Args, Descr)
#define DEF_STATIC_CONSTR_WITH_PRM(Name, Args, Params, Descr) \
    DEF_LOGICAL_FUNC_CONSTR_WITH_PRM(Name, Args, Params, Descr)


} // namespace mp

#endif // BASIC_CONSTR_H
