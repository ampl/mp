#ifndef CONSTR_NL_H
#define CONSTR_NL_H

#include <cmath>
#include <cassert>

#include "mp/error.h"
#include "mp/flat/constr_std.h"

namespace mp {

/// Class NLConstraint.
/// Algebraic range constraint with a linear part
/// and an expression term: `lb <= a'x + expr <= ub`.
/// LinConRange is a member to avoid overloading
/// when deriving from an existing constraint type.
class NLConstraint
    : public BasicConstraint, public NumericFunctionalConstraintTraits {
public:
  /// Constraint type name
  static const char* GetTypeName() {
    return "NLConstraint";
  }

  /// Constructor.
  /// @param linexpr: linear part
  /// @param expr: result variable of the expression part
  /// @param rng: {lb, ub}
  /// @param fSort=true: whether to sort linear terms
  NLConstraint(
      const LinTerms& lt, int expr, AlgConRange rng,
      bool fSort=true)
      : lcr_(lt, rng, fSort), expr_(expr) { }

  /// Get the main algebraic constraint
  const LinConRange& GetMainCon() const { return lcr_; }

  /// Has expression term?
  bool HasExpr() const { return expr_>=0; }

  /// Expression index.
  /// @note ModelAPI should call self.HasExpression()
  ///   and self.GetExpression()
  ///   to obtain the expression term.
  int ExprIndex() const { assert(HasExpr()); return expr_; }

  /// Compute violation.
  template <class VarInfo>
  Violation
  ComputeViolation(const VarInfo& x, bool logical=false) const {
    double bd = lcr_.GetBody().ComputeValue(x);
    if (HasExpr())
      bd += x[ExprIndex()];     // Add expr value. Assume it's precomputed
    if (!logical) {
      if (lcr_.lb() > bd)
        return {lcr_.lb() - bd, lcr_.lb()};
      if (bd > lcr_.ub())
        return {bd - lcr_.ub(), lcr_.ub()};
      return
          {std::max( // negative. Same for strict cmp?
               lcr_.lb() - bd, bd - lcr_.ub()),
           0.0};
    }
    return {double(!lcr_.is_valid(bd)), 1.0};
  }


private:
  LinConRange lcr_;
  int expr_ {-1};
};


/// Specialize
inline void VisitArguments(const NLConstraint& nlc,
                           std::function<void (int) > argv) {
  VisitArguments(nlc.GetMainCon(), argv);
  VisitArguments(VarArray1{nlc.ExprIndex()}, argv);
}


/// Export to JSON
inline void WriteJSON(JSONW jw,
                      const NLConstraint& nlc) {
  WriteJSON(jw["main_alg_con"], nlc.GetMainCon());
  if (nlc.HasExpr())
    WriteJSON(jw["expr_index"], nlc.ExprIndex());
}

/// Write RhsCon without name.
template <class Writer, class Names>
inline void WriteModelItem(Writer& wrt,
                           const NLConstraint& nlc,
                           Names& vnam) {
  wrt << "NLExprIndex: " << vnam.at(nlc.ExprIndex()) << " IN: ";
  WriteModelItem(wrt, nlc.GetMainCon(), vnam);
}


/// Algebraic expression explicifier.
/// Syntax sugar for the assignment: var <=/==/>= expr.
/// Can have special meaning in certain solvers.
/// Sense: equality (0), >= (1), <= (-1).
/// Can be implemented as == for all senses (e.g., GRBaddgenconstrNL),
/// the inequalities can be used to preserve convexity.
/// This is a static constraint.
template <int sense>
class NLBaseAssign
    : public BasicConstraint, public NumericFunctionalConstraintTraits {
public:
  /// Constraint type name
  static const char* GetTypeName() {
    if (0==sense) return "NLAssignEQ";
    if (-1==sense) return "NLAssignLE";
    if (1==sense) return "NLAssignGE";
    MP_RAISE("NLBaseAssign: unknown sense");
  }

  /// Construct
  NLBaseAssign(int b) : bvar_(b) { }

  /// Get var
  int GetVar() const { return bvar_; }

  /// Imitate an algebraic con
  int size() const { return 1; }
  /// Imitate
  double coef(int i) const { assert(!i); return -1.0; }
  /// Imitate
  int var(int i) const { assert(!i); return GetVar(); }

  /// Throw - should not be used
  VarArray1 GetArguments() const { MP_RAISE("No marking for NL items"); }

  // Compute violation... Should be 0

private:
  int bvar_ {-1};
};


/// Typedef NLAssignEQ
using NLAssignEQ = NLBaseAssign<0>;
/// Typedef NLAssignLE
using NLAssignLE = NLBaseAssign<-1>;
/// Typedef NLAssignGE
using NLAssignGE = NLBaseAssign<1>;


/// Specialize.
/// Don't mark arguments because it was marked 'used'
/// when extracting it from NLConstraint
template <int sense>
inline void VisitArguments(const NLBaseAssign<sense>& nlba,
                           std::function<void (int) > argv) {
  // VisitArguments(VarArray1{nlba.GetVar()}, argv);
}


/// Write a Reification
template <int sense>
inline void WriteJSON(JSONW jw,
                      const NLBaseAssign<sense>& reif) {
  jw["var_explicit_assign"] = reif.GetVar();
  jw["sense"] = sense;
}

/// Write RhsCon without name.
template <class Writer, int sense, class Names>
inline void WriteModelItem(Writer& wrt,
                           const NLBaseAssign<sense>& nlr,
                           Names& vnam) {
  wrt << "EXPLICIT ASSIGN var: " << vnam.at(nlr.GetVar());
}


/// NLComplementarity.
/// Complementarity constraint where the expression
/// has a linear and a non-linear part
class NLComplementarity
    : public BasicConstraint, public NumericFunctionalConstraintTraits {
public:
  /// Constraint type name
  static const char* GetTypeName() {
    return "NLComplementarity";
  }

  /// Constructor.
  /// @param expr_lin: linear part
  /// @param expr_nonlin_var: result variable of the expression part
  /// @param cvar: complementing variable
  NLComplementarity(AffineExpr expr_lin, int expr_nonlin_var, int cvar) :
      ccl_{std::move(expr_lin), cvar}, expr_(expr_nonlin_var) { }

  /// Get the affine part of the expression
  /// @todo Should the constant here be 0?
  const AffineExpr& GetLinearPart() const
  { return ccl_.GetExpression(); }

  /// Has expression term?
  bool HasExpr() const { return expr_>=0; }

  /// Expression index.
  /// @note ModelAPI should call self.HasExpression()
  ///   and self.GetExpression()
  ///   to obtain the expression term.
  int ExprIndex() const { assert(HasExpr()); return expr_; }

  /// the complementing variable
  int GetCVar() const { return ccl_.GetVariable(); }

  /// Compute violation
  template <class VarInfo>
  Violation ComputeViolation(const VarInfo& x) const {
    auto ve = GetLinearPart().ComputeValue(x);
    if (HasExpr())
      ve += x[ExprIndex()];   // Add expr value. Assume it's precomputed
    if (x.is_at_lb(GetCVar()))
      return {-ve, 0.0};
    else if (x.is_at_ub(GetCVar()))
      return {ve, 0.0};
    return {std::fabs(ve), 0.0};
  }


private:
  ComplementarityLinear ccl_;
  int expr_ {-1};
};


/// Specialize
inline void VisitArguments(const NLComplementarity& nlcc,
                           std::function<void (int) > argv) {
  VisitArguments(nlcc.GetLinearPart(), argv);
  if (nlcc.HasExpr())
    VisitArguments(VarArray1{nlcc.ExprIndex()}, argv);
  VisitArguments(VarArray1{nlcc.GetCVar()}, argv);
}


/// Export to JSON
inline void WriteJSON(JSONW jw,
                      const NLComplementarity& nlcc) {
  WriteJSON(jw["lin_part"], nlcc.GetLinearPart());
  if (nlcc.HasExpr())
    WriteJSON(jw["expr_index"], nlcc.ExprIndex());
  WriteJSON(jw["cvar"], nlcc.GetCVar());
}

/// Write RhsCon without name.
template <class Writer, class Names>
inline void WriteModelItem(Writer& wrt,
                           const NLComplementarity& nlcc,
                           Names& vnam) {
  if (nlcc.HasExpr())
    wrt << "NLCCExprIndex: " << vnam.at(nlcc.ExprIndex()) << " IN: ";
  WriteModelItem(wrt, nlcc.GetLinearPart(), vnam);
  wrt << "NLCCVar: " << vnam.at(nlcc.GetCVar());
}



/// NL logical constraint:
/// expr(resvar) <==> const (true or false)
class NLLogical
    : public BasicConstraint {
public:
  /// Constraint type name
  static const char* GetTypeName() {
    return "NLLogical";
  }

  /// Is logical?
  static bool IsLogical() { return true; }

  /// Construct from the result variable
  NLLogical(int rv, bool val) : resvar_(rv), value_(val)
  { assert(rv>=0); }

  /// Is const == true?
  bool IsTrue() const { return value_; }

  /// Get resvar
  int GetCapturedResultVar() const { return resvar_; }

  /// Throw - should not be used
  VarArray1 GetArguments() const { MP_RAISE("No marking for NL items"); }

  /// Compute violation
  template <class VarInfo>
  Violation
  ComputeViolation(const VarInfo& x) const {
    return
        {std::fabs(x[GetCapturedResultVar()] - value_), double(value_)};
  }

private:
  int resvar_ {-1};
  bool value_ {};
};


/// Specialize
inline void VisitArguments(const NLLogical& ,
                           std::function<void (int) > ) {
  // Do nothing because logical expressions are marked as used
}


/// Write an NLLogical
inline void WriteJSON(JSONW jw,
                      const NLLogical& nll) {
  jw["resvar"] = nll.GetCapturedResultVar();
  jw["value"] = nll.IsTrue();
}

/// Write RhsCon without name.
template <class Writer, class Names>
inline void WriteModelItem(Writer& wrt,
                           const NLLogical& nllc,
                           Names& vnam) {
  wrt << "NLLogicalExprIndex: "
      << vnam.at(nllc.GetCapturedResultVar());
  wrt << "NLLogicalExprValue: "
      << nllc.IsTrue();
}


/// Logical expression explicifier.
/// Syntax sugar for reification: b==1 <==> expr(b)==1.
/// Sense: equivalence (0), impl(-1), rimpl (1).
/// This is a static constraint.
template <int sense>
class NLBaseReif
    : public BasicConstraint {
public:
  /// Constraint type name
  static const char* GetTypeName() {
    if (0==sense) return "NLReifEquiv";
    if (-1==sense) return "NLReifImpl";
    if (1==sense) return "NLReifRimpl";
    MP_RAISE("NLReif: unknown sense");
  }

  /// Is logical?
  static bool IsLogical() { return true; }

  /// Construct
  NLBaseReif(int b) : bvar_(b) { }

  /// Get bvar
  int GetBVar() const { return bvar_; }

  /// Throw - should not be used
  VarArray1 GetArguments() const { MP_RAISE("No marking for NL items"); }

  // Compute violation... Should be 0

private:
  int bvar_ {-1};
};


/// Typedef NLReifEquiv
using NLReifEquiv = NLBaseReif<0>;
/// Typedef NLReifImpl
using NLReifImpl = NLBaseReif<-1>;
/// Typedef NLReifRImpl
using NLReifRimpl = NLBaseReif<1>;


/// Specialize
template <int sense>
inline void VisitArguments(const NLBaseReif<sense>& ,
                           std::function<void (int) > ) {
  // Again, do nothing
  // VisitArguments(VarArray1{nlbr.GetBVar()}, argv);
}



/// Write a Reification
template <int sense>
inline void WriteJSON(JSONW jw,
                      const NLBaseReif<sense>& reif) {
  jw["var_explicit_reif"] = reif.GetBVar();
  jw["sense"] = sense;
}

/// Write RhsCon without name.
template <class Writer, int sense, class Names>
inline void WriteModelItem(Writer& wrt,
                           const NLBaseReif<sense>& nlr,
                           Names& vnam) {
  wrt << "EXPLICIT REIF var: " << vnam.at(nlr.GetBVar());
}

}  // namespace mp

#endif // CONSTR_NL_H
