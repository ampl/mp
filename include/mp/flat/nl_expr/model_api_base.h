/*
 Basic expression-based model API definitions.

 Copyright (C) 2024 AMPL Optimization Inc.

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization Inc disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Gleb Belov <gleb@ampl.com>
*/
#ifndef MODEL_API_BASE_H
#define MODEL_API_BASE_H

#include <deque>
#include <utility>

#include "mp/flat/model_api_base.h"
#include "mp/flat/nl_expr/constr_nl.h"

namespace mp {

/// ModelAPIs handling expression trees should derive from
/// BasicExprModelAPI<Impl, Expr>,
/// with Impl the final implementation class (CRTP)
/// and Expr default-constructible.
template <class Impl, class ExprType=void*>
class BasicExprModelAPI
    : public BasicFlatModelAPI {
public:
  using Expr = ExprType;
  /// Placeholder for GetTypeName()
  static const char* GetTypeName()    { return "BasicExprModelAPI"; }

/// A ModelAPI accepting NL trees should declare this.
///
/// - NotAccepted: not compiled
/// - AcceptedButNotRecommended: compiled but off by default (option acc:_expr)
/// - Recommended: on by default
#define ACCEPT_EXPRESSION_INTERFACE(val) \
  static constexpr ExpressionAcceptanceLevel \
  ExpressionInterfaceAcceptanceLevel() { return ExpressionAcceptanceLevel::val; }

  /// Whether accepts NLObjective.
  /// SCIP does not.
  static int AcceptsNLObj() { return 1; }

/// Reuse inherited names
  USE_BASE_CONSTRAINT_HANDLERS(BasicFlatModelAPI)


  /// Placeholder for AddExpression<>()
  template <class Expression>
  Expr AddExpression(const Expression& e) {
    MP_RAISE(
        std::string("Not handling expression type '") +
        e.GetTypeName() +
        "'. Provide a handler or a converter method");
  }

/// Derived backends have to tell C++ to use default handlers if they are needed
/// when they overload AddExpression(), due to C++ name hiding
#define USE_BASE_EXPRESSION_HANDLERS(BaseBackend) \
  using BaseBackend::AddExpression;


  /// NL model item accessors

  /// @brief Is expression logical?
  template <class FlatCon>
  static bool IsLogical(const ExprWrapper<FlatCon>& expr)
  { return expr.GetFlatConstraint().IsLogical(); }

  /// Get num linear terms
  int GetLinSize(const NLConstraint& nlc) const {
    return nlc.GetMainCon().size();
  }

  /// Get linear coef \a i.
  double GetLinCoef(const NLConstraint& nlc, int i) const {
    return nlc.GetMainCon().coef(i);
  }

  /// Get pointer to linear part's coef array
  const double* GetLinCoefs(const NLConstraint& nlc) const
  { return nlc.GetMainCon().pcoefs(); }

  /// Get linear part var \a i.
  int GetLinVar(const NLConstraint& nlc, int i) const {
    assert(IsVarProper(nlc.GetMainCon().var(i)));
    return nlc.GetMainCon().var(i);
  }

  /// Get pointer to linear part's var array
  const int* GetLinVars(const NLConstraint& nlc) const
  { return nlc.GetMainCon().pvars(); }

  /// Get the expression term of an \a NLObjective.
  /// @note Can return the dummy expression
  ///   via Impl::GetZeroExpression().
  ExprType GetExpression(const NLObjective& nlo) {
    const auto i_expr = nlo.HasExpr()
                            ? nlo.ExprIndex() : -1;
    if (i_expr<0)
      return MPD( GetZeroExpression() );
    return GetInitExpression(i_expr);       // could be explicified
  }

  /// Get the expression term of an \a NLConstraint.
  /// @note Can return the dummy expression
  ///   via Impl::GetZeroExpression().
  ExprType GetExpression(const NLConstraint& nlc) {
    const auto i_expr = nlc.HasExpr()
                            ? nlc.ExprIndex() : -1;
    if (i_expr<0)
      return MPD( GetZeroExpression() );
    return GetInitExpression(i_expr);       // could be explicified
  }

  /// Get NLConstraint's lower bound
  double GetLower(const NLConstraint& nlc) const {
    return nlc.GetMainCon().lb();
  }

  /// Get NLConstraint's upper bound
  double GetUpper(const NLConstraint& nlc) const {
    return nlc.GetMainCon().ub();
  }

  /// Get NLConstraint's name
  const char* GetName(const NLConstraint& nl) const {
    return nl.GetName();
  }

  /// NLAssign constraint name.
  template <int sense>
  const char* GetName(const NLBaseAssign<sense>& nll) {
    return nll.GetName();
  }

  /// Get the expression term of an \a NLBaseAssign,
  /// i.e., NLAssignLE, NLAssignEQ, NLAssignGE.
  /// @note For each such 'explicification' constraint,
  /// the expression
  /// can be normally accessed only once. To access them
  /// repeatedly, call ResetIniExprRetrievedFlags()
  /// for every repetition.
  template <int sense>
  ExprType GetExpression(const NLBaseAssign<sense>& nll) {
    assert( nll.GetVar()>=0 );
    return GetPureInitExpression(nll.GetVar());
  }

  /// Get the variable of an \a NLBaseAssign.
  template <int sense>
  int GetVariable(const NLBaseAssign<sense>& nll) {
    assert( nll.GetVar()>=0 );
    return nll.GetVar();
  }

  /// Get the expression term of an \a NLLogical.
  /// @note The expression
  /// can be normally accessed only once. To access them
  /// repeatedly, call ResetIniExprRetrievedFlags()
  /// for every repetition.
  ExprType GetExpression(const NLLogical& nll) {
    assert( nll.GetCapturedResultVar()>=0 );
    return GetPureInitExpression(nll.GetCapturedResultVar());
  }

  /// Get the value of NLLogical (true/false).
  bool GetValue(const NLLogical& nll) {
    assert( nll.GetCapturedResultVar()>=0 );
    return nll.IsTrue();
  }

  /// Get the expression term of an \a NLBaseReif,
  /// i.e., NLReifImpl, NLReifEquiv, NLReifRimpl.
  /// @note For each such 'explicification' constraint,
  /// the expression
  /// can be normally accessed only once. To access them
  /// repeatedly, call ResetIniExprRetrievedFlags()
  /// for every repetition.
  template <int sense>
  ExprType GetExpression(const NLBaseReif<sense>& nll) {
    assert( nll.GetBVar()>=0 );
    return GetPureInitExpression(nll.GetBVar());
  }

  /// Get the variable of an \a NLReification.
  template <int sense>
  int GetVariable(const NLBaseReif<sense>& nll) {
    assert( nll.GetBVar()>=0 );
    return nll.GetBVar();
  }

  /// Get the expression term of an \a NLComplementarity.
  ExprType GetExpression(const NLComplementarity& nlcc) {
    const auto i_expr = nlcc.HasExpr()
                            ? nlcc.ExprIndex() : -1;
    if (i_expr<0)
      return MPD( GetZeroExpression() );
    return GetInitExpression(i_expr);       // could be explicified
  }

  /// GetLinSize(le)
  int GetLinSize(const NLAffineExpression& le) const
  { return le.GetFlatConstraint().GetAffineExpr().size(); }
  /// GetLinCoef(le, i)
  double GetLinCoef(const NLAffineExpression& le, int i) const
  { return le.GetFlatConstraint().GetAffineExpr().coef(i); }
  /// GetLinTerm(le, i)
  Expr GetLinTerm(const NLAffineExpression& le, int i)
  { return GetInitExpression(le.GetFlatConstraint().GetAffineExpr().var(i)); }
  /// GetQuadSize() for NLAffineExpression: convenience method
  int GetQuadSize(const NLAffineExpression& ) const { return 0; }
  /// GetConstTerm(le)
  double GetConstTerm(const NLAffineExpression& le) const
  { return le.GetFlatConstraint().GetAffineExpr().constant_term(); }

  /// GetLinSize(qe)
  int GetLinSize(const NLQuadExpression& qe) const
  { return qe.GetFlatConstraint().GetQuadExpr().GetBody().GetLinTerms().size(); }
  /// GetLinCoef(qe, i)
  double GetLinCoef(const NLQuadExpression& qe, int i) const
  { return qe.GetFlatConstraint().GetQuadExpr().GetBody().GetLinTerms().coef(i); }
  /// GetLinTerm(qe, i)
  Expr GetLinTerm(const NLQuadExpression& qe, int i)
  { return GetInitExpression(qe.GetFlatConstraint().GetQuadExpr().GetBody().GetLinTerms().var(i)); }

  /// GetQuadSize(qe)
  int GetQuadSize(const NLQuadExpression& qe) const
  { return qe.GetFlatConstraint().GetQuadExpr().GetBody().GetQPTerms().size(); }
  /// GetQuadCoef(qe, i)
  double GetQuadCoef(const NLQuadExpression& qe, int i) const
  { return qe.GetFlatConstraint().GetQuadExpr().GetBody().GetQPTerms().coef(i); }
  /// GetQuadTerm1(qe, i)
  Expr GetQuadTerm1(const NLQuadExpression& qe, int i)
  { return GetInitExpression(qe.GetFlatConstraint().GetQuadExpr().GetBody().GetQPTerms().var1(i)); }
  /// GetQuadTerm2(qe, i)
  Expr GetQuadTerm2(const NLQuadExpression& qe, int i)
  { return GetInitExpression(qe.GetFlatConstraint().GetQuadExpr().GetBody().GetQPTerms().var2(i)); }

  /// GetConstTerm(qe)
  double GetConstTerm(const NLQuadExpression& qe) const
  { return qe.GetFlatConstraint().GetQuadExpr().constant_term(); }

  /// Get number of arguments
  template <class FlatExpression>
  int GetNumArguments(const FlatExpression& fe)
  { return fe.GetFlatConstraint().GetArguments().size(); }

  /// Get argument expression [\a i]
  template <class FlatExpression>
  Expr GetArgExpression(const FlatExpression& fe, int i)
  { return GetInitExpression(fe.GetFlatConstraint().GetArguments().at(i)); }

  /// Get number of parameters
  template <class FlatExpression>
  int GetNumParameters(const FlatExpression& fe)
  { return fe.GetFlatConstraint().GetParameters().size(); }

  /// Get expression parameter [\a i]
  template <class FlatExpression>
  double GetParameter(const FlatExpression& fe, int i)
  { return fe.GetFlatConstraint().GetParameters().at(i); }


  /// Get the expression term of a \a ConditionalConstraint.
  template <class Rhs>
  ExprType GetExpression(
      const ExprWrapper< ConditionalConstraint<
          AlgebraicConstraint<LinTerms, Rhs> > >& cc) {
    return MPD( GetInitExpression(
        cc.GetFlatConstraint().
        GetConstraint().get_representing_variable()) );
  }

  /// Get the expression term of a \a ConditionalConstraint.
  template <class Rhs>
  double GetRHS(
      const ExprWrapper< ConditionalConstraint<
          AlgebraicConstraint<LinTerms, Rhs> > >& cc) {
    return cc.GetFlatConstraint().GetConstraint().rhs();
  }

  /// Init standard options
  void InitStandardOptions() {
    BasicFlatModelAPI::InitStandardOptions();
  }

  /// Placeholder for InitCustomOptions()
  void InitCustomOptions() { }


  ////////////////////// INTERNAL ////////////////////////

private:
  /// Get Expr for a result variable, if proper/explicit,
  /// or for the implicit InitExpression().
  Expr GetInitExpression(int i_expr) {
    /// FlatConverter should have marked all expressions;
    MP_ASSERT_ALWAYS(i_expr < (int)is_expr_stored_.size(),
                     "unexpected expression index");
    // if (i_expr >= (int)is_expr_stored_.size()) {
  //     is_expr_stored_.resize(int(i_expr*1.3)+1);
  //     expr_stored_.resize(int(i_expr*1.3)+1);
  //     is_init_expr_retrieved_.resize(int(i_expr*1.3)+1);
  //     is_var_proper_.resize(int(i_expr*1.3)+1, true);    // proper by default
  //   }
    if (!is_expr_stored_[i_expr]) {
      is_expr_stored_[i_expr] = true;
      if (IsVarProper(i_expr)) {
        expr_stored_[i_expr] = MPD( GetVarExpression(i_expr) );
      } else {
        get_init_expr_(i_expr, &expr_stored_[i_expr]);
      }
    }
    return expr_stored_[i_expr];  // ...............
  }

  /// Get Expr for the initializing expression of a result variable,
  /// even the variable is explicit.
  /// Solver expression for the given variable's init expression.
  /// This is called for NLConstraint
  /// and for NLReification. For them, during result explicification,
  /// the result is marked as variable,
  /// but we still need the expression.
  /// We have to trust that the init expression is accepted.
  /// @note GetInitExpression(\a i_expr) might still return
  ///   the Expr for the result variable \a i_expr.
  Expr GetPureInitExpression(int i_expr) {
    /// FlatConverter should have marked all expressions;
    assert(i_expr < (int)is_expr_stored_.size());
    assert(!is_init_expr_retrieved_[i_expr]);             // not twice explicified
    if (IsVarProper(i_expr)) {
      is_init_expr_retrieved_[i_expr] = true;             // is being explicified
      if (!is_init_expr_stored_[i_expr]) {
        is_init_expr_stored_[i_expr] = true;
        get_init_expr_(i_expr, &init_expr_stored_[i_expr]);
      }
      return init_expr_stored_[i_expr];  // ...............
    }
    return GetInitExpression(i_expr);                     // standard case
  }


protected:
  /// Custom type trait to check if a type is an instance of ExprWrapper
  template <typename T>
  struct is_flat_expr : std::false_type {};

  /// Specialize
  template <typename U>
  struct is_flat_expr<ExprWrapper<U>> : std::true_type {};

  /// Visit arguments of an item.
  /// @param FlatItem: underlying flat item
  /// @param Lambda: to be called on each argument's Expr
  template <class FlatItem, class Lambda>
  inline void VisitArguments(
      const FlatItem& expr, Lambda lambda,
      typename std::enable_if<is_flat_expr<
          typename std::decay<FlatItem>::type>::value>::type* = nullptr) {
    mp::VisitArguments(expr.GetFlatConstraint(),
                       [this,lambda](int v) {
                         lambda(GetInitExpression(v));
                       });
  }

  /// Other types
  template <class FlatItem, class Lambda>
  inline void VisitArguments(
      const FlatItem& expr, Lambda lambda,
      typename std::enable_if<!is_flat_expr<
          typename std::decay<FlatItem>::type>::value>::type* = nullptr) {
  }


public:
  /// Pass vector of var proper flags
  void PassVarProperFlags(std::vector<bool> isvp) {
    is_var_proper_ = std::move(isvp);
    is_expr_stored_.resize(is_var_proper_.size());    // allocate
    expr_stored_.resize(is_var_proper_.size());
    is_init_expr_retrieved_.resize(is_var_proper_.size());
    is_init_expr_stored_.resize(is_var_proper_.size());
    init_expr_stored_.resize(is_var_proper_.size());
  }

  /// Is var proper?
  bool IsVarProper(int i) const {
    assert(i>=0 && i<(int)is_var_proper_.size());
    return is_var_proper_[i];
  }

  /// Reset init_expr_retrieved_ flags,
  /// e.g., after manual marking.
  /// Necessary when visiting expressions several times,
  /// actually only explicified expressions.
  void ResetIniExprRetrievedFlags() {
    auto sz = is_init_expr_retrieved_.size();
    is_init_expr_retrieved_.clear();
    is_init_expr_retrieved_.resize(sz);
  }

  /// Init expr getter type
  using InitExprGetterType = std::function<void(int i_expr, void* pexpr)>;

  /// Provide init expr getter
  void PassInitExprGetter(InitExprGetterType gt)
  { get_init_expr_ = gt; }


private:
  std::vector<bool> is_var_proper_;
  /// actual Expr's of the result vars (for explicified exprds)
  /// or the expressions, otherwise
  std::vector<bool> is_expr_stored_;
  /// Expression cache
  std::deque<ExprType> expr_stored_;    // to keep iterators valid

  /// init exprs: for explicified expressions, whether
  /// the actual init expr-s retrieved.
  /// Should be cleared by ResetIniExpr...
  /// if retrieving repeatedly.
  std::vector<bool> is_init_expr_retrieved_;
  std::vector<bool> is_init_expr_stored_;
  /// Expression cache
  std::deque<ExprType> init_expr_stored_;    // to keep iterators valid
  InitExprGetterType get_init_expr_;
};

}  // namespace mp

#endif // MODEL_API_BASE_H
