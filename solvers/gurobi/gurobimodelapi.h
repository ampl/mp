#ifndef GUROBIMODELAPI_H
#define GUROBIMODELAPI_H

#include <deque>
#include "mp/utils-vec.h"

#include "mp/env.h"
#include "mp/flat/nl_expr/model_api_base.h"
#include "gurobicommon.h"

namespace mp {

/// How we reference an expression.
/// @note To create,
///   use methods MakeConstantExpr(), MakeVarExpr(), MakeExpr()
class GRB_Expr {
public:
  /// Default construct
  GRB_Expr() : kind_(-1) { }
  /// Construct a var/expr index.
  GRB_Expr(bool isvar, int index)
      : kind_(isvar), id_(index) { }
  /// Construct constant expression.
  GRB_Expr(double v) : val_(v), kind_(2) { }

  /// Check validity
  bool IsValid() const { return kind_>=0; }

  /// Is an expression?
  bool IsExpr() const { return 0==kind_; }
  /// Get expr index
  int GetExprIndex() const { assert(IsExpr()); return id_; }

  /// Is an variable?
  bool IsVar() const { return 1==kind_; }
  /// Get var index
  int GetVarIndex() const { assert(IsVar()); return id_; }

  /// Is a constant?
  bool IsConst() const { return 2==kind_; }
  /// Get constant
  double GetConst() const { assert(IsConst()); return val_; }
private:
  double val_; int kind_; int id_;
};


/// GurobiModelAPI
class GurobiModelAPI :
    public GurobiCommon,
    public EnvKeeper,
    public BasicExprModelAPI<GurobiModelAPI, GRB_Expr> {
  /// Typedef main base class
  using BaseModelAPI = BasicExprModelAPI<GurobiModelAPI, GRB_Expr>;

public:
  /// Model API name
  static const char* GetTypeName();
  /// Reserved
  static const char* GetLongName() { return nullptr; }

  /// Construct
  GurobiModelAPI(Env& e) : EnvKeeper(e) { }

  /// This is called before model is pushed to the Backend
  void InitProblemModificationPhase(const FlatModelInfo*);
  /// Chance to call GRBupdatemodel()
  void FinishProblemModificationPhase();

  /////////////////////////////////////////////////////////////////////////////
  //////////////////////////// MODELING ACCESSORS /////////////////////////////
  /////////////////////////////////////////////////////////////////////////////

  //////////////////////////// VARIABLES //////////////////////////////////////
  void AddVariables(const VarArrayDef& );

  //////////////////////////// OBJECTIVES /////////////////////////////////////
  void SetLinearObjective( int iobj, const LinearObjective& lo );
  /// Whether accepting quadratic objectives:
  /// 0 - no, 1 - convex, 2 - nonconvex
  static int AcceptsQuadObj() { return 2; }
  void SetQuadraticObjective( int iobj, const QuadraticObjective& qo );

  //////////////////////////// GENERAL CONSTRAINTS ////////////////////////////
  USE_BASE_CONSTRAINT_HANDLERS(BaseModelAPI)

  /// Gurobi does not properly support range linear constraints.
  /// Gurobi 9.5: GRBaddrangeconstr() creates a slack variable but
  /// does not account for it in the basis information, etc.

  /// LinCon(LE/EQ/GE) should have 'Recommended' for all backends
  /// and have an implementation.
  ACCEPT_CONSTRAINT(LinConLE, Recommended, CG_Linear)
  void AddConstraint(const LinConLE& lc);
  ACCEPT_CONSTRAINT(LinConEQ, Recommended, CG_Linear)
  void AddConstraint(const LinConEQ& lc);
  ACCEPT_CONSTRAINT(LinConGE, Recommended, CG_Linear)
  void AddConstraint(const LinConGE& lc);

  /// Ask if the solver accepts non-convex quadratic constraints
  /// (Gurobi needs option nonconvex=2 for solving)
  static constexpr bool AcceptsNonconvexQC() { return true; }

  /// If cvt:prod=7 (and not 5) default.
  /// Recommendation to return the opposite value as
  /// AcceptsNonconvexQC().
  static constexpr bool WantLogicalizedProd2Bin()
  { return !AcceptsNonconvexQC(); }

  /// Ask if the solver can recognize SOCP corner cases
  /// (non-std representations such as xy>=1, see tests)
  /// from quadratic representations
  static constexpr bool CanSOCPCornerCasesFromQC() { return true; }

  /// If using quadratics,
  /// QuadCon(LE/EQ/GE) should have 'Recommended'
  /// and have an implementation.
  /// QuadConRange is optional.
  ACCEPT_CONSTRAINT(QuadConLE, Recommended, CG_Quadratic)
  void AddConstraint(const QuadConLE& qc);
  ACCEPT_CONSTRAINT(QuadConEQ, Recommended, CG_Quadratic)
  void AddConstraint(const QuadConEQ& qc);
  ACCEPT_CONSTRAINT(QuadConGE, Recommended, CG_Quadratic)
  void AddConstraint(const QuadConGE& qc);

  /// Discrete general constraints
  ACCEPT_CONSTRAINT(MaxConstraint, Recommended, CG_General)
  void AddConstraint(const MaxConstraint& mc);
  ACCEPT_CONSTRAINT(MinConstraint, Recommended, CG_General)
  void AddConstraint(const MinConstraint& mc);
  ACCEPT_CONSTRAINT(AbsConstraint, Recommended, CG_General)
  void AddConstraint(const AbsConstraint& absc);
  ACCEPT_CONSTRAINT(AndConstraint, Recommended, CG_General)
  void AddConstraint(const AndConstraint& cc);
  ACCEPT_CONSTRAINT(OrConstraint, Recommended, CG_General)
  void AddConstraint(const OrConstraint& mc);
  ACCEPT_CONSTRAINT(IndicatorConstraintLinLE, Recommended, CG_General)
  void AddConstraint(const IndicatorConstraintLinLE& mc);
  ACCEPT_CONSTRAINT(IndicatorConstraintLinEQ, Recommended, CG_General)
  void AddConstraint(const IndicatorConstraintLinEQ& mc);
  ACCEPT_CONSTRAINT(IndicatorConstraintLinGE, Recommended, CG_General)
  void AddConstraint(const IndicatorConstraintLinGE& mc);
  ACCEPT_CONSTRAINT(PLConstraint, Recommended, CG_General)
  void AddConstraint(const PLConstraint& cc);

  /// Gurobi SOS1/2
  ACCEPT_CONSTRAINT(SOS1Constraint, Recommended, CG_SOS)
  void AddConstraint(const SOS1Constraint& cc);
  ACCEPT_CONSTRAINT(SOS2Constraint, Recommended, CG_SOS)
  void AddConstraint(const SOS2Constraint& cc);

  /// Gurobi nonlinear generals
  ACCEPT_CONSTRAINT(ExpConstraint, Recommended, CG_General)
  void AddConstraint(const ExpConstraint& cc);
  ACCEPT_CONSTRAINT(ExpAConstraint, Recommended, CG_General)
  void AddConstraint(const ExpAConstraint& cc);
  ACCEPT_CONSTRAINT(LogConstraint, Recommended, CG_General)
  void AddConstraint(const LogConstraint& cc);
  ACCEPT_CONSTRAINT(LogAConstraint, Recommended, CG_General)
  void AddConstraint(const LogAConstraint& cc);
  ACCEPT_CONSTRAINT(PowConstExpConstraint, Recommended, CG_General)
  void AddConstraint(const PowConstExpConstraint& cc);
  ACCEPT_CONSTRAINT(SinConstraint, Recommended, CG_General)
  void AddConstraint(const SinConstraint& cc);
  ACCEPT_CONSTRAINT(CosConstraint, Recommended, CG_General) // y = cos(x)
  void AddConstraint(const CosConstraint& cc);  // GRBaddgenconstrCos(x, y);
  ACCEPT_CONSTRAINT(TanConstraint, Recommended, CG_General)
  void AddConstraint(const TanConstraint& cc);

  /// Init GurobiModelAPI driver options
  void InitCustomOptions();

  //////////////////////////////////////////////////////////////////////////
  /// Expressions
#ifdef GRB_OPCODE_CONSTANT

  /// Whether accepts NLObjective
  static int AcceptsNLObj() { return 0; }

  //////////////////////////// EXPRESSION TREES ////////////////////////////
  /// Handle expression trees: inherit basic API
  USE_BASE_EXPRESSION_HANDLERS(BaseModelAPI)

  /// Overall switch
  ACCEPT_EXPRESSION_INTERFACE(Recommended)

  /// GetVarExpression(\a i): expression representing variable 0<=i<n_var.
  /// Only called for 'nonlinear' variables.
  Expr GetVarExpression(int i) { return MakeVarExpr(i); }

  /// GetZeroExpr(): constant 0.0 expression.
  /// Can be used to represent empty expression in an NLConstraint.
  Expr GetZeroExpression() { return MakeEmptyExpr(); }

  /// Gurobi 12 has no classical NL range constraint
  ACCEPT_CONSTRAINT(NLConstraint, NotAccepted, CG_Algebraic)

  /// NLAssignEQ: algebraic expression expicifier.
  /// Meaning: var == expr.
  /// @note Accessors: GetName(), GetExpression(nle), GetVariable(nle).
  ACCEPT_CONSTRAINT(NLAssignEQ, Recommended, CG_General)
  void AddConstraint(const NLAssignEQ& nle);
  /// NLAssignLE: algebraic expression expicifier in positive context.
  /// Meaning: var <= expr.
  /// @note Accessors: GetName(), GetExpression(nle), GetVariable(nle).
  ACCEPT_CONSTRAINT(NLAssignLE, Recommended, CG_General)
  void AddConstraint(const NLAssignLE& nle);
  /// NLAssignGE: algebraic expression expicifier in negative context.
  /// Meaning: var >= expr.
  /// @note Accessors: GetName(), GetExpression(nle), GetVariable(nle).
  ACCEPT_CONSTRAINT(NLAssignGE, Recommended, CG_General)
  void AddConstraint(const NLAssignGE& nle);

  /// @brief Accept NLAffineExpr.
  /// @note Use accessors, not methods;
  /// - GetLinSize(le), GetLinCoef(le, i), GetLinTerm(le, i);
  ///   GetConstTerm(le).
  ACCEPT_EXPRESSION(NLAffineExpression, Recommended);
  Expr AddExpression(const NLAffineExpression& le);

  /// Accept NLQuadExpr.
  /// @note Use accessors, not methods;
  /// - GetLinSize(le), GetLinCoef(le, i), GetLinTerm(le, i);
  ///   GetQuadSize(le), GetQuadCoef(le, i),
  ///   GetQuadTerm1(le, i), GetQuadTerm2(le, i);
  ///   GetConstTerm(le).
  ACCEPT_EXPRESSION(NLQuadExpression, Recommended);
  Expr AddExpression(const NLQuadExpression& le);

  /// Each expression can be accepted as a proper expression,
  /// or as a flat functional constraint var <=/==/>= expr
  /// (in this case, with variables as arguments).
  /// The equality/inequality type of the flat constraint is
  /// determied by GetContext().
  ///
  /// @note Use accessor: GetArgExpression(ee, 0)
  /// - don't use ...Expression's methods.
  ///
  /// Similar for other expression types.


  ACCEPT_EXPRESSION(ExpExpression, Recommended)
  Expr AddExpression(const ExpExpression& );
  ACCEPT_EXPRESSION(ExpAExpression, Recommended)
  Expr AddExpression(const ExpAExpression& );
  ACCEPT_EXPRESSION(LogExpression, Recommended)
  Expr AddExpression(const LogExpression& );
  ACCEPT_EXPRESSION(LogAExpression, Recommended)
  Expr AddExpression(const LogAExpression& );
  /// @note Use accessor: GetParameter(pe, 0)
  ///   - don't use PowConstExpExpression's methods.
  ACCEPT_EXPRESSION(PowConstExpExpression, Recommended)
  Expr AddExpression(const PowConstExpExpression& );
  ACCEPT_EXPRESSION(SinExpression, Recommended)
  Expr AddExpression(const SinExpression& );
  ACCEPT_EXPRESSION(CosExpression, Recommended)
  Expr AddExpression(const CosExpression& );
  ACCEPT_EXPRESSION(TanExpression, Recommended)
  Expr AddExpression(const TanExpression& );

#ifdef GRB_OPCODE_TANH
  ACCEPT_EXPRESSION(TanhExpression, Recommended)
  Expr AddExpression(const TanhExpression& );
#endif

  ACCEPT_EXPRESSION(DivExpression, Recommended)
  Expr AddExpression(const DivExpression& );

public:
  /// Formula
  class Formula {
  public:
    /// Construct with 1 entry
    Formula(int oc=-1, double d=-1.0, int p=-1)
        : opcode_(1, oc), data_(1, d), parent_(1, p) { }
    /// Length
    int size() const { assert(is_length_ok()); return opcode_.size(); }
    /// Opcodes
    const int* opcodes() const { return opcode_.data(); }
    /// Data
    const double* data() const { return data_.data(); }
    /// Opcodes
    const int* parents() const { return parent_.data(); }
  protected:
    friend class GurobiModelAPI;
    /// Are the array lengths equal?
    bool is_length_ok() const {
      return opcode_.size()==data_.size() && parent_.size()==data_.size();
    }
    /// Append a formula
    void Append(const Formula& frm);
  private:
    SmallVec<int, 16> opcode_;
    SmallVec<double, 16> data_;
    SmallVec<int, 16> parent_;
    // std::vector<int> opcode_;
    // std::vector<double> data_;
    // std::vector<int> parent_;
  };

protected:
  /// Create and store a formula
  template <class MPExpr>
  Expr CreateFormula(const MPExpr& mpexpr, int opcode);

  /// start formula with opcode
  static Formula StartFormula(int opcode)
  { return {opcode, -1.0, -1}; }

  /// Append argument subformula
  void AppendArgument(Formula& frm, Expr expr);

  /// Get formula from ID
  /// @note the reference should be used immediately and not stored;
  /// when calling another GetFormula(),
  /// the previous reference may become invalid.
  const Formula& GetFormula(Expr expr) const;

  /// Make a constant expression.
  static Expr MakeConstantExpr(double v) { return {v}; }

  /// Make an empty expression.
  static Expr MakeEmptyExpr() { return MakeConstantExpr(0.0); }

  /// Make an expression representing variable \a v.
  static Expr MakeVarExpr(int v) { return {true, v}; }

  /// Make a proper expression with index \a i.
  static Expr MakeExpr(int i) { return {false, i}; }

  template <class MPExpr>
  void AppendLinAndConstTerms(Formula& , const MPExpr& );

  template <class MPExpr>
  void AppendQuadTerms(Formula& , const MPExpr& );

private:
  /// Store every subformula - @todo only repeated ones
  /// @note deque: for references to stay valid
  ///   while we add subformulas
  std::deque<Formula> formulas_;
  mutable Formula formula_tmp_;

#endif  // GRB_OPCODE_CONSTANT
//////////////////////////////////////////////////////////////////

protected:
  /// First objective's sense
  void NoteGurobiMainObjSense(obj::Type s);
  obj::Type GetGurobiMainObjSense() const;


private:
  bool has_quadratic_obj_ = false;
  /// The sense of the main objective
  obj::Type main_obj_sense_;

  /// To zero out last objective
  std::vector<int> obj_ind_save_;

  /// These options are stored in the class as variables
  /// for direct access
  struct Options {
  } storedOptions_;
};

} // namespace mp

#endif // GUROBIMODELAPI_H
