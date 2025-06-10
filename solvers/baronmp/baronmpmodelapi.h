#ifndef BARONMPMODELAPI_H
#define BARONMPMODELAPI_H

#include <iostream>
#include <memory>
#include <vector>
#include <variant>


#include <functional>

#include "mp/utils-vec.h"
#include "mp/env.h"
#include "baronmpcommon.h"
#include "mp/flat/nl_expr/model_api_base.h"

namespace mp {


    enum class Opcode {
    CONSTANT,
    VARIABLE,
    ADD,
    SUBTRACT,
    MULTIPLY,
    DIVIDE,

    POW,
    LOGA,
    ABS,

    EXP,
    LOG
  };

    // Array of strings corresponding to the Opcode enum
    const std::string opcodeStrings[] = {
        "CONSTANT",
        "VARIABLE",
        "+",
        "-",
        "*",
        "/",

        "^",
        "loga",
        "abs",

        "exp",
        "log"
    };
  
  class BaronmpModelAPI;
  // Define a Node (VExpr) for the expression tree
  class VExpr {
  public:
    static std::function<std::string(int)> varName;
    void reset() {
      children.clear();
    }
    Opcode opcode;  // Stores the type of the operation or expression
    std::variant<double, int> value;  // Either a constant (double) or a variable index (int)
    std::vector<VExpr> children;  // Can store multiple children for n-ary operations
    VExpr* parent = nullptr;
    
    VExpr() : opcode(Opcode::CONSTANT), value(-1) {}

    VExpr(Opcode op) : opcode(op), value(-1) {
      assert((op != Opcode::CONSTANT) && (op != Opcode::VARIABLE));
    }
    VExpr(Opcode op, VExpr arg) : opcode(op), value(-1) {
      AddArgument(arg);
    }
    private:
    VExpr(Opcode op, std::variant<double, int> val) : opcode(op), value(val) {}
    public:
    static VExpr makeVariable(int index) {
      return VExpr(Opcode::VARIABLE, index);
    }
    static VExpr makeConstant(double value) {
      return VExpr(Opcode::CONSTANT, value);
    }
    // Add an argument (child node) to this expression
    void AddArgument(VExpr child) {
      child.parent = this;
      children.push_back(child);
    }
    void append(fmt::MemoryWriter &w, bool endl=true) const {
            assert(opcode >= Opcode::CONSTANT && opcode <= Opcode::LOG);
      if (opcode >= Opcode::EXP)
        w << opcodeStrings[(int)opcode];
      if (!children.empty()) {
        w << "(";
        for (size_t i = 0; i < children.size(); ++i) {
          children[i].append(w);
          if (i < children.size() - 1) {
            switch (opcode) {
            case Opcode::ADD: w<< " + "; break;
            case Opcode::SUBTRACT: w << " - "; break;
            case Opcode::MULTIPLY: w << " * "; break;
            case Opcode::DIVIDE: w << " / "; break;
            case Opcode::POW: w << "^"; break;
            case Opcode::LOGA: w << ","; break;
            default: break;
            }
          }
        }
          w << ")";
      }
      else {
        if (opcode == Opcode::CONSTANT) {
          w << std::get<double>(value);
        }
        else if (opcode == Opcode::VARIABLE) {
          w << varName(std::get<int>(value));
        }

      }
        if (endl && (parent == nullptr))
          w << "\n";
    }
    // Print the expression in-order
    std::string ToString(bool endl=true) const {
      
          fmt::MemoryWriter w;
          append(w, endl);
          return w.str();
    }
    
  };




/// BaronmpModelAPI.
/// @note For expression tree solvers,
///   see existing drivers, such as scipmp, gurobi and mp2nl.
class BaronmpModelAPI :
    public BaronmpCommon, public EnvKeeper,
  public BasicExprModelAPI<BaronmpModelAPI, VExpr> {
  /// Typedef main base class
  using BaseModelAPI = BasicExprModelAPI<BaronmpModelAPI, VExpr>;

  int n_unamed_constraint = 0;
  // Store the constraints to add to the file later
  // Can improve memory usage by writing two separate files then appending
  std::vector< std::string> cons;
  std::string obj;
  fmt::MemoryWriter vars_buffer;
  
  void addLinear(const std::string& name,
          double lhs, double rhs,
    size_t nvars, const int* vars, const double* coeffs);
  template <typename Args, typename Params, typename NumOrLogic, typename Id>
  void addFunctionalConstraint(const std::string& func,
    const CustomFunctionalConstraint<Args, Params, NumOrLogic, Id>& c);
  void addQuadratic(const std::string& name,
    double lhs, double rhs, const mp::LinTerms& lt,
    const mp::QuadTerms &qt);
  template <int SENSE> void addTopLevel(const NLBaseAssign<SENSE>& c);

  enum VTYPE {
    BIN = 0,
    INT = 1,
    POS = 2,
    FREE = 3
  };
public:
  std::string varName(int index);
  std::string createConName(const std::string & name);

  /// Construct
  BaronmpModelAPI(Env& e) : EnvKeeper(e) { }

  /// Class name
  static const char* GetTypeName() { return "BaronmpModelAPI"; }

  /// If any driver options added from here
  void InitCustomOptions() { }

  /// Called before problem modification.
  /// @param fmi: current problem information.
  /// @note this is called before each phase of model modification
  ///   which can happen during iterative solving.
  void InitProblemModificationPhase(const FlatModelInfo* fmi);
  /// After
  void FinishProblemModificationPhase();

  void AddVariables(const std::vector<int>& indices,
    fmt::StringRef prefix, fmt::StringRef header,
    const char* const* names=nullptr);
  void AddBounds(const double* bounds,
    const std::vector<VTYPE>& vtypes, bool lower);
  void AddVariables(const VarArrayDef& );

  void SetLinearObjective( int iobj, const LinearObjective& lo );
  /// Whether accepting quadratic objectives:
  /// 0 - no, 1 - convex, 2 - nonconvex
  static int AcceptsQuadObj() { return 2; }
  void SetQuadraticObjective(int iobj, const QuadraticObjective& qo);

  /// Whether accepts NLObjective
  static int AcceptsNLObj() { return 1; }
  void SetNLObjective(int, const NLObjective&);

  //////////////////////////// GENERAL CONSTRAINTS ////////////////////////////
  /// Handle flat constraints: inherit basic API
  USE_BASE_CONSTRAINT_HANDLERS(BaseModelAPI)

    //////////////////////////// EXPRESSION TREES ////////////////////////////
    /// Hanlde expression trees: inherit basic API
    USE_BASE_EXPRESSION_HANDLERS(BaseModelAPI)

    /// ACCEPT_EXPRESSION_INTERFACE():
    /// 'NotAccepted' or
    /// 'AcceptedButNotRecommended' would resort to flat constraints uniformly
    /// which outline each expression with an auxiliary variable:
    /// e.g., ExpConstraint(a, b), meaning a = exp(b), with a, b variables.
    /// But the latter option would enable the solver option acc:_expr,
    /// which, when set to 1, switches on the expression trees.
    ///
    /// See also per-expression and per-constraint type switches
    /// (ACCEPT_EXPRESSION and ACCEPT_CONSTRAINT.)
    /// When both are enabled for certain function
    /// (e.g., CosExpression and CosConstraint),
    /// solver option acc:cos allows switching between them.
    ///
    /// Expression trees are submitted to the solver via
    /// the high-level constraints
    /// NLConstraint, NLAssignLE, NLAssignEQ, NLAssignGE,
    /// NLComplementarity,
    /// NLLogical, NLEquivalence, NLImpl, NLRimpl, and NLObjective.
    ACCEPT_EXPRESSION_INTERFACE(Recommended);

  ACCEPT_CONSTRAINT(LinConRange, Recommended, CG_Linear)
  void AddConstraint(const LinConRange& lc);
  ACCEPT_CONSTRAINT(LinConLE, Recommended, CG_Linear)
  void AddConstraint(const LinConLE& lc);
  ACCEPT_CONSTRAINT(LinConEQ, Recommended, CG_Linear)
  void AddConstraint(const LinConEQ& lc);
  ACCEPT_CONSTRAINT(LinConGE, Recommended, CG_Linear)
  void AddConstraint(const LinConGE& lc);

  /// Ask if the solver accepts non-convex quadratic constraints
  static constexpr bool AcceptsNonconvexQC() { return true; }

  /// If cvt:prod=7 (and not 5) default.
  /// Recommendation to return the opposite value as
  /// AcceptsNonconvexQC().
  static constexpr bool WantLogicalizedProd2Bin()
  { return !AcceptsNonconvexQC(); }

  /// QuadConRange is optional.
  ACCEPT_CONSTRAINT(QuadConRange, Recommended, CG_Quadratic)
  void AddConstraint(const QuadConRange& qc);

  /// If using quadratics,
  /// QuadCon(LE/EQ/GE) should have 'Recommended'
  /// and have an implementation.
  ACCEPT_CONSTRAINT(QuadConLE, Recommended, CG_Quadratic)
  void AddConstraint(const QuadConLE& qc);
  ACCEPT_CONSTRAINT(QuadConEQ, Recommended, CG_Quadratic)
  void AddConstraint(const QuadConEQ& qc);
  ACCEPT_CONSTRAINT(QuadConGE, Recommended, CG_Quadratic)
  void AddConstraint(const QuadConGE& qc);

  /// SOCP flags
  ///////////////////////////////////////////////////////////////////////
  /// Ask if the solver can mix conic quadratic
  /// (entered via dedicated API) and direct quadratic constraints
  static constexpr bool CanMixConicQCAndQC() { return false; }

  /// Ask if the solver can recognize SOCP corner cases
  /// (non-std representations such as xy>=1, see tests)
  /// from quadratic representations
  static constexpr bool CanSOCPCornerCasesFromQC() { return false; }


  /// Some non linear constraints.
  /// See constr_std.h for more.
  //ACCEPT_CONSTRAINT(AbsConstraint, Recommended, CG_General)
 //   void AddConstraint(const AbsConstraint& absc);
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

  /// GetVarExpression(\a i): expression representing variable 0<=i<n_var.
  /// Only called for 'nonlinear' variables.
  Expr GetVarExpression(int i) { return  VExpr::makeVariable(i); }
  /// Can be used to represent empty expression in an NLConstraint.
  Expr GetZeroExpression() { return  VExpr::makeConstant(0); }

  /// Top level general NL constraint
  ACCEPT_CONSTRAINT(NLConstraint, Recommended, CG_Algebraic)
  void AddConstraint(const NLConstraint& nlc);
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
  /// The uequality/inqeuality type of the flat constraint is
  /// determied by GetContext().
  ///
  /// @note Use accessor: GetArgExpression(ee, 0)
  /// - don't AbsExpression's methods.
  /// @note Use accessor: GetParameter(pe, 0)
  /////   - don't use PowConstExpExpression's methods.
  /// Similar for other expression types.


  ACCEPT_EXPRESSION(DivExpression, Recommended)
  Expr AddExpression(const DivExpression&);
  ACCEPT_EXPRESSION(ExpAExpression, Recommended)
  Expr AddExpression(const ExpAExpression&);
  ACCEPT_EXPRESSION(ExpExpression, Recommended)
  Expr AddExpression(const ExpExpression&);
  ACCEPT_EXPRESSION(LogExpression, Recommended)
  Expr AddExpression(const LogExpression&);
  ACCEPT_EXPRESSION(LogAExpression, Recommended)
    Expr AddExpression(const LogAExpression&);

  ACCEPT_EXPRESSION(PowConstExpExpression, Recommended)
  Expr AddExpression(const PowConstExpExpression&);

  
protected:


  template <class MPExpr>
  void AppendLinAndConstTerms(Expr&, const MPExpr&);

  template <class MPExpr>
  void AppendQuadTerms(Expr&, const MPExpr&);

  private:
    void AddObjectiveHeader(fmt::MemoryWriter& w, int sense);
};

} // namespace mp

#endif // BARONMPMODELAPI_H
