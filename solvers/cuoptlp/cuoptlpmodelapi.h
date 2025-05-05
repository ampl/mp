#ifndef CUOPTLPMODELAPI_H
#define CUOPTLPMODELAPI_H

#include "mp/env.h"
#include "cuoptlpcommon.h"
// Include the following if you don't need expressions support
// #include "mp/flat/model_api_base.h"
// else include this:
#include "mp/flat/nl_expr/model_api_base.h"

namespace mp {


/// CuoptlpModelAPI.
/// @note For expression tree solvers,
///       inherit from public `BasicExprModelAPI<ModelAPIclass, expressionclass>`
///       (see below),
///       if you don't need expressions inherit from `BasicFlatModelAPI`
class CuoptlpModelAPI :
    public CuoptlpCommon, public EnvKeeper,
    public BasicFlatModelAPI
{

  /// If you don't need expression support, use this typedef
  // using BaseModelAPI = BasicFlatModelAPI;
  /// else, typedef main base class for expressions support
  using BaseModelAPI = BasicFlatModelAPI;


public:
  /// Construct
  CuoptlpModelAPI(Env& e) : EnvKeeper(e) { }

  /// Class name
  static const char* GetTypeName() { return "CuoptlpModelAPI"; }

  /// If any driver options added from here
  void InitCustomOptions() { }

  /// Called before problem modification.
  /// @param fmi: current problem information.
  /// @note this is called before each phase of model modification
  ///   which can happen during iterative solving.
  void InitProblemModificationPhase(const FlatModelInfo* fmi);
  /// After
  void FinishProblemModificationPhase();

  /// TODO Implement adding variables
  void AddVariables(const VarArrayDef& );
  /// TODO Implement setting (also changing) a linear (part of the) objective
  void SetLinearObjective( int iobj, const LinearObjective& lo );
  /// Whether accepting quadratic objectives:
  /// 0 - no, 1 - convex, 2 - nonconvex
  static int AcceptsQuadObj() { return 0; }
  /// TODO Implement setting (also changing) a quadratic objective
  void SetQuadraticObjective(int iobj, const QuadraticObjective& qo);

  //////////////////////////// GENERAL CONSTRAINTS ////////////////////////////
  USE_BASE_CONSTRAINT_HANDLERS(BaseModelAPI)

  /// TODO For each suppoted constraint type, add the ACCEPT_CONSTRAINT macro
  /// and the relative AddConstraint function.
  /// Below some typical constraint handlers of a MIP solver.
  /// Further constraint types which could be handled natively by some solvers:
  /// - IndicatorConstraint(Lin/Quad)(LE/EQ/GE)
  /// - Multidirectional indicators Cond(Lin/Quad)Con(LT/LE/EQ/GE/GT), where
  ///   the implication direction (</==/>) depends in the context
  /// - Complementarity
  /// - Logical, counting, piecewise-linear constraints.
  /// See \a constr_std.h and other drivers.


  /// The linear range constraint, if fully supported with basis info etc.
  ACCEPT_CONSTRAINT(LinConRange, Recommended, CG_Linear)
  void AddConstraint(const LinConRange& lc);

  /// LinCon(LE/EQ/GE) should have 'Recommended' for all backends
  /// and have an implementation,
  /// or a conversion rule is needed in a derived FlatConverter
  ACCEPT_CONSTRAINT(LinConLE, Recommended, CG_Linear)
  void AddConstraint(const LinConLE& lc);
  ACCEPT_CONSTRAINT(LinConEQ, Recommended, CG_Linear)
  void AddConstraint(const LinConEQ& lc);
  ACCEPT_CONSTRAINT(LinConGE, Recommended, CG_Linear)
  void AddConstraint(const LinConGE& lc);

  private:
    // Helper function to fetch variable names from the solver API stub,
    // used to print the model variable names
    std::function<std::string_view(int)> GetVarName;

    // Helper function to append an affine expression (linear + constant terms)
    // coming from MP (nla) to an expression in terms of this driver (ff)
    //template <class MPExpr>
    //void AppendLinAndConstTerms(Expr& ff, const MPExpr& nla);

    //template <class MPExpr>
    //void AppendQuadTerms(Expr&, const MPExpr&);
};

} // namespace mp

#endif // CUOPTLPMODELAPI_H
