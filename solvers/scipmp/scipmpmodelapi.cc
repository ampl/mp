#include "scipmpmodelapi.h"


namespace mp {

void ScipModelAPI::InitProblemModificationPhase(const FlatModelInfo* flat_model_info) {
  if (getPROBDATA()->nlinconss || getPROBDATA()->nlconss.size())
    SCIP_CCALL( SCIPfreeTransform(getSCIP()) );                  // allow model update
  // Allocate storage if needed:
  int n_linear_cons = flat_model_info->GetNumberOfConstraintsOfGroup(CG_Linear);
  getPROBDATA()->nlinconss = n_linear_cons;
  getPROBDATA()->linconss.resize(n_linear_cons);
  int n_nl_cons = flat_model_info->GetNumberOfConstraintsOfGroup(CG_Nonlinear);
  getPROBDATA()->nlconss.reserve(n_nl_cons);
}

void ScipModelAPI::AddVariables(const VarArrayDef& v) {
  assert( SCIPgetStage(getSCIP()) == SCIP_STAGE_PROBLEM );

  getPROBDATA()->nvars = v.size();
  SCIP_CCALL( SCIPallocBlockMemoryArray(getSCIP(), &getPROBDATA()->vars, getPROBDATA()->nvars) );
  getPROBDATA()->var_exprs.resize(v.size());

  for (int i = 0; i < v.size(); i++) {
    SCIP_VAR* var = NULL;
    double lb = std::isinf(v.plb()[i]) ? MinusInfinity() : v.plb()[i];
    double ub = std::isinf(v.pub()[i]) ? Infinity() : v.pub()[i];
    double objcoef = 0.0;
    SCIP_VARTYPE vartype;
    if (v.ptype()[i] == var::Type::INTEGER && lb == 0.0 && ub == 1.0)
      vartype = SCIP_VARTYPE_BINARY;
    else if (v.ptype()[i] == var::Type::INTEGER)
      vartype = SCIP_VARTYPE_INTEGER;
    else
      vartype = SCIP_VARTYPE_CONTINUOUS;
    if (v.pnames() != NULL)
      SCIP_CCALL( SCIPcreateVarBasic(getSCIP(), &var, v.pnames()[i], lb, ub, objcoef, vartype) );
    else
      SCIP_CCALL( SCIPcreateVarBasic(getSCIP(), &var, NULL, lb, ub, objcoef, vartype) );
    SCIP_CCALL( SCIPaddVar(getSCIP(), var) );
    getPROBDATA()->vars[i] = var;
  }
}

void ScipModelAPI::SetLinearObjective( int iobj, const LinearObjective& lo ) {
  if (iobj<1) {
    SCIP_CCALL( SCIPsetObjsense(getSCIP(),
                    obj::Type::MAX==lo.obj_sense() ? SCIP_OBJSENSE_MAXIMIZE : SCIP_OBJSENSE_MINIMIZE) );
    SCIP_VAR** vars = getPROBDATA()->vars;
    for (int i = 0; i < getPROBDATA()->nvars; i++) {
      SCIP_CCALL( SCIPchgVarObj(getSCIP(), vars[i], 0.0) );          // zero out
    }
    // @todo If becomes relevant: zero out QP/nonlinear obj terms
    for (int i = 0; i < lo.num_terms(); i++) {
      SCIP_CCALL( SCIPchgVarObj(getSCIP(), vars[lo.vars()[i]], lo.coefs()[i]) );
    }
  } 
  else {
    throw std::runtime_error("Multiple linear objectives not supported");
  }
}


void ScipModelAPI::SetQuadraticObjective(int iobj, const QuadraticObjective& qo) {
  throw std::runtime_error("Quadratic objective not supported");
}


SCIP_EXPR* ScipModelAPI::GetVarExpression(int i) {
  assert(getPROBDATA()->var_exprs.size()>i);
  if (!getPROBDATA()->var_exprs[i])
    SCIP_CCALL( SCIPcreateExprVar(
        getSCIP(), &getPROBDATA()->var_exprs[i], getPROBDATA()->vars[i], NULL, NULL) );
  return getPROBDATA()->var_exprs[i];
}

SCIP_EXPR* ScipModelAPI::GetZeroExpression() {
  if (!getPROBDATA()->dummyexpr)
    SCIP_CCALL( SCIPcreateExprValue(getSCIP(), &getPROBDATA()->dummyexpr, 0.0, NULL, NULL) );
  return getPROBDATA()->dummyexpr;
}

void ScipModelAPI::linearHelper(const int* pvars, const double* pcoefs, const size_t size, const char* name, const double lb, const double ub) {
  SCIP_VAR** vars = NULL;
  SCIP_CCALL( SCIPallocBufferArray(getSCIP(), &vars, size) );
  for (size_t i = 0; i < size; i++)
    vars[i] = getPROBDATA()->vars[pvars[i]];

  SCIP_CONS* cons;
  SCIP_CCALL( SCIPcreateConsBasicLinear(getSCIP(), &cons, name, (int)size, vars, (double*)pcoefs, lb, ub) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), cons) );
  getPROBDATA()->linconss[getPROBDATA()->i] = cons;
  getPROBDATA()->i++;

  SCIPfreeBufferArray(getSCIP(), &vars);
}

void ScipModelAPI::AddConstraint(const LinConRange& lc) {
  linearHelper(lc.pvars(), lc.pcoefs(), lc.size(), lc.GetName(), lc.lb(), lc.ub());
}

void ScipModelAPI::AddConstraint(const LinConLE& lc) {
  linearHelper(lc.pvars(), lc.pcoefs(), lc.size(), lc.GetName(), MinusInfinity(), lc.ub());
}

void ScipModelAPI::AddConstraint(const LinConEQ& lc) {
  linearHelper(lc.pvars(), lc.pcoefs(), lc.size(), lc.GetName(), lc.ub(), lc.ub());
}

void ScipModelAPI::AddConstraint(const LinConGE& lc) {
  linearHelper(lc.pvars(), lc.pcoefs(), lc.size(), lc.GetName(), lc.lb(), Infinity());
}


SCIP_EXPR* ScipModelAPI::AddExpression(const AbsExpression &abse) {
  SCIP_EXPR* result_expr;
  SCIP_CCALL( SCIPcreateExprAbs(getSCIP(), &result_expr,
                               GetArgExpression(abse, 0), NULL, NULL) );
  getPROBDATA()->exprs.push_back(result_expr);
  return result_expr;
}


void ScipModelAPI::AddConstraint(const AbsConstraint &absc)  {
  SCIP_VAR* x = getPROBDATA()->vars[absc.GetArguments()[0]];
  SCIP_VAR* res = getPROBDATA()->vars[absc.GetResultVar()];

  assert(x != NULL);
  assert(res != NULL);

  SCIP_EXPR* terms[2];
  SCIP_EXPR* xexpr;
  SCIP_EXPR* sumexpr;
  SCIP_Real coefs[2] = {1.0, -1.0};

  SCIP_CCALL( SCIPcreateExprVar(getSCIP(), &xexpr, x, NULL, NULL) );
  SCIP_CCALL( SCIPcreateExprVar(getSCIP(), &terms[0], res, NULL, NULL) );
  SCIP_CCALL( SCIPcreateExprAbs(getSCIP(), &terms[1], xexpr, NULL, NULL) );
  SCIP_CCALL( SCIPcreateExprSum(getSCIP(), &sumexpr, 2, terms, coefs, 0.0, NULL, NULL) );

  SCIP_CONS* cons;
  SCIP_CCALL( SCIPcreateConsBasicNonlinear(getSCIP(), &cons, absc.GetName(), sumexpr, 0.0, 0.0) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), cons) );
  SCIP_CCALL( SCIPreleaseCons(getSCIP(), &cons) );

  SCIP_CCALL( SCIPreleaseExpr(getSCIP(), &sumexpr) );
  SCIP_CCALL( SCIPreleaseExpr(getSCIP(), &terms[1]) );
  SCIP_CCALL( SCIPreleaseExpr(getSCIP(), &terms[0]) );
  SCIP_CCALL( SCIPreleaseExpr(getSCIP(), &xexpr) );
}

void ScipModelAPI::AddConstraint(const AndConstraint &cc)  {
  SCIP_VAR** vars = NULL;
  SCIP_CCALL( SCIPallocBufferArray(getSCIP(), &vars, cc.GetArguments().size()) );
	SCIP_Bool infeas = false;
	SCIP_CCALL( SCIPchgVarType(getSCIP(), getPROBDATA()->vars[cc.GetResultVar()],
														 SCIP_VARTYPE_BINARY, &infeas) );
	for (size_t i = 0; i < cc.GetArguments().size(); i++) {
    vars[i] = getPROBDATA()->vars[cc.GetArguments()[i]];
		SCIP_CCALL( SCIPchgVarType(getSCIP(), vars[i], SCIP_VARTYPE_BINARY, &infeas) );
	}

  SCIP_CONS* cons;
  SCIP_CCALL( SCIPcreateConsBasicAnd(getSCIP(), &cons, cc.GetName(), getPROBDATA()->vars[cc.GetResultVar()],
    (int)cc.GetArguments().size(), vars) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), cons) );
  SCIP_CCALL( SCIPreleaseCons(getSCIP(), &cons) );

  SCIPfreeBufferArray(getSCIP(), &vars);
}

void ScipModelAPI::AddConstraint(const OrConstraint &dc)  {
  SCIP_VAR** vars = NULL;
  SCIP_CCALL( SCIPallocBufferArray(getSCIP(), &vars, dc.GetArguments().size()) );
	SCIP_Bool infeas = false;
	SCIP_CCALL( SCIPchgVarType(getSCIP(), getPROBDATA()->vars[dc.GetResultVar()],
														 SCIP_VARTYPE_BINARY, &infeas) );
	for (size_t i = 0; i < dc.GetArguments().size(); i++) {
    vars[i] = getPROBDATA()->vars[dc.GetArguments()[i]];
		SCIP_CCALL( SCIPchgVarType(getSCIP(), vars[i], SCIP_VARTYPE_BINARY, &infeas) );
  }

  SCIP_CONS* cons;
  SCIP_CCALL( SCIPcreateConsBasicOr(getSCIP(), &cons, dc.GetName(), getPROBDATA()->vars[dc.GetResultVar()],
    (int)dc.GetArguments().size(), vars) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), cons) );
  SCIP_CCALL( SCIPreleaseCons(getSCIP(), &cons) );

  SCIPfreeBufferArray(getSCIP(), &vars);
}


void ScipModelAPI::helpIndicatorLin(const int* pvars, const double* pcoefs, const size_t size, const char* name, const double rhs, const int binary_var, const int binary_value, SCIP_Bool lessthanineq) {
  SCIP_VAR** vars = NULL;
  SCIP_CCALL( SCIPallocBufferArray(getSCIP(), &vars, size) );
  for (size_t i = 0; i < size; i++)
    vars[i] = getPROBDATA()->vars[pvars[i]];
  SCIP_CONS* cons;
  SCIP_CCALL( SCIPcreateConsIndicatorGeneric(getSCIP(), &cons, name, getPROBDATA()->vars[binary_var],
    (int)size, vars, (double*)pcoefs, rhs, binary_value, lessthanineq, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), cons) );
  SCIP_CCALL( SCIPreleaseCons(getSCIP(), &cons) );

  SCIPfreeBufferArray(getSCIP(), &vars);
}

void ScipModelAPI::AddConstraint(const IndicatorConstraintLinLE &ic)  {
  helpIndicatorLin(ic.get_constraint().pvars(), ic.get_constraint().pcoefs(), ic.get_constraint().size(), ic.GetName(), ic.get_constraint().rhs(), ic.get_binary_var(), ic.get_binary_value(), TRUE);
}

void ScipModelAPI::AddConstraint(const IndicatorConstraintLinEQ &ic)  {
  SCIP_VAR** vars = NULL;
  SCIP_CCALL( SCIPallocBufferArray(getSCIP(), &vars, ic.get_constraint().size()) );
  for (size_t i = 0; i < ic.get_constraint().size(); i++) {
    vars[i] = getPROBDATA()->vars[ic.get_constraint().pvars()[i]];
  }
  SCIP_CONS* cons1;
  SCIP_CONS* cons2;
  SCIP_CCALL( SCIPcreateConsIndicatorGeneric(getSCIP(), &cons1, ic.GetName(),getPROBDATA()->vars[ic.get_binary_var()],
    (int)ic.get_constraint().size(), vars, (double*)ic.get_constraint().pcoefs(), ic.get_constraint().rhs(),
    ic.get_binary_value(), TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
  SCIP_CCALL( SCIPcreateConsIndicatorGeneric(getSCIP(), &cons2, ic.GetName(),getPROBDATA()->vars[ic.get_binary_var()],
    (int)ic.get_constraint().size(), vars, (double*)ic.get_constraint().pcoefs(), ic.get_constraint().rhs(),
    ic.get_binary_value(), FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), cons1) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), cons2) );
  SCIP_CCALL( SCIPreleaseCons(getSCIP(), &cons1) );
  SCIP_CCALL( SCIPreleaseCons(getSCIP(), &cons2) );

  SCIPfreeBufferArray(getSCIP(), &vars);   
}

void ScipModelAPI::AddConstraint(const IndicatorConstraintLinGE &ic)  {
  helpIndicatorLin(ic.get_constraint().pvars(), ic.get_constraint().pcoefs(), ic.get_constraint().size(), ic.GetName(), ic.get_constraint().rhs(), ic.get_binary_var(), ic.get_binary_value(), FALSE);
}


void ScipModelAPI::helpQuad(const char* name, const int* lt_pvars, const double* lt_pcoefs, const size_t lt_size, const int* qt_pvars1, const int* qt_pvars2, const double* qt_pcoefs, const size_t qt_size, const double lb, const double ub) {
  // allocation of linterm vars
  SCIP_VAR** linvars = NULL;
  SCIP_CCALL( SCIPallocBufferArray(getSCIP(), &linvars, lt_size) );
  for (size_t i = 0; i < lt_size; i++)
    linvars[i] = getPROBDATA()->vars[lt_pvars[i]];

  // allocation of qterm1 vars
  SCIP_VAR** qvars1 = NULL;
  SCIP_CCALL( SCIPallocBufferArray(getSCIP(), &qvars1, qt_size) );
  for (size_t i = 0; i < qt_size; i++)
    qvars1[i] = getPROBDATA()->vars[qt_pvars1[i]];
  
  // allocation of qterm2 vars
  SCIP_VAR** qvars2 = NULL;
  SCIP_CCALL( SCIPallocBufferArray(getSCIP(), &qvars2, qt_size) );
  for (size_t i = 0; i < qt_size; i++)
    qvars2[i] = getPROBDATA()->vars[qt_pvars2[i]];
  

  SCIP_CONS* cons;
  SCIP_CCALL( SCIPcreateConsBasicQuadraticNonlinear(getSCIP(), &cons, name, (int)lt_size, linvars, (double*)lt_pcoefs, (int)qt_size, qvars1, qvars2, (double*)qt_pcoefs, lb, ub) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), cons) );
  SCIP_CCALL( SCIPreleaseCons(getSCIP(), &cons) );

  SCIPfreeBufferArray(getSCIP(), &linvars);
  SCIPfreeBufferArray(getSCIP(), &qvars1);
  SCIPfreeBufferArray(getSCIP(), &qvars2);
}

void ScipModelAPI::AddConstraint(const QuadConRange& qc) {
  const auto& lt = qc.GetLinTerms();
  const auto& qt = qc.GetQPTerms();

  helpQuad(qc.GetName(), lt.pvars(), lt.pcoefs(), lt.size(), qt.pvars1(), qt.pvars2(), qt.pcoefs(), qt.size(), qc.lb(), qc.ub());
}

void ScipModelAPI::AddConstraint( const QuadConLE& qc ) {
  const auto& lt = qc.GetLinTerms();
  const auto& qt = qc.GetQPTerms();

  helpQuad(qc.GetName(), lt.pvars(), lt.pcoefs(), lt.size(),
           qt.pvars1(), qt.pvars2(), qt.pcoefs(), qt.size(), MinusInfinity(), qc.rhs());
}

void ScipModelAPI::AddConstraint( const QuadConEQ& qc ) {
  const auto& lt = qc.GetLinTerms();
  const auto& qt = qc.GetQPTerms();

  helpQuad(qc.GetName(), lt.pvars(), lt.pcoefs(), lt.size(),
           qt.pvars1(), qt.pvars2(), qt.pcoefs(), qt.size(), qc.rhs(), qc.rhs());
}

void ScipModelAPI::AddConstraint( const QuadConGE& qc ) {
  const auto& lt = qc.GetLinTerms();
  const auto& qt = qc.GetQPTerms();

  helpQuad(qc.GetName(), lt.pvars(), lt.pcoefs(), lt.size(),
           qt.pvars1(), qt.pvars2(), qt.pcoefs(), qt.size(), qc.lb(), Infinity());
}

/// To access information from an NLConstraint,
/// use the following accessors (don't use methods of NLConstraint itself):
/// - GetLinSize(nlc), GetLinCoef(nlc, i), GetLinVar(nlc, i),
///   GetExpression(nlc), GetLower(nlc), GetUpper(nlc).
///
/// Implementation follows partly reader_nl.cc from SCIP.
void ScipModelAPI::AddConstraint( const NLConstraint& nlc ) {
  getPROBDATA()->nlconss.push_back(nullptr);

  SCIP_CCALL( SCIPcreateConsBasicNonlinear(
      getSCIP(), &getPROBDATA()->nlconss.back(), nlc.GetName(),
      GetExpression(nlc), GetLower(nlc), GetUpper(nlc)) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), getPROBDATA()->nlconss.back()) );
  for( int i = 0; i < GetLinSize(nlc); ++i )
    SCIP_CCALL( SCIPaddLinearVarNonlinear(
        getSCIP(), getPROBDATA()->nlconss.back(),
        getPROBDATA()->vars[ GetLinVar(nlc, i) ], GetLinCoef(nlc, i)) );
}

void ScipModelAPI::AddConstraint( const NLAssignEQ& nlae ) {
  getPROBDATA()->nlconss.push_back(nullptr);

  SCIP_CCALL( SCIPcreateConsBasicNonlinear(
      getSCIP(), &getPROBDATA()->nlconss.back(), nlae.GetName(),
      GetExpression(nlae), 0.0, 0.0) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), getPROBDATA()->nlconss.back()) );
  SCIP_CCALL( SCIPaddLinearVarNonlinear(
      getSCIP(), getPROBDATA()->nlconss.back(),
      getPROBDATA()->vars[ GetVariable(nlae) ], -1.0) );
}
void ScipModelAPI::AddConstraint( const NLAssignLE& nlae ) {
  getPROBDATA()->nlconss.push_back(nullptr);

  SCIP_CCALL( SCIPcreateConsBasicNonlinear(
      getSCIP(), &getPROBDATA()->nlconss.back(), nlae.GetName(),
      GetExpression(nlae), 0.0, Infinity()) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), getPROBDATA()->nlconss.back()) );
  SCIP_CCALL( SCIPaddLinearVarNonlinear(
      getSCIP(), getPROBDATA()->nlconss.back(),
      getPROBDATA()->vars[ GetVariable(nlae) ], -1.0) );
}
void ScipModelAPI::AddConstraint( const NLAssignGE& nlae ) {
  getPROBDATA()->nlconss.push_back(nullptr);

  SCIP_CCALL( SCIPcreateConsBasicNonlinear(
      getSCIP(), &getPROBDATA()->nlconss.back(), nlae.GetName(),
      GetExpression(nlae), MinusInfinity(), 0.0) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), getPROBDATA()->nlconss.back()) );
  SCIP_CCALL( SCIPaddLinearVarNonlinear(
      getSCIP(), getPROBDATA()->nlconss.back(),
      getPROBDATA()->vars[ GetVariable(nlae) ], -1.0) );
}

void ScipModelAPI::AddConstraint( const NLLogical& nll ) {
  getPROBDATA()->nlconss.push_back(nullptr);

  SCIP_CCALL( SCIPcreateConsBasicNonlinear(
      getSCIP(), &getPROBDATA()->nlconss.back(), nll.GetName(),
      GetExpression(nll), 1.0, 1.0) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), getPROBDATA()->nlconss.back()) );
}

void ScipModelAPI::AddConstraint( const NLReifEquiv& nll ) {
  getPROBDATA()->nlconss.push_back(nullptr);

  SCIP_CCALL( SCIPcreateConsBasicNonlinear(
      getSCIP(), &getPROBDATA()->nlconss.back(), nll.GetName(),
      GetExpression(nll), 0.0, 0.0) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), getPROBDATA()->nlconss.back()) );
  SCIP_CCALL( SCIPaddLinearVarNonlinear(
      getSCIP(), getPROBDATA()->nlconss.back(),
      getPROBDATA()->vars[ GetVariable(nll) ], -1.0) );
}
void ScipModelAPI::AddConstraint( const NLReifImpl& nll ) {
  getPROBDATA()->nlconss.push_back(nullptr);

  SCIP_CCALL( SCIPcreateConsBasicNonlinear(
      getSCIP(), &getPROBDATA()->nlconss.back(), nll.GetName(),
      GetExpression(nll), 0.0, Infinity()) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), getPROBDATA()->nlconss.back()) );
  SCIP_CCALL( SCIPaddLinearVarNonlinear(
      getSCIP(), getPROBDATA()->nlconss.back(),
      getPROBDATA()->vars[ GetVariable(nll) ], -1.0) );
}
void ScipModelAPI::AddConstraint( const NLReifRimpl& nll ) {
  getPROBDATA()->nlconss.push_back(nullptr);

  SCIP_CCALL( SCIPcreateConsBasicNonlinear(
      getSCIP(), &getPROBDATA()->nlconss.back(), nll.GetName(),
      GetExpression(nll), MinusInfinity(), 0.0) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), getPROBDATA()->nlconss.back()) );
  SCIP_CCALL( SCIPaddLinearVarNonlinear(
      getSCIP(), getPROBDATA()->nlconss.back(),
      getPROBDATA()->vars[ GetVariable(nll) ], -1.0) );
}

SCIP_EXPR* ScipModelAPI::AddExpression(const NLAffineExpression &le) {
  std::vector<SCIP_EXPR*> terms(GetLinSize(le));
  std::vector<double> coefs(GetLinSize(le));

  for (size_t i=0; i<terms.size(); ++i) {
    terms[i] = GetLinTerm(le, i);
    coefs[i] = GetLinCoef(le, i);
  }
  SCIP_EXPR* sumexpr;

  SCIP_CCALL( SCIPcreateExprSum(getSCIP(), &sumexpr,
                               terms.size(), terms.data(), coefs.data(),
                               GetConstTerm(le), NULL, NULL) );
  getPROBDATA()->exprs.push_back(sumexpr);
  return sumexpr;
}

SCIP_EXPR* ScipModelAPI::AddExpression(const NLQuadExpression &qe) {
  std::vector<SCIP_EXPR*> terms(GetLinSize(qe) + GetQuadSize(qe));
  std::vector<double> coefs(GetLinSize(qe) + GetQuadSize(qe));

  for (int i=0; i<GetLinSize(qe); ++i) {
    terms[i] = GetLinTerm(qe, i);
    coefs[i] = GetLinCoef(qe, i);
  }
  for (int i=0; i<GetQuadSize(qe); ++i) {
    int i_quad = i + GetLinSize(qe);
    SCIP_EXPR* mul_terms[] =
        { GetQuadTerm1(qe, i), GetQuadTerm2(qe, i) };
    SCIP_CCALL( SCIPcreateExprProduct(getSCIP(), &terms[i_quad],
                                     2, mul_terms, GetQuadCoef(qe, i), NULL, NULL) );
    coefs[i_quad] = 1.0;
  }
  SCIP_EXPR* sumexpr;

  SCIP_CCALL( SCIPcreateExprSum(getSCIP(), &sumexpr,
                               terms.size(), terms.data(), coefs.data(),
                               GetConstTerm(qe), NULL, NULL) );
  for (int i=0; i<GetQuadSize(qe); ++i) {
    int i_quad = i + GetLinSize(qe);
    SCIP_CCALL( SCIPreleaseExpr(getSCIP(), &terms[i_quad]) );
  }
  getPROBDATA()->exprs.push_back(sumexpr);
  return sumexpr;
}


void ScipModelAPI::AddConstraint( const QuadraticConeConstraint& qc ) {
  const auto& arg = qc.GetArguments();
  SCIP_VAR** vars = NULL;
  SCIP_CCALL( SCIPallocBufferArray(getSCIP(), &vars, arg.size()) );
  for (size_t i = 0; i < arg.size(); i++) {
    vars[i] = getPROBDATA()->vars[arg[i]];
  }

  std::vector<double> param;
  for (size_t i = 1; i < qc.GetParameters().size(); i++) {
    param.push_back(std::sqrt(qc.GetParameters()[i]));
  }

  SCIP_CONS* cons;
  SCIP_CCALL( SCIPcreateConsBasicSOCNonlinear(getSCIP(), &cons, qc.GetName(),
                                             (int)arg.size()-1, vars+1, (SCIP_Real*)param.data(), NULL,
                                             0.0, vars[0], qc.GetParameters().data()[0], 0.0) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), cons) );
  SCIP_CCALL( SCIPreleaseCons(getSCIP(), &cons) );

  SCIPfreeBufferArray(getSCIP(), &vars);
}

void ScipModelAPI::AddConstraint(const SOS1Constraint& sos) {
  SCIP_VAR** vars = NULL;
  double* weights = NULL;
  SCIP_CCALL( SCIPallocBufferArray(getSCIP(), &vars, sos.size()) );
  SCIP_CCALL( SCIPallocBufferArray(getSCIP(), &weights, sos.size()) );
  for (int i = 0; i < sos.size(); i++) {
    vars[i] = getPROBDATA()->vars[sos.get_vars().data()[i]];
    weights[i] = sos.get_weights().data()[i];
  }

  SCIP_CONS* cons = NULL;
  SCIP_CCALL( SCIPcreateConsBasicSOS1(getSCIP(), &cons, sos.GetName(),
                                     sos.size(), vars, weights) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), cons) );
  SCIP_CCALL( SCIPreleaseCons(getSCIP(), &cons) );

  SCIPfreeBufferArray(getSCIP(), &vars);
  SCIPfreeBufferArray(getSCIP(), &weights);
}

void ScipModelAPI::AddConstraint(const SOS2Constraint& sos) {
  SCIP_VAR** vars = NULL;
  double* weights = NULL;
  SCIP_CCALL( SCIPallocBufferArray(getSCIP(), &vars, sos.size()) );
  SCIP_CCALL( SCIPallocBufferArray(getSCIP(), &weights, sos.size()) );
  for (int i = 0; i < sos.size(); i++) {
    vars[i] = getPROBDATA()->vars[sos.get_vars().data()[i]];
    weights[i] = sos.get_weights().data()[i];
  }

  SCIP_CONS* cons = NULL;
  SCIP_CCALL( SCIPcreateConsBasicSOS2(getSCIP(), &cons, sos.GetName(),
                                     sos.size(), vars, weights) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), cons) );
  SCIP_CCALL( SCIPreleaseCons(getSCIP(), &cons) );

  SCIPfreeBufferArray(getSCIP(), &vars);
  SCIPfreeBufferArray(getSCIP(), &weights);
}

SCIP_EXPR* ScipModelAPI::AddExpression(const ExpExpression &ee) {
  SCIP_EXPR* result_expr;
  SCIP_CCALL( SCIPcreateExprExp(getSCIP(), &result_expr,
                               GetArgExpression(ee, 0), NULL, NULL) );
  getPROBDATA()->exprs.push_back(result_expr);
  return result_expr;
}

void ScipModelAPI::AddConstraint(const ExpConstraint &cc)  {
  SCIP_VAR* x = getPROBDATA()->vars[cc.GetArguments()[0]];
  SCIP_VAR* res = getPROBDATA()->vars[cc.GetResultVar()];

  assert(x != NULL);
  assert(res != NULL);

  SCIP_EXPR* terms[2];
  SCIP_EXPR* xexpr;
  SCIP_EXPR* sumexpr;
  SCIP_Real coefs[2] = {1.0, -1.0};

  SCIP_CCALL( SCIPcreateExprVar(getSCIP(), &xexpr, x, NULL, NULL) );
  SCIP_CCALL( SCIPcreateExprVar(getSCIP(), &terms[0], res, NULL, NULL) );
  SCIP_CCALL( SCIPcreateExprExp(getSCIP(), &terms[1], xexpr, NULL, NULL) );
  SCIP_CCALL( SCIPcreateExprSum(getSCIP(), &sumexpr, 2, terms, coefs, 0.0, NULL, NULL) );

  SCIP_CONS* cons;
  SCIP_CCALL( SCIPcreateConsBasicNonlinear(getSCIP(), &cons, cc.GetName(), sumexpr, 0.0, 0.0) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), cons) );
  SCIP_CCALL( SCIPreleaseCons(getSCIP(), &cons) );

  SCIP_CCALL( SCIPreleaseExpr(getSCIP(), &sumexpr) );
  SCIP_CCALL( SCIPreleaseExpr(getSCIP(), &terms[1]) );
  SCIP_CCALL( SCIPreleaseExpr(getSCIP(), &terms[0]) );
  SCIP_CCALL( SCIPreleaseExpr(getSCIP(), &xexpr) );
}

SCIP_EXPR* ScipModelAPI::AddExpression(const LogExpression &ee) {
  SCIP_EXPR* result_expr;
  SCIP_CCALL( SCIPcreateExprLog(getSCIP(), &result_expr,
                               GetArgExpression(ee, 0), NULL, NULL) );
  getPROBDATA()->exprs.push_back(result_expr);
  return result_expr;
}
SCIP_EXPR* ScipModelAPI::AddExpression(const LogAExpression& ee) {
  SCIP_EXPR* log_expr, * mult_expr;
  SCIP_CCALL(SCIPcreateExprLog(getSCIP(), &log_expr,
    GetArgExpression(ee, 0), NULL, NULL));
  
  SCIP_CCALL(SCIPcreateExprProduct(getSCIP(), &mult_expr,
    1, &log_expr, 1 / std::log(GetParameter(ee, 0)), NULL, NULL));

  getPROBDATA()->exprs.push_back(log_expr);
  getPROBDATA()->exprs.push_back(mult_expr);
  return mult_expr;
}

void ScipModelAPI::AddConstraint(const LogAConstraint& cc) {
  SCIP_VAR* x = getPROBDATA()->vars[cc.GetArguments()[0]];
  SCIP_VAR* res = getPROBDATA()->vars[cc.GetResultVar()];

  assert(x != NULL);
  assert(res != NULL);

  SCIP_EXPR* xexpr;      // x
  SCIP_EXPR* logexpr; // ln(x)
  SCIP_EXPR* terms[2];   // [0] = y   [1] = ln(x)/ln(a) 
  SCIP_EXPR* sumexpr;    // y - ln(x)/ln(a)

  SCIP_Real coefs[2] = { 1.0, -1.0 };

  SCIP_CCALL(SCIPcreateExprVar(getSCIP(), &xexpr, x, NULL, NULL));
  SCIP_CCALL(SCIPcreateExprVar(getSCIP(), &terms[0], res, NULL, NULL));
  SCIP_CCALL(SCIPcreateExprLog(getSCIP(), &logexpr, xexpr, NULL, NULL));
  SCIP_CCALL(SCIPcreateExprProduct(getSCIP(), &terms[1], 1, &logexpr, std::log(cc.GetParameters()[0]), NULL, NULL));
  SCIP_CCALL(SCIPcreateExprSum(getSCIP(), &sumexpr, 2, terms, coefs, 0.0, NULL, NULL));

  SCIP_CONS* cons;
  SCIP_CCALL(SCIPcreateConsBasicNonlinear(getSCIP(), &cons, cc.GetName(), sumexpr, 0.0, 0.0));
  SCIP_CCALL(SCIPaddCons(getSCIP(), cons));
  SCIP_CCALL(SCIPreleaseCons(getSCIP(), &cons));

  SCIP_CCALL(SCIPreleaseExpr(getSCIP(), &sumexpr));
  SCIP_CCALL(SCIPreleaseExpr(getSCIP(), &terms[1]));
  SCIP_CCALL(SCIPreleaseExpr(getSCIP(), &terms[0]));
  SCIP_CCALL(SCIPreleaseExpr(getSCIP(), &xexpr));
}
void ScipModelAPI::AddConstraint(const LogConstraint &cc)  {
  SCIP_VAR* x = getPROBDATA()->vars[cc.GetArguments()[0]];
  SCIP_VAR* res = getPROBDATA()->vars[cc.GetResultVar()];

  assert(x != NULL);
  assert(res != NULL);

  SCIP_EXPR* terms[2];
  SCIP_EXPR* xexpr;
  SCIP_EXPR* sumexpr;
  SCIP_Real coefs[2] = {1.0, -1.0};

  SCIP_CCALL( SCIPcreateExprVar(getSCIP(), &xexpr, x, NULL, NULL) );
  SCIP_CCALL( SCIPcreateExprVar(getSCIP(), &terms[0], res, NULL, NULL) );
  SCIP_CCALL( SCIPcreateExprLog(getSCIP(), &terms[1], xexpr, NULL, NULL) );
  SCIP_CCALL( SCIPcreateExprSum(getSCIP(), &sumexpr, 2, terms, coefs, 0.0, NULL, NULL) );

  SCIP_CONS* cons;
  SCIP_CCALL( SCIPcreateConsBasicNonlinear(getSCIP(), &cons, cc.GetName(), sumexpr, 0.0, 0.0) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), cons) );
  SCIP_CCALL( SCIPreleaseCons(getSCIP(), &cons) );

  SCIP_CCALL( SCIPreleaseExpr(getSCIP(), &sumexpr) );
  SCIP_CCALL( SCIPreleaseExpr(getSCIP(), &terms[1]) );
  SCIP_CCALL( SCIPreleaseExpr(getSCIP(), &terms[0]) );
  SCIP_CCALL( SCIPreleaseExpr(getSCIP(), &xexpr) );
}

SCIP_EXPR* ScipModelAPI::AddExpression(const PowConstExpExpression &ee) {
  SCIP_EXPR* result_expr;
  SCIP_CCALL( SCIPcreateExprPow(getSCIP(), &result_expr,
                               GetArgExpression(ee, 0),
                               GetParameter(ee, 0),
                               NULL, NULL) );
  getPROBDATA()->exprs.push_back(result_expr);
  return result_expr;
}

void ScipModelAPI::AddConstraint(const PowConstExpConstraint &cc)  {
  SCIP_VAR* x = getPROBDATA()->vars[cc.GetArguments()[0]];
  SCIP_VAR* res = getPROBDATA()->vars[cc.GetResultVar()];

  assert(x != NULL);
  assert(res != NULL);

  SCIP_EXPR* terms[2];
  SCIP_EXPR* xexpr;
  SCIP_EXPR* sumexpr;
  SCIP_Real coefs[2] = {1.0, -1.0};

  SCIP_CCALL( SCIPcreateExprVar(getSCIP(), &xexpr, x, NULL, NULL) );
  SCIP_CCALL( SCIPcreateExprVar(getSCIP(), &terms[0], res, NULL, NULL) );
  SCIP_CCALL( SCIPcreateExprPow(getSCIP(), &terms[1], xexpr, cc.GetParameters()[0], NULL, NULL) );
  SCIP_CCALL( SCIPcreateExprSum(getSCIP(), &sumexpr, 2, terms, coefs, 0.0, NULL, NULL) );

  SCIP_CONS* cons;
  SCIP_CCALL( SCIPcreateConsBasicNonlinear(getSCIP(), &cons, cc.GetName(), sumexpr, 0.0, 0.0) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), cons) );
  SCIP_CCALL( SCIPreleaseCons(getSCIP(), &cons) );

  SCIP_CCALL( SCIPreleaseExpr(getSCIP(), &sumexpr) );
  SCIP_CCALL( SCIPreleaseExpr(getSCIP(), &terms[1]) );
  SCIP_CCALL( SCIPreleaseExpr(getSCIP(), &terms[0]) );
  SCIP_CCALL( SCIPreleaseExpr(getSCIP(), &xexpr) );
}

SCIP_EXPR* ScipModelAPI::AddExpression(const SinExpression &ee) {
  SCIP_EXPR* result_expr;
  SCIP_CCALL( SCIPcreateExprSin(getSCIP(), &result_expr,
                               GetArgExpression(ee, 0), NULL, NULL) );
  getPROBDATA()->exprs.push_back(result_expr);
  return result_expr;
}

void ScipModelAPI::AddConstraint(const SinConstraint &cc)  {
  SCIP_VAR* x = getPROBDATA()->vars[cc.GetArguments()[0]];
  SCIP_VAR* res = getPROBDATA()->vars[cc.GetResultVar()];

  assert(x != NULL);
  assert(res != NULL);

  SCIP_EXPR* terms[2];
  SCIP_EXPR* xexpr;
  SCIP_EXPR* sumexpr;
  SCIP_Real coefs[2] = {1.0, -1.0};

  SCIP_CCALL( SCIPcreateExprVar(getSCIP(), &xexpr, x, NULL, NULL) );
  SCIP_CCALL( SCIPcreateExprVar(getSCIP(), &terms[0], res, NULL, NULL) );
  SCIP_CCALL( SCIPcreateExprSin(getSCIP(), &terms[1], xexpr, NULL, NULL) );
  SCIP_CCALL( SCIPcreateExprSum(getSCIP(), &sumexpr, 2, terms, coefs, 0.0, NULL, NULL) );

  SCIP_CONS* cons;
  SCIP_CCALL( SCIPcreateConsBasicNonlinear(getSCIP(), &cons, cc.GetName(), sumexpr, 0.0, 0.0) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), cons) );
  SCIP_CCALL( SCIPreleaseCons(getSCIP(), &cons) );

  SCIP_CCALL( SCIPreleaseExpr(getSCIP(), &sumexpr) );
  SCIP_CCALL( SCIPreleaseExpr(getSCIP(), &terms[1]) );
  SCIP_CCALL( SCIPreleaseExpr(getSCIP(), &terms[0]) );
  SCIP_CCALL( SCIPreleaseExpr(getSCIP(), &xexpr) );
}

SCIP_EXPR* ScipModelAPI::AddExpression(const CosExpression &ee) {
  SCIP_EXPR* result_expr;
  SCIP_CCALL( SCIPcreateExprCos(getSCIP(), &result_expr,
                               GetArgExpression(ee, 0), NULL, NULL) );
  getPROBDATA()->exprs.push_back(result_expr);
  return result_expr;
}

void ScipModelAPI::AddConstraint(const CosConstraint &cc)  {
  SCIP_VAR* x = getPROBDATA()->vars[cc.GetArguments()[0]];
  SCIP_VAR* res = getPROBDATA()->vars[cc.GetResultVar()];

  assert(x != NULL);
  assert(res != NULL);

  SCIP_EXPR* terms[2];
  SCIP_EXPR* xexpr;
  SCIP_EXPR* sumexpr;
  SCIP_Real coefs[2] = {1.0, -1.0};

  SCIP_CCALL( SCIPcreateExprVar(getSCIP(), &xexpr, x, NULL, NULL) );
  SCIP_CCALL( SCIPcreateExprVar(getSCIP(), &terms[0], res, NULL, NULL) );
  SCIP_CCALL( SCIPcreateExprCos(getSCIP(), &terms[1], xexpr, NULL, NULL) );
  SCIP_CCALL( SCIPcreateExprSum(getSCIP(), &sumexpr, 2, terms, coefs, 0.0, NULL, NULL) );

  SCIP_CONS* cons;
  SCIP_CCALL( SCIPcreateConsBasicNonlinear(getSCIP(), &cons, cc.GetName(), sumexpr, 0.0, 0.0) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), cons) );
  SCIP_CCALL( SCIPreleaseCons(getSCIP(), &cons) );

  SCIP_CCALL( SCIPreleaseExpr(getSCIP(), &sumexpr) );
  SCIP_CCALL( SCIPreleaseExpr(getSCIP(), &terms[1]) );
  SCIP_CCALL( SCIPreleaseExpr(getSCIP(), &terms[0]) );
  SCIP_CCALL( SCIPreleaseExpr(getSCIP(), &xexpr) );
}

void ScipModelAPI::FinishProblemModificationPhase() {
}


} // namespace mp
