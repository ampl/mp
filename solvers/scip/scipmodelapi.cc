#include "scipmodelapi.h"

namespace mp {

void ScipModelAPI::InitProblemModificationPhase(
    const FlatModelInfo*) {
  // Allocate storage if needed:
  // auto n_linear_cons =
  //   flat_model_info->GetNumberOfConstraintsOfGroup(CG_LINEAR);
  // preallocate_linear_cons( n_linear_cons );
}

void ScipModelAPI::AddVariables(const VarArrayDef& v) {
  assert( SCIPgetStage(getSCIP()) == SCIP_STAGE_PROBLEM );

  getPROBDATA()->nvars = v.size();
  SCIP_CCALL( SCIPallocBlockMemoryArray(getSCIP(), &getPROBDATA()->vars, getPROBDATA()->nvars) );

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
    //const char* name = v.pnames()[i];
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
    for (int i = 0; i < lo.num_terms(); i++) {
      SCIP_CCALL( SCIPchgVarObj(getSCIP(), vars[lo.vars()[i]], lo.coefs()[i]) );
    }
  } 
  else {
    throw std::runtime_error("Multiple linear objectives not supported");
  }
}


void ScipModelAPI::SetQuadraticObjective(int iobj, const QuadraticObjective& qo) {
  if (1 > iobj) {
    fmt::format("Setting first quadratic objective\n");
    SetLinearObjective(iobj, qo);                         // add the linear part
    const auto& qt = qo.GetQPTerms();
    fmt::format("Quadratic part is made of {} terms\n", qt.size());

    // Typical implementation
    //SCIP_CCALL(SCIP_SetQuadObj(lp(), qt.size(),
    //  (int*)qt.pvars1(), (int*)qt.pvars2(),
    //  (double*)qt.pcoefs()));
  }
  else {
    throw std::runtime_error("Multiple quadratic objectives not supported");
  }
}


void ScipModelAPI::linearHelper(const int* pvars, const double* pcoefs, const size_t size, const char* name, const double lb, const double ub) {
  SCIP_VAR** vars = NULL;
  SCIP_CCALL( SCIPallocBufferArray(getSCIP(), &vars, size) );
  for (size_t i = 0; i < size; i++)
    vars[i] = getPROBDATA()->vars[pvars[i]];
  SCIP_CONS* cons;
  SCIP_CCALL( SCIPcreateConsBasicLinear(getSCIP(), &cons, name, (int)size, vars, (double*)pcoefs, lb, ub) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), cons) );
  SCIP_CCALL( SCIPreleaseCons(getSCIP(), &cons) );

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


void ScipModelAPI::AddConstraint(const AndConstraint &cc)  {
  SCIP_VAR** vars = NULL;
  SCIP_CCALL( SCIPallocBufferArray(getSCIP(), &vars, cc.GetArguments().size()) );
  for (size_t i = 0; i < cc.GetArguments().size(); i++) {
    vars[i] = getPROBDATA()->vars[cc.GetArguments()[i]];
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
  for (size_t i = 0; i < dc.GetArguments().size(); i++) {
    vars[i] = getPROBDATA()->vars[dc.GetArguments()[i]];
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

  helpQuad(qc.GetName(), lt.pvars(), lt.pcoefs(), lt.size(), qt.pvars1(), qt.pvars2(), qt.pcoefs(), qt.size(), MinusInfinity(), qc.rhs());
}

void ScipModelAPI::AddConstraint( const QuadConEQ& qc ) {
  const auto& lt = qc.GetLinTerms();
  const auto& qt = qc.GetQPTerms();

  helpQuad(qc.GetName(), lt.pvars(), lt.pcoefs(), lt.size(), qt.pvars1(), qt.pvars2(), qt.pcoefs(), qt.size(), qc.rhs(), qc.rhs());
}

void ScipModelAPI::AddConstraint( const QuadConGE& qc ) {
  const auto& lt = qc.GetLinTerms();
  const auto& qt = qc.GetQPTerms();

  helpQuad(qc.GetName(), lt.pvars(), lt.pcoefs(), lt.size(), qt.pvars1(), qt.pvars2(), qt.pcoefs(), qt.size(), qc.lb(), Infinity());
}


/*
void ScipModelAPI::AddConstraint( const QuadraticConeConstraint& qc ) {
  const auto& arg = qc.GetArguments();
  SCIP_VAR** vars = NULL;
  SCIP_CCALL( SCIPallocBufferArray(getSCIP(), &vars, arg.size()) );
  for (size_t i = 0; i < arg.size(); i++) {
    vars[i] = getPROBDATA()->vars[arg[i]];
  }

  SCIP_CONS* cons;
  SCIP_CCALL( SCIPcreateConsBasicSOCNonlinear(getSCIP(), &cons, qc.GetName(), (int)arg.size(), vars, qc.GetParameters().data(), NULL, 0.0,) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), cons) );
  SCIP_CCALL( SCIPreleaseCons(getSCIP(), &cons) );

  SCIPfreeBufferArrayNull(getSCIP(), &vars, arg.size());
}*/

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
  SCIP_CCALL( SCIPcreateConsBasicSOS1(getSCIP(), &cons, sos.GetName(), sos.size(), vars, weights) );
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
  SCIP_CCALL( SCIPcreateConsBasicSOS2(getSCIP(), &cons, sos.GetName(), sos.size(), vars, weights) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), cons) );
  SCIP_CCALL( SCIPreleaseCons(getSCIP(), &cons) );

  SCIPfreeBufferArray(getSCIP(), &vars);
  SCIPfreeBufferArray(getSCIP(), &weights);
}


void ScipModelAPI::FinishProblemModificationPhase() {
}


} // namespace mp
