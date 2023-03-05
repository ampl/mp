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

void ScipModelAPI::AddConstraint(const LinConRange& lc) {
  SCIP_VAR** vars = NULL;
  SCIP_CCALL( SCIPallocBlockMemoryArray(getSCIP(), &vars, lc.size()) );
  for (size_t i = 0; i < lc.size(); i++) {
    vars[i] = getPROBDATA()->vars[lc.pvars()[i]];
  }
  SCIP_CONS* cons;
  SCIP_CCALL( SCIPcreateConsBasicLinear(getSCIP(), &cons, lc.GetName(), (int)lc.size(), vars, (double*)lc.pcoefs(), lc.lb(), lc.ub()) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), cons) );
  SCIP_CCALL( SCIPreleaseCons(getSCIP(), &cons) );

  SCIPfreeBlockMemoryArrayNull(getSCIP(), &vars, lc.size());
}
void ScipModelAPI::AddConstraint(const LinConLE& lc) {
  SCIP_VAR** vars = NULL;
  SCIP_CCALL( SCIPallocBlockMemoryArray(getSCIP(), &vars, lc.size()) );
  for (size_t i = 0; i < lc.size(); i++) {
    vars[i] = getPROBDATA()->vars[lc.pvars()[i]];
  }
  SCIP_CONS* cons;
  SCIP_CCALL( SCIPcreateConsBasicLinear(getSCIP(), &cons, lc.GetName(), (int)lc.size(), vars, (double*)lc.pcoefs(), MinusInfinity(), lc.rhs()) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), cons) );
  SCIP_CCALL( SCIPreleaseCons(getSCIP(), &cons) );

  SCIPfreeBlockMemoryArrayNull(getSCIP(), &vars, lc.size());
}
void ScipModelAPI::AddConstraint(const LinConEQ& lc) {
  SCIP_VAR** vars = NULL;
  SCIP_CCALL( SCIPallocBlockMemoryArray(getSCIP(), &vars, lc.size()) );
  for (size_t i = 0; i < lc.size(); i++) {
    vars[i] = getPROBDATA()->vars[lc.pvars()[i]];
  }
  SCIP_CONS* cons;
  SCIP_CCALL( SCIPcreateConsBasicLinear(getSCIP(), &cons, lc.GetName(), (int)lc.size(), vars, (double*)lc.pcoefs(), lc.rhs(), lc.rhs()) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), cons) );
  SCIP_CCALL( SCIPreleaseCons(getSCIP(), &cons) );

  SCIPfreeBlockMemoryArrayNull(getSCIP(), &vars, lc.size());
}
void ScipModelAPI::AddConstraint(const LinConGE& lc) {
  SCIP_VAR** vars = NULL;
  SCIP_CCALL( SCIPallocBlockMemoryArray(getSCIP(), &vars, lc.size()) );
  for (size_t i = 0; i < lc.size(); i++) {
    vars[i] = getPROBDATA()->vars[lc.pvars()[i]];
  }
  SCIP_CONS* cons;
  SCIP_CCALL( SCIPcreateConsBasicLinear(getSCIP(), &cons, lc.GetName(), (int)lc.size(), vars, (double*)lc.pcoefs(), lc.lb(), Infinity()) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), cons) );
  SCIP_CCALL( SCIPreleaseCons(getSCIP(), &cons) );

  SCIPfreeBlockMemoryArrayNull(getSCIP(), &vars, lc.size());
}

void ScipModelAPI::AddConstraint(const AndConstraint &cc)  {
  SCIP_VAR** vars = NULL;
  SCIP_CCALL( SCIPallocBlockMemoryArray(getSCIP(), &vars, cc.GetArguments().size()) );
  for (size_t i = 0; i < cc.GetArguments().size(); i++) {
    vars[i] = getPROBDATA()->vars[cc.GetArguments()[i]];
  }

  SCIP_CONS* cons;
  SCIP_CCALL( SCIPcreateConsBasicAnd(getSCIP(), &cons, cc.GetName(), getPROBDATA()->vars[cc.GetResultVar()],
    (int)cc.GetArguments().size(), vars) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), cons) );
  SCIP_CCALL( SCIPreleaseCons(getSCIP(), &cons) );

  SCIPfreeBlockMemoryArrayNull(getSCIP(), &vars, cc.GetArguments().size());
}

void ScipModelAPI::AddConstraint(const OrConstraint &dc)  {
  SCIP_VAR** vars = NULL;
  SCIP_CCALL( SCIPallocBlockMemoryArray(getSCIP(), &vars, dc.GetArguments().size()) );
  for (size_t i = 0; i < dc.GetArguments().size(); i++) {
    vars[i] = getPROBDATA()->vars[dc.GetArguments()[i]];
  }

  SCIP_CONS* cons;
  SCIP_CCALL( SCIPcreateConsBasicOr(getSCIP(), &cons, dc.GetName(), getPROBDATA()->vars[dc.GetResultVar()],
    (int)dc.GetArguments().size(), vars) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), cons) );
  SCIP_CCALL( SCIPreleaseCons(getSCIP(), &cons) );

  SCIPfreeBlockMemoryArrayNull(getSCIP(), &vars, dc.GetArguments().size());
}

void ScipModelAPI::AddConstraint(const IndicatorConstraintLinLE &ic)  {
  SCIP_VAR** vars = NULL;
  SCIP_CCALL( SCIPallocBlockMemoryArray(getSCIP(), &vars, ic.get_constraint().size()) );
  for (size_t i = 0; i < ic.get_constraint().size(); i++) {
    vars[i] = getPROBDATA()->vars[ic.get_constraint().pvars()[i]];
  }
  SCIP_CONS* cons;
  SCIP_CCALL( SCIPcreateConsIndicatorGeneric(getSCIP(), &cons, ic.GetName(), getPROBDATA()->vars[ic.get_binary_var()],
    (int)ic.get_constraint().size(), vars, (double*)ic.get_constraint().pcoefs(), ic.get_constraint().rhs(),
    ic.get_binary_value(), TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), cons) );
  SCIP_CCALL( SCIPreleaseCons(getSCIP(), &cons) );

  SCIPfreeBlockMemoryArrayNull(getSCIP(), &vars, ic.get_constraint().size());              
}
void ScipModelAPI::AddConstraint(const IndicatorConstraintLinEQ &ic)  {
  SCIP_VAR** vars = NULL;
  SCIP_CCALL( SCIPallocBlockMemoryArray(getSCIP(), &vars, ic.get_constraint().size()) );
  for (size_t i = 0; i < ic.get_constraint().size(); i++) {
    vars[i] = getPROBDATA()->vars[ic.get_constraint().pvars()[i]];
  }
  SCIP_CONS* cons1;
  SCIP_CONS* cons2;
  SCIP_CCALL( SCIPcreateConsIndicatorGeneric(getSCIP(), &cons1, ic.GetName(),getPROBDATA()->vars[ic.get_binary_var()],
    (int)ic.get_constraint().size(), vars, (double*)ic.get_constraint().pcoefs(), ic.get_constraint().rhs(),
    ic.get_binary_value(), TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
  SCIP_CCALL( SCIPcreateConsIndicatorGeneric(getSCIP(), &cons2, ic.GetName(),getPROBDATA()->vars[ic.get_binary_var()],
    (int)ic.get_constraint().size(), vars, (double*)ic.get_constraint().pcoefs(), ic.get_constraint().rhs(),
    ic.get_binary_value(), FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), cons1) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), cons2) );
  SCIP_CCALL( SCIPreleaseCons(getSCIP(), &cons1) );
  SCIP_CCALL( SCIPreleaseCons(getSCIP(), &cons2) );

  SCIPfreeBlockMemoryArrayNull(getSCIP(), &vars, ic.get_constraint().size());   
}
void ScipModelAPI::AddConstraint(const IndicatorConstraintLinGE &ic)  {
  SCIP_VAR** vars = NULL;
  SCIP_CCALL( SCIPallocBlockMemoryArray(getSCIP(), &vars, ic.get_constraint().size()) );
  for (size_t i = 0; i < ic.get_constraint().size(); i++) {
    vars[i] = getPROBDATA()->vars[ic.get_constraint().pvars()[i]];
  }
  SCIP_CONS* cons;
  SCIP_CCALL( SCIPcreateConsIndicatorGeneric(getSCIP(), &cons, ic.GetName(),getPROBDATA()->vars[ic.get_binary_var()],
    (int)ic.get_constraint().size(), vars, (double*)ic.get_constraint().pcoefs(), ic.get_constraint().rhs(),
    ic.get_binary_value(), FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), cons) );
  SCIP_CCALL( SCIPreleaseCons(getSCIP(), &cons) );

  SCIPfreeBlockMemoryArrayNull(getSCIP(), &vars, ic.get_constraint().size());
}

void ScipModelAPI::AddConstraint(const QuadConRange& qc) {
  // allocation of linterm vars
  const auto& lt = qc.GetLinTerms();
  SCIP_VAR** linvars = NULL;
  SCIP_CCALL( SCIPallocBlockMemoryArray(getSCIP(), &linvars, lt.size()) );
  for (size_t i = 0; i < lt.size(); i++) {
    linvars[i] = getPROBDATA()->vars[lt.pvars()[i]];
  }

  const auto& qt = qc.GetQPTerms();
  // allocation of qterm1 vars
  SCIP_VAR** qvars1 = NULL;
  SCIP_CCALL( SCIPallocBlockMemoryArray(getSCIP(), &qvars1, qt.size()) );
  for (size_t i = 0; i < qt.size(); i++) {
    qvars1[i] = getPROBDATA()->vars[qt.pvars1()[i]];
  }

  // allocation of qterm2 vars
  SCIP_VAR** qvars2 = NULL;
  SCIP_CCALL( SCIPallocBlockMemoryArray(getSCIP(), &qvars2, qt.size()) );
  for (size_t i = 0; i < qt.size(); i++) {
    qvars2[i] = getPROBDATA()->vars[qt.pvars2()[i]];
  }

  SCIP_CONS* cons;
  SCIP_CCALL( SCIPcreateConsBasicQuadraticNonlinear(getSCIP(), &cons, qc.GetName(), (int)lt.size(), linvars, (double*)lt.pcoefs(), (int)qt.size(), qvars1, qvars2, (double*)qt.pcoefs(), qc.lb(), qc.ub()) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), cons) );
  SCIP_CCALL( SCIPreleaseCons(getSCIP(), &cons) );

  SCIPfreeBlockMemoryArrayNull(getSCIP(), &linvars, lt.size());
  SCIPfreeBlockMemoryArrayNull(getSCIP(), &qvars1, qt.size());
  SCIPfreeBlockMemoryArrayNull(getSCIP(), &qvars2, qt.size());
}

void ScipModelAPI::AddConstraint( const QuadConLE& qc ) {
  // allocation of linterm vars
  const auto& lt = qc.GetLinTerms();
  SCIP_VAR** linvars = NULL;
  SCIP_CCALL( SCIPallocBlockMemoryArray(getSCIP(), &linvars, lt.size()) );
  for (size_t i = 0; i < lt.size(); i++) {
    linvars[i] = getPROBDATA()->vars[lt.pvars()[i]];
  }

  const auto& qt = qc.GetQPTerms();
  // allocation of qterm1 vars
  SCIP_VAR** qvars1 = NULL;
  SCIP_CCALL( SCIPallocBlockMemoryArray(getSCIP(), &qvars1, qt.size()) );
  for (size_t i = 0; i < qt.size(); i++) {
    qvars1[i] = getPROBDATA()->vars[qt.pvars1()[i]];
  }

  // allocation of qterm2 vars
  SCIP_VAR** qvars2 = NULL;
  SCIP_CCALL( SCIPallocBlockMemoryArray(getSCIP(), &qvars2, qt.size()) );
  for (size_t i = 0; i < qt.size(); i++) {
    qvars2[i] = getPROBDATA()->vars[qt.pvars2()[i]];
  }

  SCIP_CONS* cons;
  SCIP_CCALL( SCIPcreateConsBasicQuadraticNonlinear(getSCIP(), &cons, qc.GetName(), (int)lt.size(), linvars, (double*)lt.pcoefs(), (int)qt.size(), qvars1, qvars2, (double*)qt.pcoefs(), MinusInfinity(), qc.rhs()) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), cons) );
  SCIP_CCALL( SCIPreleaseCons(getSCIP(), &cons) );

  SCIPfreeBlockMemoryArrayNull(getSCIP(), &linvars, lt.size());
  SCIPfreeBlockMemoryArrayNull(getSCIP(), &qvars1, qt.size());
  SCIPfreeBlockMemoryArrayNull(getSCIP(), &qvars2, qt.size());
}

void ScipModelAPI::AddConstraint( const QuadConEQ& qc ) {
  // allocation of linterm vars
  const auto& lt = qc.GetLinTerms();
  SCIP_VAR** linvars = NULL;
  SCIP_CCALL( SCIPallocBlockMemoryArray(getSCIP(), &linvars, lt.size()) );
  for (size_t i = 0; i < lt.size(); i++) {
    linvars[i] = getPROBDATA()->vars[lt.pvars()[i]];
  }

  const auto& qt = qc.GetQPTerms();
  // allocation of qterm1 vars
  SCIP_VAR** qvars1 = NULL;
  SCIP_CCALL( SCIPallocBlockMemoryArray(getSCIP(), &qvars1, qt.size()) );
  for (size_t i = 0; i < qt.size(); i++) {
    qvars1[i] = getPROBDATA()->vars[qt.pvars1()[i]];
  }

  // allocation of qterm2 vars
  SCIP_VAR** qvars2 = NULL;
  SCIP_CCALL( SCIPallocBlockMemoryArray(getSCIP(), &qvars2, qt.size()) );
  for (size_t i = 0; i < qt.size(); i++) {
    qvars2[i] = getPROBDATA()->vars[qt.pvars2()[i]];
  }

  SCIP_CONS* cons;
  SCIP_CCALL( SCIPcreateConsBasicQuadraticNonlinear(getSCIP(), &cons, qc.GetName(), (int)lt.size(), linvars, (double*)lt.pcoefs(), (int)qt.size(), qvars1, qvars2, (double*)qt.pcoefs(), qc.rhs(), qc.rhs()) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), cons) );
  SCIP_CCALL( SCIPreleaseCons(getSCIP(), &cons) );

  SCIPfreeBlockMemoryArrayNull(getSCIP(), &linvars, lt.size());
  SCIPfreeBlockMemoryArrayNull(getSCIP(), &qvars1, qt.size());
  SCIPfreeBlockMemoryArrayNull(getSCIP(), &qvars2, qt.size());
}

void ScipModelAPI::AddConstraint( const QuadConGE& qc ) {
  // allocation of linterm vars
  const auto& lt = qc.GetLinTerms();
  SCIP_VAR** linvars = NULL;
  SCIP_CCALL( SCIPallocBlockMemoryArray(getSCIP(), &linvars, lt.size()) );
  for (size_t i = 0; i < lt.size(); i++) {
    linvars[i] = getPROBDATA()->vars[lt.pvars()[i]];
  }

  const auto& qt = qc.GetQPTerms();
  // allocation of qterm1 vars
  SCIP_VAR** qvars1 = NULL;
  SCIP_CCALL( SCIPallocBlockMemoryArray(getSCIP(), &qvars1, qt.size()) );
  for (size_t i = 0; i < qt.size(); i++) {
    qvars1[i] = getPROBDATA()->vars[qt.pvars1()[i]];
  }

  // allocation of qterm2 vars
  SCIP_VAR** qvars2 = NULL;
  SCIP_CCALL( SCIPallocBlockMemoryArray(getSCIP(), &qvars2, qt.size()) );
  for (size_t i = 0; i < qt.size(); i++) {
    qvars2[i] = getPROBDATA()->vars[qt.pvars2()[i]];
  }

  SCIP_CONS* cons;
  SCIP_CCALL( SCIPcreateConsBasicQuadraticNonlinear(getSCIP(), &cons, qc.GetName(), (int)lt.size(), linvars, (double*)lt.pcoefs(), (int)qt.size(), qvars1, qvars2, (double*)qt.pcoefs(), qc.lb(), Infinity()) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), cons) );
  SCIP_CCALL( SCIPreleaseCons(getSCIP(), &cons) );

  SCIPfreeBlockMemoryArrayNull(getSCIP(), &linvars, lt.size());
  SCIPfreeBlockMemoryArrayNull(getSCIP(), &qvars1, qt.size());
  SCIPfreeBlockMemoryArrayNull(getSCIP(), &qvars2, qt.size());
}

void ScipModelAPI::AddConstraint(const SOS1Constraint& sos) {
  SCIP_VAR** vars = NULL;
  double* weights = NULL;
  SCIP_CCALL( SCIPallocBlockMemoryArray(getSCIP(), &vars, sos.size()) );
  SCIP_CCALL( SCIPallocBlockMemoryArray(getSCIP(), &weights, sos.size()) );
  for (int i = 0; i < sos.size(); i++) {
    vars[i] = getPROBDATA()->vars[sos.get_vars().data()[i]];
    weights[i] = sos.get_weights().data()[i];
  }

  SCIP_CONS* cons = NULL;
  SCIP_CCALL( SCIPcreateConsBasicSOS1(getSCIP(), &cons, sos.GetName(), sos.size(), vars, weights) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), cons) );
  SCIP_CCALL( SCIPreleaseCons(getSCIP(), &cons) );

  SCIPfreeBlockMemoryArrayNull(getSCIP(), &vars, sos.size());
  SCIPfreeBlockMemoryArrayNull(getSCIP(), &weights, sos.size());
}

void ScipModelAPI::AddConstraint(const SOS2Constraint& sos) {
  SCIP_VAR** vars = NULL;
  double* weights = NULL;
  SCIP_CCALL( SCIPallocBlockMemoryArray(getSCIP(), &vars, sos.size()) );
  SCIP_CCALL( SCIPallocBlockMemoryArray(getSCIP(), &weights, sos.size()) );
  for (int i = 0; i < sos.size(); i++) {
    vars[i] = getPROBDATA()->vars[sos.get_vars().data()[i]];
    weights[i] = sos.get_weights().data()[i];
  }

  SCIP_CONS* cons = NULL;
  SCIP_CCALL( SCIPcreateConsBasicSOS2(getSCIP(), &cons, sos.GetName(), sos.size(), vars, weights) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), cons) );
  SCIP_CCALL( SCIPreleaseCons(getSCIP(), &cons) );

  SCIPfreeBlockMemoryArrayNull(getSCIP(), &vars, sos.size());
  SCIPfreeBlockMemoryArrayNull(getSCIP(), &weights, sos.size());
}


void ScipModelAPI::FinishProblemModificationPhase() {
}


} // namespace mp
