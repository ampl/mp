#include "gcgmpmodelapi.h"


namespace mp {

void GcgModelAPI::InitProblemModificationPhase(const FlatModelInfo* flat_model_info) {
  // Allocate storage if needed:
  int n_linear_cons = flat_model_info->GetNumberOfConstraintsOfGroup(CG_Linear);
  getPROBDATA()->nlinconss = n_linear_cons;
  getPROBDATA()->i = 0;
  GCG_CCALL( SCIPallocBlockMemoryArray(getSCIP(), &getPROBDATA()->linconss, getPROBDATA()->nlinconss) );
}

void GcgModelAPI::AddVariables(const VarArrayDef& v) {
  assert( SCIPgetStage(getSCIP()) == SCIP_STAGE_PROBLEM );

  getPROBDATA()->nvars = v.size();
  GCG_CCALL( SCIPallocBlockMemoryArray(getSCIP(), &getPROBDATA()->vars, getPROBDATA()->nvars) );

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
      GCG_CCALL( SCIPcreateVarBasic(getSCIP(), &var, v.pnames()[i], lb, ub, objcoef, vartype) );
    else
      GCG_CCALL( SCIPcreateVarBasic(getSCIP(), &var, NULL, lb, ub, objcoef, vartype) );
    GCG_CCALL( SCIPaddVar(getSCIP(), var) );
    getPROBDATA()->vars[i] = var;
  }
}

void GcgModelAPI::SetLinearObjective( int iobj, const LinearObjective& lo ) {
  if (iobj<1) {
    GCG_CCALL( SCIPsetObjsense(getSCIP(), 
                    obj::Type::MAX==lo.obj_sense() ? SCIP_OBJSENSE_MAXIMIZE : SCIP_OBJSENSE_MINIMIZE) );
    SCIP_VAR** vars = getPROBDATA()->vars;
    for (int i = 0; i < lo.num_terms(); i++) {
      GCG_CCALL( SCIPchgVarObj(getSCIP(), vars[lo.vars()[i]], lo.coefs()[i]) );
    }
  } 
  else {
    throw std::runtime_error("Multiple linear objectives not supported");
  }
}


void GcgModelAPI::SetQuadraticObjective(int iobj, const QuadraticObjective& qo) {
  if (1 > iobj) {
    throw std::runtime_error("Quadratic objectives not supported");
  }
  else {
    throw std::runtime_error("Multiple quadratic objectives not supported");
  }
}


void GcgModelAPI::linearHelper(const int* pvars, const double* pcoefs, const size_t size, const char* name, const double lb, const double ub) {
  SCIP_VAR** vars = NULL;
  GCG_CCALL( SCIPallocBufferArray(getSCIP(), &vars, size) );
  for (size_t i = 0; i < size; i++)
    vars[i] = getPROBDATA()->vars[pvars[i]];

  SCIP_CONS* cons;
  GCG_CCALL( SCIPcreateConsBasicLinear(getSCIP(), &cons, name, (int)size, vars, (double*)pcoefs, lb, ub) );
  GCG_CCALL( SCIPaddCons(getSCIP(), cons) );
  getPROBDATA()->linconss[getPROBDATA()->i] = cons;
  getPROBDATA()->i++;

  SCIPfreeBufferArray(getSCIP(), &vars);
}
void GcgModelAPI::AddConstraint(const LinConRange& lc) {
  linearHelper(lc.pvars(), lc.pcoefs(), lc.size(), lc.GetName(), lc.lb(), lc.ub());
}
void GcgModelAPI::AddConstraint(const LinConLE& lc) {
  linearHelper(lc.pvars(), lc.pcoefs(), lc.size(), lc.GetName(), MinusInfinity(), lc.ub());
}
void GcgModelAPI::AddConstraint(const LinConEQ& lc) {
  linearHelper(lc.pvars(), lc.pcoefs(), lc.size(), lc.GetName(), lc.ub(), lc.ub());
}
void GcgModelAPI::AddConstraint(const LinConGE& lc) {
linearHelper(lc.pvars(), lc.pcoefs(), lc.size(), lc.GetName(), lc.lb(), Infinity());
}

void GcgModelAPI::FinishProblemModificationPhase() {
}


} // namespace mp
