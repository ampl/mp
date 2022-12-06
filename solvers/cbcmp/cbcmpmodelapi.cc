#include "cbcmpmodelapi.h"

#include "mp/model-mgr-with-std-pb.h"
#include "mp/flat/redef/MIP/converter_mip.h"
#include "mp/flat/model_api_connect.h"

namespace mp {

std::unique_ptr<BasicModelManager>
CreateCbcmpModelMgr(CbcmpCommon& cc, Env& e,
                     pre::BasicValuePresolver*& pPre) {
  return CreateModelMgrWithFlatConverter<
      CbcmpModelAPI, MIPFlatConverter >(cc, e, pPre);
}


void CbcmpModelAPI::InitProblemModificationPhase(
    const FlatModelInfo*) {
}

void CbcmpModelAPI::AddVariables(const VarArrayDef& v) {
  for (int i = 0; i < v.size(); i++) {
    Cbc_addCol(lp(), "V", v.plb()[i], v.pub()[i],
      0, var::Type::INTEGER == v.ptype()[i], 0, NULL, NULL);
  }
}

void CbcmpModelAPI::SetLinearObjective( int iobj, const LinearObjective& lo ) {
  if (iobj<1) {
    for (int i = 0; i < lo.num_terms(); i++)
      Cbc_setObjCoeff(lp(), lo.vars()[i], lo.coefs()[i]);
    Cbc_setObjSense(lp(), lo.obj_sense()== obj::MAX ? -1 : 1);
  } else {
    throw std::runtime_error("Multiple objectives not supported");
  }
}

void CbcmpModelAPI::AddConstraint(const LinConLE& lc) {
  Cbc_addRow(lp(), lc.GetName(), lc.size(), lc.pvars(),
    lc.pcoefs(), 'L', lc.rhs());
}
void CbcmpModelAPI::AddConstraint(const LinConEQ& lc) {
  Cbc_addRow(lp(), lc.GetName(), lc.size(), lc.pvars(),
    lc.pcoefs(), 'E', lc.rhs());
}
void CbcmpModelAPI::AddConstraint(const LinConGE& lc) {
  Cbc_addRow(lp(), lc.GetName(), lc.size(), lc.pvars(),
    lc.pcoefs(), 'G', lc.rhs());
}

void CbcmpModelAPI::AddConstraint(const SOS1Constraint& sos) {
  int starts[] = { 0,  (int)sos.get_vars().size() };
  Cbc_addSOS(lp(), 1, starts,
    (int*)sos.get_vars().data(), (double*)sos.get_weights().data(), 1);
}

void CbcmpModelAPI::AddConstraint(const SOS2Constraint& sos) {
  int starts[] = { 0,  (int)sos.get_vars().size() };
  Cbc_addSOS(lp(), 1, starts,
    (int*)sos.get_vars().data(), (double*)sos.get_weights().data(), 2);
}

void CbcmpModelAPI::FinishProblemModificationPhase() {
}


} // namespace mp
