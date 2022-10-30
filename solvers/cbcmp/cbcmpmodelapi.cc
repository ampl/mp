#include "cbcmpmodelapi.h"

#include "mp/model-mgr-with-std-pb.h"
#include "mp/flat/redef/MIP/converter_mip.h"
#include "mp/flat/model_api_connect.h"

namespace mp {

/// Defining the function in ...modelapi.cc
/// for recompilation speed
std::unique_ptr<BasicModelManager>
CreateCbcmpModelMgr(CbcmpCommon& cc, Env& e,
                     pre::BasicValuePresolver*& pPre) {
  return CreateModelMgrWithFlatConverter<
      CbcmpModelAPI, MIPFlatConverter >(cc, e, pPre);
}


void CbcmpModelAPI::InitProblemModificationPhase(
    const FlatModelInfo*) {
  // Allocate storage if needed:
  // auto n_linear_cons =
  //   flat_model_info->GetNumberOfConstraintsOfGroup(CG_LINEAR);
  // preallocate_linear_cons( n_linear_cons );
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
//    TODO If we support mutiple objectives, pass them to the solver
    fmt::format("Setting {}-th linear objective: {} terms.\n", iobj, lo.num_terms());
  }
}


void CbcmpModelAPI::SetQuadraticObjective(int iobj, const QuadraticObjective& qo) {
  if (1 > iobj) {
    fmt::format("Setting first quadratic objective\n");
    SetLinearObjective(iobj, qo);                         // add the linear part
    const auto& qt = qo.GetQPTerms();
    fmt::format("Quadratic part is made of {} terms\n", qt.size());

    // Typical implementation
    //CBCMP_CCALL(CBCMP_SetQuadObj(lp(), qt.size(),
    //  (int*)qt.pvars1(), (int*)qt.pvars2(),
    //  (double*)qt.pcoefs()));
  }
  else {
    throw std::runtime_error("Multiple quadratic objectives not supported");
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

void CbcmpModelAPI::AddConstraint(const IndicatorConstraintLinLE &ic)  {
  fmt::print("Adding indicator constraint {}\n", ic.GetName());
  /*CBCMP_CCALL(CBCMP_AddIndicator(lp(),
    ic.get_binary_var(), ic.get_binary_value(),
    (int)ic.get_constraint().size(),
    ic.get_constraint().pvars(),
    ic.get_constraint().pcoefs(),
    CBCMP_LESS_EQUAL,
    ic.get_constraint().rhs()));*/
                               
}
void CbcmpModelAPI::AddConstraint(const IndicatorConstraintLinEQ &ic)  {
  fmt::print("Adding indicator constraint {}\n", ic.GetName());
  /*CBCMP_CCALL(CBCMP_AddIndicator(lp(),
    ic.get_binary_var(), ic.get_binary_value(),
    (int)ic.get_constraint().size(),
    ic.get_constraint().pvars(),
    ic.get_constraint().pcoefs(),
    CBCMP_EQUAL,
    ic.get_constraint().rhs()));*/
}
void CbcmpModelAPI::AddConstraint(const IndicatorConstraintLinGE &ic)  {
  fmt::print("Adding indicator constraint {}\n", ic.GetName());
  /*CBCMP_CCALL(CBCMP_AddIndicator(lp(),
    ic.get_binary_var(), ic.get_binary_value(),
    (int)ic.get_constraint().size(),
    ic.get_constraint().pvars(),
    ic.get_constraint().pcoefs(),
    CBCMP_GREATER_EQUAL,
    ic.get_constraint().rhs()));*/

}

void CbcmpModelAPI::AddConstraint(const SOS1Constraint& sos) {
  fmt::print("Adding SOS1 constraint {}\n", sos.GetName());
/*  int type = CBCMP_SOS_TYPE1;
  int beg = 0;
  const int size = sos.size();
  CBCMP_CCALL(CBCMP_AddSOSs(lp(), 1, &type, &beg,
    &size, (int*)sos.get_vars().data(),
    (double*)sos.get_weights().data())); */
}

void CbcmpModelAPI::AddConstraint(const SOS2Constraint& sos) {
  fmt::print("Adding SOS1 constraint {}\n", sos.GetName());
  /*int type = CBCMP_SOS_TYPE2;
  int beg = 0;
  const int size = sos.size();
  CBCMP_CCALL(CBCMP_AddSOSs(lp(), 1, &type, &beg,
    &size, (int*)sos.get_vars().data(),
    (double*)sos.get_weights().data()));*/
}


void CbcmpModelAPI::FinishProblemModificationPhase() {
}


} // namespace mp
