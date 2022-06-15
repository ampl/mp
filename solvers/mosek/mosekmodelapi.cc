#include "mosekmodelapi.h"

#include "mp/model-mgr-with-std-pb.h"
#include "mp/flat/redef/MIP/converter_mip.h"
#include "mp/flat/model_api_connect.h"

namespace mp {

/// Defining the function in ...modelapi.cc
/// for recompilation speed
std::unique_ptr<BasicModelManager>
CreateMosekModelMgr(MosekCommon& cc, Env& e,
                     pre::BasicValuePresolver*& pPre) {
  return CreateModelMgrWithFlatConverter<
      MosekModelAPI, MIPFlatConverter >(cc, e, pPre);
}


void MosekModelAPI::InitProblemModificationPhase(const FlatModelInfo* info) { 
}

void MosekModelAPI::AddVariables(const VarArrayDef& v) {
  MOSEK_CCALL(MSK_appendvars(lp(), v.size()));
  for (size_t i = v.size(); i--; ) {
    MSKboundkeye k;
    bool freelb = v.plb()[i] == -std::numeric_limits<double>::infinity();
    bool freeub = v.pub()[i] == std::numeric_limits<double>::infinity();
    if (freelb && freeub)
      k = MSK_BK_FR;
    else if (freelb == freeub)
    {
      if (v.plb()[i] == v.pub()[i])
        k = MSK_BK_FX;
      else
        k = MSK_BK_RA;
    } else {
      if (freelb)
        k = MSK_BK_UP;
      else
        k = MSK_BK_LO;
    }
    MOSEK_CCALL(MSK_putvarbound(lp(), i, k, v.plb()[i], v.pub()[i]));
    MOSEK_CCALL(MSK_putvartype(lp(), i,
      var::Type::CONTINUOUS == v.ptype()[i] ? MSK_VAR_TYPE_CONT :
      MSK_VAR_TYPE_INT));
   }
  if (v.pnames()) 
    for (size_t i = v.size(); i--; )
      MOSEK_CCALL(MSK_putvarname(lp(), i, v.pnames()[i]));
}

void MosekModelAPI::SetLinearObjective( int iobj, const LinearObjective& lo ) {
  if (iobj < 1) {
    MOSEK_CCALL(MSK_putobjsense(lp(),
      obj::Type::MAX == lo.obj_sense() ? MSK_OBJECTIVE_SENSE_MAXIMIZE : MSK_OBJECTIVE_SENSE_MINIMIZE));
    for (size_t i = 0; i < lo.num_terms(); i++) {
      MOSEK_CCALL(MSK_putcj(lp(), lo.vars()[i], lo.coefs()[i]));
    }
  }
  else {
//    TODO If we support mutiple objectives, pass them to the solver
    fmt::format("Setting {}-th linear objective: {} terms.\n", iobj, lo.num_terms());
  }
}


void MosekModelAPI::SetQuadraticObjective(int iobj, const QuadraticObjective& qo) {
  if (1 > iobj) {
    fmt::format("Setting first quadratic objective\n");
    SetLinearObjective(iobj, qo);                         // add the linear part
    const auto& qt = qo.GetQPTerms();
    fmt::format("Quadratic part is made of {} terms\n", qt.size());
    // TODO
  }
  else {
    throw std::runtime_error("Multiple quadratic objectives not supported");
  }
}
void addLinearConstraint(MSKtask_t lp, size_t size, MSKboundkey_enum key, double lb, double ub,
  const int* vindex, const double* values, const char* name) {

  MOSEK_CCALL(MSK_appendcons(lp, 1));
  int nc;
  MSK_getnumcon(lp, &nc);
  if (lb == -std::numeric_limits<double>::infinity())
    lb = -MSK_DPAR_DATA_TOL_BOUND_INF;
  if (ub == std::numeric_limits<double>::infinity())
    ub = MSK_DPAR_DATA_TOL_BOUND_INF;

  MOSEK_CCALL(MSK_putarow(lp, --nc, size, vindex, values));
  MOSEK_CCALL(MSK_putconbound(lp, nc, key, lb, ub));
  MOSEK_CCALL(MSK_putconname(lp, nc, name));

}

void MosekModelAPI::AddConstraint(const LinConRange& lc) {
  addLinearConstraint(lp(), lc.size(), MSK_BK_RA, lc.lb(), lc.ub(),
    lc.pvars(), lc.pcoefs(), lc.name());
}
void MosekModelAPI::AddConstraint(const LinConLE& lc) {
  addLinearConstraint(lp(), lc.size(), MSK_BK_UP, lc.lb(), lc.ub(),
    lc.pvars(), lc.pcoefs(), lc.name());
}
void MosekModelAPI::AddConstraint(const LinConEQ& lc) {
  addLinearConstraint(lp(), lc.size(), MSK_BK_FX, lc.lb(), lc.ub(),
    lc.pvars(), lc.pcoefs(), lc.name());
}
void MosekModelAPI::AddConstraint(const LinConGE& lc) {
  addLinearConstraint(lp(), lc.size(), MSK_BK_LO, lc.lb(), lc.ub(),
    lc.pvars(), lc.pcoefs(), lc.name());
}

void MosekModelAPI::AddConstraint(const IndicatorConstraintLinLE &ic)  {
  fmt::print("Adding indicator constraint {}\n", ic.GetTypeName());
  // TODO
}
void MosekModelAPI::AddConstraint(const IndicatorConstraintLinEQ &ic)  {
  fmt::print("Adding indicator constraint {}\n", ic.GetTypeName());
  // TODO
}
void MosekModelAPI::AddConstraint(const IndicatorConstraintLinGE &ic)  {
  fmt::print("Adding indicator constraint {}\n", ic.GetTypeName());
  // TODO
}

void MosekModelAPI::AddConstraint(const QuadConRange& qc) {
  fmt::print("Adding quadratic constraint {}\n", qc.GetTypeName());
  // TODO
  /*
  const auto& lt = qc.GetLinTerms();
  const auto& qt = qc.GetQPTerms();
  GRB_CALL( GRBaddqrangeconstr(model(), lt.size(), (int*)lt.pvars(), (double*)lt.pcoefs(),
                          qt.size(), (int*)qt.pvars1(), (int*)qt.pvars2(),
                          (double*)qt.pcoefs(), qc.lb(), qc.ub(), NULL) );
  */
}

void MosekModelAPI::AddConstraint( const QuadConLE& qc ) {
  fmt::print("Adding quadratic constraint {}\n", qc.GetTypeName());
  // TODO
  /*
  const auto& lt = qc.GetLinTerms();
  const auto& qt = qc.GetQPTerms();
  GRB_CALL( GRBaddqconstr(model(), lt.size(), (int*)lt.pvars(), (double*)lt.pcoefs(),
                          qt.size(), (int*)qt.pvars1(), (int*)qt.pvars2(),
                          (double*)qt.pcoefs(), GRB_LESS_EQUAL, qc.rhs(), NULL) );
                          */
}

void MosekModelAPI::AddConstraint( const QuadConEQ& qc ) {
  fmt::print("Adding quadratic constraint {}\n", qc.GetTypeName());
  // TODO
  /*
  const auto& lt = qc.GetLinTerms();
  const auto& qt = qc.GetQPTerms();
  GRB_CALL( GRBaddqconstr(model(), lt.size(), (int*)lt.pvars(), (double*)lt.pcoefs(),
                          qt.size(), (int*)qt.pvars1(), (int*)qt.pvars2(),
                          (double*)qt.pcoefs(), GRB_EQUAL, qc.rhs(), NULL) );
                          */
}

void MosekModelAPI::AddConstraint( const QuadConGE& qc ) {
  fmt::print("Adding quadratic constraint {}\n", qc.GetTypeName());
  // TODO
  /*
  const auto& lt = qc.GetLinTerms();
  const auto& qt = qc.GetQPTerms();
  GRB_CALL( GRBaddqconstr(model(), lt.size(), (int*)lt.pvars(), (double*)lt.pcoefs(),
                          qt.size(), (int*)qt.pvars1(), (int*)qt.pvars2(),
                          (double*)qt.pcoefs(), GRB_GREATER_EQUAL, qc.rhs(), NULL) );
                          */
}

void MosekModelAPI::AddConstraint(const SOS1Constraint& sos) {
  fmt::print("Adding SOS1 constraint {}\n", sos.GetTypeName());
  // TODO
/*  int type = MOSEK_SOS_TYPE1;
  int beg = 0;
  const int size = sos.size();
  MOSEK_CCALL(MOSEK_AddSOSs(lp(), 1, &type, &beg,
    &size, (int*)sos.get_vars().data(),
    (double*)sos.get_weights().data())); */
}

void MosekModelAPI::AddConstraint(const SOS2Constraint& sos) {
  fmt::print("Adding SOS1 constraint {}\n", sos.GetTypeName());
  // TODO
  /*int type = MOSEK_SOS_TYPE2;
  int beg = 0;
  const int size = sos.size();
  MOSEK_CCALL(MOSEK_AddSOSs(lp(), 1, &type, &beg,
    &size, (int*)sos.get_vars().data(),
    (double*)sos.get_weights().data()));*/
}


void MosekModelAPI::FinishProblemModificationPhase() {
}


} // namespace mp
