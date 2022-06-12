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


void MosekModelAPI::InitProblemModificationPhase() { }

void MosekModelAPI::AddVariables(const VarArrayDef& v) {
  // TODO Add variables using solver API; typically, first convert the variable
  // type to the appropriate solver-defined type, then add them as a bunch
  int nInteger = 0, nContinuous = 0;
  for (auto i = v.size(); i--; )
    if (var::Type::CONTINUOUS == v.ptype()[i])
      nContinuous++;
    else
      nInteger++;
  fmt::print("Adding {} continuous and {} integer variables.\n", nContinuous, nInteger);

  /* Typical implementation
  for (size_t i=v.size(); i--; )
    vtypes[i] = var::Type::CONTINUOUS == v.ptype()[i] ?
          MOSEK_CONTINUOUS : MOSEK_INTEGER;
  MOSEK_CCALL(MOSEK_AddCols(lp(), (int)v.size(), NULL, NULL,
    NULL, NULL, NULL, vtypes.data(), v.plb(), v.pub(),  NULL)); */
}

void MosekModelAPI::SetLinearObjective( int iobj, const LinearObjective& lo ) {
  if (iobj<1) {
    fmt::format("Setting first linear objective: {} terms.\n", lo.num_terms());
    /*
    MOSEK_CCALL(MOSEK_SetObjSense(lp(), 
                    obj::Type::MAX==lo.obj_sense() ? MOSEK_MAXIMIZE : MOSEK_MINIMIZE) );
    MOSEK_CCALL(MOSEK_SetColObj(lp(), lo.num_terms(),
                           lo.vars().data(), lo.coefs().data()) ); */
  } else {
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

    // Typical implementation
    //MOSEK_CCALL(MOSEK_SetQuadObj(lp(), qt.size(), (int*)qt.pvars1(), (int*)qt.pvars2(),
    //  (double*)qt.pcoefs()));
  }
  else {
    throw std::runtime_error("Multiple quadratic objectives not supported");
  }
}

void MosekModelAPI::AddConstraint(const LinConRange& lc) {
  fmt::print("Adding range linear constraint {}\n", lc.GetTypeName());
  fmt::print("{} <=", lc.lb());
  for (size_t i = 0; i < lc.size(); i++)
  {
    fmt::print(" {}", lc.pcoefs()[i] >= 0 ? '+' : '-');
    if (std::fabs(lc.pcoefs()[i]) != 1.0)
      fmt::print("{}*", lc.pcoefs()[i]);
    fmt::print("x{}", lc.pvars()[i]);
  }
  fmt::print(" <= {}\n", lc.ub());
//  MOSEK_CCALL(MOSEK_AddRow(lp(), lc.size(), lc.pvars(), lc.pcoefs(), 
 //   NULL, lc.lb(), lc.ub(), NULL));
}
void MosekModelAPI::AddConstraint(const LinConLE& lc) {
  fmt::print("Adding <= linear constraint {}\n", lc.GetTypeName());
 // char sense = MOSEK_LESS_EQUAL;
 // MOSEK_CCALL(MOSEK_AddRow(lp(), lc.size(), lc.pvars(), lc.pcoefs(),
  //  sense, lc.rhs(), 0, NULL));
}
void MosekModelAPI::AddConstraint(const LinConEQ& lc) {
  fmt::print("Adding == linear constraint {}\n", lc.GetTypeName());
//  char sense = MOSEK_EQUAL;
// MOSEK_CCALL(MOSEK_AddRow(lp(), lc.size(), lc.pvars(), lc.pcoefs(),
//   sense, lc.rhs(), 0, NULL));
}
void MosekModelAPI::AddConstraint(const LinConGE& lc) {
  fmt::print("Adding >= linear constraint {}\n", lc.GetTypeName());
  //char sense = MOSEK_GREATER_EQUAL;
  //MOSEK_CCALL(MOSEK_AddRow(lp(), lc.size(), lc.pvars(), lc.pcoefs(),
  //  sense, lc.rhs(), 0, NULL));
}

void MosekModelAPI::AddConstraint(const IndicatorConstraintLinLE &ic)  {
  fmt::print("Adding indicator constraint {}\n", ic.GetTypeName());
  /*MOSEK_CCALL(MOSEK_AddIndicator(lp(),
    ic.get_binary_var(), ic.get_binary_value(),
    (int)ic.get_constraint().size(),
    ic.get_constraint().pvars(),
    ic.get_constraint().pcoefs(),
    MOSEK_LESS_EQUAL,
    ic.get_constraint().rhs()));*/
                               
}
void MosekModelAPI::AddConstraint(const IndicatorConstraintLinEQ &ic)  {
  fmt::print("Adding indicator constraint {}\n", ic.GetTypeName());
  /*MOSEK_CCALL(MOSEK_AddIndicator(lp(),
    ic.get_binary_var(), ic.get_binary_value(),
    (int)ic.get_constraint().size(),
    ic.get_constraint().pvars(),
    ic.get_constraint().pcoefs(),
    MOSEK_EQUAL,
    ic.get_constraint().rhs()));*/
}
void MosekModelAPI::AddConstraint(const IndicatorConstraintLinGE &ic)  {
  fmt::print("Adding indicator constraint {}\n", ic.GetTypeName());
  /*MOSEK_CCALL(MOSEK_AddIndicator(lp(),
    ic.get_binary_var(), ic.get_binary_value(),
    (int)ic.get_constraint().size(),
    ic.get_constraint().pvars(),
    ic.get_constraint().pcoefs(),
    MOSEK_GREATER_EQUAL,
    ic.get_constraint().rhs()));*/

}

void MosekModelAPI::AddConstraint(const QuadConRange& qc) {
  fmt::print("Adding quadratic constraint {}\n", qc.GetTypeName());
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
/*  int type = MOSEK_SOS_TYPE1;
  int beg = 0;
  const int size = sos.size();
  MOSEK_CCALL(MOSEK_AddSOSs(lp(), 1, &type, &beg,
    &size, (int*)sos.get_vars().data(),
    (double*)sos.get_weights().data())); */
}

void MosekModelAPI::AddConstraint(const SOS2Constraint& sos) {
  fmt::print("Adding SOS1 constraint {}\n", sos.GetTypeName());
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
