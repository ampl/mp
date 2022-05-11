#include "visitormodelapi.h"

#include "mp/model-mgr-with-std-pb.h"
#include "mp/flat/problem_flattener.h"
#include "mp/flat/redef/MIP/converter_mip.h"
#include "mp/flat/model_api_connect.h"

namespace mp {

/// Defining the function in ...modelapi.cc
/// for recompilation speed
std::unique_ptr<BasicModelManager>
CreateVisitorModelMgr(VisitorCommon& cc, Env& e,
                     pre::BasicPresolver*& pPre) {
  return CreateModelMgrWithFlatConverter<
      VisitorModelAPI, MIPFlatConverter >(cc, e, pPre);
}


void VisitorModelAPI::InitProblemModificationPhase() { }

void VisitorModelAPI::AddVariables(const VarArrayDef& v) {
  // TODO Add variables using solver API; typically, first convert the variable
  // type to the appropriate solver-defined type, then add them as a bunch
  int nInteger = 0, nContinuous = 0;
  for (size_t i = v.size(); i--; )
    if (var::Type::CONTINUOUS == v.ptype()[i])
      nContinuous++;
    else
      nInteger++;
  fmt::print("Adding {} continuous and {} integer variables.\n", nContinuous, nInteger);

  /* Typical implementation
  for (size_t i=v.size(); i--; )
    vtypes[i] = var::Type::CONTINUOUS == v.ptype()[i] ?
          VISITOR_CONTINUOUS : VISITOR_INTEGER;
  VISITOR_CCALL(VISITOR_AddCols(lp(), (int)v.size(), NULL, NULL,
    NULL, NULL, NULL, vtypes.data(), v.plb(), v.pub(),  NULL)); */
}

void VisitorModelAPI::SetLinearObjective( int iobj, const LinearObjective& lo ) {
  if (iobj<1) {
    fmt::format("Setting first linear objective: {} terms.\n", lo.num_terms());
    /*
    VISITOR_CCALL(VISITOR_SetObjSense(lp(), 
                    obj::Type::MAX==lo.obj_sense() ? VISITOR_MAXIMIZE : VISITOR_MINIMIZE) );
    VISITOR_CCALL(VISITOR_SetColObj(lp(), lo.num_terms(),
                           lo.vars().data(), lo.coefs().data()) ); */
  } else {
//    TODO If we support mutiple objectives, pass them to the solver
    fmt::format("Setting {}-th linear objective: {} terms.\n", iobj, lo.num_terms());
  }
}


void VisitorModelAPI::SetQuadraticObjective(int iobj, const QuadraticObjective& qo) {
  if (1 > iobj) {
    fmt::format("Setting first quadratic objective\n");
    SetLinearObjective(iobj, qo);                         // add the linear part
    const auto& qt = qo.GetQPTerms();
    fmt::format("Quadratic part is made of {} terms\n", qt.size());

    // Typical implementation
    //VISITOR_CCALL(VISITOR_SetQuadObj(lp(), qt.size(), (int*)qt.pvars1(), (int*)qt.pvars2(),
    //  (double*)qt.pcoefs()));
  }
  else {
    throw std::runtime_error("Multiple quadratic objectives not supported");
  }
}

void VisitorModelAPI::AddConstraint(const LinConRange& lc) {
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
//  VISITOR_CCALL(VISITOR_AddRow(lp(), lc.size(), lc.pvars(), lc.pcoefs(), 
 //   NULL, lc.lb(), lc.ub(), NULL));
}
void VisitorModelAPI::AddConstraint(const LinConLE& lc) {
  fmt::print("Adding <= linear constraint {}\n", lc.GetTypeName());
 // char sense = VISITOR_LESS_EQUAL;
 // VISITOR_CCALL(VISITOR_AddRow(lp(), lc.size(), lc.pvars(), lc.pcoefs(),
  //  sense, lc.rhs(), 0, NULL));
}
void VisitorModelAPI::AddConstraint(const LinConEQ& lc) {
  fmt::print("Adding == linear constraint {}\n", lc.GetTypeName());
//  char sense = VISITOR_EQUAL;
// VISITOR_CCALL(VISITOR_AddRow(lp(), lc.size(), lc.pvars(), lc.pcoefs(),
//   sense, lc.rhs(), 0, NULL));
}
void VisitorModelAPI::AddConstraint(const LinConGE& lc) {
  fmt::print("Adding >= linear constraint {}\n", lc.GetTypeName());
  //char sense = VISITOR_GREATER_EQUAL;
  //VISITOR_CCALL(VISITOR_AddRow(lp(), lc.size(), lc.pvars(), lc.pcoefs(),
  //  sense, lc.rhs(), 0, NULL));
}

void VisitorModelAPI::AddConstraint(const IndicatorConstraintLinLE &ic)  {
  fmt::print("Adding indicator constraint {}\n", ic.GetTypeName());
  /*VISITOR_CCALL(VISITOR_AddIndicator(lp(),
    ic.get_binary_var(), ic.get_binary_value(),
    (int)ic.get_constraint().size(),
    ic.get_constraint().pvars(),
    ic.get_constraint().pcoefs(),
    VISITOR_LESS_EQUAL,
    ic.get_constraint().rhs()));*/
                               
}
void VisitorModelAPI::AddConstraint(const IndicatorConstraintLinEQ &ic)  {
  fmt::print("Adding indicator constraint {}\n", ic.GetTypeName());
  /*VISITOR_CCALL(VISITOR_AddIndicator(lp(),
    ic.get_binary_var(), ic.get_binary_value(),
    (int)ic.get_constraint().size(),
    ic.get_constraint().pvars(),
    ic.get_constraint().pcoefs(),
    VISITOR_EQUAL,
    ic.get_constraint().rhs()));*/
}
void VisitorModelAPI::AddConstraint(const IndicatorConstraintLinGE &ic)  {
  fmt::print("Adding indicator constraint {}\n", ic.GetTypeName());
  /*VISITOR_CCALL(VISITOR_AddIndicator(lp(),
    ic.get_binary_var(), ic.get_binary_value(),
    (int)ic.get_constraint().size(),
    ic.get_constraint().pvars(),
    ic.get_constraint().pcoefs(),
    VISITOR_GREATER_EQUAL,
    ic.get_constraint().rhs()));*/

}

void VisitorModelAPI::AddConstraint(const QuadConRange& qc) {
  fmt::print("Adding quadratic constraint {}\n", qc.GetTypeName());
  /*
  const auto& lt = qc.GetLinTerms();
  const auto& qt = qc.GetQPTerms();
  GRB_CALL( GRBaddqrangeconstr(model(), lt.size(), (int*)lt.pvars(), (double*)lt.pcoefs(),
                          qt.size(), (int*)qt.pvars1(), (int*)qt.pvars2(),
                          (double*)qt.pcoefs(), qc.lb(), qc.ub(), NULL) );
  */
}

void VisitorModelAPI::AddConstraint( const QuadConLE& qc ) {
  fmt::print("Adding quadratic constraint {}\n", qc.GetTypeName());
  /*
  const auto& lt = qc.GetLinTerms();
  const auto& qt = qc.GetQPTerms();
  GRB_CALL( GRBaddqconstr(model(), lt.size(), (int*)lt.pvars(), (double*)lt.pcoefs(),
                          qt.size(), (int*)qt.pvars1(), (int*)qt.pvars2(),
                          (double*)qt.pcoefs(), GRB_LESS_EQUAL, qc.rhs(), NULL) );
                          */
}

void VisitorModelAPI::AddConstraint( const QuadConEQ& qc ) {
  fmt::print("Adding quadratic constraint {}\n", qc.GetTypeName());
  /*
  const auto& lt = qc.GetLinTerms();
  const auto& qt = qc.GetQPTerms();
  GRB_CALL( GRBaddqconstr(model(), lt.size(), (int*)lt.pvars(), (double*)lt.pcoefs(),
                          qt.size(), (int*)qt.pvars1(), (int*)qt.pvars2(),
                          (double*)qt.pcoefs(), GRB_EQUAL, qc.rhs(), NULL) );
                          */
}

void VisitorModelAPI::AddConstraint( const QuadConGE& qc ) {
  fmt::print("Adding quadratic constraint {}\n", qc.GetTypeName());
  /*
  const auto& lt = qc.GetLinTerms();
  const auto& qt = qc.GetQPTerms();
  GRB_CALL( GRBaddqconstr(model(), lt.size(), (int*)lt.pvars(), (double*)lt.pcoefs(),
                          qt.size(), (int*)qt.pvars1(), (int*)qt.pvars2(),
                          (double*)qt.pcoefs(), GRB_GREATER_EQUAL, qc.rhs(), NULL) );
                          */
}

void VisitorModelAPI::AddConstraint(const SOS1Constraint& sos) {
  fmt::print("Adding SOS1 constraint {}\n", sos.GetTypeName());
/*  int type = VISITOR_SOS_TYPE1;
  int beg = 0;
  const int size = sos.size();
  VISITOR_CCALL(VISITOR_AddSOSs(lp(), 1, &type, &beg,
    &size, (int*)sos.get_vars().data(),
    (double*)sos.get_weights().data())); */
}

void VisitorModelAPI::AddConstraint(const SOS2Constraint& sos) {
  fmt::print("Adding SOS1 constraint {}\n", sos.GetTypeName());
  /*int type = VISITOR_SOS_TYPE2;
  int beg = 0;
  const int size = sos.size();
  VISITOR_CCALL(VISITOR_AddSOSs(lp(), 1, &type, &beg,
    &size, (int*)sos.get_vars().data(),
    (double*)sos.get_weights().data()));*/
}


void VisitorModelAPI::FinishProblemModificationPhase() {
}


} // namespace mp
