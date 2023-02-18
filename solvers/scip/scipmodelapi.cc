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
  // TODO Add variables using solver API; typically,
  // first convert the variable
  // type to the appropriate solver-defined type,
  // then add them as a bunch
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
          SCIP_CONTINUOUS : SCIP_INTEGER;
  SCIP_CCALL(SCIP_AddCols(lp(), (int)v.size(), NULL, NULL,
    NULL, NULL, NULL, vtypes.data(), v.plb(), v.pub(), v.pnames())); */
}

void ScipModelAPI::SetLinearObjective( int iobj, const LinearObjective& lo ) {
  if (iobj<1) {
    fmt::format("Setting first linear objective \"{}\": {} terms.\n",
                lo.name(), lo.num_terms());
    /*
    SCIP_CCALL(SCIP_SetObjSense(lp(), 
                    obj::Type::MAX==lo.obj_sense() ? SCIP_MAXIMIZE : SCIP_MINIMIZE) );
    SCIP_CCALL(SCIP_SetColObj(lp(), lo.num_terms(),
                           lo.vars().data(), lo.coefs().data()) ); */
  } else {
//    TODO If we support mutiple objectives, pass them to the solver
    fmt::format("Setting {}-th linear objective: {} terms.\n", iobj, lo.num_terms());
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
  fmt::print("Adding range linear constraint {}\n", lc.name());
  fmt::print("{} <=", lc.lb());
  for (size_t i = 0; i < lc.size(); i++)
  {
    fmt::print(" {}", lc.pcoefs()[i] >= 0 ? '+' : '-');
    if (std::fabs(lc.pcoefs()[i]) != 1.0)
      fmt::print("{}*", lc.pcoefs()[i]);
    fmt::print("x{}", lc.pvars()[i]);
  }
  fmt::print(" <= {}\n", lc.ub());
//  SCIP_CCALL(SCIP_AddRow(lp(), lc.size(), lc.pvars(), lc.pcoefs(), 
 //   NULL, lc.lb(), lc.ub(), lc.name()));
}
void ScipModelAPI::AddConstraint(const LinConLE& lc) {
  fmt::print("Adding <= linear constraint {}\n", lc.GetName());
 // char sense = SCIP_LESS_EQUAL;
 // SCIP_CCALL(SCIP_AddRow(lp(), lc.size(), lc.pvars(), lc.pcoefs(),
  //  sense, lc.rhs(), 0, NULL));
}
void ScipModelAPI::AddConstraint(const LinConEQ& lc) {
  fmt::print("Adding == linear constraint {}\n", lc.GetName());
//  char sense = SCIP_EQUAL;
// SCIP_CCALL(SCIP_AddRow(lp(), lc.size(), lc.pvars(), lc.pcoefs(),
//   sense, lc.rhs(), 0, NULL));
}
void ScipModelAPI::AddConstraint(const LinConGE& lc) {
  fmt::print("Adding >= linear constraint {}\n", lc.GetName());
  //char sense = SCIP_GREATER_EQUAL;
  //SCIP_CCALL(SCIP_AddRow(lp(), lc.size(), lc.pvars(), lc.pcoefs(),
  //  sense, lc.rhs(), 0, NULL));
}

void ScipModelAPI::AddConstraint(const IndicatorConstraintLinLE &ic)  {
  fmt::print("Adding indicator constraint {}\n", ic.GetName());
  /*SCIP_CCALL(SCIP_AddIndicator(lp(),
    ic.get_binary_var(), ic.get_binary_value(),
    (int)ic.get_constraint().size(),
    ic.get_constraint().pvars(),
    ic.get_constraint().pcoefs(),
    SCIP_LESS_EQUAL,
    ic.get_constraint().rhs()));*/
                               
}
void ScipModelAPI::AddConstraint(const IndicatorConstraintLinEQ &ic)  {
  fmt::print("Adding indicator constraint {}\n", ic.GetName());
  /*SCIP_CCALL(SCIP_AddIndicator(lp(),
    ic.get_binary_var(), ic.get_binary_value(),
    (int)ic.get_constraint().size(),
    ic.get_constraint().pvars(),
    ic.get_constraint().pcoefs(),
    SCIP_EQUAL,
    ic.get_constraint().rhs()));*/
}
void ScipModelAPI::AddConstraint(const IndicatorConstraintLinGE &ic)  {
  fmt::print("Adding indicator constraint {}\n", ic.GetName());
  /*SCIP_CCALL(SCIP_AddIndicator(lp(),
    ic.get_binary_var(), ic.get_binary_value(),
    (int)ic.get_constraint().size(),
    ic.get_constraint().pvars(),
    ic.get_constraint().pcoefs(),
    SCIP_GREATER_EQUAL,
    ic.get_constraint().rhs()));*/

}

void ScipModelAPI::AddConstraint(const QuadConRange& qc) {
  fmt::print("Adding quadratic constraint {}\n", qc.GetName());
  /*
  const auto& lt = qc.GetLinTerms();
  const auto& qt = qc.GetQPTerms();
  GRB_CALL( GRBaddqrangeconstr(model(), lt.size(), (int*)lt.pvars(), (double*)lt.pcoefs(),
                          qt.size(), (int*)qt.pvars1(), (int*)qt.pvars2(),
                          (double*)qt.pcoefs(), qc.lb(), qc.ub(), NULL) );
  */
}

void ScipModelAPI::AddConstraint( const QuadConLE& qc ) {
  fmt::print("Adding quadratic constraint {}\n", qc.GetName());
  /*
  const auto& lt = qc.GetLinTerms();
  const auto& qt = qc.GetQPTerms();
  GRB_CALL( GRBaddqconstr(model(), lt.size(), (int*)lt.pvars(), (double*)lt.pcoefs(),
                          qt.size(), (int*)qt.pvars1(), (int*)qt.pvars2(),
                          (double*)qt.pcoefs(), GRB_LESS_EQUAL, qc.rhs(), NULL) );
                          */
}

void ScipModelAPI::AddConstraint( const QuadConEQ& qc ) {
  fmt::print("Adding quadratic constraint {}\n", qc.GetName());
  /*
  const auto& lt = qc.GetLinTerms();
  const auto& qt = qc.GetQPTerms();
  GRB_CALL( GRBaddqconstr(model(), lt.size(), (int*)lt.pvars(), (double*)lt.pcoefs(),
                          qt.size(), (int*)qt.pvars1(), (int*)qt.pvars2(),
                          (double*)qt.pcoefs(), GRB_EQUAL, qc.rhs(), NULL) );
                          */
}

void ScipModelAPI::AddConstraint( const QuadConGE& qc ) {
  fmt::print("Adding quadratic constraint {}\n", qc.GetName());
  /*
  const auto& lt = qc.GetLinTerms();
  const auto& qt = qc.GetQPTerms();
  GRB_CALL( GRBaddqconstr(model(), lt.size(), (int*)lt.pvars(), (double*)lt.pcoefs(),
                          qt.size(), (int*)qt.pvars1(), (int*)qt.pvars2(),
                          (double*)qt.pcoefs(), GRB_GREATER_EQUAL, qc.rhs(), NULL) );
                          */
}

void ScipModelAPI::AddConstraint(const SOS1Constraint& sos) {
  fmt::print("Adding SOS1 constraint {}\n", sos.GetName());
/*  int type = SCIP_SOS_TYPE1;
  int beg = 0;
  const int size = sos.size();
  SCIP_CCALL(SCIP_AddSOSs(lp(), 1, &type, &beg,
    &size, (int*)sos.get_vars().data(),
    (double*)sos.get_weights().data())); */
}

void ScipModelAPI::AddConstraint(const SOS2Constraint& sos) {
  fmt::print("Adding SOS1 constraint {}\n", sos.GetName());
  /*int type = SCIP_SOS_TYPE2;
  int beg = 0;
  const int size = sos.size();
  SCIP_CCALL(SCIP_AddSOSs(lp(), 1, &type, &beg,
    &size, (int*)sos.get_vars().data(),
    (double*)sos.get_weights().data()));*/
}


void ScipModelAPI::FinishProblemModificationPhase() {
}


} // namespace mp
