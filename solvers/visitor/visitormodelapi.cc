#include "visitormodelapi.h"


namespace mp {

void VisitorModelAPI::InitProblemModificationPhase(
    const FlatModelInfo*) {
  // Allocate storage if needed:
  // auto n_linear_cons =
  //   flat_model_info->GetNumberOfConstraintsOfGroup(CG_LINEAR);
  // preallocate_linear_cons( n_linear_cons );
}

void VisitorModelAPI::AddVariables(const VarArrayDef& v) {
  // TODO Add variables using solver API; typically,
  // first convert the variable
  // type to the appropriate solver-defined type,
  // then add them as a bunch
  lp()->allocateVars(v.size());
  int nInteger = 0, nContinuous = 0;
  for (auto i = v.size(); i--; )
    lp()->setVariable(i, v.ptype()[i] == var::Type::INTEGER);
  
  fmt::print("Adding {} continuous and {} integer variables.\n",
    lp()->getNumVars(false), lp()->getNumVars(true));
  
  /* Typical implementation
  for (size_t i=v.size(); i--; )
    vtypes[i] = var::Type::CONTINUOUS == v.ptype()[i] ?
          VISITOR_CONTINUOUS : VISITOR_INTEGER;
  VISITOR_CCALL(VISITOR_AddCols(lp(), (int)v.size(), NULL, NULL,
    NULL, NULL, NULL, vtypes.data(), v.plb(), v.pub(), v.pnames())); */
}

void VisitorModelAPI::SetLinearObjective( int iobj, const LinearObjective& lo ) {
  lp()->addObjective();
  if (iobj<1) {

    fmt::print("Setting first linear objective \"{}\": {} terms.\n",
                lo.name(), lo.num_terms());
    /*
    VISITOR_CCALL(VISITOR_SetObjSense(lp(), 
                    obj::Type::MAX==lo.obj_sense() ? VISITOR_MAXIMIZE : VISITOR_MINIMIZE) );
    VISITOR_CCALL(VISITOR_SetColObj(lp(), lo.num_terms(),
                           lo.vars().data(), lo.coefs().data()) ); */
  } else {
//    TODO If we support mutiple objectives, pass them to the solver
    fmt::print("Setting {}-th linear objective: {} terms.\n", iobj, lo.num_terms());
  }
}


void VisitorModelAPI::SetQuadraticObjective(int iobj, const QuadraticObjective& qo) {
  lp()->addQuadTerms();
  if (1 > iobj) {
    fmt::format("Setting first quadratic objective\n");
    SetLinearObjective(iobj, qo);                         // add the linear part
    const auto& qt = qo.GetQPTerms();
    fmt::format("Quadratic part is made of {} terms\n", qt.size());

    // Typical implementation
    //VISITOR_CCALL(VISITOR_SetQuadObj(lp(), qt.size(),
    //  (int*)qt.pvars1(), (int*)qt.pvars2(),
    //  (double*)qt.pcoefs()));
  }
  else {
    throw std::runtime_error("Multiple quadratic objectives not supported");
  }
}

void VisitorModelAPI::AddConstraint(const LinConRange& lc) {
  fmt::print("Adding range linear constraint \"{}\"\n", lc.name());
  fmt::print("{} <=", lc.lb());
  for (size_t i = 0; i < lc.size(); i++)
  {
    fmt::print(" {}", lc.pcoefs()[i] >= 0 ? '+' : '-');
    if (std::fabs(lc.pcoefs()[i]) != 1.0)
      fmt::print("{}*", lc.pcoefs()[i]);
    fmt::print("x{}", lc.pvars()[i]);
  }
  fmt::print(" <= \"{}\"\n", lc.ub());
  lp()->addEntity(Solver::CONS_LIN);
//  VISITOR_CCALL(VISITOR_AddRow(lp(), lc.size(), lc.pvars(), lc.pcoefs(), 
 //   NULL, lc.lb(), lc.ub(), lc.name()));
}
void VisitorModelAPI::AddConstraint(const LinConLE& lc) {
  fmt::print("Adding <= linear constraint \"{}\"\n", lc.GetName());
  lp()->addEntity(Solver::CONS_LIN);
 // char sense = VISITOR_LESS_EQUAL;
 // VISITOR_CCALL(VISITOR_AddRow(lp(), lc.size(), lc.pvars(), lc.pcoefs(),
  //  sense, lc.rhs(), 0, NULL));
}
void VisitorModelAPI::AddConstraint(const LinConEQ& lc) {
  fmt::print("Adding == linear constraint \"{}\"\n", lc.GetName());
  lp()->addEntity(Solver::CONS_LIN);
//  char sense = VISITOR_EQUAL;
// VISITOR_CCALL(VISITOR_AddRow(lp(), lc.size(), lc.pvars(), lc.pcoefs(),
//   sense, lc.rhs(), 0, NULL));
}
void VisitorModelAPI::AddConstraint(const LinConGE& lc) {
  fmt::print("Adding >= linear constraint \"{}\"\n", lc.GetName());
  lp()->addEntity(Solver::CONS_LIN);
  //char sense = VISITOR_GREATER_EQUAL;
  //VISITOR_CCALL(VISITOR_AddRow(lp(), lc.size(), lc.pvars(), lc.pcoefs(),
  //  sense, lc.rhs(), 0, NULL));
}

void VisitorModelAPI::AddConstraint(const IndicatorConstraintLinLE &ic)  {
  fmt::print("Adding indicator constraint \"{}\"\n", ic.GetName());
  lp()->addEntity(Solver::CONS_INDIC);
  /*VISITOR_CCALL(VISITOR_AddIndicator(lp(),
    ic.get_binary_var(), ic.get_binary_value(),
    (int)ic.get_constraint().size(),
    ic.get_constraint().pvars(),
    ic.get_constraint().pcoefs(),
    VISITOR_LESS_EQUAL,
    ic.get_constraint().rhs()));*/
                               
}
void VisitorModelAPI::AddConstraint(const IndicatorConstraintLinEQ &ic)  {
  fmt::print("Adding indicator constraint \"{}\"\n", ic.GetName());
  lp()->addEntity(Solver::CONS_INDIC);
  /*VISITOR_CCALL(VISITOR_AddIndicator(lp(),
    ic.get_binary_var(), ic.get_binary_value(),
    (int)ic.get_constraint().size(),
    ic.get_constraint().pvars(),
    ic.get_constraint().pcoefs(),
    VISITOR_EQUAL,
    ic.get_constraint().rhs()));*/
}
void VisitorModelAPI::AddConstraint(const IndicatorConstraintLinGE &ic)  {
  fmt::print("Adding indicator constraint \"{}\"\n", ic.GetName());
  lp()->addEntity(Solver::CONS_INDIC);
  /*VISITOR_CCALL(VISITOR_AddIndicator(lp(),
    ic.get_binary_var(), ic.get_binary_value(),
    (int)ic.get_constraint().size(),
    ic.get_constraint().pvars(),
    ic.get_constraint().pcoefs(),
    VISITOR_GREATER_EQUAL,
    ic.get_constraint().rhs()));*/

}

void VisitorModelAPI::AddConstraint(const QuadConRange& qc) {
  fmt::print("Adding quadratic constraint \"{}\"\n", qc.GetName());
  lp()->addEntity(Solver::CONS_QUAD);
  /*
  const auto& lt = qc.GetLinTerms();
  const auto& qt = qc.GetQPTerms();
  VISITOR_CCALL( VISITOR_addqrangeconstr(model(), lt.size(), (int*)lt.pvars(), (double*)lt.pcoefs(),
                          qt.size(), (int*)qt.pvars1(), (int*)qt.pvars2(),
                          (double*)qt.pcoefs(), qc.lb(), qc.ub(), NULL) );
  */
}

void VisitorModelAPI::AddConstraint( const QuadConLE& qc ) {
  fmt::print("Adding quadratic constraint \"{}\"\n", qc.GetName());
  lp()->addEntity(Solver::CONS_QUAD);
  /*
  const auto& lt = qc.GetLinTerms();
  const auto& qt = qc.GetQPTerms();
  VISITOR_CCALL( VISITOR_addqconstr(model(), lt.size(), (int*)lt.pvars(), (double*)lt.pcoefs(),
                          qt.size(), (int*)qt.pvars1(), (int*)qt.pvars2(),
                          (double*)qt.pcoefs(), VISITOR__LESS_EQUAL, qc.rhs(), NULL) );
                          */
}

void VisitorModelAPI::AddConstraint( const QuadConEQ& qc ) {
  fmt::print("Adding quadratic constraint \"{}\"\n", qc.GetName());
  lp()->addEntity(Solver::CONS_QUAD);
  /*
  const auto& lt = qc.GetLinTerms();
  const auto& qt = qc.GetQPTerms();
  VISITOR_CCALL( VISITOR_addqconstr(model(), lt.size(), (int*)lt.pvars(), (double*)lt.pcoefs(),
                          qt.size(), (int*)qt.pvars1(), (int*)qt.pvars2(),
                          (double*)qt.pcoefs(), VISITOR__EQUAL, qc.rhs(), NULL) );
                          */
}

void VisitorModelAPI::AddConstraint( const QuadConGE& qc ) {
  fmt::print("Adding quadratic constraint \"{}\"\n", qc.GetName());
  lp()->addEntity(Solver::CONS_QUAD);
  /*
  const auto& lt = qc.GetLinTerms();
  const auto& qt = qc.GetQPTerms();
  VISITOR_CCALL( VISITOR_addqconstr(model(), lt.size(), (int*)lt.pvars(), (double*)lt.pcoefs(),
                          qt.size(), (int*)qt.pvars1(), (int*)qt.pvars2(),
                          (double*)qt.pcoefs(), VISITOR__GREATER_EQUAL, qc.rhs(), NULL) );
                          */
}

void VisitorModelAPI::AddConstraint( const QuadraticConeConstraint& qc ) {
  fmt::print("Adding quadratic cone constraint \"{}\"\n", qc.GetName());
  lp()->addEntity(Solver::CONS_QUAD_CONE);
}

void VisitorModelAPI::AddConstraint(
    const RotatedQuadraticConeConstraint& qc ) {
  fmt::print("Adding rotated quadratic cone constraint \"{}\"\n", qc.GetName());
  lp()->addEntity(Solver::CONS_QUAD_CONE_ROTATED);
}

void VisitorModelAPI::AddConstraint(const SOS1Constraint& sos) {
  fmt::print("Adding SOS1 constraint \"{}\"\n", sos.GetName());
  lp()->addEntity(Solver::CONS_SOS);
/*  int type = VISITOR_SOS_TYPE1;
  int beg = 0;
  const int size = sos.size();
  VISITOR_CCALL(VISITOR_AddSOSs(lp(), 1, &type, &beg,
    &size, (int*)sos.get_vars().data(),
    (double*)sos.get_weights().data())); */
}

void VisitorModelAPI::AddConstraint(const SOS2Constraint& sos) {
  fmt::print("Adding SOS2 constraint \"{}\"\n", sos.GetName());
  lp()->addEntity(Solver::CONS_SOS);
  /*int type = VISITOR_SOS_TYPE2;
  int beg = 0;
  const int size = sos.size();
  VISITOR_CCALL(VISITOR_AddSOSs(lp(), 1, &type, &beg,
    &size, (int*)sos.get_vars().data(),
    (double*)sos.get_weights().data()));*/
}



void VisitorModelAPI::AddConstraint(const MaxConstraint& mc) {
  fmt::print("Adding Max constraint \"{}\"\n", mc.GetName());
  lp()->addEntity(Solver::CONS_MAX);
  /*
  const auto& args = mc.GetArguments();
  VISITOR_CCALL(VISITOR_addgenconstrMax(model(), mc.name(),
    mc.GetResultVar(),
    (int)args.size(), args.data(),
    MinusInfinity()));
    */
}

void VisitorModelAPI::AddConstraint(const MinConstraint& mc) {
  fmt::print("Adding Min constraint \"{}\"\n", mc.GetName());
  lp()->addEntity(Solver::CONS_MIN);
  /*
  VISITOR_CCALL(VISITOR_addgenconstrMin(model(), mc.name(),
    mc.GetResultVar(),
    (int)args.size(), args.data(),
    Infinity()));
    */
}

void VisitorModelAPI::AddConstraint(const AbsConstraint& absc) {
  fmt::print("Adding Abs constraint \"{}\"\n", absc.GetName());
  lp()->addEntity(Solver::CONS_ABS);
  /*
  const auto& args = absc.GetArguments();
  VISITOR_CCALL(VISITOR_addgenconstrAbs(model(), absc.name(),
    absc.GetResultVar(), args[0]));*/
}

void VisitorModelAPI::AddConstraint(const AndConstraint& cc) {
  fmt::print("Adding And constraint \"{}\"\n", cc.GetName());
  lp()->addEntity(Solver::CONS_AND);
  /*
  const auto& args = cc.GetArguments();
  VISITOR_CCALL(VISITOR_addgenconstrAnd(model(), cc.name(),
    cc.GetResultVar(),
    (int)args.size(), args.data()));*/
}

void VisitorModelAPI::AddConstraint(const OrConstraint& dc) {
  fmt::print("Adding Or constraint \"{}\"\n", dc.GetName());
  lp()->addEntity(Solver::CONS_OR);
}


void VisitorModelAPI::AddConstraint(const ExpConstraint& cc) {
  fmt::print("Adding Exp constraint \"{}\"\n", cc.GetName());
  lp()->addEntity(Solver::CONS_EXP);
}

void VisitorModelAPI::AddConstraint(const ExpAConstraint& cc) {
  fmt::print("Adding ExpA constraint \"{}\"\n", cc.GetName());
  lp()->addEntity(Solver::CONS_EXPA);
}

void VisitorModelAPI::AddConstraint(const LogConstraint& cc) {
  fmt::print("Adding Log constraint \"{}\"\n", cc.GetName());
  lp()->addEntity(Solver::CONS_LOG);
  /*VISITOR_CCALL(VISITOR_addgenconstrLog(model(), cc.name(),
    cc.GetArguments()[0], cc.GetResultVar(), ""));*/
}

void VisitorModelAPI::AddConstraint(const LogAConstraint& cc) {
  fmt::print("Adding LogA constraint \"{}\"\n", cc.GetName());
  lp()->addEntity(Solver::CONS_LOGA);
  /*
  VISITOR_CCALL(VISITOR_addgenconstrLogA(model(), cc.name(),
    cc.GetArguments()[0], cc.GetResultVar(), cc.GetParameters()[0], ""));*/
}

void VisitorModelAPI::AddConstraint(const PowConstraint& cc) {
  fmt::print("Adding Pow constraint \"{}\"\n", cc.GetName());
  lp()->addEntity(Solver::CONS_POW);
  /*
  VISITOR_CCALL(VISITOR_addgenconstrPow(model(), cc.name(),
    cc.GetArguments()[0], cc.GetResultVar(), cc.GetParameters()[0], ""));*/
}

void VisitorModelAPI::AddConstraint(const SinConstraint& cc) {
  fmt::print("Adding Sin constraint \"{}\"\n", cc.GetName());
  lp()->addEntity(Solver::CONS_SIN);
  /*VISITOR_CCALL(VISITOR_addgenconstrSin(model(), cc.name(),
    cc.GetArguments()[0], cc.GetResultVar(), ""));*/
}

void VisitorModelAPI::AddConstraint(const CosConstraint& cc) {
  fmt::print("Adding Cos constraint \"{}\"\n", cc.GetName());
  lp()->addEntity(Solver::CONS_COS);
  /*
  VISITOR_CCALL(VISITOR_addgenconstrCos(model(), cc.name(),
    cc.GetArguments()[0], cc.GetResultVar(), ""));*/
}

void VisitorModelAPI::AddConstraint(const TanConstraint& cc) {
  fmt::print("Adding Tan constraint \"{}\"\n", cc.GetName());
  lp()->addEntity(Solver::CONS_TAN);
  
  /*VISITOR_CCALL(VISITOR_addgenconstrTan(model(), cc.name(),
    cc.GetArguments()[0], cc.GetResultVar(), ""));*/
}

void VisitorModelAPI::AddConstraint(const PLConstraint& plc) {
  fmt::print("Adding PL constraint \"{}\"\n", plc.GetName());
  lp()->addEntity(Solver::CONS_PL);
  /*
  const auto& plp = plc.GetParameters().GetPLPoints();
  VISITOR_CCALL(VISITOR_addgenconstrPWL(model(), plc.name(),
    plc.GetArguments()[0], plc.GetResultVar(),
    plp.x_.size(), plp.x_.data(), plp.y_.data()));*/
}


void VisitorModelAPI::FinishProblemModificationPhase() {
}


} // namespace mp
