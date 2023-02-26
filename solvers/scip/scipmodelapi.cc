#include "scipmodelapi.h"

#include "scip/cons_linear.h"

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
    SCIP_VARTYPE vartype = var::Type::INTEGER == v.ptype()[i] ? SCIP_VARTYPE_INTEGER : SCIP_VARTYPE_CONTINUOUS; //implement bin and impint
    //const char* name = v.pnames()[i];
    SCIP_CCALL( SCIPcreateVarBasic(getSCIP(), &var, NULL, lb, ub, objcoef, vartype) );
    SCIP_CCALL( SCIPaddVar(getSCIP(), var) );
    getPROBDATA()->vars[i] = var;
    //SCIP_CCALL( SCIPreleaseVar(getSCIP(), &var) );
    //vars_.push_back(var);
  }
}

void ScipModelAPI::SetLinearObjective( int iobj, const LinearObjective& lo ) {
  if (iobj<1) {
    //fmt::format("Setting linear objective \"{}\": {} terms.\n",
    //            lo.name(), lo.num_terms());
    
    SCIP_CCALL( SCIPsetObjsense(getSCIP(), 
                    obj::Type::MAX==lo.obj_sense() ? SCIP_OBJSENSE_MAXIMIZE : SCIP_OBJSENSE_MINIMIZE) );
    SCIP_VAR** vars = getPROBDATA()->vars;
    for (int i = 0; i < lo.num_terms(); i++) {
      SCIP_CCALL( SCIPchgVarObj(getSCIP(), vars[i], lo.coefs()[i]) );
    }
  } else {
    fmt::format("Multiple linear objectives not supported");
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
  //fmt::print("Adding range linear constraint {}\n", lc.name());
  //fmt::print("{} <=", lc.lb());
  //for (size_t i = 0; i < lc.size(); i++) {
  //  fmt::print(" {}", lc.pcoefs()[i] >= 0 ? '+' : '-');
  //  if (std::fabs(lc.pcoefs()[i]) != 1.0)
  //    fmt::print("{}*", lc.pcoefs()[i]);
  //  fmt::print("x{}", lc.pvars()[i]);
  //}
  //fmt::print(" <= {}\n", lc.ub());

  SCIP_VAR** vars = NULL;
  SCIP_CCALL( SCIPallocBlockMemoryArray(getSCIP(), &vars, lc.size()) );
  for (size_t i = 0; i < lc.size(); i++) {
    vars[i] = getPROBDATA()->vars[lc.pvars()[i]];
  }
  SCIP_CONS* cons;
  SCIP_CCALL( SCIPcreateConsBasicLinear(getSCIP(), &cons, lc.GetName(), (int)lc.size(), vars, (double*)lc.pcoefs(), lc.lb(), lc.ub()) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), cons) );
  SCIP_CCALL( SCIPreleaseCons(getSCIP(), &cons) );
//  SCIP_CCALL(SCIP_AddRow(lp(), lc.size(), lc.pvars(), lc.pcoefs(), 
 //   NULL, lc.lb(), lc.ub(), lc.name()));

  SCIPfreeBlockMemoryArrayNull(getSCIP(), &vars, lc.size());
}
void ScipModelAPI::AddConstraint(const LinConLE& lc) {
  //fmt::print("Adding <= linear constraint {}\n", lc.GetName());

  SCIP_VAR** vars = NULL;
  SCIP_CCALL( SCIPallocBlockMemoryArray(getSCIP(), &vars, lc.size()) );
  for (size_t i = 0; i < lc.size(); i++) {
    vars[i] = getPROBDATA()->vars[lc.pvars()[i]];
  }
  SCIP_CONS* cons;
  SCIP_CCALL( SCIPcreateConsBasicLinear(getSCIP(), &cons, lc.GetName(), (int)lc.size(), vars, (double*)lc.pcoefs(), MinusInfinity(), lc.rhs()) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), cons) );
  SCIP_CCALL( SCIPreleaseCons(getSCIP(), &cons) );

//  SCIP_CCALL(SCIP_AddRow(lp(), lc.size(), lc.pvars(), lc.pcoefs(), 
 //   NULL, lc.lb(), lc.ub(), lc.name()));

  SCIPfreeBlockMemoryArrayNull(getSCIP(), &vars, lc.size());
 // char sense = SCIP_LESS_EQUAL;
 // SCIP_CCALL(SCIP_AddRow(lp(), lc.size(), lc.pvars(), lc.pcoefs(),
  //  sense, lc.rhs(), 0, NULL));
}
void ScipModelAPI::AddConstraint(const LinConEQ& lc) {
  //fmt::print("Adding == linear constraint {}\n", lc.GetName());

  SCIP_VAR** vars = NULL;
  SCIP_CCALL( SCIPallocBlockMemoryArray(getSCIP(), &vars, lc.size()) );
  for (size_t i = 0; i < lc.size(); i++) {
    vars[i] = getPROBDATA()->vars[lc.pvars()[i]];
  }
  SCIP_CONS* cons;
  SCIP_CCALL( SCIPcreateConsBasicLinear(getSCIP(), &cons, lc.GetName(), (int)lc.size(), vars, (double*)lc.pcoefs(), lc.rhs(), lc.rhs()) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), cons) );
  SCIP_CCALL( SCIPreleaseCons(getSCIP(), &cons) );

//  SCIP_CCALL(SCIP_AddRow(lp(), lc.size(), lc.pvars(), lc.pcoefs(), 
 //   NULL, lc.lb(), lc.ub(), lc.name()));

  SCIPfreeBlockMemoryArrayNull(getSCIP(), &vars, lc.size());
//  char sense = SCIP_EQUAL;
// SCIP_CCALL(SCIP_AddRow(lp(), lc.size(), lc.pvars(), lc.pcoefs(),
//   sense, lc.rhs(), 0, NULL));
}
void ScipModelAPI::AddConstraint(const LinConGE& lc) {
  //fmt::print("Adding >= linear constraint {}\n", lc.GetName());

  SCIP_VAR** vars = NULL;
  SCIP_CCALL( SCIPallocBlockMemoryArray(getSCIP(), &vars, lc.size()) );
  for (size_t i = 0; i < lc.size(); i++) {
    vars[i] = getPROBDATA()->vars[lc.pvars()[i]];
  }
  SCIP_CONS* cons;
  SCIP_CCALL( SCIPcreateConsBasicLinear(getSCIP(), &cons, lc.GetName(), (int)lc.size(), vars, (double*)lc.pcoefs(), lc.lb(), Infinity()) );
  SCIP_CCALL( SCIPaddCons(getSCIP(), cons) );
  SCIP_CCALL( SCIPreleaseCons(getSCIP(), &cons) );

//  SCIP_CCALL(SCIP_AddRow(lp(), lc.size(), lc.pvars(), lc.pcoefs(), 
 //   NULL, lc.lb(), lc.ub(), lc.name()));

  SCIPfreeBlockMemoryArrayNull(getSCIP(), &vars, lc.size());
  //char sense = SCIP_GREATER_EQUAL;
  //SCIP_CCALL(SCIP_AddRow(lp(), lc.size(), lc.pvars(), lc.pcoefs(),
  //  sense, lc.rhs(), 0, NULL));
}

//void ScipModelAPI::AddConstraint(const IndicatorConstraintLinLE &ic)  {
//  fmt::print("Adding indicator constraint {}\n", ic.GetName());
  /*SCIP_CCALL(SCIP_AddIndicator(lp(),
    ic.get_binary_var(), ic.get_binary_value(),
    (int)ic.get_constraint().size(),
    ic.get_constraint().pvars(),
    ic.get_constraint().pcoefs(),
    SCIP_LESS_EQUAL,
    ic.get_constraint().rhs()));*/
                               
//}
//void ScipModelAPI::AddConstraint(const IndicatorConstraintLinEQ &ic)  {
//  fmt::print("Adding indicator constraint {}\n", ic.GetName());
  /*SCIP_CCALL(SCIP_AddIndicator(lp(),
    ic.get_binary_var(), ic.get_binary_value(),
    (int)ic.get_constraint().size(),
    ic.get_constraint().pvars(),
    ic.get_constraint().pcoefs(),
    SCIP_EQUAL,
    ic.get_constraint().rhs()));*/
//}
//void ScipModelAPI::AddConstraint(const IndicatorConstraintLinGE &ic)  {
//  fmt::print("Adding indicator constraint {}\n", ic.GetName());
  /*SCIP_CCALL(SCIP_AddIndicator(lp(),
    ic.get_binary_var(), ic.get_binary_value(),
    (int)ic.get_constraint().size(),
    ic.get_constraint().pvars(),
    ic.get_constraint().pcoefs(),
    SCIP_GREATER_EQUAL,
    ic.get_constraint().rhs()));*/

//}

//void ScipModelAPI::AddConstraint(const QuadConRange& qc) {
//  fmt::print("Adding quadratic constraint {}\n", qc.GetName());
  /*
  const auto& lt = qc.GetLinTerms();
  const auto& qt = qc.GetQPTerms();
  GRB_CALL( GRBaddqrangeconstr(model(), lt.size(), (int*)lt.pvars(), (double*)lt.pcoefs(),
                          qt.size(), (int*)qt.pvars1(), (int*)qt.pvars2(),
                          (double*)qt.pcoefs(), qc.lb(), qc.ub(), NULL) );
  */
//}

//void ScipModelAPI::AddConstraint( const QuadConLE& qc ) {
//  fmt::print("Adding quadratic constraint {}\n", qc.GetName());
  /*
  const auto& lt = qc.GetLinTerms();
  const auto& qt = qc.GetQPTerms();
  GRB_CALL( GRBaddqconstr(model(), lt.size(), (int*)lt.pvars(), (double*)lt.pcoefs(),
                          qt.size(), (int*)qt.pvars1(), (int*)qt.pvars2(),
                          (double*)qt.pcoefs(), GRB_LESS_EQUAL, qc.rhs(), NULL) );
                          */
//}

//void ScipModelAPI::AddConstraint( const QuadConEQ& qc ) {
//  fmt::print("Adding quadratic constraint {}\n", qc.GetName());
  /*
  const auto& lt = qc.GetLinTerms();
  const auto& qt = qc.GetQPTerms();
  GRB_CALL( GRBaddqconstr(model(), lt.size(), (int*)lt.pvars(), (double*)lt.pcoefs(),
                          qt.size(), (int*)qt.pvars1(), (int*)qt.pvars2(),
                          (double*)qt.pcoefs(), GRB_EQUAL, qc.rhs(), NULL) );
                          */
//}

//void ScipModelAPI::AddConstraint( const QuadConGE& qc ) {
//  fmt::print("Adding quadratic constraint {}\n", qc.GetName());
  /*
  const auto& lt = qc.GetLinTerms();
  const auto& qt = qc.GetQPTerms();
  GRB_CALL( GRBaddqconstr(model(), lt.size(), (int*)lt.pvars(), (double*)lt.pcoefs(),
                          qt.size(), (int*)qt.pvars1(), (int*)qt.pvars2(),
                          (double*)qt.pcoefs(), GRB_GREATER_EQUAL, qc.rhs(), NULL) );
                          */
//}

//void ScipModelAPI::AddConstraint(const SOS1Constraint& sos) {
//  fmt::print("Adding SOS1 constraint {}\n", sos.GetName());
/*  int type = SCIP_SOS_TYPE1;
  int beg = 0;
  const int size = sos.size();
  SCIP_CCALL(SCIP_AddSOSs(lp(), 1, &type, &beg,
    &size, (int*)sos.get_vars().data(),
    (double*)sos.get_weights().data())); */
//}

//void ScipModelAPI::AddConstraint(const SOS2Constraint& sos) {
//  fmt::print("Adding SOS1 constraint {}\n", sos.GetName());
  /*int type = SCIP_SOS_TYPE2;
  int beg = 0;
  const int size = sos.size();
  SCIP_CCALL(SCIP_AddSOSs(lp(), 1, &type, &beg,
    &size, (int*)sos.get_vars().data(),
    (double*)sos.get_weights().data()));*/
//}


void ScipModelAPI::FinishProblemModificationPhase() {
}


} // namespace mp
