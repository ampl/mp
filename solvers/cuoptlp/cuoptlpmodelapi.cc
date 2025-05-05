#include "cuoptlpmodelapi.h"


namespace mp {

void CuoptlpModelAPI::InitProblemModificationPhase(
    const FlatModelInfo*) {
  auto varname = std::bind(&Solver::SolverModel::var_name, this->lp(),
    std::placeholders::_1);
  GetVarName = varname;
  // Allocate storage if needed:
  // auto n_linear_cons =
  //   flat_model_info->GetNumberOfConstraintsOfGroup(CG_LINEAR);
  // reallocate_linear_cons( n_linear_cons );
}

void CuoptlpModelAPI::AddVariables(const VarArrayDef& v) {
  // TODO Add variables using solver API; typically,
  // first convert the MP variable type to the appropriate solver-defined type,
  // then add them 

  // Preallocate in solver mem
  lp()->allocateVars(v.size());
  // Assign one by one
  for (auto i = 0; i<v.size();  i++)
    lp()->AddVariable(v.ptype()[i],
      v.plb()[i], v.pub()[i]), 
      v.pnames() == NULL ? nullptr : v.pnames()[i];
  
  fmt::print("Added {} continuous, {} integer and {} binary variables.\n",
    lp()->getNumVars(Solver::VarType::CONTINUOUS), 
    lp()->getNumVars(Solver::VarType::INTEGER),
    lp()->getNumVars(Solver::VarType::BINARY));
  
  /* Typical implementation when passing all the arrays
  for (size_t i=v.size(); i--; )
    vtypes[i] = var::Type::CONTINUOUS == v.ptype()[i] ?
          CUOPTLP_CONTINUOUS : CUOPTLP_INTEGER;
  CUOPTLP_CCALL(CUOPTLP_AddCols(lp(), (int)v.size(), NULL, NULL,
    NULL, NULL, NULL, vtypes.data(), v.plb(), v.pub(), v.pnames())); */
}


void PrintCoefficient(double value, bool first) {
  if (value > 0) {
    if (!first) fmt::print("+");
    // Omit '1*' for positive coefficient 1
    if (value != 1)
      fmt::print("{}*", value);
  }
  else { // Coefficient is negative
    // If the coefficient is -1, just print ' - ', otherwise print the coefficient
    if (value == -1) fmt::print("-");
    else  fmt::print("{}*", value);
  }
}
void PrintQuadTerm(double coeff, const char* v1, const char* v2 = nullptr, bool first = false)
{
  PrintCoefficient(coeff, first);
  if (v2)
    fmt::print("{}*{}", v1, v2);
  else
    fmt::print("{}^2", v1);
}
void PrintQuadBody(std::function<std::string_view(int)> varName,
  int size, const int* vars1, const int* vars2, const double* values)
{
  for (int i = 0; i < size; i++) {
    PrintQuadTerm(values[i], varName(vars1[i]).data(),
      vars1[i] != vars2[i] ? varName(vars2[i]).data() : nullptr,
      i == 0);
  }
}

void PrintConsName(std::string_view name) {
  fmt::print("{}: ", name.data());
}
void PrintLhs(double lhs, double rhs) {
  bool hasLhs = lhs != -std::numeric_limits<double>::infinity();
  bool hasRhs = rhs != std::numeric_limits<double>::infinity();
  if ((hasLhs && hasRhs) && (lhs != rhs))
    fmt::print("{} <= ", lhs);

}
void PrintRhs(double lhs, double rhs) {
  bool hasLhs = lhs != -std::numeric_limits<double>::infinity();
  bool hasRhs = rhs != std::numeric_limits<double>::infinity();
  if (hasRhs)
    fmt::print("{}{}", ((lhs != rhs) ? "<=" : "=="), rhs);
  else if (hasLhs) // range
    fmt::print(" >= {}", lhs);
  fmt::print("\n");
}
void PrintLinearBody(std::function<std::string_view(int)> varName,
  int size, const int* vars, const double* values, bool first=true) {
  for (auto i = 0; i < size; i++) {
    PrintCoefficient(values[i], (i == 0) && first);
    fmt::print(varName(vars[i]).data());
   }
}

void PrintLinearConstraint(std::string_view name, std::function<std::string_view(int)> varName,
  int size, const int* vars, const double* values, double lhs, double rhs) {
  PrintConsName(name);
  PrintLhs(lhs, rhs);
  PrintLinearBody(varName, size, vars, values);
  PrintRhs(lhs, rhs);
}
void PrintLinearObjective(std::string_view name, std::function<std::string_view(int)> varName,
  int size, const int* vars, const double* values, bool maximize) {
  PrintConsName(name);
  std::string direction = maximize ? "maximize" : "minimize";
  fmt::print("{} ", direction);
  PrintLinearBody(varName, size, vars, values);
  fmt::print("\n");
}
void PrintQuadraticObjective(std::string_view name, std::function<std::string_view(int)> varName,
  int size, const int* vars, const double* values, bool maximize,
  int sizequad, const int* vars1, const int* vars2, const double* valuesquad) {
  PrintConsName(name);
  std::string direction = maximize ? "maximize" : "minimize";
  fmt::print("{} ", direction);
  PrintQuadBody(varName, sizequad, vars1, vars2, valuesquad);
  PrintLinearBody(varName, size, vars, values, false);
  fmt::print("\n");
}
void PrintQuadraticConstraint(std::string_view name, std::function<std::string_view(int)> varName,
  int size, const int* vars, const double* values,
  double lhs, double rhs,
  int sizequad, const int* vars1, const int* vars2, const double* valuesquad)
{
  PrintConsName(name);
  PrintLhs(lhs, rhs);
  PrintQuadBody(varName, sizequad, vars1, vars2, valuesquad);
  PrintLinearBody(varName, size, vars, values, false);
  PrintRhs(lhs, rhs);
}

void PrintIndicator(std::string_view name,
  std::function<std::string_view(int)> varName,
  int binaryVar, int binaryValue, int size, const int* vars, const double* values,
  double lhs, double rhs) {
  PrintConsName(name);
  // print condition
  fmt::print("{}=={} ==> ", varName(binaryVar).data(), binaryValue);
  PrintLhs(lhs, rhs);
  PrintLinearBody(varName, size, vars, values);
  PrintRhs(lhs, rhs);

}
void PrintFunctionalConstraintNoParam(std::string_view name, std::function<std::string_view(int)> varName,
  int resvar, int nargs, const int* argvars, std::string_view func) {
  PrintConsName(name);
  fmt::print("{} = {}({}", varName(resvar).data(), func.data(), varName(argvars[0]).data());
  for (int i = 1; i < nargs; i++)
    fmt::print(",{}", varName(argvars[i]).data());
  fmt::print(")\n");
}
void PrintExpConstraint(std::string_view name, std::function<std::string_view(int)> varName,
  int resvar, int npar, double param) {
  PrintConsName(name);
  std::string exponent = npar == 1 ? fmt::format("{}", param) : "e";
  fmt::print("{} = {}^{}\n", varName(resvar).data(), exponent);
}

void PrintPowConstraint(std::string_view name, const std::string& res, const std::string& base, const std::string& exponent) {
  PrintConsName(name);
  fmt::print("{} = {}^{}\n", res, base, exponent);
}


void CuoptlpModelAPI::SetLinearObjective( int iobj, const LinearObjective& lo ) {
  std::string name = lp()->AddObjective(lo.name(), Solver::OBJ_LIN);
  if (lp()->GetVerbosity() < 1)
    return;
  if (iobj<1) {
    PrintLinearObjective(name, GetVarName, lo.num_terms(), lo.vars().data(), lo.coefs().data(),
      lo.obj_sense() == mp::obj::MAX);
    /*
    CUOPTLP_CCALL(CUOPTLP_SetObjSense(lp(), 
                    obj::Type::MAX==lo.obj_sense() ? CUOPTLP_MAXIMIZE : CUOPTLP_MINIMIZE) );
    // This should set the objective exactly as given,
    // even when changing from a previous objective.
    CUOPTLP_CCALL(CUOPTLP_SetColObj(lp(), lo.num_terms(),
                           lo.vars().data(), lo.coefs().data()) ); */
  } else {
//    TODO If we support mutiple objectives, pass them to the solver
    fmt::print("Setting {}-th linear objective\n");
    PrintLinearObjective(name, GetVarName, lo.num_terms(), lo.vars().data(), lo.coefs().data(),
      lo.obj_sense() == mp::obj::MAX);
  }
}


void CuoptlpModelAPI::SetQuadraticObjective(int iobj, const QuadraticObjective& qo) {
  std::string name = lp()->AddObjective(qo.name(), Solver::OBJ_QUAD);
  if (lp()->GetVerbosity() < 1) 
    return;

  const auto &q = qo.GetQPTerms();
  const auto &l = qo.GetLinTerms();
  if (1 > iobj) {
    PrintQuadraticObjective(name, GetVarName, l.size(), l.pvars(), l.pcoefs(),
      qo.obj_sense() == mp::obj::MAX,
      q.size(), q.pvars1(), q.pvars2(), q.pcoefs());

    // Typical implementation
    //CUOPTLP_CCALL(CUOPTLP_SetQuadObj(lp(), qt.size(),
    //  (int*)qt.pvars1(), (int*)qt.pvars2(),
    //  (double*)qt.pcoefs()));
  }
  else {
    fmt::print("Setting {}-th objective\n");
    PrintQuadraticObjective(name, GetVarName, l.size(), l.pvars(), l.pcoefs(),
      qo.obj_sense() == mp::obj::MAX,
      q.size(), q.pvars1(), q.pvars2(), q.pcoefs());

  }
}



void CuoptlpModelAPI::AddConstraint(const LinConRange& lc) {
  std::string name = lp()->AddConstraintFlat(Solver::ConsType::CONS_LIN, lc.name());  
  if (lp()->GetVerbosity() < 1)
    return;
  PrintLinearConstraint(name, GetVarName,
    lc.size(), lc.pvars(), lc.pcoefs(), lc.lb(), lc.ub());
  
//  CUOPTLP_CCALL(CUOPTLP_AddRow(lp(), lc.size(), lc.pvars(), lc.pcoefs(), 
 //   NULL, lc.lb(), lc.ub(), lc.name()));
}
void CuoptlpModelAPI::AddConstraint(const LinConLE& lc) {
  std::string name = lp()->AddConstraintFlat(Solver::ConsType::CONS_LIN, lc.name());
  if (lp()->GetVerbosity() < 1)
    return;
  PrintLinearConstraint(name, GetVarName,
    lc.size(), lc.pvars(), lc.pcoefs(), MinusInfinity(), lc.rhs());

}
void CuoptlpModelAPI::AddConstraint(const LinConEQ& lc) {
  std::string name = lp()->AddConstraintFlat(Solver::ConsType::CONS_LIN, lc.name());
  if (lp()->GetVerbosity() < 1)
    return;
  PrintLinearConstraint(name, GetVarName,
    lc.size(), lc.pvars(), lc.pcoefs(), lc.rhs(), lc.rhs());
  
}
void CuoptlpModelAPI::AddConstraint(const LinConGE& lc) {
  std::string name = lp()->AddConstraintFlat(Solver::ConsType::CONS_LIN, lc.name());
  if (lp()->GetVerbosity() < 1)
    return;
  PrintLinearConstraint(name, GetVarName,
    lc.size(), lc.pvars(), lc.pcoefs(), lc.rhs(), Infinity());
  
}



void CuoptlpModelAPI::AddConstraint(const IndicatorConstraintLinLE &ic)  {
  std::string name = lp()->AddConstraintFlat(Solver::CONS_INDIC, ic.name());
  if (lp()->GetVerbosity() < 1)
    return;
  auto c = ic.get_constraint();
  PrintIndicator(name, GetVarName,
    ic.get_binary_var(), ic.get_binary_value(),
    c.size(), c.pvars(), c.pcoefs(), MinusInfinity(),
    c.rhs());
 // lp()->addEntity(Solver::CONS_INDIC);
  /*CUOPTLP_CCALL(CUOPTLP_AddIndicator(lp(),
    ic.get_binary_var(), ic.get_binary_value(),
    (int)ic.get_constraint().size(),
    ic.get_constraint().pvars(),
    ic.get_constraint().pcoefs(),
    CUOPTLP_LESS_EQUAL,
    ic.get_constraint().rhs()));*/
                               
}
void CuoptlpModelAPI::AddConstraint(const IndicatorConstraintLinEQ &ic)  {
  std::string name = lp()->AddConstraintFlat(Solver::CONS_INDIC, ic.name());
  if (lp()->GetVerbosity() < 1)
    return;
  auto c = ic.get_constraint();
  PrintIndicator(name, GetVarName,
    ic.get_binary_var(), ic.get_binary_value(),
    c.size(), c.pvars(), c.pcoefs(), c.rhs(),
    c.rhs());
}
void CuoptlpModelAPI::AddConstraint(const IndicatorConstraintLinGE &ic)  {
  std::string name = lp()->AddConstraintFlat(Solver::CONS_INDIC, ic.name());
  if (lp()->GetVerbosity() < 1)
    return;
  auto c = ic.get_constraint();
  PrintIndicator(name, GetVarName,
    ic.get_binary_var(), ic.get_binary_value(),
    c.size(), c.pvars(), c.pcoefs(), c.rhs(),
    Infinity());

}

void CuoptlpModelAPI::AddConstraint(const QuadConRange& qc) {
  std::string name = lp()->AddConstraintFlat(Solver::CONS_QUAD, qc.name());
  if (lp()->GetVerbosity() < 1)
    return;
  const auto& q = qc.GetQPTerms();
  const auto& l = qc.GetLinTerms();
  PrintQuadraticConstraint(name, GetVarName,
    l.size(), l.pvars(), l.pcoefs(), qc.lb(), qc.ub(),
    q.size(), q.pvars1(), q.pvars2(), q.pcoefs());
}

void CuoptlpModelAPI::AddConstraint( const QuadConLE& qc ) {
  std::string name = lp()->AddConstraintFlat(Solver::CONS_QUAD, qc.name());
  if (lp()->GetVerbosity() < 1)
    return;
  const auto& q = qc.GetQPTerms();
  const auto& l = qc.GetLinTerms();
  PrintQuadraticConstraint(name, GetVarName,
    l.size(), l.pvars(), l.pcoefs(), MinusInfinity(), qc.rhs(),
    q.size(), q.pvars1(), q.pvars2(), q.pcoefs());
}

void CuoptlpModelAPI::AddConstraint( const QuadConEQ& qc ) {
  std::string name = lp()->AddConstraintFlat(Solver::CONS_QUAD, qc.name());
  if (lp()->GetVerbosity() < 1)
    return;
  const auto& q = qc.GetQPTerms();
  const auto& l = qc.GetLinTerms();
  PrintQuadraticConstraint(name, GetVarName,
    l.size(), l.pvars(), l.pcoefs(), qc.rhs(), qc.rhs(),
    q.size(), q.pvars1(), q.pvars2(), q.pcoefs());
}

void CuoptlpModelAPI::AddConstraint( const QuadConGE& qc ) {
  std::string name = lp()->AddConstraintFlat(Solver::CONS_QUAD, qc.name());
  if (lp()->GetVerbosity() < 1)
    return;
  const auto& q = qc.GetQPTerms();
  const auto& l = qc.GetLinTerms();
  PrintQuadraticConstraint(name, GetVarName,
    l.size(), l.pvars(), l.pcoefs(), qc.rhs(), Infinity(),
    q.size(), q.pvars1(), q.pvars2(), q.pcoefs());
}

void CuoptlpModelAPI::AddConstraint( const QuadraticConeConstraint& qc ) {
  std::string name = lp()->AddConstraintFlat(Solver::CONS_QUAD_CONE);
  if (lp()->GetVerbosity() < 1)
    return;
  fmt::print("{}: ", name);
  PrintCoefficient(qc.GetParameters()[0], true);
  fmt::print(GetVarName(qc.GetArguments()[0]).data());
  fmt::print(">=sqrt(");

  for (int i = 1; i < qc.GetArguments().size(); i++) {
    PrintCoefficient(qc.GetParameters()[i], i == 0);
    fmt::print("{}^2", GetVarName(qc.GetArguments()[i]).data());
  }
  fmt::print(")\n");
  
}

void CuoptlpModelAPI::AddConstraint(
  const RotatedQuadraticConeConstraint& qc) {
  std::string name = lp()->AddConstraintFlat(Solver::CONS_QUAD_CONE_ROTATED);
  if (lp()->GetVerbosity() < 1)
    return;
  fmt::print("{}: ", name);
  PrintCoefficient(qc.GetParameters()[0]* qc.GetParameters()[1], true);
  fmt::print("{}*{} >=", GetVarName(qc.GetArguments()[0]).data(), GetVarName(qc.GetArguments()[1]).data());
  for (int i = 2; i < qc.GetArguments().size(); i++) {
    PrintCoefficient(qc.GetParameters()[i], i == 0);
    fmt::print("{}^2", GetVarName(qc.GetArguments()[i]).data());
  }
  fmt::print("\n");
}

void CuoptlpModelAPI::AddConstraint(
  const ExponentialConeConstraint& qc) {
  std::string name = lp()->AddConstraintFlat(Solver::CONS_QUAD_CONE_EXP);
  if (lp()->GetVerbosity() < 1)
    return;
  fmt::print("{}: ", name);
  PrintCoefficient(qc.GetParameters()[0] * qc.GetParameters()[1], true);
  fmt::print("{} >=", GetVarName(qc.GetArguments()[0]).data());
  PrintCoefficient(qc.GetParameters()[1], true);
  fmt::print(GetVarName(qc.GetArguments()[1]).data());
  fmt::print("*exp({}/{})\n", GetVarName(qc.GetArguments()[2]).data(), GetVarName(qc.GetArguments()[1]).data());
}

void CuoptlpModelAPI::AddConstraint(const SOS1Constraint& sos) {
  std::string name = lp()->AddConstraintFlat(Solver::CONS_SOS);
  if (lp()->GetVerbosity() < 1)
    return;
  fmt::print("{} ", name);
  fmt::print("({}):", "SOS1");
  PrintLinearBody(GetVarName, sos.size(), sos.get_vars().data(), sos.get_weights().data());
  fmt::print("\n");
}

void CuoptlpModelAPI::AddConstraint(const SOS2Constraint& sos) {
  std::string name = lp()->AddConstraintFlat(Solver::CONS_SOS);
  if (lp()->GetVerbosity() < 1)
    return;
  fmt::print("{} ", name);
  fmt::print("({}):", "SOS2");
  PrintLinearBody(GetVarName, sos.size(), sos.get_vars().data(), sos.get_weights().data());
  fmt::print("\n");
}

void CuoptlpModelAPI::AddConstraint(const MaxConstraint& mc) {
  std::string name = lp()->AddConstraintFlat(Solver::ConsType::CONS_MAX, mc.name());
  if (lp()->GetVerbosity() < 1)
    return;
  PrintFunctionalConstraintNoParam(name, GetVarName, mc.GetResultVar(),
    mc.GetArguments().size(), mc.GetArguments().data(), "max");
  /*
  const auto& args = mc.GetArguments();
  CUOPTLP_CCALL(CUOPTLP_addgenconstrMax(model(), mc.name(),
    mc.GetResultVar(),
    (int)args.size(), args.data(),
    MinusInfinity()));
    */
}

void CuoptlpModelAPI::AddConstraint(const MinConstraint& mc) {
  std::string name = lp()->AddConstraintFlat(Solver::ConsType::CONS_MAX, mc.name());
  if (lp()->GetVerbosity() < 1)
    return;
  PrintFunctionalConstraintNoParam(name, GetVarName, mc.GetResultVar(),
    mc.GetArguments().size(), mc.GetArguments().data(), "min");
}

void CuoptlpModelAPI::AddConstraint(const AbsConstraint& absc) {
  assert(absc.GetArguments().size() == 0);
  std::string name = lp()->AddConstraintFlat(Solver::ConsType::CONS_ABS, absc.name());
  if (lp()->GetVerbosity() < 1)
    return;
  PrintFunctionalConstraintNoParam(name, GetVarName, absc.GetResultVar(),
    absc.GetArguments().size(), absc.GetArguments().data(), "abs");
  /*
  const auto& args = absc.GetArguments();
  CUOPTLP_CCALL(CUOPTLP_addgenconstrAbs(model(), absc.name(),
    absc.GetResultVar(), args[0]));*/
}

void CuoptlpModelAPI::AddConstraint(const AndConstraint& cc) {
  std::string name = lp()->AddConstraintFlat(Solver::ConsType::CONS_AND, cc.name());
  if (lp()->GetVerbosity() < 1)
    return;
  PrintFunctionalConstraintNoParam(name, GetVarName, cc.GetResultVar(),
    cc.GetArguments().size(), cc.GetArguments().data(), "and");
}

void CuoptlpModelAPI::AddConstraint(const OrConstraint& dc) {
  std::string name = lp()->AddConstraintFlat(Solver::ConsType::CONS_OR, dc.name());
  if (lp()->GetVerbosity() < 1)
    return;
  PrintFunctionalConstraintNoParam(name, GetVarName, dc.GetResultVar(),
    dc.GetArguments().size(), dc.GetArguments().data(), "or");
}

void CuoptlpModelAPI::AddConstraint(const ExpConstraint& cc) {
  std::string name = lp()->AddConstraintFlat(Solver::ConsType::CONS_EXP, cc.name());
  if (lp()->GetVerbosity() < 1)
    return;
  PrintExpConstraint(name, GetVarName, cc.GetResultVar(), 0, 0);
}

void CuoptlpModelAPI::AddConstraint(const ExpAConstraint& cc) {
  std::string name = lp()->AddConstraintFlat(Solver::ConsType::CONS_EXP, cc.name());
  if (lp()->GetVerbosity() < 1)
    return;
  PrintExpConstraint(name, GetVarName, cc.GetResultVar(), 1, cc.GetParameters()[0]);
}

void CuoptlpModelAPI::AddConstraint(const LogConstraint& cc) {
  std::string name = lp()->AddConstraintFlat(Solver::ConsType::CONS_LOG, cc.name());
  if (lp()->GetVerbosity() < 1)
    return;
  PrintFunctionalConstraintNoParam(name, GetVarName, cc.GetResultVar(),
    cc.GetArguments().size(), cc.GetArguments().data(), "log");
};


void CuoptlpModelAPI::AddConstraint(const LogAConstraint& cc) {
  std::string name = lp()->AddConstraintFlat(Solver::ConsType::CONS_LOGA, cc.name());
  if (lp()->GetVerbosity() < 1)
    return;
  PrintConsName(name);
  fmt::print("{} = loga_{}({})\n", GetVarName(cc.GetResultVar()).data(), cc.GetParameters()[0],
    GetVarName(cc.GetArguments()[0]).data());
}


void CuoptlpModelAPI::AddConstraint(const PowConstExpConstraint& cc) {
  std::string name = lp()->AddConstraintFlat(Solver::ConsType::CONS_POW, cc.name());
  if (lp()->GetVerbosity() < 1)
    return;
  std::string base = std::string(GetVarName(cc.GetArguments()[0]));
  std::string exp = fmt::format("{}", cc.GetParameters()[0]);
  PrintPowConstraint(name, std::string(GetVarName(cc.GetResultVar())),
    base, exp);
}

void CuoptlpModelAPI::AddConstraint(const PowConstraint& cc) {
  std::string name = lp()->AddConstraintFlat(Solver::ConsType::CONS_POW, cc.name());
  if (lp()->GetVerbosity() < 1)
    return;
  std::string base = std::string(GetVarName(cc.GetArguments()[0]));
  std::string exp = std::string(GetVarName(cc.GetArguments()[1]));
  PrintPowConstraint(name, std::string(GetVarName(cc.GetResultVar())),
    base, exp);
}

void CuoptlpModelAPI::AddConstraint(const SinConstraint& cc) {
  std::string name = lp()->AddConstraintFlat(Solver::ConsType::CONS_SIN, cc.name());
  if (lp()->GetVerbosity() < 1)
    return;
  PrintFunctionalConstraintNoParam(name, GetVarName, cc.GetResultVar(),
    cc.GetArguments().size(), cc.GetArguments().data(), "sin");
}

void CuoptlpModelAPI::AddConstraint(const CosConstraint& cc) {
  std::string name = lp()->AddConstraintFlat(Solver::ConsType::CONS_COS, cc.name());
  if (lp()->GetVerbosity() < 1)
    return;
  PrintFunctionalConstraintNoParam(name, GetVarName, cc.GetResultVar(),
    cc.GetArguments().size(), cc.GetArguments().data(), "cos");
}

void CuoptlpModelAPI::AddConstraint(const TanConstraint& cc) {
  std::string name = lp()->AddConstraintFlat(Solver::ConsType::CONS_TAN, cc.name());
  if (lp()->GetVerbosity() < 1)
    return;
  PrintFunctionalConstraintNoParam(name, GetVarName, cc.GetResultVar(),
    cc.GetArguments().size(), cc.GetArguments().data(), "tan");
}

void CuoptlpModelAPI::AddConstraint(const PLConstraint& plc) {
  std::string name = lp()->AddConstraintFlat(Solver::ConsType::CONS_PL, plc.name());
  if (lp()->GetVerbosity() < 1)
    return;
  fmt::print("{} = <", GetVarName(plc.GetResultVar()).data());
  const auto& plp = plc.GetParameters().GetPLPoints();
  // Print slopes point
  for (int i = 0; i < plp.x_.size(); i++)
    fmt::print("<{},{}>", plp.x_.data()[i], plp.y_.data()[i]);
  // Print variable
  fmt::print("> {}", GetVarName(plc.GetArguments()[0]).data());

 // lp()->addEntity(Solver::CONS_PL);
  /*
  const auto& plp = plc.GetParameters().GetPLPoints();
  CUOPTLP_CCALL(CUOPTLP_addgenconstrPWL(model(), plc.name(),
    plc.GetArguments()[0], plc.GetResultVar(),
    plp.x_.size(), plp.x_.data(), plp.y_.data()));*/
}

template <class MPExpr>
void CuoptlpModelAPI::AppendLinAndConstTerms(
  Expr& ff, const MPExpr& nla) {
  if (double ct = GetConstTerm(nla))
    ff.AddArgument(Expr::makeConstant(ct));

  for (int i = 0; i < GetLinSize(nla); ++i) {
    if (double coef = GetLinCoef(nla, i)) {
      if (1.0 != coef) {
        auto f1 = Expr(Expr::Opcode::MULTIPLY);
        f1.AddArgument(Expr::makeConstant(coef));
        f1.AddArgument(GetLinTerm(nla, i));
        ff.AddArgument(f1);
      }
      else {
        ff.AddArgument(GetLinTerm(nla, i));
      }
    }
  }
}
template <class MPExpr>
void CuoptlpModelAPI::AppendQuadTerms(
  Expr& ff, const MPExpr& nlq) {
  for (int i = 0; i < GetQuadSize(nlq); ++i) {
    if (double coef = GetQuadCoef(nlq, i)) {
      auto f1 = Expr(Expr::Opcode::MULTIPLY);
      f1.AddArgument(GetQuadTerm1(nlq, i));
      f1.AddArgument(GetQuadTerm2(nlq, i));
      if (1.0 != coef) {
        f1.AddArgument(Expr::makeConstant(coef));
      }
      ff.AddArgument(f1);
    }
  }
}

Solver::VExpr CuoptlpModelAPI::AddExpression(const NLAffineExpression& nae) {
  Expr tree;
  tree.opcode = Solver::VExpr::Opcode::ADD;
  AppendLinAndConstTerms(tree, nae);
  return tree;
}

Solver::VExpr CuoptlpModelAPI::AddExpression(const NLQuadExpression& nlq) {
  Expr tree;
  tree.opcode = Expr::Opcode::ADD;
  AppendLinAndConstTerms(tree, nlq);
  AppendQuadTerms(tree, nlq);
  return tree;

}

void CuoptlpModelAPI::SetNLObjective(int i, const NLObjective& nlo) {
  const auto& exp = GetExpression(nlo);
  lp()->AddNLObjective(nlo.name(), exp, nlo.obj_sense() == mp::obj::MAX);
}







void CuoptlpModelAPI::AddConstraint(const NLConstraint& nlc) {
  const auto& nlexp = GetExpression(nlc);
  double lhs = GetLower(nlc), rhs = GetUpper(nlc);
  // Extract the linear part (similar to AppendLinAndConstantTerms but 
  // for NLConstraint
  Expr linearexp = Expr(Expr::Opcode::ADD);
  for (int i = 0; i < GetLinSize(nlc); ++i)
  {
    // coeff*var
    Expr term = Expr(Expr::Opcode::MULTIPLY);
    term.AddArgument(Expr::makeConstant(GetLinCoef(nlc, i)));
    term.AddArgument(Expr::makeVariable(GetLinVar(nlc, i)));
    linearexp.AddArgument(term);
  }
  // sum the nonlinear part
  linearexp.AddArgument(nlexp);
  lp()->AddConstraintNL(linearexp, lhs, rhs, nlc.GetName());
}

void  CuoptlpModelAPI::AddConstraint(const NLLogical& nll) {
  std::string name = lp()->AddConstraintNLAssign(nll.GetResultVar(),
    Solver::NLLogical, GetExpression(nll), nll.name());
}


void  CuoptlpModelAPI::AddConstraint(const NLReifEquiv& nle) {
  std::string name = lp()->AddConstraintNLAssign(GetVariable(nle),
    Solver::NLReifEquiv, GetExpression(nle), nle.name());

}

void  CuoptlpModelAPI::AddConstraint(const NLReifImpl& nle) {
    std::string name = lp()->AddConstraintNLAssign(GetVariable(nle),
    Solver::NLReifImpl, GetExpression(nle), nle.name());
}

void  CuoptlpModelAPI::AddConstraint(const NLReifRimpl& nle) {
std::string name = lp()->AddConstraintNLAssign(GetVariable(nle),
    Solver::NLReifRimpl, GetExpression(nle), nle.name());
}


void CuoptlpModelAPI::AddConstraint(const NLAssignEQ& nle) {
  std::string name = lp()->AddConstraintNLAssign(GetVariable(nle),
    Solver::AssignEQ, GetExpression(nle), nle.name());
}

/// NLAssignLE: algebraic expression expicifier in positive context.
/// Meaning: var <= expr.
void CuoptlpModelAPI::AddConstraint(const NLAssignLE& nle) {
  std::string name = lp()->AddConstraintNLAssign(GetVariable(nle),
    Solver::AssignLE, GetExpression(nle), nle.name());
}

/// NLAssignGE: algebraic expression expicifier in negative context.
/// Meaning: var >= expr.
void CuoptlpModelAPI::AddConstraint(const NLAssignGE& nle) {

  std::string name = lp()->AddConstraintNLAssign(GetVariable(nle),
    Solver::AssignGE, GetExpression(nle), nle.name());
}

/// res =log(exp)
Solver::VExpr  CuoptlpModelAPI::AddExpression(const LogExpression& e) {
   return Expr(Expr::Opcode::LOG, GetArgExpression(e, 0));
}
/// res = log_a(exp) where a is the parameter
Solver::VExpr  CuoptlpModelAPI::AddExpression(const LogAExpression& e) {
  Expr loga=Expr(Expr::Opcode::LOGA, GetArgExpression(e, 0));
  loga.AddArgument(Expr::makeConstant(GetParameter(e, 0)));
  return loga;
}
/// res = e^expr
Solver::VExpr  CuoptlpModelAPI::AddExpression(const ExpExpression& expr) {
  
  return Expr(Expr::Opcode::EXP, GetArgExpression(expr, 0));
}
/// res = a^expr
Solver::VExpr  CuoptlpModelAPI::AddExpression(const ExpAExpression& e) {
  
  Expr res = Expr(Expr::Opcode::POW);
  res.AddArgument(Expr::makeConstant(GetParameter(e, 0)));
  res.AddArgument(GetArgExpression(e, 0));
  return res;
}
/// res = x^y
Solver::VExpr  CuoptlpModelAPI::AddExpression(const PowExpression& e) {
  
  auto expr = Expr(Expr::Opcode::POW, GetArgExpression(e, 0));
  expr.AddArgument(GetArgExpression(e, 1));
  return expr;
}
/// res = x^const
Solver::VExpr CuoptlpModelAPI::AddExpression(const PowConstExpExpression& e) {
  auto expr = Expr(Expr::Opcode::POW, GetArgExpression(e, 0));
  expr.AddArgument(Solver::VExpr::makeConstant(GetParameter(e, 0)));
  return expr;
}
/// r = sin(e)
Solver::VExpr CuoptlpModelAPI::AddExpression(const SinExpression& e) {
  return Expr(Expr::Opcode::SIN, GetArgExpression(e, 0));
}
// r = cos(e)
Solver::VExpr CuoptlpModelAPI::AddExpression(const CosExpression& e) {
  return Expr(Expr::Opcode::COS, GetArgExpression(e, 0));
}
/// r=e1/e2
Solver::VExpr CuoptlpModelAPI::AddExpression(const DivExpression& e) {
  
  auto res = Expr(Expr::Opcode::DIVIDE);
  res.AddArgument(GetArgExpression(e, 0));
  res.AddArgument(GetArgExpression(e, 1));
  return res;
}

Solver::VExpr CuoptlpModelAPI::AddExpression(const TanExpression& e) {
  return Expr(Expr::Opcode::TAN, GetArgExpression(e, 0));
}

/// r = asin(v)
Solver::VExpr CuoptlpModelAPI::AddExpression(const AsinExpression& e) {
  return Expr(Expr::Opcode::ASIN, GetArgExpression(e, 0));
}

/// r = acos(v)
Solver::VExpr CuoptlpModelAPI::AddExpression(const AcosExpression& e) {
  return Expr(Expr::Opcode::ACOS, GetArgExpression(e, 0));
}

/// r = atan(v)
Solver::VExpr CuoptlpModelAPI::AddExpression(const AtanExpression& e) {
  return Expr(Expr::Opcode::ATAN, GetArgExpression(e, 0));
}

/// r = sinh(v)
Solver::VExpr CuoptlpModelAPI::AddExpression(const SinhExpression& e) {
  return Expr(Expr::Opcode::SINH, GetArgExpression(e, 0));
}

/// r = cosh(v)
Solver::VExpr CuoptlpModelAPI::AddExpression(const CoshExpression& e) {
  return Expr(Expr::Opcode::COSH, GetArgExpression(e, 0));
}

/// r = tanh(v)
Solver::VExpr CuoptlpModelAPI::AddExpression(const TanhExpression& e) {
  return Expr(Expr::Opcode::TANH, GetArgExpression(e, 0));
}

/// r = asinh(v)
Solver::VExpr CuoptlpModelAPI::AddExpression(const AsinhExpression& e) {
  return Expr(Expr::Opcode::ASINH, GetArgExpression(e, 0));
}

/// r = acosh(v)
Solver::VExpr CuoptlpModelAPI::AddExpression(const AcoshExpression& e) {
  return Expr(Expr::Opcode::ACOSH, GetArgExpression(e, 0));
}

/// r = atanh(v)
Solver::VExpr CuoptlpModelAPI::AddExpression(const AtanhExpression& e) {
  return Expr(Expr::Opcode::ATANH, GetArgExpression(e, 0));
}

void CuoptlpModelAPI::FinishProblemModificationPhase() {
}


} // namespace mp
