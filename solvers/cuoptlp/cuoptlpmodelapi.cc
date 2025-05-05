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


void CuoptlpModelAPI::FinishProblemModificationPhase() {
}


} // namespace mp
