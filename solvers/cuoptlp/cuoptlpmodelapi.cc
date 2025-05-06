#include "cuoptlpmodelapi.h"


namespace mp {

void CuoptlpModelAPI::InitProblemModificationPhase(
    const FlatModelInfo* flat_model_info) {
  //auto varname = std::bind(&Solver::SolverModel::var_name, this->lp(),
  //  std::placeholders::_1);
  //GetVarName = varname;
  // Allocate storage if needed:
  // auto n_linear_cons =
  //   flat_model_info->GetNumberOfConstraintsOfGroup(CG_LINEAR);
  // reallocate_linear_cons( n_linear_cons );
  fmt::print("Initializing problem modification phase\n");
  //fmt::print("Number of variables: {}\n", flat_model_info->GetNumberOfVariables());
  fmt::print("Number of constraints: {}\n", flat_model_info->GetNumberOfConstraintsOfGroup(mp::ConstraintGroup::CG_Linear));

  if (lp_ == nullptr) {
    return;
  }

  lp_->num_constraints = 0;
  lp_->num_variables = 0;
  lp_->constraint_matrix_row_offsets.clear();
  lp_->constraint_matrix_coefficients.clear();
  lp_->constraint_matrix_column_indices.clear();
  lp_->constraint_sense.clear();
  lp_->rhs.clear();
  lp_->lower_bounds.clear();
  lp_->upper_bounds.clear();
  lp_->variable_types.clear();
  lp_->objective_coefficients.clear();
  lp_->objective_sense = CUOPT_MINIMIZE;
  lp_->objective_offset = 0.0;
}

void CuoptlpModelAPI::AddVariables(const VarArrayDef& v) {
  // TODO Add variables using solver API; typically,
  // first convert the MP variable type to the appropriate solver-defined type,
  // then add them
  fmt::print("Adding {} variables\n", v.size());
  lp_->lower_bounds.resize(v.size());
  lp_->upper_bounds.resize(v.size());
  lp_->variable_types.resize(v.size());
  lp_->num_variables = v.size();
  for (auto i = 0; i<v.size();  i++) {
    lp_->lower_bounds[i] = v.plb()[i];
    lp_->upper_bounds[i] = v.pub()[i];
    lp_->variable_types[i] = v.ptype()[i] == var::Type::CONTINUOUS ? CUOPT_CONTINUOUS : CUOPT_INTEGER;
    fmt::print("Adding variable of type {} lower bound {} upper bound {}\n", v.ptype()[i], v.plb()[i], v.pub()[i]);
  }


  if (v.pnames() != nullptr) {
    for (auto i = 0; i<v.size();  i++) {
      fmt::print("Adding variable of name {}\n", v.pnames()[i]);
    }
  } else {
    fmt::print("No variable names provided\n");
  }

  // Preallocate in solver mem
  //lp()->allocateVars(v.size());
  // Assign one by one
  //for (auto i = 0; i<v.size();  i++)
  //  lp()->AddVariable(v.ptype()[i],
  //    v.plb()[i], v.pub()[i]),
  //    v.pnames() == NULL ? nullptr : v.pnames()[i];

  //fmt::print("Added {} continuous, {} integer and {} binary variables.\n",
  //  lp()->getNumVars(Solver::VarType::CONTINUOUS),
  //  lp()->getNumVars(Solver::VarType::INTEGER),
  //  lp()->getNumVars(Solver::VarType::BINARY));

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


void CuoptlpModelAPI::SetLinearObjective( int iobj, const LinearObjective& lo ) {
  fmt::print("Setting linear objective\n");

  lp_->objective_sense = lo.obj_sense() == mp::obj::MAX ? CUOPT_MAXIMIZE : CUOPT_MINIMIZE;
  lp_->objective_coefficients.resize(lp_->variable_types.size());

  for (auto k = 0; k < lo.num_terms(); k++) {
    lp_->objective_coefficients[lo.vars()[k]] = lo.coefs()[k];
    fmt::print("Setting objective coefficient {} for variable {}\n", lo.coefs()[k], lo.vars()[k]);
  }


  //std::string name = lp()->AddObjective(lo.name(), Solver::OBJ_LIN);
  //if (lp()->GetVerbosity() < 1)
  //  return;
  //if (iobj<1) {
  //  PrintLinearObjective(name, GetVarName, lo.num_terms(), lo.vars().data(), lo.coefs().data(),
  //    lo.obj_sense() == mp::obj::MAX);
    /*
    CUOPTLP_CCALL(CUOPTLP_SetObjSense(lp(),
                    obj::Type::MAX==lo.obj_sense() ? CUOPTLP_MAXIMIZE : CUOPTLP_MINIMIZE) );
    // This should set the objective exactly as given,
    // even when changing from a previous objective.
    CUOPTLP_CCALL(CUOPTLP_SetColObj(lp(), lo.num_terms(),
                           lo.vars().data(), lo.coefs().data()) ); */
  //} else {
//    TODO If we support mutiple objectives, pass them to the solver
    //fmt::print("Setting {}-th linear objective\n");
    //PrintLinearObjective(name, GetVarName, lo.num_terms(), lo.vars().data(), lo.coefs().data(),
    //  lo.obj_sense() == mp::obj::MAX);
  //}
}


#if 0
void CuoptlpModelAPI::AddConstraint(const LinConRange& lc) {
  fmt::print("Adding linear constraint with range\n");
  //std::string name = lp()->AddConstraintFlat(Solver::ConsType::CONS_LIN, lc.name());
  //if (lp()->GetVerbosity() < 1)
  //  return;
  //PrintLinearConstraint(lc.name(), GetVarName,
  //  lc.size(), lc.pvars(), lc.pcoefs(), lc.lb(), lc.ub());

  if (constraint_matrix_row_offsets_.size() == 0) {
    constraint_matrix_row_offsets_.push_back(0);
  }

  for (auto i = 0; i < lc.size(); i++) {
    fmt::print("Adding constraint coefficient {} for variable {}\n", lc.pcoefs()[i], lc.pvars()[i]);
    constraint_matrix_coefficients_.push_back(lc.pcoefs()[i]);
    constraint_matrix_column_indices_.push_back(lc.pvars()[i]);
    nnz_++;
  }
  constraint_matrix_row_offsets_.push_back(nnz_);
  constraint_lower_bounds_.push_back(lc.lb());
  constraint_upper_bounds_.push_back(lc.ub());
  num_constraints_++;
}
#endif


void CuoptlpModelAPI::cuOptAddConstraint(size_t num_coefficients, const double* coefficients, const int* variables, char sense, double rhs) {
  fmt::print("Adding constraint with {} coefficients\n", num_coefficients);
  fmt::print("Adding constraint with sense {}\n", sense);
  fmt::print("Adding constraint with rhs {}\n", rhs);

  if (lp_->constraint_matrix_row_offsets.size() == 0) {
    lp_->constraint_matrix_row_offsets.push_back(0);
  }

  for (auto k = 0; k < num_coefficients; k++) {
    lp_->constraint_matrix_coefficients.push_back(coefficients[k]);
    lp_->constraint_matrix_column_indices.push_back(variables[k]);
    lp_->nnz++;
  }
  lp_->constraint_matrix_row_offsets.push_back(lp_->nnz);
  lp_->constraint_sense.push_back(sense);
  lp_->rhs.push_back(rhs);
  lp_->num_constraints++;
}

void CuoptlpModelAPI::AddConstraint(const LinConLE& lc) {
  fmt::print("Adding linear constraint with less than or equal to\n");
  cuOptAddConstraint(lc.size(), lc.pcoefs(), lc.pvars(), CUOPT_LESS_THAN, lc.rhs());
  //std::string name = lp()->AddConstraintFlat(Solver::ConsType::CONS_LIN, lc.name());
  //if (lp()->GetVerbosity() < 1)
  //  return;
  //PrintLinearConstraint(name, GetVarName,
  //  lc.size(), lc.pvars(), lc.pcoefs(), MinusInfinity(), lc.rhs());

}
void CuoptlpModelAPI::AddConstraint(const LinConEQ& lc) {
  fmt::print("Adding linear constraint with equal to\n");
  cuOptAddConstraint(lc.size(), lc.pcoefs(), lc.pvars(), CUOPT_EQUAL, lc.rhs());
  //std::string name = lp()->AddConstraintFlat(Solver::ConsType::CONS_LIN, lc.name());
  //if (lp()->GetVerbosity() < 1)
  //  return;
  //PrintLinearConstraint(name, GetVarName,
  //  lc.size(), lc.pvars(), lc.pcoefs(), lc.rhs(), lc.rhs());

}
void CuoptlpModelAPI::AddConstraint(const LinConGE& lc) {
  fmt::print("Adding linear constraint with greater than or equal to\n");
  cuOptAddConstraint(lc.size(), lc.pcoefs(), lc.pvars(), CUOPT_GREATER_THAN, lc.rhs());
  //std::string name = lp()->AddConstraintFlat(Solver::ConsType::CONS_LIN, lc.name());
  //if (lp()->GetVerbosity() < 1)
  //  return;
  //PrintLinearConstraint(name, GetVarName,
  //  lc.size(), lc.pvars(), lc.pcoefs(), lc.rhs(), Infinity());

}


void CuoptlpModelAPI::FinishProblemModificationPhase() {
  fmt::print("Finishing problem modification phase\n");
  fmt::print("Number of constraints: {}\n", lp_->num_constraints);
  fmt::print("Number of variables: {}\n", lp_->num_variables);
  fmt::print("Number of non-zeros: {}\n", lp_->nnz);


  cuopt_int_t status = cuOptCreateProblem(
    lp_->num_constraints,
    lp_->num_variables,
    lp_->objective_sense,
    lp_->objective_offset,
    lp_->objective_coefficients.data(),
    lp_->constraint_matrix_row_offsets.data(),
    lp_->constraint_matrix_column_indices.data(),
    lp_->constraint_matrix_coefficients.data(),
    lp_->constraint_sense.data(),
    lp_->rhs.data(),
    lp_->lower_bounds.data(),
    lp_->upper_bounds.data(),
    lp_->variable_types.data(),
    &lp_->problem
  );
  if (status != CUOPT_SUCCESS) {
    throw std::runtime_error(fmt::format("Error creating problem: {}", status));
  }
}


} // namespace mp
