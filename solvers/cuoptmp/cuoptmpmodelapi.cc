#include "cuoptmpmodelapi.h"


namespace mp {

void CuoptmpModelAPI::InitProblemModificationPhase(
    const FlatModelInfo* flat_model_info) {
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

void CuoptmpModelAPI::AddVariables(const VarArrayDef& v) {
  lp_->lower_bounds.resize(v.size());
  lp_->upper_bounds.resize(v.size());
  lp_->variable_types.resize(v.size());
  lp_->num_variables = v.size();
  for (auto i = 0; i<v.size();  i++) {
    lp_->lower_bounds[i] = v.plb()[i];
    lp_->upper_bounds[i] = v.pub()[i];
    lp_->variable_types[i] = v.ptype()[i] == var::Type::CONTINUOUS ? CUOPT_CONTINUOUS : CUOPT_INTEGER;
  }
}


void CuoptmpModelAPI::SetLinearObjective( int iobj, const LinearObjective& lo ) {
  lp_->objective_sense = lo.obj_sense() == mp::obj::MAX ? CUOPT_MAXIMIZE : CUOPT_MINIMIZE;
  lp_->objective_coefficients.clear();
  lp_->objective_coefficients.resize(lp_->variable_types.size());

  for (auto k = 0; k < lo.num_terms(); k++) {
    lp_->objective_coefficients[lo.vars()[k]] = lo.coefs()[k];
  }
}

void CuoptmpModelAPI::cuOptAddConstraint(size_t num_coefficients, const double* coefficients, const int* variables, char sense, double rhs) {
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

void CuoptmpModelAPI::AddConstraint(const LinConLE& lc) {
  cuOptAddConstraint(lc.size(), lc.pcoefs(), lc.pvars(), CUOPT_LESS_THAN, lc.rhs());
}
void CuoptmpModelAPI::AddConstraint(const LinConEQ& lc) {
  cuOptAddConstraint(lc.size(), lc.pcoefs(), lc.pvars(), CUOPT_EQUAL, lc.rhs());
}

void CuoptmpModelAPI::AddConstraint(const LinConGE& lc) {
  cuOptAddConstraint(lc.size(), lc.pcoefs(), lc.pvars(), CUOPT_GREATER_THAN, lc.rhs());
}


void CuoptmpModelAPI::FinishProblemModificationPhase() {
  if (lp_->num_constraints == 0) {
    // Add a dummy constraint to make the problem non-empty
    // 0.0 * x[0] == 0.0
    lp_->num_constraints = 1;
    lp_->constraint_matrix_row_offsets.push_back(0);
    lp_->constraint_matrix_row_offsets.push_back(1);
    lp_->constraint_matrix_coefficients.push_back(0.0);
    lp_->constraint_matrix_column_indices.push_back(0);
    lp_->constraint_sense.push_back(CUOPT_EQUAL);
    lp_->rhs.push_back(0.0);
    lp_->nnz = 1;
  }

  if (lp_->objective_coefficients.size() == 0)
  {
    lp_->objective_coefficients.resize(lp_->num_variables, 0.0);
  }

  for (int j = 0; j < lp_->num_variables; j++) {
    if (lp_->lower_bounds[j] < -1e20) {
      lp_->lower_bounds[j] = -INFINITY;
    }
    if (lp_->upper_bounds[j] > 1e20) {
      lp_->upper_bounds[j] = INFINITY;
    }
  }


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
