#include "cuoptmodelapi.h"


namespace mp {

void CuoptModelAPI::InitProblemModificationPhase(const FlatModelInfo*) {
  json* prob = get_json_prob();

  (*prob)["variable_types"] = json::array();

  (*prob)["variable_bounds"]["lower_bounds"] = json::array();
  (*prob)["variable_bounds"]["upper_bounds"] = json::array();

  (*prob)["constraint_bounds"]["lower_bounds"] = json::array();
  (*prob)["constraint_bounds"]["upper_bounds"] = json::array();

  (*prob)["csr_constraint_matrix"]["offsets"] = json::array();
  (*prob)["csr_constraint_matrix"]["offsets"].push_back(0);
}

void CuoptModelAPI::AddVariables(const VarArrayDef& v) {
  json* prob = get_json_prob();

  nvars_ = v.size();

  for (int i = 0; i < v.size(); i++) {
    if (std::isinf(v.plb()[i])) {
      (*prob)["variable_bounds"]["lower_bounds"].push_back("ninf");
    } else {
      (*prob)["variable_bounds"]["lower_bounds"].push_back(v.plb()[i]);
    }

    if (std::isinf(v.pub()[i])) {
      (*prob)["variable_bounds"]["upper_bounds"].push_back("inf");
    } else {
      (*prob)["variable_bounds"]["upper_bounds"].push_back(v.pub()[i]);
    }

    if (v.ptype()[i] == var::Type::INTEGER) {
      (*prob)["variable_types"].push_back("I");
      SetIsMIP(true);
    } else {
      (*prob)["variable_types"].push_back("C");
    }
  }
}

void CuoptModelAPI::SetLinearObjective( int iobj, const LinearObjective& lo ) {
  json* prob = get_json_prob();

  if (lo.obj_sense() == mp::obj::MAX) {
    (*prob)["maximize"] = true;
  } else {
    (*prob)["maximize"] = false;
  }

  (*prob)["objective_data"]["coefficients"] = std::vector<double>(nvars_, 0.0);
  for (std::size_t i = 0; i < lo.vars().size(); i++) {
    (*prob)["objective_data"]["coefficients"][lo.vars()[i]]= lo.coefs()[i];
  }
}


void CuoptModelAPI::AddConstraint(const LinConRange& lc) {
  json* prob = get_json_prob();

  if (std::isinf(lc.lb())) {
    (*prob)["constraint_bounds"]["lower_bounds"].push_back("ninf");
  } else {
    (*prob)["constraint_bounds"]["lower_bounds"].push_back(lc.lb());
  }
  if (std::isinf(lc.ub())) {
    (*prob)["constraint_bounds"]["upper_bounds"].push_back("inf");
  } else {
    (*prob)["constraint_bounds"]["upper_bounds"].push_back(lc.ub());
  }

  for (std::size_t i = 0; i < lc.size(); i++) {
    (*prob)["csr_constraint_matrix"]["values"].push_back(lc.pcoefs()[i]);
    (*prob)["csr_constraint_matrix"]["indices"].push_back(lc.pvars()[i]);
  }
  (*prob)["csr_constraint_matrix"]["offsets"].push_back(
    (*prob)["csr_constraint_matrix"]["values"].size());
}

void CuoptModelAPI::AddConstraint(const LinConLE& lc) {
  json* prob = get_json_prob();

  (*prob)["constraint_bounds"]["lower_bounds"].push_back("ninf");
  (*prob)["constraint_bounds"]["upper_bounds"].push_back(lc.rhs());

  for (std::size_t i = 0; i < lc.size(); i++) {
    (*prob)["csr_constraint_matrix"]["values"].push_back(lc.pcoefs()[i]);
    (*prob)["csr_constraint_matrix"]["indices"].push_back(lc.pvars()[i]);
  }
  (*prob)["csr_constraint_matrix"]["offsets"].push_back(
    (*prob)["csr_constraint_matrix"]["values"].size());
}

void CuoptModelAPI::AddConstraint(const LinConEQ& lc) {
  json* prob = get_json_prob();

  (*prob)["constraint_bounds"]["lower_bounds"].push_back(lc.rhs());
  (*prob)["constraint_bounds"]["upper_bounds"].push_back(lc.rhs());

  for (std::size_t i = 0; i < lc.size(); i++) {
    (*prob)["csr_constraint_matrix"]["values"].push_back(lc.pcoefs()[i]);
    (*prob)["csr_constraint_matrix"]["indices"].push_back(lc.pvars()[i]);
  }
  (*prob)["csr_constraint_matrix"]["offsets"].push_back(
    (*prob)["csr_constraint_matrix"]["values"].size());
}

void CuoptModelAPI::AddConstraint(const LinConGE& lc) {
  json* prob = get_json_prob();

  (*prob)["constraint_bounds"]["lower_bounds"].push_back(lc.rhs());
  (*prob)["constraint_bounds"]["upper_bounds"].push_back("inf");

  for (std::size_t i = 0; i < lc.size(); i++) {
    (*prob)["csr_constraint_matrix"]["values"].push_back(lc.pcoefs()[i]);
    (*prob)["csr_constraint_matrix"]["indices"].push_back(lc.pvars()[i]);
  }
  (*prob)["csr_constraint_matrix"]["offsets"].push_back(
    (*prob)["csr_constraint_matrix"]["values"].size());
}

void CuoptModelAPI::SetQuadraticObjective(int iobj, const QuadraticObjective& qo) {
  throw std::runtime_error("Quadratic objective not supported");
}

void CuoptModelAPI::FinishProblemModificationPhase() { 
}

} // namespace mp
