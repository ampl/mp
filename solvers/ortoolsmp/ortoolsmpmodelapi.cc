#include "ortoolsmpmodelapi.h"


namespace mp {

void OrtoolsModelAPI::InitProblemModificationPhase(
    const FlatModelInfo*) { }

void OrtoolsModelAPI::AddVariables(const VarArrayDef& v) {
  for (int i = 0; i < v.size(); ++i) {
    lp()->MakeVar(v.plb()[i], v.pub()[i], 
      v.ptype()[i] == var::Type::INTEGER, "");
  }
}

void OrtoolsModelAPI::SetLinearObjective( int iobj, const LinearObjective& lo ) {
  if (iobj<1) {
    operations_research::MPObjective* const objective = lp()->MutableObjective();
    for (int j = 0; j < lo.num_terms(); ++j) {
      objective->SetCoefficient(lp()->variable(lo.vars()[j]), lo.coefs()[j]);
    }
    if (lo.obj_sense() == obj::Type::MAX)
      objective->SetMaximization();
    else
      objective->SetMinimization();
  } else {
    throw std::runtime_error("Multiple objectives not supported");
  }
}


void OrtoolsModelAPI::SetQuadraticObjective(int iobj, const QuadraticObjective& qo) {
  throw std::runtime_error("Quadratic objective not supported");
}

void OrtoolsModelAPI::AddConstraint(const LinConRange& lc) {
  auto c = lp()->MakeRowConstraint(lc.lb(), lc.ub());
  for (size_t i = 0; i < lc.size(); i++)
    c->SetCoefficient(lp()->variable(lc.var(i)), lc.coef(i));
}

void OrtoolsModelAPI::AddConstraint(const LinConLE& lc) {
  auto c = lp()->MakeRowConstraint(lc.lb(), lc.ub());
  for (size_t i = 0; i < lc.size(); i++)
    c->SetCoefficient(lp()->variable(lc.var(i)), lc.coef(i));
}
void OrtoolsModelAPI::AddConstraint(const LinConEQ& lc) {
  auto c = lp()->MakeRowConstraint(lc.lb(), lc.ub());
  for (size_t i = 0; i < lc.size(); i++)
    c->SetCoefficient(lp()->variable(lc.var(i)), lc.coef(i));
}
void OrtoolsModelAPI::AddConstraint(const LinConGE& lc) {
  auto c = lp()->MakeRowConstraint(lc.lb(), lc.ub());
  for (size_t i = 0; i < lc.size(); i++)
    c->SetCoefficient(lp()->variable(lc.var(i)), lc.coef(i));
}



void OrtoolsModelAPI::FinishProblemModificationPhase() {
}


} // namespace mp
