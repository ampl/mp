#include <cmath>

#include "highsmpmodelapi.h"


namespace mp {

void HighsModelAPI::InitProblemModificationPhase(const FlatModelInfo*) { }

void HighsModelAPI::AddVariables(const VarArrayDef& v) {
  std::vector<int> intIndices;
  std::vector<double> costs(v.size(), 0);
  std::vector<double> lbs(v.size());
  std::vector<double> ubs(v.size());
  for (int i = 0; i < v.size(); i++) {
    if (var::Type::INTEGER == v.ptype()[i]) intIndices.push_back(i);
    lbs[i] = std::isinf(v.plb()[i]) ? MinusInfinity() : v.plb()[i];
    ubs[i] = std::isinf(v.pub()[i]) ?  Infinity() : v.pub()[i];

  }
  HIGHS_CCALL(Highs_addCols(lp(), v.size(), costs.data(), lbs.data(), ubs.data(), 0, NULL, NULL, NULL));
  if (intIndices.size() > 0) {
    std::vector<int> types(intIndices.size(), 1); // TODO get the 1 from solver API?
    HIGHS_CCALL(Highs_changeColsIntegralityBySet(lp(), intIndices.size(),
      intIndices.data(), types.data()));
  }
}

void HighsModelAPI::SetLinearObjective( int iobj, const LinearObjective& lo ) {
  if (iobj < 1) {
  // Highs_changeObjectiveOffset(highs, offset); TODO offset?
    HIGHS_CCALL(Highs_changeColsCostBySet(lp(), lo.num_terms(), 
      lo.vars().data(), lo.coefs().data()));
    HIGHS_CCALL(Highs_changeObjectiveSense(lp(), 
                    obj::Type::MAX==lo.obj_sense() ? 
      kHighsObjSenseMaximize : kHighsObjSenseMinimize) );
  } else {
      throw std::runtime_error("HighS does not support multiple objectives.");
  }
}


void HighsModelAPI::SetQuadraticObjective(int iobj, const QuadraticObjective& qo) {
  if (1 > iobj) {
    SetLinearObjective(iobj, qo);
    const auto& qt = qo.GetQPTerms();
    std::vector<int> startCols(NumVars());
    std::vector<double> coeffs(qt.size());
    int currentCol = 0, newCol = 0;
    startCols[0] = 0;
    // Convert to Highs Hessian upper triangular format
    for (size_t i = 0; i < qt.size(); i++)
    {
      newCol = qt.vars1()[i];
      coeffs[i] = (qt.vars2()[i] == newCol) ? qt.coefs()[i]*2 : qt.coefs()[i];
      if (newCol == currentCol)
        continue;
      else if (newCol < currentCol)
        throw std::runtime_error("Check order of quadratic hessian");
      else // vars1 > currentCol
      {
        startCols[newCol] = i;
        currentCol = newCol;
      }
    }
    HIGHS_CCALL(Highs_passHessian(lp(), NumVars(), qt.size(), kHighsHessianFormatTriangular,
      startCols.data(), qt.pvars2(), coeffs.data()));
  }
  else {
    throw std::runtime_error("Multiple quadratic objectives not supported");
  }
}

void HighsModelAPI::AddConstraint(const LinConRange& lc) {
  AccConstraints.add(lc);
}
void HighsModelAPI::AddConstraint(const LinConLE& lc) {
  AccConstraints.add(lc);
}
void HighsModelAPI::AddConstraint(const LinConEQ& lc) {
  AccConstraints.add(lc);
}
void HighsModelAPI::AddConstraint(const LinConGE& lc) {
  AccConstraints.add(lc);
}

void HighsModelAPI::FinishProblemModificationPhase() {
  HIGHS_CCALL(Highs_addRows(lp(),
    AccConstraints.lb.size(),
    AccConstraints.lb.data(),
    AccConstraints.ub.data(),
    AccConstraints.coeffs.size(),
    AccConstraints.starts.data(),
    AccConstraints.indices.data(),
    AccConstraints.coeffs.data()));
}
} // namespace mp
