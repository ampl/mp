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
  if (v.pnames())
    for (int i = 0; i < v.size(); i++)
      HIGHS_CCALL(Highs_passColName(lp(), i, v.pnames()[i]));
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
    std::vector<int> index;
    std::vector<double> coeffs;
    index.reserve(NumVars());
    coeffs.reserve(NumVars());
    // Convert to Highs Hessian upper triangular format.
    // As of Highs 1.5.3, we need a sparse element for every
    // column, so we fill 0's where needed.
    size_t q=0;              // the index in qt
    for (size_t j=0; j<startCols.size(); ++j) {
      assert(q>=qt.size() || j<=(size_t)qt.var1(q)); // qt sorted
      if (q<qt.size() && j==(size_t)qt.var1(q)) {
        startCols[j] = index.size();
        for ( ; q<qt.size() && (size_t)qt.var1(q) == j; ++q) {
          assert(j <= (size_t)qt.var2(q));      // upper triangular
          index.push_back(qt.var2(q));
          coeffs.push_back(
                (j == (size_t)qt.var2(q) ? 2.0 : 1.0) * qt.coef(q));
        }
      } else {
        startCols[j] = index.size();
        index.push_back(j);
        coeffs.push_back(0.0);        // empty diagonal element
      }
    }
    HIGHS_CCALL(Highs_passHessian(lp(), NumVars(), index.size(),
                                  kHighsHessianFormatTriangular,
      startCols.data(), index.data(), coeffs.data()));
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
