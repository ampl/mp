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
  accObjectives().setNumVars(v.size());
}

void HighsModelAPI::SetLinearObjective( int iobj, const LinearObjective& lo ) {
  // Note:
  // In case of native multiobjectives, HiGHS wants priorities and other meta info
  // passed when adding the objectives, that is available in Backend::InputExtras.
  // So we accumulate the objectives here and set them in HiGHS in Backend::InputExtras.
  accObjectives().add(lo.vars(), lo.coefs(), lo.obj_sense()==obj::Type::MAX);
}


void HighsModelAPI::SetQuadraticObjective(int iobj, const QuadraticObjective& qo) {
  if (1 > iobj) {
    auto lo= qo.GetLinTerms();
    // Note that the linear part is only stored (see SetLinearObjective)
    accObjectives().add(lo.vars(), lo.coefs(), qo.obj_sense() == obj::Type::MAX);
    

    const auto& qt = qo.GetQPTerms();
    std::vector<int> startCols(NumVars());
    std::vector<double> coeffs(qt.size());
    // Convert to Highs Hessian upper triangular format.
    size_t q=0;              // the index in qt
    for (size_t j=0; j<startCols.size(); ++j) {
      assert(q>=qt.size() || j<=(size_t)qt.var1(q)); // qt sorted
      if (q<qt.size() && j==(size_t)qt.var1(q)) {
        startCols[j] = q;
        for ( ; q<qt.size() && (size_t)qt.var1(q) == j; ++q) {
          assert(j <= (size_t)qt.var2(q));      // upper triangular
          coeffs[q] =
                (j == (size_t)qt.var2(q) ? 2.0 : 1.0) * qt.coef(q);
        }
      } else {
        startCols[j] = q;
      }
    }
    HIGHS_CCALL(Highs_passHessian(lp(), NumVars(), qt.size(),
                                  kHighsHessianFormatTriangular,
      startCols.data(), qt.pvars2(), coeffs.data()));
  }
  else {
    throw std::runtime_error("Multiple quadratic objectives not supported natively, try using multi-objective\nemulator by setting option multiobj=2");
  }
}

void HighsModelAPI::AddConstraint(const LinConRange& lc) {
  acc_constraints_.add(lc);
}
void HighsModelAPI::AddConstraint(const LinConLE& lc) {
  acc_constraints_.add(lc);
}
void HighsModelAPI::AddConstraint(const LinConEQ& lc) {
  acc_constraints_.add(lc);
}
void HighsModelAPI::AddConstraint(const LinConGE& lc) {
  acc_constraints_.add(lc);
}

void HighsModelAPI::FinishProblemModificationPhase() {
  HIGHS_CCALL(Highs_addRows(lp(),
    acc_constraints_.lb.size(),
    acc_constraints_.lb.data(),
    acc_constraints_.ub.data(),
    acc_constraints_.coeffs.size(),
    acc_constraints_.starts.data(),
    acc_constraints_.indices.data(),
    acc_constraints_.coeffs.data()));
    // reinitialize accumulator for model modification 
    acc_constraints_ = AccConstraints();     

  // If in multiobjective simulator, set the objective each time
  if (accObjectives().hadEmulatedMultiObj()) {
    accObjectives().setAllInHighs(lp());
    accObjectives().clear();
  }
}

} // namespace mp
