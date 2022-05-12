#include <cmath>

#include "highsdirectmodelapi.h"

#include "mp/model-mgr-with-std-pb.h"
#include "mp/flat/problem_flattener.h"
#include "mp/flat/redef/MIP/converter_mip.h"
#include "mp/flat/model_api_connect.h"

namespace mp {

/// Defining the function in ...modelapi.cc
/// for recompilation speed
std::unique_ptr<BasicModelManager>
CreateHighsModelMgr(HighsCommon& cc, Env& e,
                     pre::BasicPresolver*& pPre) {
  return CreateModelMgrWithFlatConverter<
      HighsModelAPI, MIPFlatConverter >(cc, e, pPre);
}


void HighsModelAPI::InitProblemModificationPhase() { }

void HighsModelAPI::AddVariables(const VarArrayDef& v) {
  std::vector<int> intIndices;
  std::vector<double> costs(v.size(), 0);
  std::vector<double> lbs(v.size());
  std::vector<double> ubs(v.size());
  for (size_t i = 0; i < v.size(); i++) {
    if (var::Type::INTEGER == v.ptype()[i]) intIndices.push_back(i);
    lbs[i] = std::isinf(v.plb()[i]) ? MinusInfinity() : v.plb()[i];
    ubs[i] = std::isinf(v.pub()[i]) ?  Infinity() : v.pub()[i];

  }
  HIGHS_CCALL(Highs_addCols(lp(), v.size(), costs.data(), lbs.data(), ubs.data(), 0, NULL, NULL, NULL));
  std::vector<int> types(intIndices.size(), 1); // TODO get the 1 from solver API?
  HIGHS_CCALL(Highs_changeColsIntegralityBySet(lp(), intIndices.size(),
    intIndices.data(), types.data()));
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
//    TODO If we support mutiple objectives, pass them to the solver
    fmt::format("Setting {}-th linear objective: {} terms.\n", iobj, lo.num_terms());
  }
}


void HighsModelAPI::SetQuadraticObjective(int iobj, const QuadraticObjective& qo) {
  if (1 > iobj) {
    SetLinearObjective(iobj, qo);                         // add the linear part
    const auto& qt = qo.GetQPTerms();
    int nv = NumVars();
    std::vector<int> startCols(nv);
    int currentCol = 0;
    startCols[0] = 0;
    int newCol;

    // Convert to Highs Hessian format
    for (int i = 0; i < qt.size(); i++)
    {
      newCol = qt.vars1()[i];
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
      startCols.data(), qt.pvars2(), qt.pcoefs()));
  }
  else {
    throw std::runtime_error("Multiple quadratic objectives not supported");
  }
}

void HighsModelAPI::AddConstraint(const LinConRange& lc) {
  HIGHS_CCALL(Highs_addRow(lp(), lc.lb(), lc.ub(), lc.size(), lc.pvars(), lc.pcoefs()));
}
void HighsModelAPI::AddConstraint(const LinConLE& lc) {
  HIGHS_CCALL(Highs_addRow(lp(), lc.lb(), lc.ub(), lc.size(), lc.pvars(), lc.pcoefs()));
}
void HighsModelAPI::AddConstraint(const LinConEQ& lc) {
  HIGHS_CCALL(Highs_addRow(lp(), lc.lb(), lc.ub(), lc.size(), lc.pvars(), lc.pcoefs()));
}
void HighsModelAPI::AddConstraint(const LinConGE& lc) {
  HIGHS_CCALL(Highs_addRow(lp(), lc.lb(), lc.ub(), lc.size(), lc.pvars(), lc.pcoefs()));
}


void HighsModelAPI::FinishProblemModificationPhase() {
  Highs_writeModel(lp(), "d:/highsout.lp");
}



} // namespace mp
