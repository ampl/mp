#include "coptmodelapi.h"

#include "mp/model-mgr-with-std-pb.h"
#include "mp/flat/model_flattener.h"
#include "mp/flat/redef/MIP/converter_mip.h"


namespace mp {


std::unique_ptr<BasicModelManager>
CreateCoptModelMgr(CoptCommon& cc, Env& e,
                     pre::BasicPresolver*& pPre) {
  using CoptFlatCvt = FlatCvtImpl<MIPFlatConverter, CoptModelAPI>;
  using CoptProblemFlattener = mp::ModelFltImpl<
    mp::ModelFlattener, mp::Problem, CoptFlatCvt>;
  auto pcvt = new CoptProblemFlattener(e);
  auto res = CreateModelManagerWithStdBuilder(
        std::unique_ptr< BasicConverter<mp::Problem> >{ pcvt } );
  pcvt->GetFlatCvt().GetBasicBackend().set_other_copt(&cc);
  cc.set_other_copt(
        &pcvt->GetFlatCvt().GetBasicBackend());
  pPre = &pcvt->GetFlatCvt().GetPresolver();
  return res;
}


void CoptModelAPI::InitProblemModificationPhase() { }

void CoptModelAPI::AddVariables(const VarArrayDef& v) {
  std::vector<char> vtypes(v.size());
  for (size_t i=v.size(); i--; )
    vtypes[i] = var::Type::CONTINUOUS==v.ptype()[i] ?
          COPT_CONTINUOUS : COPT_INTEGER;
  COPT_CCALL(COPT_AddCols(lp(), (int)v.size(), NULL, NULL,
    NULL, NULL, NULL, vtypes.data(), v.plb(), v.pub(),  NULL));
}

void CoptModelAPI::SetLinearObjective( int iobj, const LinearObjective& lo ) {
  if (1>iobj) {
    COPT_CCALL(COPT_SetObjSense(lp(), 
                    obj::Type::MAX==lo.obj_sense() ? COPT_MAXIMIZE : COPT_MINIMIZE) );
    COPT_CCALL(COPT_SetColObj(lp(), lo.num_terms(),
                           lo.vars().data(), lo.coefs().data()) );
  } else {
//    TODO
  }
}
void CoptModelAPI::AddConstraint(const RangeLinCon& lc) {
  /*char sense = COPT_EQUAL;                     // good to initialize things
  double rhs = lc.lb();
  if (lc.lb()==lc.ub())
    sense = COPT_EQUAL;
  else {                                // Let solver deal with lb>~ub etc.
    if (lc.lb()>MinusInfinity()) {
      sense = COPT_GREATER_EQUAL;
    }
    if (lc.ub()<Infinity()) {
      if (COPT_GREATER_EQUAL == sense)
        sense = 0;
      else {
        sense = COPT_LESS_EQUAL;
        rhs = lc.ub();
      }
    }
  }
  */
  COPT_CCALL(COPT_AddRow(lp(), lc.size(), lc.pvars(), lc.pcoefs(), 
    NULL, lc.lb(), lc.ub(), NULL));
  /*
  if ('R'==sense) {
    int indices = NumLinCons()-1;
    double range = lc.ub()-lc.lb();
    COPT_CCALL( CPXchgrngval (env(), lp(), 1, &indices, &range) );
  }
  */
}
// TODO Check if the following make sense, they are not reccomended
void CoptModelAPI::AddConstraint(const LinConLE& lc) {
  char sense = COPT_LESS_EQUAL;
  COPT_CCALL(COPT_AddRow(lp(), lc.size(), lc.pvars(), lc.pcoefs(),
    sense, lc.rhs(), 0, NULL));
}
void CoptModelAPI::AddConstraint(const LinConEQ& lc) {
  char sense = COPT_EQUAL;
  COPT_CCALL(COPT_AddRow(lp(), lc.size(), lc.pvars(), lc.pcoefs(),
    sense, lc.rhs(), 0, NULL));
}
void CoptModelAPI::AddConstraint(const LinConGE& lc) {
  char sense = COPT_GREATER_EQUAL;
  COPT_CCALL(COPT_AddRow(lp(), lc.size(), lc.pvars(), lc.pcoefs(),
    sense, lc.rhs(), 0, NULL));
}

// int COPT_AddIndicator(copt_prob *prob, int binColIdx, int binColVal,
//int nRowMatCnt, const int* rowMatIdx, const double* rowMatElem, char
//cRowSense, double dRowBound)

void CoptModelAPI::AddConstraint(const IndicatorConstraintLinLE &ic)  {
  /*
  COPT_CCALL(COPT_AddIndicator(lp(),
                               ic.get_binary_var(), !ic.get_binary_value(),
                               (int)ic.get_constraint().size(),
                               ic.get_constraint().rhs(), 'L',
                               ic.get_constraint().pvars(),
                               ic.get_constraint().pcoefs(), NULL) );
                               */
}
void CoptModelAPI::AddConstraint(const IndicatorConstraintLinEQ &ic)  {
  /*
  COPT_CCALL( CPXaddindconstr (env(), lp(),
                               ic.get_binary_var(), !ic.get_binary_value(),
                               (int)ic.get_constraint().size(),
                               ic.get_constraint().rhs(), 'E',
                               ic.get_constraint().pvars(),
                               ic.get_constraint().pcoefs(), NULL) );
                               */
}

void CoptModelAPI::AddConstraint(const QuadraticConstraint& qc) {



}

void CoptModelAPI::FinishProblemModificationPhase() {
}


} // namespace mp
