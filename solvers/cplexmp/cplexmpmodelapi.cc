#include "cplexmpmodelapi.h"


namespace mp {

void CplexModelAPI::InitProblemModificationPhase(const FlatModelInfo*) { }

void CplexModelAPI::AddVariables(const VarArrayDef& v) {
  int nint = 0;
  std::vector<char> vtypes(v.size());
  for (size_t i = v.size(); i--; )
    if (var::Type::CONTINUOUS == v.ptype()[i])
      vtypes[i] = CPX_CONTINUOUS;
    else {
      nint++;
      vtypes[i] = (v.plb()[i] == 0) && (v.pub()[i] == 1) ? CPX_BINARY : CPX_INTEGER;
    }
  if (nint > 0)
    CPLEX_CALL(CPXnewcols(env(), lp(), (int)v.size(), nullptr,
      v.plb(), v.pub(), vtypes.data(), const_cast<char**>(v.pnames())));
  else
    CPLEX_CALL(CPXnewcols(env(), lp(), (int)v.size(), nullptr,
      v.plb(), v.pub(), nullptr , const_cast<char**>(v.pnames())));
}

void CplexModelAPI::NoteCPLEXMainObjSense(obj::Type s) { main_obj_sense_ = s; }
obj::Type CplexModelAPI::GetCPLEXMainObjSense() const { return main_obj_sense_; }

void CplexModelAPI::SetLinearObjective( int iobj, const LinearObjective& lo ) {
  char ERROR[512];

  if (1>iobj) {
    CPLEX_CALL( CPXchgobjsen (env(), lp(),
                    obj::Type::MAX==lo.obj_sense() ? CPX_MAX : CPX_MIN) );
    NoteCPLEXMainObjSense(lo.obj_sense());
    CPLEX_CALL( CPXchgobj (env(), lp(), lo.num_terms(),
                           lo.vars().data(), lo.coefs().data()) );
  } else {
    CPXsetnumobjs(env(), lp(), iobj+1);
    auto status = CPXmultiobjsetobj(env(), lp(), iobj,
      lo.num_terms(),
      lo.vars().data(), lo.coefs().data(), CPX_NO_OFFSET_CHANGE,
      lo.obj_sense() == GetCPLEXMainObjSense() ? CPX_NO_WEIGHT_CHANGE : -1.0,
      CPX_NO_PRIORITY_CHANGE, CPX_NO_ABSTOL_CHANGE, CPX_NO_RELTOL_CHANGE,
      lo.name());
    
    if (status)
    {
      CPXgeterrorstring(env(), status, ERROR);
      printf(ERROR);
    }

  }
}
void CplexModelAPI::SetQuadraticObjective(int iobj, const QuadraticObjective& qo) {
  if (1 > iobj) {
    SetLinearObjective(iobj, qo);
    
    const auto& qt = qo.GetQPTerms();
    std::vector<int> qmatbeg(NumVars() + 1), qmatcnt(NumVars());
    int nnondiagonal = 0;
    std::vector<std::map<int, double>> acc(NumVars());
    for (int i = 0; i < qt.size(); i++) {
      if(qt.var1(i)==qt.var2(i))
        acc[qt.var1(i)][qt.var2(i)] = qt.coef(i)*2;
      else {
        nnondiagonal++; // Need to add the corresponding entries
        acc[qt.var1(i)][qt.var2(i)] = qt.coef(i) ;
        acc[qt.var2(i)][qt.var1(i)] = qt.coef(i) ;
      }
    }

    std::vector<int> qmatind(qt.size()+nnondiagonal);
    std::vector<double> qmatval(qt.size() + nnondiagonal);
    int currentbeg = 0;
    for (int i = 0; i < NumVars(); i++)
    {
      qmatbeg[i] = currentbeg;
      qmatcnt[i] = acc[i].size();
      for (auto c : acc[i])
      {
        qmatind[currentbeg] = c.first;
        qmatval[currentbeg] = c.second;
        currentbeg++;
      }
    }
    qmatbeg[NumVars()] = currentbeg;
    CPLEX_CALL((CPXcopyquad(env(), lp(),
      qmatbeg.data(), qmatcnt.data(), qmatind.data(), qmatval.data())));
  }
  else {
    throw std::runtime_error("Multiple quadratic objectives not supported");
  }
}

void CplexModelAPI::AddConstraint(const LinConRange& lc) {
  char sense = 'E';                     // good to initialize things
  double rhs = lc.lb();
  if (lc.lb()==lc.ub())
    sense = 'E';
  else {                                // Let solver deal with lb>~ub etc.
    if (lc.lb()>MinusInfinity()) {
      sense = 'G';
    }
    if (lc.ub()<Infinity()) {
      if ('G'==sense)
        sense = 'R';
      else {
        sense = 'L';
        rhs = lc.ub();
      }
    }
  }
  int rmatbeg[] = { 0 };
  CPLEX_CALL( CPXaddrows (env(), lp(), 0, 1, lc.size(), &rhs,
                          &sense, rmatbeg, lc.pvars(), lc.pcoefs(),
                          NULL, NULL) );
  if ('R'==sense) {
    int indices = NumLinCons()-1;
    double range = lc.ub()-lc.lb();
    CPLEX_CALL( CPXchgrngval (env(), lp(), 1, &indices, &range) );
  }
}
void CplexModelAPI::AddConstraint(const LinConLE& lc) {
  char sense = 'L';                     // good to initialize things
  double rhs = lc.rhs();
  int rmatbeg[] = { 0 };
  CPLEX_CALL( CPXaddrows (env(), lp(), 0, 1, lc.size(), &rhs,
                          &sense, rmatbeg, lc.pvars(), lc.pcoefs(),
                          NULL, NULL) );
}
void CplexModelAPI::AddConstraint(const LinConEQ& lc) {
  char sense = 'E';                     // good to initialize things
  double rhs = lc.rhs();
  int rmatbeg[] = { 0 };
  CPLEX_CALL( CPXaddrows (env(), lp(), 0, 1, lc.size(), &rhs,
                          &sense, rmatbeg, lc.pvars(), lc.pcoefs(),
                          NULL, NULL) );
}
void CplexModelAPI::AddConstraint(const LinConGE& lc) {
  char sense = 'G';                     // good to initialize things
  double rhs = lc.rhs();
  int rmatbeg[] = { 0 };
  CPLEX_CALL( CPXaddrows (env(), lp(), 0, 1, lc.size(), &rhs,
                          &sense, rmatbeg, lc.pvars(), lc.pcoefs(),
                          NULL, NULL) );
}
#define addqc(qc, sense)\
const auto& lt = qc.GetLinTerms();\
const auto& qt = qc.GetQPTerms();\
CPLEX_CALL(CPXaddqconstr(env(), lp(), lt.size(), qt.size(), qc.rhs(),\
  sense, (int*)lt.pvars(), (double*)lt.pcoefs(),\
  (int*)qt.pvars1(), (int*)qt.pvars2(), (double*)qt.pcoefs(), qc.name()));

void CplexModelAPI::AddConstraint(const QuadConLE& qc) {
  addqc(qc, 'L');
}

void CplexModelAPI::AddConstraint(const QuadConEQ& qc) {
  addqc(qc, 'L');
  CPLEX_CALL(CPXaddqconstr(env(), lp(), lt.size(), qt.size(), qc.rhs(),
    'G', (int*)lt.pvars(), (double*)lt.pcoefs(),
    (int*)qt.pvars1(), (int*)qt.pvars2(), (double*)qt.pcoefs(), qc.name()));
}

void CplexModelAPI::AddConstraint(const QuadConGE& qc) {
  addqc(qc, 'G');
}

#define addsos(cc, sostype)\
  int beg = 0;\
  char type = sostype;\
  CPLEX_CALL(CPXaddsos(env(), lp(), 1, cc.size(), &type, &beg,\
    (int*)cc.get_vars().data(), (double*)cc.get_weights().data(),\
    NULL) );\

void CplexModelAPI::AddConstraint(const SOS1Constraint& cc) {
  addsos(cc, CPX_TYPE_SOS1);
}

void CplexModelAPI::AddConstraint(const SOS2Constraint& cc) {
  addsos(cc, CPX_TYPE_SOS2);
}
#define addindic(ic, direction)\
  CPLEX_CALL(CPXaddindconstr(env(), lp(),\
    ic.get_binary_var(), !ic.get_binary_value(),\
    (int)ic.get_constraint().size(),\
    ic.get_constraint().rhs(), direction,\
    ic.get_constraint().pvars(),\
    ic.get_constraint().pcoefs(), NULL));

void CplexModelAPI::AddConstraint(const IndicatorConstraintLinLE &ic)  {
  addindic(ic, 'L');
}
void CplexModelAPI::AddConstraint(const IndicatorConstraintLinEQ &ic)  {
  addindic(ic, 'F');
}
void CplexModelAPI::AddConstraint(const IndicatorConstraintLinGE &ic)  {
  addindic(ic, 'G');
}

void CplexModelAPI::AddConstraint(const PLConstraint& plc) {
  const auto& plp = plc.GetParameters().GetPLPoints();
  CPLEX_CALL( CPXaddpwl(env(), lp(),
              plc.GetResultVar(), plc.GetArguments()[0],
              plp.PreSlope(), plp.PostSlope(),
              plp.x_.size(), plp.x_.data(), plp.y_.data(), NULL) );
}

void CplexModelAPI::FinishProblemModificationPhase() {
}


} // namespace mp
