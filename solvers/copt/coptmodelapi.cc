#include "coptmodelapi.h"

namespace mp {

void CoptModelAPI::InitProblemModificationPhase(const FlatModelInfo*) { }

void CoptModelAPI::AddVariables(const VarArrayDef& v) {
  std::vector<char> vtypes(v.size());
  for (size_t i=v.size(); i--; )
    vtypes[i] = var::Type::CONTINUOUS==v.ptype()[i] ?
          COPT_CONTINUOUS : COPT_INTEGER;
  COPT_CCALL(COPT_AddCols(lp(), (int)v.size(), NULL, NULL,
    NULL, NULL, NULL, vtypes.data(), v.plb(), v.pub(),  v.pnames()));
}

void CoptModelAPI::SetLinearObjective( int iobj, const LinearObjective& lo ) {
  if (iobj<1) {
    COPT_CCALL(COPT_SetObjSense(lp(), 
                    obj::Type::MAX==lo.obj_sense() ? COPT_MAXIMIZE : COPT_MINIMIZE) );
    double zero_out = 0.0;
    for (int i=NumVars(); i--; )
      COPT_CCALL(COPT_SetColObj(lp(), 1, &i, &zero_out));
    COPT_CCALL(COPT_SetColObj(lp(), lo.num_terms(),
                           lo.vars().data(), lo.coefs().data()) );
  } else {
//    TODO
  }
}


void CoptModelAPI::SetQuadraticObjective(int iobj, const QuadraticObjective& qo) {
  if (1 > iobj) {
    SetLinearObjective(iobj, qo);                         // add the linear part
    const auto& qt = qo.GetQPTerms();
    COPT_CCALL(COPT_SetQuadObj(lp(), qt.size(), (int*)qt.pvars1(), (int*)qt.pvars2(),
      (double*)qt.pcoefs()));
  }
  else {
    throw std::runtime_error("Multiple quadratic objectives not supported");
  }
}

void CoptModelAPI::AddConstraint(const LinConRange& lc) {
	COPT_CCALL(COPT_AddRow(lp(), lc.size(),
												 lc.pvars(), lc.pcoefs(), 0,
												 lc.lb() < -COPT_INFINITY ? -COPT_INFINITY: lc.lb(),
												 lc.ub() > COPT_INFINITY ? COPT_INFINITY : lc.ub(),
												 lc.name()));
}

void CoptModelAPI::AddConstraint(const LinConLE& lc) {
  char sense = COPT_LESS_EQUAL;
  COPT_CCALL(COPT_AddRow(lp(), lc.size(), lc.pvars(), lc.pcoefs(),
		sense, lc.rhs(), 0, lc.name()));
}
void CoptModelAPI::AddConstraint(const LinConEQ& lc) {
  char sense = COPT_EQUAL;
  COPT_CCALL(COPT_AddRow(lp(), lc.size(), lc.pvars(), lc.pcoefs(),
		sense, lc.rhs(), 0, lc.name()));
}
void CoptModelAPI::AddConstraint(const LinConGE& lc) {
  char sense = COPT_GREATER_EQUAL;
  COPT_CCALL(COPT_AddRow(lp(), lc.size(), lc.pvars(), lc.pcoefs(),
		sense, lc.rhs(), 0, lc.name()));
}

void CoptModelAPI::AddConstraint(const QuadConLE& qc) {
  const auto& lt = qc.GetLinTerms();
  const auto& qt = qc.GetQPTerms();
  COPT_CCALL(COPT_AddQConstr(lp(), lt.size(), (int*)lt.pvars(), (double*)lt.pcoefs(),
                             qt.size(), (int*)qt.pvars1(), (int*)qt.pvars2(),
														 (double*)qt.pcoefs(), COPT_LESS_EQUAL, qc.rhs(), qc.name()));
}

void CoptModelAPI::AddConstraint(const QuadConEQ& qc) {
  const auto& lt = qc.GetLinTerms();
  const auto& qt = qc.GetQPTerms();
  COPT_CCALL(COPT_AddQConstr(lp(), lt.size(), (int*)lt.pvars(), (double*)lt.pcoefs(),
                             qt.size(), (int*)qt.pvars1(), (int*)qt.pvars2(),
														 (double*)qt.pcoefs(), COPT_EQUAL, qc.rhs(), qc.name()));
}

void CoptModelAPI::AddConstraint(const QuadConGE& qc) {
  const auto& lt = qc.GetLinTerms();
  const auto& qt = qc.GetQPTerms();
  COPT_CCALL(COPT_AddQConstr(lp(), lt.size(), (int*)lt.pvars(), (double*)lt.pcoefs(),
                             qt.size(), (int*)qt.pvars1(), (int*)qt.pvars2(),
														 (double*)qt.pcoefs(), COPT_GREATER_EQUAL, qc.rhs(), qc.name()));
}



void CoptModelAPI::AddConstraint(const IndicatorConstraintLinLE &ic)  {
  COPT_CCALL(COPT_AddIndicator(lp(),
    ic.get_binary_var(), ic.get_binary_value(),
    (int)ic.get_constraint().size(),
    ic.get_constraint().pvars(),
    ic.get_constraint().pcoefs(),
    COPT_LESS_EQUAL,
    ic.get_constraint().rhs()));
                               
}
void CoptModelAPI::AddConstraint(const IndicatorConstraintLinEQ &ic)  {
  COPT_CCALL(COPT_AddIndicator(lp(),
    ic.get_binary_var(), ic.get_binary_value(),
    (int)ic.get_constraint().size(),
    ic.get_constraint().pvars(),
    ic.get_constraint().pcoefs(),
    COPT_EQUAL,
    ic.get_constraint().rhs()));
}
void CoptModelAPI::AddConstraint(const IndicatorConstraintLinGE &ic)  {
  COPT_CCALL(COPT_AddIndicator(lp(),
    ic.get_binary_var(), ic.get_binary_value(),
    (int)ic.get_constraint().size(),
    ic.get_constraint().pvars(),
    ic.get_constraint().pcoefs(),
    COPT_GREATER_EQUAL,
    ic.get_constraint().rhs()));
}

void CoptModelAPI::AddConstraint(const SOS1Constraint& sos) {
  int type = COPT_SOS_TYPE1;
  int beg = 0;
  const int size = sos.size();
  COPT_CCALL(COPT_AddSOSs(lp(), 1, &type, &beg,
    &size, (int*)sos.get_vars().data(),
    (double*)sos.get_weights().data()));
}

void CoptModelAPI::AddConstraint(const SOS2Constraint& sos) {
  int type = COPT_SOS_TYPE2;
  int beg = 0;
  const int size = sos.size();
  COPT_CCALL(COPT_AddSOSs(lp(), 1, &type, &beg,
    &size, (int*)sos.get_vars().data(),
    (double*)sos.get_weights().data()));
}


void CoptModelAPI::FinishProblemModificationPhase() {
}

template <class MPExpr>
void CoptModelAPI::AppendLinAndConstTerms(Expr& exp, const MPExpr& ae) {
    double ct = GetConstTerm(ae);
    int size = GetLinSize(ae);
    if (ct)
    {
        // TODO - what about a constant?
        throw std::runtime_error("Constants are not supported yet");
    }
    
    for (int i = 0; i < size; ++i) {
        auto index = GetLinTerm(ae, i).tokens()[0];
        auto coeff = GetLinCoef(ae, i);
        exp.addLinear(index, coeff);
    }
   
}

NLParams CoptModelAPI::AddExpression(const NLAffineExpression& ae) {
    NLParams exp;


    double ct = GetConstTerm(ae);
    bool hasConst = ct != 0;
    int size = GetLinSize(ae);
  
    if (size > 1 || hasConst)
    {
        exp.addOp(COPT_NL_SUM);
        exp.addVar(size + hasConst);
    }
    for (int i = 0; i < size; ++i) {
        auto index = GetLinTerm(ae, i);
        auto coeff = GetLinCoef(ae, i);
        if (coeff != 1.0) {
            exp.addOp(COPT_NL_MULT);
            exp.addConstant(coeff);
        }
        exp.addMembers(index);
    }
    if (hasConst)
        exp.addConstant(ct);
    return exp;
}

NLParams CoptModelAPI::AddExpression(const NLQuadExpression& qe) {
    NLParams quad;
    AppendLinAndConstTerms(quad, qe);
    int size = GetQuadSize(qe);
  
    if (size > 1)
    {
        quad.addOp(COPT_NL_SUM);
        quad.addVar(size);
    }
    for (int i = 0; i < GetQuadSize(qe); ++i) {
        if (double coef = GetQuadCoef(qe, i)) {
            if (1.0 != coef) {
                quad.addOp(COPT_NL_MULT);
                quad.addConstant(coef);
            }
            quad.addOp(COPT_NL_MULT);
            quad.addMembers(GetQuadTerm1(qe, i));
            quad.addMembers(GetQuadTerm2(qe, i));
        }
    }
    return quad;
}
NLParams CoptModelAPI::AddExpression(const SinExpression& e) {
    return CreateExpressionOneArg(e, COPT_NL_SIN);
}
NLParams CoptModelAPI::AddExpression(const CosExpression& e) {
    return CreateExpressionOneArg(e, COPT_NL_COS);
}
NLParams CoptModelAPI::AddExpression(const TanExpression& e) {
    return CreateExpressionOneArg(e, COPT_NL_TAN);
}


void CoptModelAPI::AddConstraint(const NLConstraint& nl) {
    auto exp = GetExpression(nl);
    double lhs = GetLower(nl), rhs = GetUpper(nl);
    char type;
    double range;
    double* prange = nullptr;
   
    lhs = lhs < -COPT_INFINITY ? -COPT_INFINITY : lhs;
    rhs = rhs > COPT_INFINITY ? COPT_INFINITY : rhs;

    if (GetLinSize(nl) > 0) {
        exp.reserveLinear(GetLinSize(nl));
        for (int i = 0; i < GetLinSize(nl); ++i)
        {
            exp.addLinear(GetLinVar(nl, i), GetLinCoef(nl, i));
            // coeff*var

        }
    }
    COPT_CCALL(COPT_AddNLConstr(lp(), exp.nTokens(), exp.nTokenElements(),
        exp.tokens(), exp.tokenElements(), exp.nLinear(), exp.linearIndices(),
        exp.linearCoeffs(), 0, lhs, rhs, 0));
}
/*
* sin2 is
c0: x3<=0.5
c1: -inf <= ((-1 * x2) + sin(((2 * x1)))) <= 0
c2: x3 >= sin(x1)
o0: maximize x1+x2
*/
void CoptModelAPI::AddGlobalConstraint(const NLParams& exp, char sense) {
    COPT_CCALL(COPT_AddNLConstr(lp(), exp.nTokens(), exp.nTokenElements(),
        exp.tokens(), exp.tokenElements(), exp.nLinear(), exp.linearIndices(),
        exp.linearCoeffs(), sense, 0, 0, "name"));
   

}
void CoptModelAPI::AddConstraint(const NLAssignEQ& neq) {

    NLParams params;
    params.addOp(COPT_NL_SUM);
    params.addVar(2); // number of items 
    params.addOp(COPT_NL_NEG);
    params.addVar(GetVariable(neq));
    params.addMembers(GetExpression(neq));

    AddGlobalConstraint(params, COPT_EQUAL);

}
void CoptModelAPI::AddConstraint(const NLAssignGE& nge) {
    NLParams params;
    params.addOp(COPT_NL_SUM);
    params.addVar(2); // number of items 
    params.addOp(COPT_NL_NEG);
    params.addVar(GetVariable(nge));
    params.addMembers(GetExpression(nge));

    AddGlobalConstraint(params, COPT_LESS_EQUAL);

}
void CoptModelAPI::AddConstraint(const NLAssignLE& nle) {
    NLParams params;
    params.addOp(COPT_NL_SUM);
    params.addVar(2); // number of items 
    params.addOp(COPT_NL_NEG);
    params.addVar(GetVariable(nle));
    params.addMembers(GetExpression(nle));

    AddGlobalConstraint(params, COPT_GREATER_EQUAL);
}
NLParams CoptModelAPI::AddExpression(const DivExpression& e) {
    NLParams exp;
    exp.addOp(COPT_NL_DIV);
    exp.addMembers(GetArgExpression(e, 0));
    exp.addMembers(GetArgExpression(e, 1));
    return exp;
}



NLParams CoptModelAPI::AddExpression(const LogExpression& e) {
    return CreateExpressionOneArg(e, COPT_NL_LOG);
}

NLParams CoptModelAPI::AddExpression(const LogAExpression& e) {
    NLParams exp;
    exp.addOp(COPT_NL_DIV);
    exp.addOp(COPT_NL_LOG);
    auto ex = GetArgExpression(e, 0);
    exp.addMembers(ex);
    exp.addOp(COPT_NL_LOG);
    auto par = GetParameter(e, 0);
    exp.addConstant(par);
    return exp;

}
NLParams CoptModelAPI::AddExpression(const ExpExpression& e) {
    return CreateExpressionOneArg(e, COPT_NL_EXP);
}

NLParams CoptModelAPI::AddExpression(const ExpAExpression&e) {
    // base ^ x = e ^ (x * ln(b)) )
    NLParams arg;
    arg.addOp(COPT_NL_EXP);
    arg.addOp(COPT_NL_MULT);
    arg.addMembers(GetArgExpression(e, 0));
    arg.addOp(COPT_NL_LOG);
    arg.addConstant(GetParameter(e, 0));
    return arg;
}



NLParams CoptModelAPI::AddExpression(const PowConstExpExpression& e) {
    NLParams arg;
    arg.addOp(COPT_NL_POW);
    arg.addMembers(GetArgExpression(e, 0));
    arg.addConstant(GetParameter(e, 0));
    return arg;
}

} // namespace mp
