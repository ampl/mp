#include "coptmodelapi.h"

namespace mp {

void CoptModelAPI::InitProblemModificationPhase(const FlatModelInfo* fi) {
  num_lin_obj_ = fi->GetObjInfo()[0];
}

void CoptModelAPI::AddVariables(const VarArrayDef& v) {
  std::vector<char> vtypes(v.size());
  for (size_t i=v.size(); i--; )
    vtypes[i] = var::Type::CONTINUOUS==v.ptype()[i] ?
          COPT_CONTINUOUS : COPT_INTEGER;
  COPT_CCALL(COPT_AddCols(lp(), (int)v.size(), NULL, NULL,
    NULL, NULL, NULL, vtypes.data(), v.plb(), v.pub(),  v.pnames()));
}

void CoptModelAPI::SetLinearObjective( int iobj, const LinearObjective& lo ) {
  int sense = obj::Type::MAX == lo.obj_sense() ? COPT_MAXIMIZE : COPT_MINIMIZE;
  senses()[iobj] = sense;

  if ((iobj<1) && (num_lin_obj_ <= 1)) {
    COPT_CCALL(COPT_SetObjSense(lp(), sense ));
    double zero_out = 0.0;
    for (int i=NumVars(); i--; )
      COPT_CCALL(COPT_SetColObj(lp(), 1, &i, &zero_out));
    COPT_CCALL(COPT_SetColObj(lp(), lo.num_terms(),
                              lo.vars().data(), lo.coefs().data()) );
    
  } else {
    COPT_CCALL(COPT_MultiObjSetColObj(lp(), iobj, lo.num_terms(),
                                      lo.vars().data(), lo.coefs().data()));
    COPT_CCALL(COPT_MultiObjSetObjSense(lp(), iobj, senses()[0]));
    COPT_CCALL(COPT_MultiObjSetObjConst(lp(), iobj, 0.0)); // Seems necessary for MultiObjParam in 8.0.3.
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


void CoptModelAPI::AddGenericCone(int conetype,
    const int* vars,
    const double* coeffs,
    int dim,
    const char* name) {

    if (conetype ==COPT_EXPCONE_PRIMAL)
    {
        if (dim != 3)
            throw std::runtime_error(
                "ExponentialConeConstraint must have exactly 3 entries");
    }

    int coneType[1] = { conetype };
    int coneBeg[1] = { 0 };
    int coneCnt[1] = { dim };

    if (std::all_of(coeffs, coeffs+dim, [](double x) { // check if all are 1.0
        return std::abs(x - 1.0) <= 1e-12;
        })) {

		// If all coefficients are one, use the shortcut functions:
        switch(conetype) {
            case COPT_CONE_QUAD: // quadratic cone
                COPT_CCALL(COPT_AddCones(lp(), 1, coneType, coneBeg, coneCnt, vars));
                break;
            case COPT_CONE_RQUAD: 
                COPT_CCALL(COPT_AddCones(lp(), 1, coneType, coneBeg, coneCnt, vars));
                break;
            case COPT_EXPCONE_PRIMAL:
                COPT_CCALL(COPT_AddExpCones(lp(), 1, coneType, vars));
                break;
            default:
                throw std::runtime_error("Unsupported cone type");
		}
    }
    // Otherwise use affine cones
    // General affine case.
    std::vector<int> rowMatBeg(dim);
    std::vector<int> rowMatCount(dim, 1);
    std::vector<double> rowMatConst(dim, 0.0);

    for (int i = 0; i < dim; i++)
        rowMatBeg[i] = i;

    int ret = COPT_AddAffineCone(lp(), conetype,
        dim,  // dimension
        0, NULL, // unused
        NULL, NULL, NULL, NULL, // PSD data
        rowMatBeg.data(),
        rowMatCount.data(),
        vars, coeffs,
        rowMatConst.data(),
        name);
    if (ret != 0) throw std::runtime_error("COPT_AddAffineCone failed");



}

void CoptModelAPI::AddConstraint(
    const QuadraticConeConstraint& qc) {
    const auto& vars = qc.GetArguments();
    const auto& coeffs = qc.GetParameters();
    AddGenericCone(COPT_CONE_QUAD, vars.data(), coeffs.data(),
        static_cast<int>(vars.size()), qc.GetName());
}


void CoptModelAPI::AddConstraint(
    const RotatedQuadraticConeConstraint& qc) {
    const auto& vars = qc.GetArguments();
    const auto& coeffs = qc.GetParameters();
    AddGenericCone(COPT_CONE_RQUAD, vars.data(), coeffs.data(),
        static_cast<int>(vars.size()), qc.GetName());
}


void CoptModelAPI::AddConstraint(const ExponentialConeConstraint& ec) {
        const auto& vars = ec.GetArguments();
        const auto& coeffs = ec.GetParameters();
        AddGenericCone(COPT_EXPCONE_PRIMAL, vars.data(), coeffs.data(),
            static_cast<int>(vars.size()), ec.GetName());
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

NLParams CoptModelAPI::AddExpression(const AbsExpression& e) {
    return CreateExpressionOneArg(e, COPT_NL_ABS);
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


NLParams CoptModelAPI::AddExpression(const AsinExpression& e) {
    return CreateExpressionOneArg(e, COPT_NL_ASIN);
}
NLParams CoptModelAPI::AddExpression(const AcosExpression& e) {
    return CreateExpressionOneArg(e, COPT_NL_ACOS);
}
NLParams CoptModelAPI::AddExpression(const AtanExpression& e) {
    return CreateExpressionOneArg(e, COPT_NL_ATAN);
}

NLParams CoptModelAPI::AddExpression(const SinhExpression& e) {
    return CreateExpressionOneArg(e, COPT_NL_SINH);
}
NLParams CoptModelAPI::AddExpression(const CoshExpression& e) {
    return CreateExpressionOneArg(e, COPT_NL_COSH);
}
NLParams CoptModelAPI::AddExpression(const TanhExpression& e) {
    return CreateExpressionOneArg(e, COPT_NL_TANH);
}
NLParams CoptModelAPI::AddExpression(const AsinhExpression& e) {
    return CreateExpressionOneArg(e, COPT_NL_ASINH);
}
NLParams CoptModelAPI::AddExpression(const AcoshExpression& e) {
    return CreateExpressionOneArg(e, COPT_NL_ACOSH);
}
NLParams CoptModelAPI::AddExpression(const AtanhExpression& e) {
    return CreateExpressionOneArg(e, COPT_NL_ATANH);
}


void CoptModelAPI::SetNLObjective(int i, const NLObjective& nlo) {
    const auto& exp = GetExpression(nlo);
    if (i == 0)
    {
        COPT_CCALL(COPT_SetObjSense(lp(),
            obj::Type::MAX == nlo.obj_sense() ? COPT_MAXIMIZE : COPT_MINIMIZE));
        COPT_CCALL(COPT_SetNLObj(lp(), exp.nTokens(), exp.nTokenElements(), exp.tokens(), exp.tokenElements()));
    }
    else {
        MP_RAISE("Multiple non-linear objectives not supported natively. Use multi-objective emulator by setting obj:multi=2");
    }
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
