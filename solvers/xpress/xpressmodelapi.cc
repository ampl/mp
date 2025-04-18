#include "xpressmodelapi.h"


namespace mp {

void XpressmpModelAPI::InitProblemModificationPhase(
    const FlatModelInfo*) {
}

std::string& myreplace(std::string& s, const std::string& from, const std::string& to)
{
  for (size_t pos = 0; (pos = s.find(from, pos)) != std::string::npos; pos += to.size())
    s.replace(pos, from.size(), to);
  return s;
}

std::string sanitizeName(std::string n, const std::string &lastName) {
  // Xpress does not like square brackets or spaces in variable names
  std::replace(n.begin(), n.end(), '[', '(');
  std::replace(n.begin(), n.end(), ']', ')');
  std::replace(n.begin(), n.end(), ' ', '_');
  myreplace(n,"\'", "-");
  myreplace(n, "\"", "--");
  // Workaround for duplicated names:
  if (n == lastName) n += "_a";
  return n;
}
void XpressmpModelAPI::AddVariables(const VarArrayDef& v) {
  std::vector<double> objs(v.size(), 0);
  std::vector<int> iii(v.size(), 0);
  XPRESSMP_CCALL(XPRSaddcols(lp(), v.size(), 0, objs.data(), iii.data(), NULL,
    NULL, v.plb(), v.pub()));
  std::string name,lastname;
  
  if (v.pnames() != NULL)
  {
    fmt::MemoryWriter w;
    for (int i = 0; i < v.size(); ++i)
    {
      name = sanitizeName(v.pnames()[i], lastname);
      w << name << '\0';
      lastname = name;
    }   
    XPRESSMP_CCALL(XPRSaddnames(lp(), 2, w.c_str(), 0, v.size() - 1));
  }
  // All variables are continuous by default, set the integer ones
  std::vector<int> intIndices;
  intIndices.reserve(v.size());
  std::vector<char> types;
  types.reserve(v.size());
  for (auto i = 0; i < v.size(); i++)
    if (v.ptype()[i] == var::Type::INTEGER) {
      intIndices.push_back(i);
      if (!v.plb()[i] && 1.0==v.pub()[i]) {
        types.push_back('B');   // Asks explicitly for indicators in 9.2.0
      }
      else {
        types.push_back('I');
      }
    }
  get_other()->numIntVars(intIndices.size());
  if (get_other()->numIntVars() > 0)
  {
    XPRESSMP_CCALL(XPRSchgcoltype(lp(), intIndices.size(),
      intIndices.data(), types.data()));
  }
}

void XpressmpModelAPI::SetLinearObjective( int iobj, const LinearObjective& lo ) {
  if (iobj<1) {
    if (lo.obj_sense() == obj::Type::MAX)
      XPRESSMP_CCALL(XPRSchgobjsense(lp(), XPRS_OBJ_MAXIMIZE));
    if (obj_ind_save_.size()) {
      std::vector<double> obj_coef_0(obj_ind_save_.size(), 0.0);
      XPRESSMP_CCALL(XPRSchgobj(lp(), obj_ind_save_.size(),
                                obj_ind_save_.data(), obj_coef_0.data()));
    }
    XPRESSMP_CCALL(XPRSchgobj(lp(), lo.num_terms(), lo.vars().data(), lo.coefs().data()));
    obj_ind_save_ = lo.vars();
  } else {
    // All objectives must have the same sense, so we will have to automatically 
    // set a conflicting objective's weight to -1 
    obj::Type mainObjSense = getDblAttr(XPRS_OBJSENSE) == -1.0 ? obj::Type::MAX : obj::Type::MIN;
    double weight = lo.obj_sense() == mainObjSense ? 1.0 : -1.0; 
    XPRESSMP_CCALL(XPRSaddobj(lp(), lo.num_terms(), lo.vars().data(), lo.coefs().data(), 0, weight));
  }
}

void XpressmpModelAPI::SetQuadraticObjective(int iobj, const QuadraticObjective& qo) {
  if (1 > iobj) {
    // fmt::format("Setting first quadratic objective\n");
    SetLinearObjective(iobj, qo);                         // add the linear part
    const auto& qt = qo.GetQPTerms();
    std::vector<double> coeffs(qt.coefs());
    for (std::size_t i = 0; i < qt.size(); i++)
      if (qt.pvars1()[i] == qt.pvars2()[i]) coeffs[i] *= 2;

    if (qobj_ind1_save_.size()) {
      assert(qobj_ind1_save_.size()==qobj_ind2_save_.size());
      std::vector<double> qobj_coef_0(qobj_ind1_save_.size(), 0.0);
      XPRESSMP_CCALL(XPRSchgmqobj(lp(), qobj_ind1_save_.size(),
                                  qobj_ind1_save_.data(), qobj_ind2_save_.data(),
                                  qobj_coef_0.data()));
    }
    // fmt::format("Quadratic part is made of {} terms\n", qt.size());
    XPRESSMP_CCALL(XPRSchgmqobj(lp(), qt.size(),
      (int*)qt.pvars1(), (int*)qt.pvars2(),
      (double*)coeffs.data()));
    qobj_ind1_save_ = qt.vars1();
    qobj_ind2_save_ = qt.vars2();
  }
  else {
    throw std::runtime_error("Multiple quadratic objectives not supported");
  }
}

void XpressmpModelAPI::AddConstraint(const LinConLE& lc) {
  char type[] = { 'L' };
  double rhs[] = { lc.rhs() };
  int start[] = { 0 };
  XPRESSMP_CCALL(XPRSaddrows(lp(), 1, lc.coefs().size(), type,
    rhs, NULL, start, lc.pvars(), lc.pcoefs()));
}
void XpressmpModelAPI::AddConstraint(const LinConEQ& lc) {
  char type[] = { 'E' };
  double rhs[] = { lc.rhs() };
  int start[] = { 0 };
  XPRESSMP_CCALL(XPRSaddrows(lp(), 1, lc.coefs().size(), type,
    rhs, NULL, start, lc.pvars(), lc.pcoefs()));
}
void XpressmpModelAPI::AddConstraint(const LinConGE& lc) {
  char type[] = { 'G' };
  double rhs[] = { lc.rhs() };
  int start[] = { 0 };
  XPRESSMP_CCALL(XPRSaddrows(lp(), 1, lc.coefs().size(), type,
    rhs, NULL, start, lc.pvars(), lc.pcoefs()));
}

// Workaround as mp does not distinguish binary variables and
// xpress does not consider "binary" an integer where lb=ub=0 or 1.
// This creates a problem if that variable is used as binary in an indicator
// constraint, so we manually override this here
#define _addIndicator_mp AddConstraint(ic.get_constraint());\
int rowindex[] = { NumLinCons() - 1 };\
int colindex[] = { ic.get_binary_var() };\
int complement[] = { ic.get_binary_value() ? 1 : -1 };\
XPRESSMP_CCALL(XPRSsetindicators(lp(), 1, rowindex, colindex, complement));

void XpressmpModelAPI::AddConstraint(const IndicatorConstraintLinLE &ic)  {
  _addIndicator_mp
}
void XpressmpModelAPI::AddConstraint(const IndicatorConstraintLinEQ &ic)  {
  _addIndicator_mp
}
void XpressmpModelAPI::AddConstraint(const IndicatorConstraintLinGE &ic)  {
  _addIndicator_mp
}

void XpressmpModelAPI::AddLinTerms(
    XPRSprob lp, const LinTerms& lt, double rhsc,  const char typec) {
  char type[] = { typec };
  double rhs[] = { rhsc };
  int start[] = { 0 };
  XPRESSMP_CCALL(XPRSaddrows(lp, 1, lt.coefs().size(), type,
    rhs, NULL, start, lt.pvars(), lt.pcoefs()));
}

#define addqp(type)    numQuadCons(numQuadCons()+1);\
const auto& lt = qc.GetLinTerms();\
AddLinTerms(lp(), lt, qc.rhs(), type); \
const auto& qt = qc.GetQPTerms();\
  int row = NumLinCons() - 1;\
  for (int i = 0; i < qt.size(); i++)\
    XPRESSMP_CCALL(XPRSchgqrowcoeff(lp(), row, \
      qt.var1(i), qt.var2(i), \
      qt.var1(i)==qt.var2(i) ? qt.coef(i) : 0.5*qt.coef(i)));

void XpressmpModelAPI::AddConstraint( const QuadConLE& qc ) {
  addqp('L')
}

void XpressmpModelAPI::AddConstraint( const QuadConEQ& qc ) {
  addqp('E');
}

void XpressmpModelAPI::AddConstraint( const QuadConGE& qc ) {
  addqp('G');
}

void XpressmpModelAPI::AddConstraint(const SOS1Constraint& sos) {
  char type[] = { '1'};
  const int beg = 0;
  const int size = sos.size();
  XPRESSMP_CCALL(XPRSaddsets(lp(), 1, sos.size(), type, &beg, (int*)sos.get_vars().data(),
    (double*)sos.get_weights().data()));
}

void XpressmpModelAPI::AddConstraint(const SOS2Constraint& sos) {
  char type[] = { '2' };
  const int beg = 0;
  const int size = sos.size();
  XPRESSMP_CCALL(XPRSaddsets(lp(), 1, sos.size(), type, &beg, (int*)sos.get_vars().data(),
    (double*)sos.get_weights().data()));
}

template <class Args, class Params, class NumOrLogic, class Id>
void XpressmpModelAPI::addGenCon(
    const CustomFunctionalConstraint<Args, Params, NumOrLogic, Id>& c,
    int xpressConType)
{
  int type[] = { xpressConType };
  int resultant[] = { c.GetResultVar() };
  int colstart[] = { 0 };
  const auto args = c.GetArguments();
  auto colindices = args.data();
  XPRESSMP_CCALL(XPRSaddgencons(lp(), 1, (int)args.size(),
    0, type, resultant, colstart, colindices, NULL, NULL));
}

void XpressmpModelAPI::AddConstraint(const AbsConstraint& c) {
  addGenCon(c, XPRS_GENCONS_ABS);
}

void XpressmpModelAPI::AddConstraint(const MaxConstraint& c) {
  addGenCon(c, XPRS_GENCONS_MAX);
}

void XpressmpModelAPI::AddConstraint(const MinConstraint& c) {
  addGenCon(c, XPRS_GENCONS_MIN);
}

void XpressmpModelAPI::AddConstraint(const AndConstraint& c) {
  addGenCon(c, XPRS_GENCONS_AND);
}

void XpressmpModelAPI::AddConstraint(const OrConstraint& c) {
  addGenCon(c, XPRS_GENCONS_OR);
}


void XpressmpModelAPI::AddConstraint(const DivConstraint& cc) {
  int v1 = cc.GetArguments()[0];
  int v2 = cc.GetArguments()[1];
  NLParams params(cc.GetResultVar());
  params.addMember(XPRS_TOK_COL, v1);
  params.addMember(XPRS_TOK_COL, v2);
  params.addMember(XPRS_TOK_OP, XPRS_OP_DIVIDE);
  params.addMember(XPRS_TOK_EOF, 0);
  AddGlobalConstraint(params, 'E');
}

void XpressmpModelAPI::AddConstraint(const PowConstExpConstraint& cc) {
  double exponent = cc.GetParameters()[0];
  if (exponent == 0.5)
    AddGlobalConstraint(cc.GetResultVar(), cc.GetArguments()[0], XPRS_IFUN_SQRT);
  else
  {
    NLParams params(cc.GetResultVar());
    params.addMember(XPRS_TOK_COL, cc.GetArguments()[0]);
    params.addMember(XPRS_TOK_CON, cc.GetParameters()[0]);
    params.addMember(XPRS_TOK_OP, XPRS_OP_EXPONENT);
    params.addMember(XPRS_TOK_EOF, 0);
    AddGlobalConstraint(params, 'E');
  }
}
  
void XpressmpModelAPI::AddConstraint(const ExpConstraint& cc) {
  AddGlobalConstraint(cc.GetResultVar(), cc.GetArguments()[0], XPRS_IFUN_EXP);
}
void XpressmpModelAPI::AddConstraint(const SinConstraint& cc) {
    AddGlobalConstraint(cc.GetResultVar(), cc.GetArguments()[0], XPRS_IFUN_SIN);
}
void XpressmpModelAPI::AddConstraint(const CosConstraint& cc) {
  AddGlobalConstraint(cc.GetResultVar(), cc.GetArguments()[0], XPRS_IFUN_COS);
}
void XpressmpModelAPI::AddConstraint(const TanConstraint& cc) {
  AddGlobalConstraint(cc.GetResultVar(), cc.GetArguments()[0], XPRS_IFUN_TAN);
}
void XpressmpModelAPI::AddConstraint(const AsinConstraint& cc) {
  AddGlobalConstraint(cc.GetResultVar(), cc.GetArguments()[0], XPRS_IFUN_ARCSIN);
}
void XpressmpModelAPI::AddConstraint(const AcosConstraint& cc) {
  AddGlobalConstraint(cc.GetResultVar(), cc.GetArguments()[0], XPRS_IFUN_ARCCOS);
}
void XpressmpModelAPI::AddConstraint(const AtanConstraint& cc) {
  AddGlobalConstraint(cc.GetResultVar(), cc.GetArguments()[0], XPRS_IFUN_ARCTAN);
}

void XpressmpModelAPI::AddConstraint(const LogConstraint& cc) {
  AddGlobalConstraint(cc.GetResultVar(), cc.GetArguments()[0], XPRS_IFUN_LN);
}
void XpressmpModelAPI::AddConstraint(const LogAConstraint& cc) {
  // log(x) / log(a) RPN:   ) x ln ) a ln /
  NLParams params(cc.GetResultVar());
  params.addMember(XPRS_TOK_RB, 0);
  params.addMember(XPRS_TOK_COL, cc.GetArguments()[0]);
  params.addMember(XPRS_TOK_IFUN, XPRS_IFUN_LN);
  params.addMember(XPRS_TOK_RB,0 );
  params.addMember(XPRS_TOK_CON, cc.GetParameters()[0]);
  params.addMember(XPRS_TOK_IFUN, XPRS_IFUN_LN);
  params.addMember(XPRS_TOK_OP, XPRS_OP_DIVIDE);
  params.addMember(XPRS_TOK_EOF, 0);
  AddGlobalConstraint(params, 'E');
}
void XpressmpModelAPI::AddConstraint(const ExpAConstraint& cc) {
  // base^x = e^(x*ln(b)) ) b  ln  x  *  exp
  NLParams params(cc.GetResultVar());
  params.addMember(XPRS_TOK_RB, 0);
  params.addMember(XPRS_TOK_COL, cc.GetArguments()[0]);
  params.addMember(XPRS_TOK_RB, 0);
  params.addMember(XPRS_TOK_CON, cc.GetParameters()[0]);
  params.addMember(XPRS_TOK_IFUN, XPRS_IFUN_LN);
  params.addMember(XPRS_TOK_OP, XPRS_OP_MULTIPLY);
  params.addMember(XPRS_TOK_IFUN, XPRS_IFUN_EXP);
  params.addMember(XPRS_TOK_EOF, 0);
  AddGlobalConstraint(params, 'E');
}

void XpressmpModelAPI::AddConstraint(const PLConstraint& plc) {
  const auto& points = plc.GetParameters().GetPLPoints();
  int x = plc.GetArguments()[0];
  int y = plc.GetResultVar();
  int size = points.x_.size();
  int starts[] = { 0, size };
  XPRESSMP_CCALL(XPRSaddpwlcons(lp(), 1, size,
    &x, &y, starts, points.x_.data(), points.y_.data()));
}

void XpressmpModelAPI::AddGlobalConstraint(int resultVar, int argumentVar, int functionId) {
  NLParams params(resultVar);
  params.addMember(XPRS_TOK_RB, 0);
  params.addMember(NLParams::var(argumentVar));
  params.addMember(NLParams::func(functionId));
  params.addMember(XPRS_TOK_EOF, 0);
  AddGlobalConstraint(params, 'E');
}

void XpressmpModelAPI::AddGlobalConstraint(const NLParams& params, char type) {
  char BUFFER[512];
  int status;

  double rhs = 0, coef = -1;
  int start = 0;
  XPRESSMP_CCALL(XPRSaddrows(lp(), 1, 1, &type, &rhs, NULL, &start, params.resultVar(), &coef));
  int rowindex = NumLinCons() - 1;
  int formulaStart[] = { 0, params.size() };
  status = XPRSnlpaddformulas(lp(), 1, &rowindex, formulaStart, true, params.types(), params.values());
  if (status) {
    XPRSgetlasterror(lp(), BUFFER);
    printf(BUFFER);
  }
}


void XpressmpModelAPI::AddConstraint(const SinhConstraint& cc) {
  // sinh(x) = 0.5*(e^x-e^-x)
  NLParams params(cc.GetResultVar());
  params.addMember(XPRS_TOK_RB, 0);
  params.addMember(XPRS_TOK_COL, cc.GetArguments()[0]);
  params.addMember(XPRS_TOK_IFUN, XPRS_IFUN_EXP);
  params.addMember(XPRS_TOK_RB, 0);
  params.addMember(XPRS_TOK_COL, cc.GetArguments()[0]);
  params.addMember(XPRS_TOK_OP, XPRS_OP_UMINUS);
  params.addMember(XPRS_TOK_IFUN, XPRS_IFUN_EXP);
  params.addMember(XPRS_TOK_OP, XPRS_OP_MINUS);
  params.addMember(XPRS_TOK_CON, 2);
  params.addMember(XPRS_TOK_OP, XPRS_OP_DIVIDE);
  params.addMember(XPRS_TOK_EOF, 0);
  AddGlobalConstraint(params, 'E');
}



void XpressmpModelAPI::AddConstraint(const CoshConstraint& cc) {
  // sinh(x) = 0.5*(e^x+e^-x)
  NLParams params(cc.GetResultVar());
  params.addMember(XPRS_TOK_RB, 0);
  params.addMember(XPRS_TOK_COL, cc.GetArguments()[0]);
  params.addMember(XPRS_TOK_IFUN, XPRS_IFUN_EXP);
  params.addMember(XPRS_TOK_RB, 0);
  params.addMember(XPRS_TOK_COL, cc.GetArguments()[0]);
  params.addMember(XPRS_TOK_OP, XPRS_OP_UMINUS);
  params.addMember(XPRS_TOK_IFUN, XPRS_IFUN_EXP);
  params.addMember(XPRS_TOK_OP, XPRS_OP_PLUS);
  params.addMember(XPRS_TOK_CON, 2);
  params.addMember(XPRS_TOK_OP, XPRS_OP_DIVIDE);
  params.addMember(XPRS_TOK_EOF, 0);
  AddGlobalConstraint(params, 'E');
}

void XpressmpModelAPI::AddConstraint(const TanhConstraint& cc) {
  int x = cc.GetArguments()[0];
  // sinh(x) = (e^x-e^-x)/(e^x+e^-x)
  NLParams params(cc.GetResultVar());

  params.addMember(XPRS_TOK_RB, 0); // e^x
  params.addMember(XPRS_TOK_COL, x);
  params.addMember(XPRS_TOK_IFUN, XPRS_IFUN_EXP);

  params.addMember(XPRS_TOK_RB, 0); // e^-x
  params.addMember(XPRS_TOK_COL, x);
  params.addMember(XPRS_TOK_OP, XPRS_OP_UMINUS);
  params.addMember(XPRS_TOK_IFUN, XPRS_IFUN_EXP);
  
  params.addMember(XPRS_TOK_OP, XPRS_OP_MINUS); // e^x-e^-x

  params.addMember(XPRS_TOK_RB, 0);   // e^x
  params.addMember(XPRS_TOK_COL, x);
  params.addMember(XPRS_TOK_IFUN, XPRS_IFUN_EXP);
  
  params.addMember(XPRS_TOK_RB, 0);   // e^-x
  params.addMember(XPRS_TOK_COL, x);
  params.addMember(XPRS_TOK_OP, XPRS_OP_UMINUS);
  params.addMember(XPRS_TOK_IFUN, XPRS_IFUN_EXP);
  params.addMember(XPRS_TOK_OP, XPRS_OP_PLUS);  // e^x +e^x

  params.addMember(XPRS_TOK_OP, XPRS_OP_DIVIDE);
  params.addMember(XPRS_TOK_EOF, 0);

  AddGlobalConstraint(params, 'E');
}

void XpressmpModelAPI::AddConstraint(const AsinhConstraint& cc) {
   // ln(x+sqrt(x^2+1))
  auto x = cc.GetArguments()[0];
  NLParams params(cc.GetResultVar());

  params.addMember(XPRS_TOK_RB, 0); // sqrt(x^2+1)

  params.addMember(XPRS_TOK_RB, 0); // x^2
  params.addMember(XPRS_TOK_COL, x);
  params.addMember(XPRS_TOK_CON, 2);
  params.addMember(XPRS_TOK_OP, XPRS_OP_EXPONENT);

  params.addMember(XPRS_TOK_CON, 1);  
  params.addMember(XPRS_TOK_OP, XPRS_OP_PLUS);
  params.addMember(XPRS_TOK_IFUN, XPRS_IFUN_SQRT);
  
  params.addMember(XPRS_TOK_COL,x);
  params.addMember(XPRS_TOK_OP, XPRS_OP_PLUS);

  params.addMember(XPRS_TOK_IFUN, XPRS_IFUN_LN);
  params.addMember(XPRS_TOK_EOF, 0);
  AddGlobalConstraint(params, 'E');
}

void XpressmpModelAPI::AddConstraint(const AcoshConstraint& cc) {
  // ln(x+sqrt(x-1)*sqrt(x+1) )

  auto x = cc.GetArguments()[0];
  NLParams params(cc.GetResultVar());

  params.addMember(XPRS_TOK_RB, 0); // ln

  params.addMember(XPRS_TOK_RB, 0); // sqrt(x-1)
  params.addMember(XPRS_TOK_COL, x);
  params.addMember(XPRS_TOK_CON, 1);
  params.addMember(XPRS_TOK_OP, XPRS_OP_MINUS);
  params.addMember(XPRS_TOK_IFUN, XPRS_IFUN_SQRT);

  params.addMember(XPRS_TOK_RB, 0); // sqrt(x+1)
  params.addMember(XPRS_TOK_COL, x);
  params.addMember(XPRS_TOK_CON, 1);
  params.addMember(XPRS_TOK_OP, XPRS_OP_PLUS);
  params.addMember(XPRS_TOK_IFUN, XPRS_IFUN_SQRT);

  params.addMember(XPRS_TOK_OP, XPRS_OP_MULTIPLY);
  
  params.addMember(XPRS_TOK_COL, x);
  params.addMember(XPRS_TOK_OP, XPRS_OP_PLUS);

  params.addMember(XPRS_TOK_IFUN, XPRS_IFUN_LN);
  params.addMember(XPRS_TOK_EOF, 0);
  AddGlobalConstraint(params, 'E');


}

void XpressmpModelAPI::AddConstraint(const AtanhConstraint& cc) {
  //y = 0.5*(ln(x+1)-ln(1-x))
  auto x = cc.GetArguments()[0];
  NLParams params(cc.GetResultVar());
  
  params.addMember(XPRS_TOK_RB, 0); // ln(x+1)
  params.addMember(XPRS_TOK_CON, 1);
  params.addMember(XPRS_TOK_COL, x);
  params.addMember(XPRS_TOK_OP, XPRS_OP_PLUS);
  params.addMember(XPRS_TOK_IFUN, XPRS_IFUN_LN);

  params.addMember(XPRS_TOK_RB, 0); // ln(1-x)
  params.addMember(XPRS_TOK_CON, 1);
  params.addMember(XPRS_TOK_COL, x);
  params.addMember(XPRS_TOK_OP, XPRS_OP_MINUS);
  params.addMember(XPRS_TOK_IFUN, XPRS_IFUN_LN);


  params.addMember(XPRS_TOK_OP, XPRS_OP_MINUS); 

  params.addMember(XPRS_TOK_CON, 0.5);
  params.addMember(XPRS_TOK_OP, XPRS_OP_MULTIPLY);

  params.addMember(XPRS_TOK_EOF, 0);
  AddGlobalConstraint(params, 'E');
}

void XpressmpModelAPI::AddConstraint(const NLConstraint& nl) {
  auto& exp = GetExpression(nl);
  exp.reverse();
  exp.addMember(XPRS_TOK_EOF, 0);
  double lhs = GetLower(nl), rhs = GetUpper(nl);
  char type;
  double range;
  double* prange = nullptr;
  if (lhs == rhs)
    type = 'E';
  else
  {
    if ((lhs != (-std::numeric_limits<double>::infinity()) && rhs != std::numeric_limits<double>::infinity())) // range
    {
      type = 'R';
      prange = &range;
      range = rhs - lhs;
    }
    else
    {
      if (lhs != (-std::numeric_limits<double>::infinity())) // GE
      {
        type = 'G';
        rhs = lhs;
      }
      else {
        type = 'L';
      }
    }
  }
  int rowindex = addLinearRow(nl, type, rhs, prange);
  char BUFFER[512];
  int status;
  int formulaStart[] = { 0, exp.size() };
  status = XPRSnlpaddformulas(lp(), 1, &rowindex, formulaStart, true, exp.types(), exp.values());
  if (status) {
    XPRSgetlasterror(lp(), BUFFER);
    printf(BUFFER);
  }
}
void XpressmpModelAPI::AddConstraint(const NLAssignEQ& neq) {
  NLParams params(GetVariable(neq));
  params.addMembers(GetExpression(neq));
  params.reverse();
  params.addMember(XPRS_TOK_EOF, 0);
  AddGlobalConstraint(params, 'E');

}
void XpressmpModelAPI::AddConstraint(const NLAssignGE& nge) {
  NLParams params(GetVariable(nge));
  params.addMembers(GetExpression(nge));
  params.reverse();
  params.addMember(XPRS_TOK_EOF, 0);

  AddGlobalConstraint(params, 'L');

}
void XpressmpModelAPI::AddConstraint(const NLAssignLE& nle) {
  NLParams params(GetVariable(nle));
  params.addMembers(GetExpression(nle));
  params.reverse();
  params.addMember(XPRS_TOK_EOF, 0);
  // Note, the sense is for a constraint of type:
  // -var + expr sense 0, so they have to be inverted
  AddGlobalConstraint(params, 'G'); 
}

template <class MPExpr> 
void XpressmpModelAPI::AppendLinAndConstTerms(Expr& exp, const MPExpr& ae) {
  double ct = GetConstTerm(ae);
  int size = GetLinSize(ae);
  // cases
  // 4   -> const(4)
  // 4+x -> PLUS | const(4) | var(y)
  // 4+x+y -> PLUS | const(4) | PLUS | var(x) | var(y)
  // x     -> var(x)
  // x+y   -> PLUS | var(x) | var(y)
  // x+y+z -> PLUS | var(x) | PLUS | var(y) | var(z)
  if ((size > 1) || ct)
    exp.addMember(NLParams::op(XPRS_OP_PLUS));
  if (ct)
  {
    exp.addMember(NLParams::constant(ct));
    // is we have to add more than one term, must put the operator now
    if(size > 1) exp.addMember(NLParams::op(XPRS_OP_PLUS));
  }

  for (int i = 0; i < size; ++i) {
    if (double coef = GetLinCoef(ae, i)) {
      if (1.0 != coef) {
        exp.addMember(NLParams::op(XPRS_OP_MULTIPLY));
        exp.addMember(NLParams::constant(coef));
        exp.addMembers(GetLinTerm(ae, i));
      }
      else {
        exp.addMembers(GetLinTerm(ae, i));
      }
    }
    if (i < size - 2)
    {
      exp.addMember(NLParams::op(XPRS_OP_PLUS));
    }
  }
}

NLParams XpressmpModelAPI::AddExpression(const NLAffineExpression& ae) {
  NLParams affine;
  AppendLinAndConstTerms(affine, ae);
  return affine;
}

NLParams XpressmpModelAPI::AddExpression(const NLQuadExpression& qe) {
  NLParams quad;
  
  if(GetLinSize(qe) > 0)
    quad.addMember(NLParams::op(XPRS_OP_PLUS));

  AppendLinAndConstTerms(quad, qe);
  int size = GetQuadSize(qe);
  if (size > 1)
    quad.addMember(NLParams::op(XPRS_OP_PLUS));
  for (int i = 0; i < GetQuadSize(qe); ++i) {
    if (double coef = GetQuadCoef(qe, i)) {
      if (1.0 != coef) {
        quad.addMember(NLParams::op(XPRS_OP_MULTIPLY));
        quad.addMember(NLParams::constant(coef));
      }
      quad.addMember(NLParams::op(XPRS_OP_MULTIPLY));
      quad.addMembers(GetQuadTerm1(qe, i));
      quad.addMembers(GetQuadTerm2(qe, i));
      if (i < size - 2)
        quad.addMember(NLParams::op(XPRS_OP_PLUS));
    }
  }
  return quad;
}
NLParams XpressmpModelAPI::AddExpression(const SinExpression& e) {
  return CreateExpressionOneArg(e, XPRS_IFUN_SIN);
}
NLParams XpressmpModelAPI::AddExpression(const CosExpression& e) {
  return CreateExpressionOneArg(e, XPRS_IFUN_COS);
}
NLParams XpressmpModelAPI::AddExpression(const TanExpression& e) {
  return CreateExpressionOneArg(e, XPRS_IFUN_TAN);
}
NLParams XpressmpModelAPI::AddExpression(const AsinExpression& e) {
  return CreateExpressionOneArg(e, XPRS_IFUN_ARCSIN);
}
NLParams XpressmpModelAPI::AddExpression(const AcosExpression& e) {
  return CreateExpressionOneArg(e, XPRS_IFUN_ARCCOS);
}
NLParams XpressmpModelAPI::AddExpression(const AtanExpression& e) {
  return CreateExpressionOneArg(e, XPRS_IFUN_ARCTAN);
}

NLParams XpressmpModelAPI::AddExpression(const LogExpression& e ) {
  return CreateExpressionOneArg(e, XPRS_IFUN_LN);
}

NLParams XpressmpModelAPI::AddExpression(const ExpExpression& e) {
  return CreateExpressionOneArg(e, XPRS_IFUN_EXP);
}

NLParams XpressmpModelAPI::AddExpression(const PowConstExpExpression& e) {
  double exponent = GetParameter(e, 0);
  
  if (exponent == 0.5)
    return CreateExpressionOneArg(e, XPRS_IFUN_SQRT);
  NLParams exp;
  exp.addMember(NLParams::op(XPRS_OP_EXPONENT));
  exp.addMember(NLParams::constant(exponent));
  exp.addMembers(GetArgExpression(e, 0));
  return exp;
}
NLParams XpressmpModelAPI::AddExpression(const DivExpression& e) {
  NLParams exp;
  exp.addMember(NLParams::op(XPRS_OP_DIVIDE));
  exp.addMembers(GetArgExpression(e, 1));
  exp.addMembers(GetArgExpression(e, 0));
  return exp;
}

NLParams XpressmpModelAPI::AddExpression(const MinExpression& e) {
  return CreateExpressionNArgs(e, XPRS_IFUN_MIN);
}
NLParams XpressmpModelAPI::AddExpression(const MaxExpression& e) {
  return CreateExpressionNArgs(e, XPRS_IFUN_MAX);
}
NLParams XpressmpModelAPI::AddExpression(const AbsExpression& e) {
  return CreateExpressionOneArg(e, XPRS_IFUN_ABS);
}


NLParams 
XpressmpModelAPI::AddExpression(const LogAExpression& e) {
  NLParams exp;
  exp.addMember(NLParams::op(XPRS_OP_DIVIDE));
  exp.addMember(NLParams::func(XPRS_IFUN_LN));
  auto ex = GetArgExpression(e, 0);
  exp.addMembers(ex);
  exp.addMember(XPRS_TOK_RB, 0);
  exp.addMember(NLParams::func(XPRS_IFUN_LN));
  auto par = GetParameter(e, 0);
  exp.addMember(NLParams::constant(par));
  exp.addMember(XPRS_TOK_RB, 0);
  return exp;
}
void XpressmpModelAPI::FinishProblemModificationPhase() {
}


} // namespace mp

