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

std::string sanitizeName(std::string n) {
  // Xpress does not like square brackets or spaces in variable names
  std::replace(n.begin(), n.end(), '[', '(');
  std::replace(n.begin(), n.end(), ']', ')');
  std::replace(n.begin(), n.end(), ' ', '_');
  myreplace(n,"\'", "-");
  myreplace(n, "\"", "--");
  return n;
}
void XpressmpModelAPI::AddVariables(const VarArrayDef& v) {
  std::vector<double> objs(v.size(), 0);
  std::vector<int> iii(v.size(), 0);
  XPRESSMP_CCALL(XPRSaddcols(lp(), v.size(), 0, objs.data(), iii.data(), NULL,
    NULL, v.plb(), v.pub()));

  if (v.pnames() != NULL)
  {
    fmt::MemoryWriter w;
    for (int i = 0; i < v.size(); ++i)
      w << sanitizeName(v.pnames()[i]) << '\0';
    XPRESSMP_CCALL(XPRSaddnames(lp(), 2, w.c_str(), 0, v.size() - 1));
  }
  // All variables are continuous by default, set the integer ones
  std::vector<int> intIndices;
  for (auto i = 0; i < v.size(); i++)
    if (v.ptype()[i] == var::Type::INTEGER)
      intIndices.push_back(i);
  get_other()->numIntVars(intIndices.size());
  if (get_other()->numIntVars() > 0)
  {
    std::vector<char> types(intIndices.size(), 'I');
    XPRESSMP_CCALL(XPRSchgcoltype(lp(), intIndices.size(),
      intIndices.data(), types.data()));
  }
}

void XpressmpModelAPI::SetLinearObjective( int iobj, const LinearObjective& lo ) {
  if (iobj<1) {
    if (lo.obj_sense() == obj::Type::MAX)
      XPRESSMP_CCALL(XPRSchgobjsense(lp(), XPRS_OBJ_MAXIMIZE));
    XPRESSMP_CCALL(XPRSchgobj(lp(), lo.num_terms(), lo.vars().data(), lo.coefs().data()));
    
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
    fmt::format("Setting first quadratic objective\n");
    SetLinearObjective(iobj, qo);                         // add the linear part
    const auto& qt = qo.GetQPTerms();
    std::vector<double> coeffs(qt.coefs());
    for (std::size_t i = 0; i < qt.size(); i++)
      if (qt.pvars1()[i] == qt.pvars2()[i]) coeffs[i] *= 2;

    fmt::format("Quadratic part is made of {} terms\n", qt.size());
    XPRESSMP_CCALL(XPRSchgmqobj(lp(), qt.size(),
      (int*)qt.pvars1(), (int*)qt.pvars2(),
      (double*)coeffs.data()));
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
char type[] ={'I'};\
XPRESSMP_CCALL(XPRSchgcoltype(lp(), 1, colindex, type));\
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

void XpressmpModelAPI::AddLinTerms(XPRSprob lp, const LinTerms& lt, double rhsc,  const char typec) {
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
    XPRESSMP_CCALL(XPRSchgqrowcoeff(lp(), row, qt.var1(i), qt.var2(i), qt.coef(i)));

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

template <class Args, class Params, class NumOrLogic, class Id> void XpressmpModelAPI::addGenCon(
  const CustomFunctionalConstraint<Args, Params, NumOrLogic, Id>& c, int xpressConType)
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

void XpressmpModelAPI::FinishProblemModificationPhase() {
}


} // namespace mp
