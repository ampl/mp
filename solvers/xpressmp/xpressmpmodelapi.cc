#include "xpressmpmodelapi.h"

#include "mp/model-mgr-with-std-pb.h"
#include "mp/flat/redef/MIP/converter_mip.h"
#include "mp/flat/model_api_connect.h"

namespace mp {

/// Defining the function in ...modelapi.cc
/// for recompilation speed
std::unique_ptr<BasicModelManager>
CreateXpressmpModelMgr(XpressmpCommon& cc, Env& e,
                     pre::BasicValuePresolver*& pPre) {
  return CreateModelMgrWithFlatConverter<
      XpressmpModelAPI, MIPFlatConverter >(cc, e, pPre);
}


void XpressmpModelAPI::InitProblemModificationPhase(
    const FlatModelInfo*) {
}

void XpressmpModelAPI::AddVariables(const VarArrayDef& v) {
  std::vector<double> objs(v.size(), 0);
  std::vector<int> iii(v.size(), 0);
  XPRESSMP_CCALL(XPRSaddcols(lp(), v.size(), 0, objs.data(), iii.data(), NULL,
    NULL, v.plb(), v.pub()));
  if(v.pnames() != NULL)
    XPRESSMP_CCALL(XPRSaddnames(lp(), 2, (const char*)v.pnames(), 0, v.size()-1));
  // All variables are continuous by default, set the integer ones
  std::vector<int> intIndices;
  for (auto i = 0; i < v.size(); i++)
    if (v.ptype()[i] == var::Type::INTEGER)
      intIndices.push_back(i);
  get_other()->numIntVars(intIndices.size());
  if (numIntVars() > 0)
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
    throw std::runtime_error("Multiple objectives not supported");
  }
}


void XpressmpModelAPI::SetQuadraticObjective(int iobj, const QuadraticObjective& qo) {
  if (1 > iobj) {
    fmt::format("Setting first quadratic objective\n");
    SetLinearObjective(iobj, qo);                         // add the linear part
    const auto& qt = qo.GetQPTerms();
    fmt::format("Quadratic part is made of {} terms\n", qt.size());
    XPRESSMP_CCALL(XPRSchgmqobj(lp(), qt.size(),
      (int*)qt.pvars1(), (int*)qt.pvars2(),
      (double*)qt.pcoefs()));
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

#define _addIndicator_mp AddConstraint(ic.get_constraint());\
int rowindex[] = { NumLinCons() - 1 };\
int colindex[] = { ic.get_binary_var() };\
int complement[] = { ic.get_binary_value() };\
XPRESSMP_CCALL(XPRSsetindicators(lp(), 1, rowindex, colindex, complement));\

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


void XpressmpModelAPI::FinishProblemModificationPhase() {
}


} // namespace mp
