#include "gurobimodelapi.h"

namespace mp {

const char* GurobiModelAPI::GetTypeName()
  { return "GurobiModelAPI"; }


//////////////////////////////////////////////////////////////////////////
/////////////////////////// Modeling interface ///////////////////////////
//////////////////////////////////////////////////////////////////////////
void GurobiModelAPI::InitProblemModificationPhase(const FlatModelInfo*) { }

///////////////////////////////////////////////////////
void GurobiModelAPI::FinishProblemModificationPhase() {
  // Update before adding statuses etc
  GRB_CALL( GRBupdatemodel(model()) );
}


void GurobiModelAPI::AddVariables(const VarArrayDef& v) {
  std::vector<char> vtypes(v.size());
  for (size_t i=v.size(); i--; )
    vtypes[i] = var::Type::CONTINUOUS==v.ptype()[i] ?
          GRB_CONTINUOUS : GRB_INTEGER;
  GRB_CALL( GRBaddvars(model(), (int)v.size(), 0,
                       NULL, NULL, NULL, NULL, // placeholders, no matrix here
                       (double*)v.plb(), (double*)v.pub(),
                       vtypes.data(), (char**)v.pnames()) );
}
void ThrowForUnsupportedObjectives()
{
  throw std::runtime_error(format_error("Multiple quadratic objectives are not supported natively by Gurobi; "
                                        "enable multi-objective emulator by setting the option obj:multi=2"));
}

void GurobiModelAPI::SetLinearObjective( int iobj, const LinearObjective& lo ) {
  if (1>iobj) {
    GrbSetIntAttr( GRB_INT_ATTR_MODELSENSE,
                  obj::Type::MAX==lo.obj_sense() ? GRB_MAXIMIZE : GRB_MINIMIZE);
    NoteGurobiMainObjSense(lo.obj_sense());
    if (obj_ind_save_.size()) {
      std::vector<double> obj_coef_0(obj_ind_save_.size(), 0.0);
      GrbSetDblAttrList( GRB_DBL_ATTR_OBJ, obj_ind_save_, obj_coef_0 );
    }
    GrbSetDblAttrList( GRB_DBL_ATTR_OBJ, lo.vars(), lo.coefs() );
    obj_ind_save_ = lo.vars();
  } else {
    if (has_quadratic_obj_) ThrowForUnsupportedObjectives();

    GRB_CALL( GRBsetobjectiven(model(), iobj, 0,           // default priority 0
                               /// Gurobi allows opposite sense by weight sign
                               lo.obj_sense()==GetGurobiMainObjSense() ? 1.0 : -1.0,
                               0.0, 0.0, lo.name(),
                               0.0, lo.num_terms(),
                               (int*)lo.vars().data(), (double*)lo.coefs().data()) );
  }
}

void GurobiModelAPI::SetQuadraticObjective(int iobj, const QuadraticObjective &qo) {
  has_quadratic_obj_ = true;
  if (1>iobj) {
    GRB_CALL( GRBdelq(model()) );                         // delete current QP terms
    SetLinearObjective(iobj, qo);                         // add the linear part
    const auto& qt = qo.GetQPTerms();
    GRB_CALL( GRBaddqpterms(model(), qt.size(),
                                (int*)qt.pvars1(), (int*)qt.pvars2(),
                            (double*)qt.pcoefs()) );
  } else {
    ThrowForUnsupportedObjectives();
  }
}

void GurobiModelAPI::NoteGurobiMainObjSense(obj::Type s) { main_obj_sense_ = s; }

obj::Type GurobiModelAPI::GetGurobiMainObjSense() const { return main_obj_sense_; }


void GurobiModelAPI::AddConstraint( const LinConLE& lc ) {
  GRB_CALL( GRBaddconstr(model(), lc.size(),
                         (int*)lc.pvars(), (double*)lc.pcoefs(),
                         GRB_LESS_EQUAL, lc.rhs(), lc.name()) );
}
void GurobiModelAPI::AddConstraint( const LinConEQ& lc ) {
  GRB_CALL( GRBaddconstr(model(), lc.size(),
                         (int*)lc.pvars(), (double*)lc.pcoefs(),
                         GRB_EQUAL, lc.rhs(), lc.name()) );
}
void GurobiModelAPI::AddConstraint( const LinConGE& lc ) {
  GRB_CALL( GRBaddconstr(model(), lc.size(),
                         (int*)lc.pvars(), (double*)lc.pcoefs(),
                         GRB_GREATER_EQUAL, lc.rhs(), lc.name()) );
}

void GurobiModelAPI::AddConstraint( const QuadConLE& qc ) {
  const auto& lt = qc.GetLinTerms();
  const auto& qt = qc.GetQPTerms();
  GRB_CALL( GRBaddqconstr(model(), lt.size(), (int*)lt.pvars(), (double*)lt.pcoefs(),
                          qt.size(), (int*)qt.pvars1(), (int*)qt.pvars2(),
                          (double*)qt.pcoefs(), GRB_LESS_EQUAL, qc.rhs(), qc.name()) );
}

void GurobiModelAPI::AddConstraint( const QuadConEQ& qc ) {
  const auto& lt = qc.GetLinTerms();
  const auto& qt = qc.GetQPTerms();
  GRB_CALL( GRBaddqconstr(model(), lt.size(), (int*)lt.pvars(), (double*)lt.pcoefs(),
                          qt.size(), (int*)qt.pvars1(), (int*)qt.pvars2(),
                          (double*)qt.pcoefs(), GRB_EQUAL, qc.rhs(), qc.name()) );
}

void GurobiModelAPI::AddConstraint( const QuadConGE& qc ) {
  const auto& lt = qc.GetLinTerms();
  const auto& qt = qc.GetQPTerms();
  GRB_CALL( GRBaddqconstr(model(), lt.size(), (int*)lt.pvars(), (double*)lt.pcoefs(),
                          qt.size(), (int*)qt.pvars1(), (int*)qt.pvars2(),
                          (double*)qt.pcoefs(), GRB_GREATER_EQUAL, qc.rhs(), qc.name()) );
}


void GurobiModelAPI::AddConstraint(const MaxConstraint &mc)  {
  const auto& args = mc.GetArguments();
  GRB_CALL( GRBaddgenconstrMax(model(), mc.name(),
                               mc.GetResultVar(),
                               (int)args.size(), args.data(),
                               MinusInfinity()) );
}

void GurobiModelAPI::AddConstraint(const MinConstraint &mc)  {
  const auto& args = mc.GetArguments();
  GRB_CALL( GRBaddgenconstrMin(model(), mc.name(),
                               mc.GetResultVar(),
                               (int)args.size(), args.data(),
                               Infinity()) );
}

void GurobiModelAPI::AddConstraint(const AbsConstraint &absc)  {
  const auto& args = absc.GetArguments();
  GRB_CALL( GRBaddgenconstrAbs(model(), absc.name(),
                               absc.GetResultVar(), args[0]) );
}

void GurobiModelAPI::AddConstraint(const AndConstraint &cc)  {
  const auto& args = cc.GetArguments();
  GRB_CALL( GRBaddgenconstrAnd(model(), cc.name(),
                               cc.GetResultVar(),
                               (int)args.size(), args.data()) );
}

void GurobiModelAPI::AddConstraint(const OrConstraint &dc)  {
  const auto& args = dc.GetArguments();
  GRB_CALL( GRBaddgenconstrOr(model(), dc.name(),
                               dc.GetResultVar(),
                               (int)args.size(), args.data()) );
}

void GurobiModelAPI::AddConstraint(const IndicatorConstraintLinLE &ic)  {
  GRB_CALL( GRBaddgenconstrIndicator(model(), ic.name(),
                               ic.get_binary_var(), ic.get_binary_value(),
                                     (int)ic.get_constraint().size(),
                               ic.get_constraint().pvars(),
                                     ic.get_constraint().pcoefs(),
                                     GRB_LESS_EQUAL,
                                     ic.get_constraint().rhs() ) );
}
void GurobiModelAPI::AddConstraint(const IndicatorConstraintLinEQ &ic)  {
  GRB_CALL( GRBaddgenconstrIndicator(model(), ic.name(),
                                     ic.get_binary_var(), ic.get_binary_value(),
                                           (int)ic.get_constraint().size(),
                                     ic.get_constraint().pvars(),
                                           ic.get_constraint().pcoefs(),
                                           GRB_EQUAL,
                                           ic.get_constraint().rhs() ) );
}
void GurobiModelAPI::AddConstraint(const IndicatorConstraintLinGE &ic)  {
  GRB_CALL( GRBaddgenconstrIndicator(model(), ic.name(),
                               ic.get_binary_var(), ic.get_binary_value(),
                                     (int)ic.get_constraint().size(),
                               ic.get_constraint().pvars(),
                                     ic.get_constraint().pcoefs(),
                                     GRB_GREATER_EQUAL,
                                     ic.get_constraint().rhs() ) );
}

//////////////////// General constraints /////////////////////
void GurobiModelAPI::AddConstraint(const SOS1Constraint &sos)  {
  int type = GRB_SOS_TYPE1;
  int beg = 0;
  GRB_CALL( GRBaddsos(model(), 1, sos.size(), &type, &beg,
              (int*)sos.get_vars().data(),
                      (double*)sos.get_weights().data()) );
}

void GurobiModelAPI::AddConstraint(const SOS2Constraint &sos)  {
  int type = GRB_SOS_TYPE2;
  int beg = 0;
  GRB_CALL( GRBaddsos(model(), 1, sos.size(), &type, &beg,
              (int*)sos.get_vars().data(),
                      (double*)sos.get_weights().data()) );
}

void GurobiModelAPI::AddConstraint(const ExpConstraint &cc)  {
  GRB_CALL( GRBaddgenconstrExp(model(), cc.name(),
              cc.GetArguments()[0], cc.GetResultVar(), "") );
}

void GurobiModelAPI::AddConstraint(const ExpAConstraint &cc)  {
  GRB_CALL( GRBaddgenconstrExpA(model(), cc.name(),
              cc.GetArguments()[0], cc.GetResultVar(), cc.GetParameters()[0], "") );
}

void GurobiModelAPI::AddConstraint(const LogConstraint &cc)  {
  GRB_CALL( GRBaddgenconstrLog(model(), cc.name(),
              cc.GetArguments()[0], cc.GetResultVar(), "") );
}

void GurobiModelAPI::AddConstraint(const LogAConstraint &cc)  {
  GRB_CALL( GRBaddgenconstrLogA(model(), cc.name(),
              cc.GetArguments()[0], cc.GetResultVar(), cc.GetParameters()[0], "") );
}

void GurobiModelAPI::AddConstraint(const PowConstExpConstraint &cc)  {
  GRB_CALL( GRBaddgenconstrPow(model(), cc.name(),
              cc.GetArguments()[0], cc.GetResultVar(), cc.GetParameters()[0], "") );
}

void GurobiModelAPI::AddConstraint(const SinConstraint &cc)  {
  GRB_CALL( GRBaddgenconstrSin(model(), cc.name(),
              cc.GetArguments()[0], cc.GetResultVar(), "") );
}

void GurobiModelAPI::AddConstraint(const CosConstraint &cc)  {
  GRB_CALL( GRBaddgenconstrCos(model(), cc.name(),
              cc.GetArguments()[0], cc.GetResultVar(), "") );
}

void GurobiModelAPI::AddConstraint(const TanConstraint &cc)  {
  GRB_CALL( GRBaddgenconstrTan(model(), cc.name(),
              cc.GetArguments()[0], cc.GetResultVar(), "") );
}

void GurobiModelAPI::AddConstraint(const PLConstraint& plc) {
  const auto& plp = plc.GetParameters().GetPLPoints();
  GRB_CALL( GRBaddgenconstrPWL(model(), plc.name(),
              plc.GetArguments()[0], plc.GetResultVar(),
              plp.x_.size(), plp.x_.data(), plp.y_.data()) );
}

void GurobiModelAPI::InitCustomOptions() { }


/////////////////////////////// EXPRESSIONS ////////////////////////////////////
#ifdef GRB_OPCODE_CONSTANT

void GurobiModelAPI::AddConstraint(const NLAssignEQ& nla) {
  const auto& frm = GetFormula(GetExpression(nla));
  GRB_CALL( GRBaddgenconstrNL(
      model(), GetName(nla), GetVariable(nla), frm.size(),
      frm.opcodes(), frm.data(), frm.parents()) );
}
void GurobiModelAPI::AddConstraint(const NLAssignLE& nla) {
  const auto& frm = GetFormula(GetExpression(nla));
  GRB_CALL( GRBaddgenconstrNL(
      model(), GetName(nla), GetVariable(nla), frm.size(),
      frm.opcodes(), frm.data(), frm.parents()) );
}
void GurobiModelAPI::AddConstraint(const NLAssignGE& nla) {
  const auto& frm = GetFormula(GetExpression(nla));
  GRB_CALL( GRBaddgenconstrNL(
      model(), GetName(nla), GetVariable(nla), frm.size(),
      frm.opcodes(), frm.data(), frm.parents()) );
}

GRB_Expr GurobiModelAPI::AddExpression(const NLAffineExpression& nla) {
  int new_pos = (int)formulas_.size();   // we may add more in GetArg...
  formulas_.push_back( StartFormula(GRB_OPCODE_PLUS) );
  AppendLinAndConstTerms(formulas_[new_pos], nla);
  return MakeExpr(new_pos);
}

GRB_Expr GurobiModelAPI::AddExpression(const NLQuadExpression& nlq) {
  int new_pos = (int)formulas_.size();   // we may add more in GetArg...
  formulas_.push_back( StartFormula(GRB_OPCODE_PLUS) );
  AppendLinAndConstTerms(formulas_[new_pos], nlq);
  AppendQuadTerms(formulas_[new_pos], nlq);
  return MakeExpr(new_pos);
}

template <class MPExpr>
void GurobiModelAPI::AppendLinAndConstTerms(
    GurobiModelAPI::Formula& ff, const MPExpr& nla) {
  if (double ct = GetConstTerm(nla))
    AppendArgument(ff, MakeConstantExpr(ct));
  for (int i=0; i<GetLinSize(nla); ++i) {
    if (double coef = GetLinCoef(nla, i)) {
      if (1.0 != coef) {
        auto f1 = StartFormula(GRB_OPCODE_MULTIPLY);
        AppendArgument(f1, MakeConstantExpr(coef));
        AppendArgument(f1, GetLinTerm(nla, i));
        ff.Append(f1);
      }  else {
        AppendArgument(ff, GetLinTerm(nla, i));
      }
    }
  }
}

template <class MPExpr>
void GurobiModelAPI::AppendQuadTerms(
    GurobiModelAPI::Formula& ff, const MPExpr& nlq) {
  for (int i=0; i<GetQuadSize(nlq); ++i) {
    if (double coef = GetQuadCoef(nlq, i)) {
      auto f1 = StartFormula(GRB_OPCODE_MULTIPLY);
      AppendArgument(f1, GetQuadTerm1(nlq, i));
      AppendArgument(f1, GetQuadTerm2(nlq, i));
      if (1.0 != coef) {
        AppendArgument(f1, MakeConstantExpr(coef));
      }
      ff.Append(f1);
    }
  }
}


GRB_Expr GurobiModelAPI::AddExpression(const ExpExpression& e) {
  return CreateFormula(e, GRB_OPCODE_EXP);
}
GRB_Expr GurobiModelAPI::AddExpression(const ExpAExpression& e) {
  int new_pos = (int)formulas_.size();   // we may add more in GetArg...
  formulas_.push_back( StartFormula(GRB_OPCODE_EXP) );
  auto f1 = StartFormula(GRB_OPCODE_MULTIPLY);
  AppendArgument(f1,
                 MakeConstantExpr(std::log(GetParameter(e, 0))));
  AppendArgument(f1, GetArgExpression(e, 0));
  formulas_[new_pos].Append(f1);
  return MakeExpr(new_pos);
}
GRB_Expr GurobiModelAPI::AddExpression(const LogExpression& e) {
  return CreateFormula(e, GRB_OPCODE_LOG);
}
GRB_Expr GurobiModelAPI::AddExpression(const LogAExpression& e) {
  int new_pos = (int)formulas_.size();   // we may add more in GetArg...
  formulas_.push_back( StartFormula(GRB_OPCODE_MULTIPLY) );
  auto a = GetParameter(e, 0);
  if (1.0==a)
    MP_RAISE("Logarithm with base 1");
  AppendArgument(formulas_[new_pos],
                 MakeConstantExpr(1.0 / std::log(a)));
  auto f1 = StartFormula(GRB_OPCODE_LOG);
  AppendArgument(f1, GetArgExpression(e, 0));
  formulas_[new_pos].Append(f1);
  return MakeExpr(new_pos);
}

GRB_Expr GurobiModelAPI::AddExpression(const PowConstExpExpression& e) {
  return CreateFormula(e, GRB_OPCODE_POW);
}
GRB_Expr GurobiModelAPI::AddExpression(const SignpowConstExpExpression& e) {
  return CreateFormula(e, GRB_OPCODE_SIGNPOW);
}
GRB_Expr GurobiModelAPI::AddExpression(const SinExpression& e) {
  return CreateFormula(e, GRB_OPCODE_SIN);
}
GRB_Expr GurobiModelAPI::AddExpression(const CosExpression& e) {
  return CreateFormula(e, GRB_OPCODE_COS);
}
GRB_Expr GurobiModelAPI::AddExpression(const TanExpression &e)
{
  return CreateFormula(e, GRB_OPCODE_TAN);
}

#ifdef GRB_OPCODE_TANH
GRB_Expr GurobiModelAPI::AddExpression(const TanhExpression &e)
{
  return CreateFormula(e, GRB_OPCODE_TANH);
}
#endif

GRB_Expr GurobiModelAPI::AddExpression(const DivExpression& e) {
  return CreateFormula(e, GRB_OPCODE_DIVIDE);
}


template <class MPExpr>
GRB_Expr GurobiModelAPI::CreateFormula(
    const MPExpr& mpexpr, int opcode) {
  int new_pos = (int)formulas_.size();   // we may add more in GetArg...
  formulas_.push_back( StartFormula(opcode) );
  for (int i=0; i<GetNumArguments(mpexpr); ++i)
    AppendArgument(formulas_[new_pos], GetArgExpression(mpexpr, i));
  for (int i=0; i<GetNumParameters(mpexpr); ++i)
    AppendArgument(formulas_[new_pos],
                   MakeConstantExpr(GetParameter(mpexpr, i)));

  return MakeExpr(new_pos);
}

void GurobiModelAPI::AppendArgument(Formula& frm, Expr expr) {
  frm.Append(GetFormula(expr));
}

const GurobiModelAPI::Formula&
GurobiModelAPI::GetFormula(GRB_Expr expr) const {
  if (expr.IsExpr())
    return formulas_.at(expr.GetExprIndex());
  if (expr.IsVar()) {
    return
        formula_tmp_ = {GRB_OPCODE_VARIABLE, (double)expr.GetVarIndex()};
  }
  assert(expr.IsConst());
  return
      formula_tmp_ = {GRB_OPCODE_CONSTANT, expr.GetConst()};
}

void GurobiModelAPI::Formula::Append(const Formula& frm) {
  assert(size());
  int sz_last = size();
  assert(frm.size());
  assert(-1==frm.parent_.front());
  assert(is_length_ok());
  if (capacity() < size()+frm.size()) {  // manually 2x size
    size_t sz_new = std::max(size()*2, size()+frm.size());
    opcode_.reserve(sz_new);
    data_.reserve(sz_new);
    parent_.reserve(sz_new);
  }
  opcode_.insert(opcode_.end(), frm.opcode_.begin(), frm.opcode_.end());
  data_.insert(data_.end(), frm.data_.begin(), frm.data_.end());
  parent_.insert(parent_.end(), frm.parent_.begin(), frm.parent_.end());
  parent_[sz_last] = 0;           // make it child of the current root
  for (int i=sz_last+1; i<size(); ++i)
    parent_[i] += sz_last;
}

#endif  // GRB_OPCODE_CONSTANT

} // namespace mp
