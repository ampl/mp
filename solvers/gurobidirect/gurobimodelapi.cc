#include "gurobimodelapi.h"
#include "mp/model-mgr-with-std-pb.h"
#include "mp/flat/expr_flattener.h"
#include "mp/flat/MIP/mp2mip.h"

namespace mp {

std::unique_ptr<BasicModelManager>
CreateGurobiModelMgr(GurobiCommon& gc, Env& e,
                     pre::BasicPresolver*& pPre) {
  using FlatCvtType = Interface<MIPFlatConverter, GurobiModelAPI>;
  using ConverterType = mp::ExprFlattenerImpl<
    mp::ExprFlattener, mp::Problem, FlatCvtType>;
  auto pcvt = new ConverterType(e);
  auto res = CreateModelManagerWithStdBuilder(
        std::unique_ptr<BasicConverter<mp::Problem> >{ pcvt } );
  pcvt->GetFlatCvt().GetBasicBackend().set_other_gurobi(&gc);
  gc.set_other_gurobi(
        &pcvt->GetFlatCvt().GetBasicBackend());
  pPre = &pcvt->GetFlatCvt().GetPresolver();
  return res;
}


const char* GurobiModelAPI::GetBackendName()
  { return "GurobiModelAPI"; }


//////////////////////////////////////////////////////////////////////////
/////////////////////////// Modeling interface ///////////////////////////
//////////////////////////////////////////////////////////////////////////
void GurobiModelAPI::InitProblemModificationPhase() { }

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
                                          vtypes.data(), NULL) );
}

void GurobiModelAPI::SetLinearObjective( int iobj, const LinearObjective& lo ) {
  if (1>iobj) {
    GrbSetIntAttr( GRB_INT_ATTR_MODELSENSE,
                  obj::Type::MAX==lo.obj_sense() ? GRB_MAXIMIZE : GRB_MINIMIZE);
    NoteGurobiMainObjSense(lo.obj_sense());
    GrbSetDblAttrList( GRB_DBL_ATTR_OBJ, lo.vars(), lo.coefs() );
  } else {
    GRB_CALL( GRBsetobjectiven(model(), iobj, 0,           // default priority 0
                               /// Gurobi allows opposite sense by weight sign
                               lo.obj_sense()==GetGurobiMainObjSense() ? 1.0 : -1.0,
                               0.0, 0.0, nullptr,
                               0.0, lo.num_terms(),
                               (int*)lo.vars().data(), (double*)lo.coefs().data()) );
  }
}

void GurobiModelAPI::SetQuadraticObjective(int iobj, const QuadraticObjective &qo) {
  if (1>iobj) {
    SetLinearObjective(iobj, qo);                         // add the linear part
    const auto& qt = qo.GetQPTerms();
    GRB_CALL( GRBaddqpterms(model(), qt.size(),
                                (int*)qt.pvars1(), (int*)qt.pvars2(),
                            (double*)qt.pcoefs()) );
  } else {
    throw std::runtime_error("Multiple quadratic objectives not supported");
  }
}

void GurobiModelAPI::NoteGurobiMainObjSense(obj::Type s) { main_obj_sense_ = s; }

obj::Type GurobiModelAPI::GetGurobiMainObjSense() const { return main_obj_sense_; }


void GurobiModelAPI::AddConstraint( const LinConLE& lc ) {
  GRB_CALL( GRBaddconstr(model(), lc.size(),
                         (int*)lc.pvars(), (double*)lc.pcoefs(),
                         GRB_LESS_EQUAL, lc.rhs(), NULL) );
}
void GurobiModelAPI::AddConstraint( const LinConEQ& lc ) {
  GRB_CALL( GRBaddconstr(model(), lc.size(),
                         (int*)lc.pvars(), (double*)lc.pcoefs(),
                         GRB_EQUAL, lc.rhs(), NULL) );
}
void GurobiModelAPI::AddConstraint( const LinConGE& lc ) {
  GRB_CALL( GRBaddconstr(model(), lc.size(),
                         (int*)lc.pvars(), (double*)lc.pcoefs(),
                         GRB_GREATER_EQUAL, lc.rhs(), NULL) );
}

void GurobiModelAPI::AddConstraint( const QuadraticConstraint& qc ) {
  const auto& qt = qc.GetQPTerms();
  if (qc.lb()==qc.ub())
    GRB_CALL( GRBaddqconstr(model(), qc.size(), (int*)qc.pvars(), (double*)qc.pcoefs(),
                            qt.size(), (int*)qt.pvars1(), (int*)qt.pvars2(),
                            (double*)qt.pcoefs(), GRB_EQUAL, qc.lb(), NULL) );
  else {            // Let solver deal with lb>~ub etc.
    if (qc.lb()>MinusInfinity()) {
      GRB_CALL( GRBaddqconstr(model(), qc.size(), (int*)qc.pvars(), (double*)qc.pcoefs(),
                              qt.size(), (int*)qt.pvars1(), (int*)qt.pvars2(),
                              (double*)qt.pcoefs(), GRB_GREATER_EQUAL, qc.lb(), NULL) );
    }
    if (qc.ub()<Infinity()) {
      GRB_CALL( GRBaddqconstr(model(), qc.size(), (int*)qc.pvars(), (double*)qc.pcoefs(),
                              qt.size(), (int*)qt.pvars1(), (int*)qt.pvars2(),
                              (double*)qt.pcoefs(), GRB_LESS_EQUAL, qc.ub(), NULL) );
    }
  }
}

void GurobiModelAPI::AddConstraint(const MaximumConstraint &mc)  {
  const auto& args = mc.GetArguments();
  GRB_CALL( GRBaddgenconstrMax(model(), NULL,
                               mc.GetResultVar(),
                               (int)args.size(), args.data(),
                               MinusInfinity()) );
}

void GurobiModelAPI::AddConstraint(const MinimumConstraint &mc)  {
  const auto& args = mc.GetArguments();
  GRB_CALL( GRBaddgenconstrMin(model(), NULL,
                               mc.GetResultVar(),
                               (int)args.size(), args.data(),
                               Infinity()) );
}

void GurobiModelAPI::AddConstraint(const AbsConstraint &absc)  {
  const auto& args = absc.GetArguments();
  GRB_CALL( GRBaddgenconstrAbs(model(), NULL,
                               absc.GetResultVar(), args[0]) );
}

void GurobiModelAPI::AddConstraint(const ConjunctionConstraint &cc)  {
  const auto& args = cc.GetArguments();
  GRB_CALL( GRBaddgenconstrAnd(model(), NULL,
                               cc.GetResultVar(),
                               (int)args.size(), args.data()) );
}

void GurobiModelAPI::AddConstraint(const DisjunctionConstraint &dc)  {
  const auto& args = dc.GetArguments();
  GRB_CALL( GRBaddgenconstrOr(model(), NULL,
                               dc.GetResultVar(),
                               (int)args.size(), args.data()) );
}

void GurobiModelAPI::AddConstraint(const IndicatorConstraintLinLE &ic)  {
  GRB_CALL( GRBaddgenconstrIndicator(model(), NULL,
                               ic.get_binary_var(), ic.get_binary_value(),
                                     (int)ic.get_constraint().size(),
                               ic.get_constraint().pvars(),
                                     ic.get_constraint().pcoefs(),
                                     GRB_LESS_EQUAL,
                                     ic.get_constraint().rhs() ) );
}
void GurobiModelAPI::AddConstraint(const IndicatorConstraintLinEQ &ic)  {
  GRB_CALL( GRBaddgenconstrIndicator(model(), NULL,
                                     ic.get_binary_var(), ic.get_binary_value(),
                                           (int)ic.get_constraint().size(),
                                     ic.get_constraint().pvars(),
                                           ic.get_constraint().pcoefs(),
                                           GRB_EQUAL,
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
  GRB_CALL( GRBaddgenconstrExp(model(), NULL,
              cc.GetArguments()[0], cc.GetResultVar(), "") );
}

void GurobiModelAPI::AddConstraint(const ExpAConstraint &cc)  {
  GRB_CALL( GRBaddgenconstrExpA(model(), NULL,
              cc.GetArguments()[0], cc.GetResultVar(), cc.GetParameters()[0], "") );
}

void GurobiModelAPI::AddConstraint(const LogConstraint &cc)  {
  GRB_CALL( GRBaddgenconstrLog(model(), NULL,
              cc.GetArguments()[0], cc.GetResultVar(), "") );
}

void GurobiModelAPI::AddConstraint(const LogAConstraint &cc)  {
  GRB_CALL( GRBaddgenconstrLogA(model(), NULL,
              cc.GetArguments()[0], cc.GetResultVar(), cc.GetParameters()[0], "") );
}

void GurobiModelAPI::AddConstraint(const PowConstraint &cc)  {
  GRB_CALL( GRBaddgenconstrPow(model(), NULL,
              cc.GetArguments()[0], cc.GetResultVar(), cc.GetParameters()[0], "") );
}

void GurobiModelAPI::AddConstraint(const SinConstraint &cc)  {
  GRB_CALL( GRBaddgenconstrSin(model(), NULL,
              cc.GetArguments()[0], cc.GetResultVar(), "") );
}

void GurobiModelAPI::AddConstraint(const CosConstraint &cc)  {
  GRB_CALL( GRBaddgenconstrCos(model(), NULL,
              cc.GetArguments()[0], cc.GetResultVar(), "") );
}

void GurobiModelAPI::AddConstraint(const TanConstraint &cc)  {
  GRB_CALL( GRBaddgenconstrTan(model(), NULL,
              cc.GetArguments()[0], cc.GetResultVar(), "") );
}

void GurobiModelAPI::AddConstraint(const PLConstraint& plc) {
  PLPoints plp(plc.GetParameters());
  GRB_CALL( GRBaddgenconstrPWL(model(), NULL,
              plc.GetArguments()[0], plc.GetResultVar(),
              plp.x_.size(), plp.x_.data(), plp.y_.data()) );
}

static const mp::OptionValueInfo values_item_acceptance[] = {
  { "0", "Not accepted natively, automatic redefinition will be attempted", 0},
  { "1", "Accepted but automatic redefinition will be used where possible", 1},
  { "2", "Accepted natively and preferred", 2}
};

void GurobiModelAPI::InitOptions() {
  /// Constraint acceptance
  GetEnv().AddStoredOption("acc:abs",
                  "Acceptance level for 'abs' expression, default 2:\n"
                  "\n.. value-table::\n"
                  "\n"
                  "It is advisable to experiment with the acceptance levels"
                  " to find best-performing settings.",
                  storedOptions_.acc_abs_, values_item_acceptance);
  GetEnv().AddStoredOption("acc:min",
                  "Acceptance level for 'min' expression, default 2.",
                  storedOptions_.acc_min_);
  GetEnv().AddStoredOption("acc:max",
                  "Acceptance level for 'max' expression, default 2.",
                  storedOptions_.acc_max_);
  GetEnv().AddStoredOption("acc:and",
                  "Acceptance level for 'and'/'forall' expressions, default 2.",
                  storedOptions_.acc_and_);
  GetEnv().AddStoredOption("acc:or",
                  "Acceptance level for 'or'/'exists' expressions, default 2.",
                  storedOptions_.acc_or_);
  GetEnv().AddStoredOption("acc:ind:le acc:ind_le",
                  "Acceptance level for 'implied-less-equal' expression, default 2.",
                  storedOptions_.acc_ind_le_);
  GetEnv().AddStoredOption("acc:ind:eq acc:ind_eq",
                  "Acceptance level for 'implied-equal' expression, default 2.",
                  storedOptions_.acc_ind_eq_);

}

} // namespace mp
