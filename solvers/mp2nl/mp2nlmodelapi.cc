#include <cstring>
#include <unordered_map>

#include "mp2nlmodelapi.h"
#include "mp/nl-solver.hpp"
#include "mp/nl-opcodes.h"
#include "mp/sol-handler.h"


namespace mp {

void MP2NLModelAPI::InitCustomOptions() {
  GetEnv().AddStoredOption("nl:stub stub nlstub",
                           "Filename stub for the NL and SOL files.",
                           storedOptions_.stub_);
  GetEnv().AddStoredOption("nl:text nltextformat nltext",
                           "0*/1: write NL file in text format.",
                           storedOptions_.nl_format_text_);
  GetEnv().AddStoredOption("nl:comments nlcomments",
                           "0*/1: add comments to the text-format NL file.",
                           storedOptions_.nl_comments_);

}

void MP2NLModelAPI::InitProblemModificationPhase(const FlatModelInfo* flat_model_info) {
  alg_con_info_.reserve(
      flat_model_info->GetNumberOfConstraintsOfGroup(CG_Algebraic));
  log_con_info_.reserve(
      flat_model_info->GetNumberOfConstraintsOfGroup(CG_Logical));
  sos_info_.reserve(
      flat_model_info->GetNumberOfConstraintsOfGroup(CG_SOS));
  hdr_is_current_ = false;
}

void MP2NLModelAPI::AddVariables(const VarArrayDef& vad) {
  var_lbs_ = {vad.plb(), (size_t)vad.size()};
  var_ubs_ = {vad.pub(), (size_t)vad.size()};
  var_types_ = {vad.ptype(), (size_t)vad.size()};
  if (vad.pnames())
    var_names_ = {vad.pnames(), (size_t)vad.size()};
  mark_data_.col_sizes_orig_.resize(var_lbs_.size());
  is_var_nlo_.resize(vad.size());
  is_var_nlc_.resize(vad.size());
}


void MP2NLModelAPI::SetLinearObjective( int iobj, const LinearObjective& lo ) {
  assert(iobj == (int)obj_info_.size()
         || (!iobj && 1==obj_info_.size()));    // replacing objective
  auto ii = MakeItemInfo(lo, StaticItemTypeID::ID_LinearObjective, false);
  if (iobj == (int)obj_info_.size())
    obj_info_.push_back(std::move(ii));
  else
    obj_info_[iobj] = std::move(ii);
}

void MP2NLModelAPI::SetQuadraticObjective(int iobj, const QuadraticObjective& qo) {
  assert(iobj == (int)obj_info_.size()
         || (!iobj && 1==obj_info_.size()));    // replacing objective
  /// @todo ?
  throw std::runtime_error("Quadratic objective not supported");
}

void MP2NLModelAPI::SetNLObjective( int iobj, const NLObjective& nlo ) {
  assert(iobj == (int)obj_info_.size()
         || (!iobj && 1==obj_info_.size()));    // replacing objective
  auto ii = MakeItemInfo(nlo, StaticItemTypeID::ID_NLObjective, false);
  if (iobj == (int)obj_info_.size())
    obj_info_.push_back(std::move(ii));
  else obj_info_[iobj] = std::move(ii);
}


void MP2NLModelAPI::AddConstraint(const LinConRange& lc) {
  alg_con_info_.push_back(MakeItemInfo(lc, StaticItemTypeID::ID_LinConRange, false));
}

void MP2NLModelAPI::AddConstraint(const LinConLE& lc)
{ alg_con_info_.push_back(MakeItemInfo(lc, StaticItemTypeID::ID_LinConLE, false)); }

void MP2NLModelAPI::AddConstraint(const LinConEQ& lc)
{ alg_con_info_.push_back(MakeItemInfo(lc, StaticItemTypeID::ID_LinConEQ, false)); }

void MP2NLModelAPI::AddConstraint(const LinConGE& lc)
{ alg_con_info_.push_back(MakeItemInfo(lc, StaticItemTypeID::ID_LinConGE, false)); }


/// To access information from an NLConstraint,
/// use the following accessors (don't use methods of NLConstraint itself):
/// - GetLinSize(nlc), GetLinCoef(nlc, i), GetLinVar(nlc, i),
///   GetExpression(nlc), GetLower(nlc), GetUpper(nlc).
///
/// Implementation follows partly reader_nl.cc from SCIP.
void MP2NLModelAPI::AddConstraint( const NLConstraint& nlc )
{ alg_con_info_.push_back(MakeItemInfo(nlc, StaticItemTypeID::ID_NLConstraint, false)); }

void MP2NLModelAPI::AddConstraint( const NLAssignEQ& nlae )
{ alg_con_info_.push_back(MakeItemInfo(nlae, StaticItemTypeID::ID_NLAssignEQ, false)); }
void MP2NLModelAPI::AddConstraint( const NLAssignLE& nlae )
{ alg_con_info_.push_back(MakeItemInfo(nlae, StaticItemTypeID::ID_NLAssignLE, false)); }
void MP2NLModelAPI::AddConstraint( const NLAssignGE& nlae )
{ alg_con_info_.push_back(MakeItemInfo(nlae, StaticItemTypeID::ID_NLAssignGE, false)); }

void MP2NLModelAPI::AddConstraint(const NLComplementarity& cc)
{ alg_con_info_.push_back(MakeItemInfo(cc, StaticItemTypeID::ID_NLComplementarity, false)); }


void MP2NLModelAPI::AddConstraint( const NLLogical& nll )
{ log_con_info_.push_back(MakeItemInfo(nll, StaticItemTypeID::ID_NLLogical, true)); }

void MP2NLModelAPI::AddConstraint( const NLReifEquiv& nll )
{ log_con_info_.push_back(MakeItemInfo(nll, StaticItemTypeID::ID_NLReifEquiv, true)); }
void MP2NLModelAPI::AddConstraint( const NLReifImpl& nll )
{ log_con_info_.push_back(MakeItemInfo(nll, StaticItemTypeID::ID_NLReifImpl, true)); }
void MP2NLModelAPI::AddConstraint( const NLReifRimpl& nll )
{ log_con_info_.push_back(MakeItemInfo(nll, StaticItemTypeID::ID_NLReifRimpl, true)); }


void MP2NLModelAPI::AddConstraint(const IndicatorConstraintLinLE &ic)
{ log_con_info_.push_back(
      MakeItemInfo(ic, StaticItemTypeID::ID_IndicatorConstraintLinLE, true)); }
void MP2NLModelAPI::AddConstraint(const IndicatorConstraintLinEQ &ic)
{ log_con_info_.push_back(
      MakeItemInfo(ic, StaticItemTypeID::ID_IndicatorConstraintLinEQ, true)); }
void MP2NLModelAPI::AddConstraint(const IndicatorConstraintLinGE &ic)
{ log_con_info_.push_back(
      MakeItemInfo(ic, StaticItemTypeID::ID_IndicatorConstraintLinGE, true)); }



void MP2NLModelAPI::AddConstraint(const SOS1Constraint& sos)
{ sos_info_.push_back(MakeItemInfo(sos, StaticItemTypeID::ID_SOS1Constraint, false)); }
void MP2NLModelAPI::AddConstraint(const SOS2Constraint& sos)
{ sos_info_.push_back(MakeItemInfo(sos, StaticItemTypeID::ID_SOS2Constraint, false)); }


template <class Expr2>
MP2NL_Expr MP2NLModelAPI::AddExpression(
    const Expr2 &expr, ExpressionTypeID eid) {
  SparsityTmp spstmp;
  // Only should visit such arguments
  // which are true expressions, i.e.,
  // not the pure-variable parts?
  // Thus, need to specialize for NLDefVar?
  // But currently seems to be fine as is,
  // because RegisterExpression() only takes true expressions.
  VisitArguments(expr, [this,&spstmp](MP2NL_Expr mp2nle){
    RegisterExpression(mp2nle);
    MergeSparsityTmp(spstmp, mp2nle);
  });
  return StoreMP2NLExprID(expr, eid, spstmp);
}

template <class Expr2>
MP2NL_Expr MP2NLModelAPI::StoreMP2NLExprID(
    const Expr2 &expr, ExpressionTypeID eid, const SparsityTmp& spars) {
  // std::printf("   Storing MP2NL_Expr[%ld]: type %s, logical = %d, @%p\n",
  //             expr_info_.size(),
  //             expr.GetFlatConstraint().GetTypeName(),
  //             IsLogical(expr), &expr);
  expr_info_.push_back(
      MakeItemInfo(expr, eid, IsLogical(expr)));
  expr_counter_.push_back(0);
  expr_sparsity_.push_back( {spars.begin(), spars.end()} );
  return MakeExprID( int(expr_info_.size()-1) );
}

void MP2NLModelAPI::RegisterExpression(MP2NL_Expr expr) {
  CountExpressionOccurrences(expr);
}

void MP2NLModelAPI::MergeSparsityTmp(
    MP2NLModelAPI::SparsityTmp& spstmp, MP2NL_Expr expr) {
  if (expr.IsVariable()) {
    spstmp.insert(expr.GetVarIndex());
  } else if (expr.IsExpression()) {
    const auto& sps4expr =  expr_sparsity_.at(expr.GetExprIndex());
    spstmp.insert(sps4expr.begin(), sps4expr.end());
  }     // That's it
}

/// Count expression depending on its kind.
/// This duplicates the counters in FlatConverter
///   -- could check equality.
void MP2NLModelAPI::CountExpressionOccurrences(MP2NL_Expr expr) {
  if (expr.IsExpression()) {
    auto index = expr.GetExprIndex();
    assert(index >=0 && index < (int)expr_counter_.size()
           && index < (int)expr_info_.size());
    ++expr_counter_[index];
  }   // @todo here also variable usage dep. on top-level item?
}


MP2NL_Expr MP2NLModelAPI::AddExpression(const NLAffineExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_NLAffine); }
MP2NL_Expr MP2NLModelAPI::AddExpression(const NLQuadExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_NLQuad); }
MP2NL_Expr MP2NLModelAPI::AddExpression(const AbsExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_Abs); }

MP2NL_Expr MP2NLModelAPI::AddExpression(const AllDiffExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_AllDiff); }
MP2NL_Expr MP2NLModelAPI::AddExpression(const NumberofConstExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_NumberofConst); }
MP2NL_Expr MP2NLModelAPI::AddExpression(const NumberofVarExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_NumberofVar); }
MP2NL_Expr MP2NLModelAPI::AddExpression(const mp::CountExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_Count); }

MP2NL_Expr MP2NLModelAPI::AddExpression(const AndExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_And); }
MP2NL_Expr MP2NLModelAPI::AddExpression(const OrExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_Or); }
MP2NL_Expr MP2NLModelAPI::AddExpression(const EquivalenceExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_Equivalence); }

MP2NL_Expr MP2NLModelAPI::AddExpression(const CondLTExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_CondLT); }
MP2NL_Expr MP2NLModelAPI::AddExpression(const CondLEExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_CondLE); }
MP2NL_Expr MP2NLModelAPI::AddExpression(const CondEQExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_CondEQ); }
MP2NL_Expr MP2NLModelAPI::AddExpression(const CondGEExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_CondGE); }
MP2NL_Expr MP2NLModelAPI::AddExpression(const CondGTExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_CondGT); }

MP2NL_Expr MP2NLModelAPI::AddExpression(const IfThenExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_IfThen); }
MP2NL_Expr MP2NLModelAPI::AddExpression(const ImplicationExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_Implication); }

MP2NL_Expr MP2NLModelAPI::AddExpression(const NotExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_Not); }
MP2NL_Expr MP2NLModelAPI::AddExpression(const MinExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_Min); }
MP2NL_Expr MP2NLModelAPI::AddExpression(const MaxExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_Max); }

MP2NL_Expr MP2NLModelAPI::AddExpression(const ExpExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_Exp); }
MP2NL_Expr MP2NLModelAPI::AddExpression(const ExpAExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_ExpA); }
MP2NL_Expr MP2NLModelAPI::AddExpression(const LogExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_Log); }
MP2NL_Expr MP2NLModelAPI::AddExpression(const LogAExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_LogA); }
MP2NL_Expr MP2NLModelAPI::AddExpression(const PowExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_Pow); }
MP2NL_Expr MP2NLModelAPI::AddExpression(const PowConstExpExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_PowConstExp); }

MP2NL_Expr MP2NLModelAPI::AddExpression(const SinExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_Sin); }
MP2NL_Expr MP2NLModelAPI::AddExpression(const CosExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_Cos); }
MP2NL_Expr MP2NLModelAPI::AddExpression(const TanExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_Tan); }
MP2NL_Expr MP2NLModelAPI::AddExpression(const AsinExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_Asin); }
MP2NL_Expr MP2NLModelAPI::AddExpression(const AcosExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_Acos); }
MP2NL_Expr MP2NLModelAPI::AddExpression(const AtanExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_Atan); }

MP2NL_Expr MP2NLModelAPI::AddExpression(const SinhExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_Sinh); }
MP2NL_Expr MP2NLModelAPI::AddExpression(const CoshExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_Cosh); }
MP2NL_Expr MP2NLModelAPI::AddExpression(const TanhExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_Tanh); }
MP2NL_Expr MP2NLModelAPI::AddExpression(const AsinhExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_Asinh); }
MP2NL_Expr MP2NLModelAPI::AddExpression(const AcoshExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_Acosh); }
MP2NL_Expr MP2NLModelAPI::AddExpression(const AtanhExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_Atanh); }

MP2NL_Expr MP2NLModelAPI::AddExpression(const DivExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_Div); }


void MP2NLModelAPI::FinishProblemModificationPhase() { }




NLHeader MP2NLModelAPI::Header() {
  if (!hdr_is_current_) {
    PrepareModel();
    hdr_ = DoMakeHeader();
    hdr_is_current_ = true;
  }
  return hdr_;
}

void MP2NLModelAPI::PrepareModel() {
  MapExpressions();
  MarkVars();
  SortVars();
  MarkAlgCons();
  SortAlgCons();
}

void MP2NLModelAPI::MapExpressions() {
  ResetObjMetaInfo();
  is_alg_con_nl_.resize(alg_con_info_.size());  // grow flags array
  for (int i=0; i<(int)obj_info_.size(); ++i) {
    MapExprTreeFromItemInfo(i, obj_info_[i], 2);
  }
  for (int i=0; i<(int)alg_con_info_.size(); ++i) {
    if (!alg_con_info_[i].IsExprTreeMarked()) {
      MapExprTreeFromItemInfo(i, alg_con_info_[i], 0);
      alg_con_info_[i].SetExprTreeMarked();
    }
  }
  for (int i=0; i<(int)log_con_info_.size(); ++i) {
    if (!log_con_info_[i].IsExprTreeMarked()) {
      MapExprTreeFromItemInfo(i, log_con_info_[i], 1);
      log_con_info_[i].SetExprTreeMarked();
    }
  }
  FinishMapExpressions();
  ResetIniExprRetrievedFlags();
}

void MP2NLModelAPI::MapExprTreeFromItemInfo(
    int i_item, const ItemInfo& info, int kind) {
  auto mp2nlexpr        // so that AddExpression() is called
      = info.GetDispatcher().GetExpression(info.GetPItem());
  RegisterExpression(mp2nlexpr);     // tree root
  if (1 != kind) {                   // alg con / obj
    MergeItemSparsity(info, kind, mp2nlexpr);
    UpdateAlgebraicMetaInfo(info, kind);  // nnz, colsizes, n ranges/eqns
  }
  MarkNLVars(i_item, mp2nlexpr, kind);            // also for logical cons
}

void MP2NLModelAPI::ResetObjMetaInfo() {
  mark_data_.n_obj_nz_ = 0;
  mark_data_.nnlo_ = 0;
  is_var_nlo_.clear();
  is_var_nlo_.resize(var_lbs_.size());
}

void MP2NLModelAPI::MergeItemSparsity(
    const ItemInfo &info, int , MP2NL_Expr expr) {
  if (expr.IsEmptyExpr()) {    // Just use existing linear part
    auto lp = GetLinPart(info);
    info.SetExtLinPart(std::move(lp));
  } else {                     // proper expr or variable
    std::unordered_map<int, double> linp_ext;
    auto lp = GetLinPart(info);
    for (auto i=lp.size(); i--; ) {
      linp_ext[lp.vars()[i]] += lp.coefs()[i];
    }
    if (expr.IsVariable()) {   // @todo should have been inlined
      linp_ext[expr.GetVarIndex()] += 1.0;      // Merge 1 var
    } else {
      const auto& spars = expr_sparsity_.at(expr.GetExprIndex());
      for (auto v: spars)
        linp_ext[v];           // Merge expression sparsity
    } // https://stackoverflow.com/questions/33325470/
      // c-is-there-a-way-to-initialize-unordered-map-using-
      // contents-of-vector-as-keys
    std::vector<int> vars;     // @todo SmallVec
    vars.reserve(linp_ext.size());
    std::vector<double> coefs;
    coefs.reserve(linp_ext.size());
    for (const auto& val: linp_ext) {
      coefs.push_back(val.second);
      vars.push_back(val.first);
    }
    info.SetExtLinPart(        // Put joint result into \a info
        {std::move(coefs), std::move(vars)} );
  }
}

inline
    MP2NLModelAPI::LinPartRefOrStorage
    GetLPRoS(const LinTerms& lt) {
  return {lt.coefs(), lt.vars()};
}

template <int sense>
inline
    MP2NLModelAPI::LinPartRefOrStorage
    GetLPRoS(const NLBaseAssign<sense>& nla) {
  std::vector<double> c{{nla.coef(0)}};    // @todo: SmallVec
  std::vector<int> v{{nla.var(0)}};
  return {std::move(c), std::move(v)};      // std::move to pass storage
}

MP2NLModelAPI::LinPartRefOrStorage
MP2NLModelAPI::GetLinPart(const ItemInfo &info) {
  const auto* pitem = info.GetPItem();
  switch (info.GetStaticTypeID()) {
  case StaticItemTypeID::ID_LinearObjective:
  case StaticItemTypeID::ID_NLObjective: {
    return GetLPRoS( ((LinearObjective*)(pitem))->GetLinTerms() );
  }
  case StaticItemTypeID::ID_LinConRange: {
    return GetLPRoS( ((LinConRange*)(pitem))->GetLinTerms() );
  }
  case StaticItemTypeID::ID_LinConLE: {
    return GetLPRoS( ((LinConLE*)(pitem))->GetLinTerms() );
  }
  case StaticItemTypeID::ID_LinConEQ: {
    return GetLPRoS( ((LinConEQ*)(pitem))->GetLinTerms() );
  }
  case StaticItemTypeID::ID_LinConGE: {
    return GetLPRoS( ((LinConGE*)(pitem))->GetLinTerms() );
  }
  case StaticItemTypeID::ID_NLConstraint: {
    return GetLPRoS(
        ((const NLConstraint*)(pitem))->GetMainCon().GetLinTerms() );
  }
  case StaticItemTypeID::ID_NLAssignLE: {
    return GetLPRoS( *((NLAssignLE*)(pitem)) );
  }
  case StaticItemTypeID::ID_NLAssignEQ: {
    return GetLPRoS( *((NLAssignEQ*)(pitem)) );
  }
  case StaticItemTypeID::ID_NLAssignGE: {
    return GetLPRoS( *((NLAssignGE*)(pitem)) );
  }
  case StaticItemTypeID::ID_NLComplementarity: {
    return GetLPRoS(                 // Do we need the \a compl_var?
        ((const NLComplementarity*)(pitem))->GetExpression().GetBody() );
  }
  default:
    MP_RAISE("Unknown objective or algebraic constraint type");
  }
}

void MP2NLModelAPI::UpdateAlgebraicMetaInfo(
    const ItemInfo &info, int kind) {
  const auto& lp_ext = info.GetExtLinPart();
  if (kind) {
    mark_data_.n_obj_nz_ += lp_ext.size();
  } else {
    mark_data_.n_con_nz_ += lp_ext.size();
    Add2ColSizes(lp_ext.vars());
    info.GetDispatcher().MarkRangeOrEqn(info.GetPItem());
  }
}

void MP2NLModelAPI::MarkNLVars(
    int i_item, MP2NL_Expr expr, int kind) {
  if (expr.IsEmptyExpr()) {
    // nothing
  } else {
    std::vector<bool>& var_nl_flags
        = kind ? is_var_nlo_ : is_var_nlc_;
    if (expr.IsVariable()) {         // we still count it as non-linear
      var_nl_flags.at(expr.GetVarIndex()) = true;
    } else {
      for (auto v: expr_sparsity_.at(expr.GetExprIndex())) {
        assert(v>=0 && v<(int)var_lbs_.size());
        var_nl_flags[v] = true;
      }
    }
    if (0==kind) {
      is_alg_con_nl_.at(i_item) = true;   // nonlinear cons
      ++mark_data_.nnlc_;
    } else {
      ++mark_data_.nnlo_;
    }
  }
}

void MP2NLModelAPI::FinishMapExpressions() { }

void MP2NLModelAPI::MarkVars() {
  mark_data_.var_prior_.clear();
  mark_data_.var_prior_.resize(var_lbs_.size());
  mark_data_.n_bin_var_lin_ = 0;
  mark_data_.n_int_var_lin_ = 0;
  mark_data_.n_cont_var_nl_con_ = 0;
  mark_data_.n_cont_var_nl_obj_ = 0;
  mark_data_.n_cont_var_nl_both_ = 0;
  mark_data_.n_int_var_nl_con_ = 0;
  mark_data_.n_int_var_nl_obj_ = 0;
  mark_data_.n_int_var_nl_both_ = 0;
  for (auto i=var_lbs_.size(); i--; ) {
    mark_data_.var_prior_[i].second = i;
    bool f_int_non_bin{}, f_bin{}, f_int_or_bin{};
    if (var::Type::INTEGER == var_types_[i]) {
      f_int_or_bin = true;
      if (!var_lbs_[i] && 1.0==var_ubs_[i]) {    // binary var
        mark_data_.var_prior_[i].first = 1;
        f_bin = true;
      } else {
        mark_data_.var_prior_[i].first = 2;
        f_int_non_bin = true;
      }
    }
    if (is_var_nlc_[i]) {
      if (is_var_nlo_[i]) {
        mark_data_.var_prior_[i].first -= 9;
        f_int_or_bin
            ? ++mark_data_.n_int_var_nl_both_
            : ++mark_data_.n_cont_var_nl_both_;
      } else {
        mark_data_.var_prior_[i].first -= 6;
        f_int_or_bin
            ? ++mark_data_.n_int_var_nl_con_
            : ++mark_data_.n_cont_var_nl_con_;
      }
    } else if (is_var_nlo_[i]) {
      mark_data_.var_prior_[i].first -= 3;
      f_int_or_bin
          ? ++mark_data_.n_int_var_nl_obj_
          : ++mark_data_.n_cont_var_nl_obj_;
    } else {               // 'linear' var
      if (f_bin)
        ++mark_data_.n_bin_var_lin_;
      else
        if (f_int_non_bin)
          ++mark_data_.n_int_var_lin_;
    }
  }
}

void MP2NLModelAPI::SortVars() {
  std::stable_sort(mark_data_.var_prior_.begin(), mark_data_.var_prior_.end());
  mark_data_.var_order_12_.resize(var_lbs_.size());
  mark_data_.var_order_21_.resize(var_lbs_.size());
  for (auto i=var_lbs_.size(); i--; ) {
    mark_data_.var_order_12_[i] = mark_data_.var_prior_[i].second;
    mark_data_.var_order_21_[mark_data_.var_prior_[i].second] = i;
  }
}

void MP2NLModelAPI::MarkAlgCons() {
  mark_data_.con_prior_.clear();
  mark_data_.con_prior_.resize(alg_con_info_.size());
  for (auto i=alg_con_info_.size(); i--; ) {
    mark_data_.con_prior_[i] = { -(int)is_alg_con_nl_[i], i };
  }
}

void MP2NLModelAPI::SortAlgCons() {
  std::stable_sort(mark_data_.con_prior_.begin(), mark_data_.con_prior_.end());
  mark_data_.con_order_12_.resize(alg_con_info_.size());
  mark_data_.con_order_21_.resize(alg_con_info_.size());
  for (auto i=alg_con_info_.size(); i--; ) {
    mark_data_.con_order_12_[i] = mark_data_.con_prior_[i].second;
    mark_data_.con_order_21_[mark_data_.con_prior_[i].second] = i;
  }
}

void MP2NLModelAPI::Add2ColSizes(ArrayRef<int> vars) {
  for (auto v: vars)
    ++mark_data_.col_sizes_orig_[v];
}


NLHeader MP2NLModelAPI::DoMakeHeader() {
  mp::NLHeader hdr;

  hdr.num_vars = var_lbs_.size();
  hdr.num_algebraic_cons = alg_con_info_.size();
  hdr.num_objs = obj_info_.size();
  hdr.num_ranges = mark_data_.n_ranges_;
  hdr.num_eqns = mark_data_.n_eqns_;
  hdr.num_logical_cons = log_con_info_.size();

  /** Total number of nonlinear constraints. */
  hdr.num_nl_cons = mark_data_.nnlc_;
  hdr.num_nl_objs = mark_data_.nnlo_;
  hdr.num_compl_conds = 0;
  hdr.num_nl_compl_conds = 0;
  hdr.num_compl_dbl_ineqs = 0;
  hdr.num_compl_vars_with_nz_lb = 0;

  /** Number of nonlinear network constraints. */
  hdr.num_nl_net_cons = 0;
  hdr.num_linear_net_cons = 0;

  /** Number of nonlinear variables in both constraints and objectives. */
  hdr.num_nl_vars_in_both
      = mark_data_.n_cont_var_nl_both_ + mark_data_.n_int_var_nl_both_;

  /**
      Number of nonlinear variables in constraints including nonlinear
      variables in both constraints and objectives.
     */
  hdr.num_nl_vars_in_cons
      = hdr.num_nl_vars_in_both
        + mark_data_.n_cont_var_nl_con_ + mark_data_.n_int_var_nl_con_;

  /**
      Number of nonlinear variables in objectives including nonlinear
      variables in both constraints and objectives.
     */
  hdr.num_nl_vars_in_objs
      = hdr.num_nl_vars_in_cons
        + mark_data_.n_cont_var_nl_obj_ + mark_data_.n_int_var_nl_obj_;

  // Miscellaneous
  // -------------

  /** Number of linear network variables (arcs). */
  hdr.num_linear_net_vars = 0;

  /** Number of functions. */
  hdr.num_funcs = 0;

  // Information about discrete variables
  // ------------------------------------

  /** Number of linear binary variables. */
  hdr.num_linear_binary_vars = mark_data_.n_bin_var_lin_;

  /** Number of linear non-binary integer variables. */
  hdr.num_linear_integer_vars = mark_data_.n_int_var_lin_;

  /**
      Number of integer nonlinear variables in both constraints and objectives.
     */
  hdr.num_nl_integer_vars_in_both = mark_data_.n_int_var_nl_both_;

  /** Number of integer nonlinear variables just in constraints. */
  hdr.num_nl_integer_vars_in_cons = mark_data_.n_int_var_nl_con_;

  /** Number of integer nonlinear variables just in objectives. */
  hdr.num_nl_integer_vars_in_objs = mark_data_.n_int_var_nl_obj_;

  // Information about nonzeros
  // --------------------------

  /** Number of nonzeros in constraints' Jacobian. */
  hdr.num_con_nonzeros = mark_data_.n_con_nz_;

  /** Number of nonzeros in all objective gradients. */
  hdr.num_obj_nonzeros = mark_data_.n_obj_nz_;

  // Information about names
  // -----------------------

  /** Length of longest con/obj name if names are available. */
  hdr.max_con_name_len = 0;    // no need to set

  /** Length of longest variable name if names are available. */
  hdr.max_var_name_len = 0;    // no need to set

  // Information about common expressions
  // ------------------------------------

  /**
      Number of common expressions that appear both in constraints
      and objectives.
     */
  hdr.num_common_exprs_in_both = 0;

  /**
      Number of common expressions that appear in multiple constraints
      and don't appear in objectives.
     */
  hdr.num_common_exprs_in_cons = 0;

  /**
      Number of common expressions that appear in multiple objectives
      and don't appear in constraints.
     */
  hdr.num_common_exprs_in_objs = 0;

  /**
      Number of common expressions that only appear in a single constraint
      and don't appear in objectives.
     */
  hdr.num_common_exprs_in_single_cons = 0;

  /**
      Number of common expressions that only appear in a single objective
      and don't appear in constraints.
     */
  hdr.num_common_exprs_in_single_objs = 0;

  hdr.prob_name = "mp2nl_model";

  hdr.format = storedOptions_.nl_format_text_
                   ? mp::NLHeader::TEXT : mp::NLHeader::BINARY;

  return hdr;
}

int MP2NLModelAPI::ObjType(int i) {
  switch (obj_info_[i].GetStaticTypeID()) {
  case StaticItemTypeID::ID_LinearObjective:
  case StaticItemTypeID::ID_NLObjective:
    return obj::MIN==((LinearObjective*)(obj_info_[i].GetPItem()))->obj_sense()
               ? 0 : 1;
  default:
    MP_RAISE("Unknown objective type");
  }
}

template <class ObjGradWriterFactory>
void MP2NLModelAPI::FeedObjGradient(int i, ObjGradWriterFactory& svwf) {
  FeedExtLinPart(obj_info_[i], svwf);
}

template <class ObjExprWriter>
void MP2NLModelAPI::FeedObjExpression(int iobj, ObjExprWriter& ew) {
  const auto& obj_info = obj_info_.at(iobj);
  switch (obj_info.GetStaticTypeID()) {
  case StaticItemTypeID::ID_LinearObjective:
    ew.NPut(0.0);
    break;
  case StaticItemTypeID::ID_NLObjective: {
    const auto& obj = *((NLObjective*)(obj_info.GetPItem()));
    if (obj.HasExpr()) {
      FeedExpr(GetExpression(obj), ew);
    } else
      ew.NPut(0.0);
  } break;
  default:
    MP_RAISE("Unknown objective type");
  }
}

template <class VarBoundsWriter>
void MP2NLModelAPI::FeedVarBounds(VarBoundsWriter& vbw) {
  for (size_t i = 0; i < var_lbs_.size(); i++) {
    auto i_old = GetOldVarIndex(i);
    vbw.WriteLbUb(var_lbs_[i_old], var_ubs_[i_old]);
  }
}

template <class ConBoundsWriter>
void MP2NLModelAPI::FeedConBounds(ConBoundsWriter& cbw) {
  for (size_t i_new=0; i_new<alg_con_info_.size(); ++i_new) {
    auto i = GetOldAlgConIndex(i_new);
    switch (alg_con_info_[i].GetStaticTypeID()) {
    case StaticItemTypeID::ID_LinConRange: {
      const auto& lcon = *((LinConRange*)(alg_con_info_[i].GetPItem()));
      cbw.WriteAlgConRange( AlgConRange{lcon.lb(), lcon.ub()} );
    } break;
    case StaticItemTypeID::ID_LinConLE: {
      const auto& lcon = *((LinConLE*)(alg_con_info_[i].GetPItem()));
      cbw.WriteAlgConRange( AlgConRange{lcon.lb(), lcon.ub()} );
    } break;
    case StaticItemTypeID::ID_LinConEQ: {
      const auto& lcon = *((LinConEQ*)(alg_con_info_[i].GetPItem()));
      cbw.WriteAlgConRange( AlgConRange{lcon.lb(), lcon.ub()} );
    } break;
    case StaticItemTypeID::ID_LinConGE: {
      const auto& lcon = *((LinConGE*)(alg_con_info_[i].GetPItem()));
      cbw.WriteAlgConRange( AlgConRange{lcon.lb(), lcon.ub()} );
    } break;
    case StaticItemTypeID::ID_NLConstraint: {
      const auto& lcon = *((NLConstraint*)(alg_con_info_[i].GetPItem()));
      cbw.WriteAlgConRange( AlgConRange{GetLower(lcon), GetUpper(lcon)} );
    } break;
    case StaticItemTypeID::ID_NLAssignLE: {  // -resvar + (expr) >= 0
      const auto& lcon = *((NLAssignLE*)(alg_con_info_[i].GetPItem()));
      cbw.WriteAlgConRange( AlgConRange{0.0, Infinity()} );
    } break;
    case StaticItemTypeID::ID_NLAssignEQ: {
      const auto& lcon = *((NLAssignEQ*)(alg_con_info_[i].GetPItem()));
      cbw.WriteAlgConRange( AlgConRange{0.0, 0.0} );
    } break;
    case StaticItemTypeID::ID_NLAssignGE: {
      const auto& lcon = *((NLAssignGE*)(alg_con_info_[i].GetPItem()));
      cbw.WriteAlgConRange( AlgConRange{MinusInfinity(), 0.0} );
    } break;
    case StaticItemTypeID::ID_NLComplementarity: {
      const auto& lcon = *((NLComplementarity*)(alg_con_info_[i].GetPItem()));
      AlgConRange bnd;
      auto j = lcon.GetVariable();
      auto ct = lcon.GetExpression().constant_term(); // @todo NLExpression?
      bnd.L = MinusInfinity();
      bnd.U = Infinity();
      bnd.k = 0;
      if (var_lbs_[j] > MinusInfinity()) {
        bnd.k = 1;
        bnd.L = -ct;         // empty expr then
      }
      if (var_ubs_[j] < Infinity()) {
        bnd.k |= 2;
        bnd.U = -ct;
      }
      if (3==bnd.k) {
        bnd.L = MinusInfinity();
        bnd.U = Infinity();  // expr should be ct.
      }
      assert(bnd.k);
      bnd.cvar = GetNewVarIndex(j);
      cbw.WriteAlgConRange( bnd );
    } break;
    default:
      MP_RAISE("Unknown algebraic constraint type");
    }
  }
}

template <class ConLinearExprWriterFactory>
void MP2NLModelAPI::FeedLinearConExpr(int i, ConLinearExprWriterFactory& svwf) {
  FeedExtLinPart(alg_con_info_[GetOldAlgConIndex(i)], svwf);
}

template <class ConLinearExprWriterFactory>
void MP2NLModelAPI::FeedExtLinPart(
    const ItemInfo& item,
    ConLinearExprWriterFactory& svwf) {
  const auto& lp_ext = item.GetExtLinPart();
  if (lp_ext.size()) {
    std::vector<int> vars_tmp = lp_ext.vars();
    for (auto j=vars_tmp.size(); j--; )
      vars_tmp[j] = GetNewVarIndex( vars_tmp[j] );        // new ordering
    LinTerms lt_srt {lp_ext.coefs(), std::move(vars_tmp)};
    lt_srt.sort_terms__leave_0s();   // Leave sparsity pattern
    assert(lt_srt.size());
    {   // Sparsity pattern can have coefs 0.0
      auto svw = svwf.MakeVectorWriter(lt_srt.size());
      for (int j=0; j<lt_srt.size(); ++j) {
        svw.Write(lt_srt.var(j), lt_srt.coef(j));
      }
    }
  }
}

template <class ConExprWriter>
void MP2NLModelAPI::FeedConExpression(
    int icon, ConExprWriter& ew) {
  if (icon < (int)alg_con_info_.size())
    FeedAlgConExpression(GetOldAlgConIndex(icon), ew);
  else
    FeedLogicalConExpression(icon - alg_con_info_.size(), ew);
}

template <class ConExprWriter>
void MP2NLModelAPI::FeedAlgConExpression(
    int icon, ConExprWriter& ew) {
  // We could use a virtual function
  // to call GetExpression() instead of switch().
  const auto& con_info = alg_con_info_.at(icon);
  const auto* pitem = con_info.GetPItem();
  switch (con_info.GetStaticTypeID()) {
  case StaticItemTypeID::ID_NLConstraint: {
    const auto& con = *((NLConstraint*)(pitem));
    if (con.HasExpr()) {
      FeedExpr(GetExpression(con), ew);
    } else
      ew.NPut(0.0);                   // linear case
  } break;
  case StaticItemTypeID::ID_NLAssignLE:
    FeedExpr(GetExpression(*((NLAssignLE*)(pitem))), ew);
    break;
  case StaticItemTypeID::ID_NLAssignEQ:
    FeedExpr(GetExpression(*((NLAssignEQ*)(pitem))), ew);
    break;
  case StaticItemTypeID::ID_NLAssignGE:
    FeedExpr(GetExpression(*((NLAssignGE*)(pitem))), ew);
    break;
  case StaticItemTypeID::ID_NLComplementarity: {
    const auto& cc = *((NLComplementarity*)(pitem));
    auto j = cc.GetVariable();
    if (var_lbs_[j]<=MinusInfinity()   // @todo proper expressions
        || var_ubs_[j]>=Infinity())
      ew.NPut(0.0);
    else
      ew.NPut(cc.GetExpression().constant_term());
  } break;
  default:
    ew.NPut(0.0);                     // must be linear constr
    break;
  case StaticItemTypeID::ID_NLLogical:
  case StaticItemTypeID::ID_NLReifImpl:
  case StaticItemTypeID::ID_NLReifEquiv:
  case StaticItemTypeID::ID_NLReifRimpl:
  case StaticItemTypeID::ID_LinearObjective:
  case StaticItemTypeID::ID_NLObjective:
  case StaticItemTypeID::ID_None:
    MP_RAISE("Unknown algebraic constraint type");
  }
}

template <class ConExprWriter>
void MP2NLModelAPI::FeedLogicalConExpression(
    int icon, ConExprWriter& ew) {
  const auto& con_info = log_con_info_.at(icon);
  const auto* pitem = con_info.GetPItem();
  switch (con_info.GetStaticTypeID()) {
  case StaticItemTypeID::ID_NLLogical:
    FeedNLLogical(*((NLLogical*)(pitem)), ew);
    break;
  case StaticItemTypeID::ID_NLReifImpl: {
    FeedReification((*((NLReifImpl*)(pitem))), ew);
  } break;
  case StaticItemTypeID::ID_NLReifEquiv:
    FeedReification((*((NLReifEquiv*)(pitem))), ew);
    break;
  case StaticItemTypeID::ID_NLReifRimpl:
    FeedReification((*((NLReifRimpl*)(pitem))), ew);
    break;
  default:
    MP_RAISE("Unknown logical constraint type");
  }
}

template <class ExprWriter>
void MP2NLModelAPI::FeedExpr(Expr expr, ExprWriter& ew) {
  if (expr.IsEmptyExpr())
    ew.NPut(0.0);
  else if (expr.IsVariable()) {       // @todo def vars
    ew.VPut( GetNewVarIndex( expr.GetVarIndex() ),
            GetVarName(expr.GetVarIndex()));
  } else
    FeedOpcode(expr, ew);
}

#define HANDLE_OPCODE_CASE_1_ARG(Expr, Opcode, fdr) \
case ExpressionTypeID::ID_ ## Expr: \
    fdr(*(const Expr ## Expression*)pitem, ew.OPut1(nl::Opcode)); \
    break;
#define HANDLE_OPCODE_CASE_2_ARG(Expr, Opcode, fdr) \
case ExpressionTypeID::ID_ ## Expr: \
    fdr(*(const Expr ## Expression*)pitem, ew.OPut2(nl::Opcode)); \
    break;
#define HANDLE_OPCODE_CASE_3_ARG(Expr, Opcode, fdr) \
case ExpressionTypeID::ID_ ## Expr: \
    fdr(*(const Expr ## Expression*)pitem, ew.OPut3(nl::Opcode)); \
    break;
#define HANDLE_OPCODE_CASE_N_ARG(Expr, Opcode, fdr) \
case ExpressionTypeID::ID_ ## Expr: \
    fdr(*(const mp::Expr ## Expression*)pitem, \
      ew.OPutN(nl::Opcode, \
        GetNumArguments(*(const mp::Expr ## Expression*)pitem) \
        + GetNumParameters(*(const mp::Expr ## Expression*)pitem))); \
    break;
#define HANDLE_OPCODE_CASE_N_OR_2_ARG(Expr, Opcode, Opcode2, fdr) \
case ExpressionTypeID::ID_ ## Expr: { \
    auto n_args = GetNumArguments(*(const Expr ## Expression*)pitem); \
    fdr(*(const Expr ## Expression*)pitem, \
      2==n_args ? ew.OPut2(nl::Opcode2) \
      : ew.OPutN(nl::Opcode, n_args)); \
} break;
#define HANDLE_RELATIONAL(Expr, Opcode) \
case ExpressionTypeID::ID_ ## Expr: \
    FeedRelational(*(const Expr ## Expression*)pitem, ew.OPut2(nl::Opcode)); \
    break;

template <class ExprWriter>
void MP2NLModelAPI::FeedOpcode(Expr expr, ExprWriter& ew) {
  assert(expr.IsExpression());
  const auto& expr_info = expr_info_.at(expr.GetExprIndex());
  const auto* pitem = expr_info.GetPItem();
  switch (expr_info.GetExprTypeID()) {

  case ExpressionTypeID::ID_NLAffine:
    FeedAlgebraic(*(const NLAffineExpression*)pitem, ew);
    break;
  case ExpressionTypeID::ID_NLQuad:
    FeedAlgebraic(*(const NLQuadExpression*)pitem, ew);
    break;

    HANDLE_OPCODE_CASE_N_ARG(AllDiff, ALLDIFF, FdArgs)
    HANDLE_OPCODE_CASE_N_ARG(NumberofConst, NUMBEROF, FdArgs_Params1st)
    HANDLE_OPCODE_CASE_N_ARG(NumberofVar, NUMBEROF, FdArgs)
    HANDLE_OPCODE_CASE_N_ARG(Count, COUNT, FdLogicArgs)
    HANDLE_OPCODE_CASE_N_OR_2_ARG(And, FORALL, AND, FdLogicArgs)
    HANDLE_OPCODE_CASE_N_OR_2_ARG(Or, EXISTS, OR, FdLogicArgs)
    HANDLE_OPCODE_CASE_2_ARG(Equivalence, IFF, FdLogicArgs)

    HANDLE_OPCODE_CASE_3_ARG(IfThen, IF, FdIfArgs)    // @todo logic?
    HANDLE_OPCODE_CASE_3_ARG(Implication, IMPLICATION, FdLogicArgs)
    HANDLE_OPCODE_CASE_1_ARG(Not, NOT, FdLogicArgs)

    HANDLE_RELATIONAL(CondLT, LT)
    HANDLE_RELATIONAL(CondLE, LE)
    HANDLE_RELATIONAL(CondEQ, EQ)
    HANDLE_RELATIONAL(CondGE, GE)
    HANDLE_RELATIONAL(CondGT, GT)

    HANDLE_OPCODE_CASE_1_ARG(Abs, ABS, FdArgs)
    HANDLE_OPCODE_CASE_N_ARG(Min, MIN, FdArgs)
    HANDLE_OPCODE_CASE_N_ARG(Max, MAX, FdArgs)

    HANDLE_OPCODE_CASE_1_ARG(Exp, EXP, FdArgs)
    HANDLE_OPCODE_CASE_2_ARG(ExpA, POW, FdArgs_Params1st)
    HANDLE_OPCODE_CASE_1_ARG(Log, LOG, FdArgs)
  case ExpressionTypeID::ID_LogA:
    FdLogA(*(const LogAExpression*)pitem, ew);
    break;

    HANDLE_OPCODE_CASE_2_ARG(Pow, POW, FdArgs)
    HANDLE_OPCODE_CASE_2_ARG(PowConstExp, POW, FdArgs)
    HANDLE_OPCODE_CASE_1_ARG(Sin, SIN, FdArgs)
    HANDLE_OPCODE_CASE_1_ARG(Cos, COS, FdArgs)
    HANDLE_OPCODE_CASE_1_ARG(Tan, TAN, FdArgs)
    HANDLE_OPCODE_CASE_1_ARG(Asin, ASIN, FdArgs)
    HANDLE_OPCODE_CASE_1_ARG(Acos, ACOS, FdArgs)
    HANDLE_OPCODE_CASE_1_ARG(Atan, ATAN, FdArgs)
    HANDLE_OPCODE_CASE_1_ARG(Sinh, SINH, FdArgs)
    HANDLE_OPCODE_CASE_1_ARG(Cosh, COSH, FdArgs)
    HANDLE_OPCODE_CASE_1_ARG(Tanh, TANH, FdArgs)
    HANDLE_OPCODE_CASE_1_ARG(Asinh, ASINH, FdArgs)
    HANDLE_OPCODE_CASE_1_ARG(Acosh, ACOSH, FdArgs)
    HANDLE_OPCODE_CASE_1_ARG(Atanh, ATANH, FdArgs)

    HANDLE_OPCODE_CASE_2_ARG(Div, DIV, FdArgs)

  default:
    MP_RAISE("MP2NL: unknown expression type");
  }
}

template <class AlgMPExpr, class ExprWriter>
void MP2NLModelAPI::FeedAlgebraic(
    const AlgMPExpr& e, ExprWriter& ew) {
  int n_args = GetLinSize(e) + GetQuadSize(e) + bool(GetConstTerm(e));
	auto write_args = [&](const auto& e, auto ew_args0) {
    for (int i=0; i<GetLinSize(e); ++i) {
      if (1.0==GetLinCoef(e, i))
        ew_args0.EPut(GetLinTerm(e, i));
      else {
        auto ew_args1 = ew_args0.OPut2(nl::MUL);
        ew_args1.NPut(GetLinCoef(e, i));
        ew_args1.EPut(GetLinTerm(e, i));
      }
    }
    if constexpr (std::is_same_v<AlgMPExpr, NLQuadExpression>) {
      for (int i=0; i<GetQuadSize(e); ++i) {
        if (1.0==GetQuadCoef(e, i)) {
          auto ew_args1 = ew_args0.OPut2(nl::MUL);
          ew_args1.EPut(GetQuadTerm1(e, i));
          ew_args1.EPut(GetQuadTerm2(e, i));
        } else {
          auto ew_args1 = ew_args0.OPut2(nl::MUL);
          ew_args1.NPut(GetQuadCoef(e, i));
          auto ew_args2 = ew_args1.OPut2(nl::MUL);
          ew_args2.EPut(GetQuadTerm1(e, i));
          ew_args2.EPut(GetQuadTerm2(e, i));
        }
      }
    }
    if (GetConstTerm(e))
      ew_args0.NPut(GetConstTerm(e));
  };
  if (!n_args)
    ew.NPut(0.0);
  else if (1==n_args)
    write_args(e, std::move(ew));
  else if (2==n_args)
    write_args(e, ew.OPut2(nl::ADD));
  else
    write_args(e, ew.OPutN(nl::SUM, n_args));
}

template <class CondMPExpr, class ExprWriter>
void MP2NLModelAPI::FeedRelational(
    const CondMPExpr& e, ExprWriter ew) {
  ew.EPut(GetExpression(e));
  ew.NPut(GetRHS(e));
}

template <class ExprWriter>
void MP2NLModelAPI::FeedNLLogical(const NLLogical& e, ExprWriter& ew) {
  if (!e.IsTrue()) {
    auto ew_args = ew.OPut1(nl::NOT);
    FeedExpr(GetExpression(e), ew_args);
  } else {
    FeedExpr(GetExpression(e), ew);
  }
}

template <int sense, class ExprWriter>
void MP2NLModelAPI::FeedReification(
    const NLBaseReif<sense>& e, ExprWriter& ew) {
  if (sense<0) {    // @todo could do some via IFF; except, e.g., AllDiff
    auto args = ew.OPut2(nl::OR);    // However this would set CTX_MIX
    auto args_lhs = args.OPut2(nl::EQ);
    args_lhs.VPut( GetNewVarIndex(GetVariable(e)) );
    args_lhs.NPut(0.0);
    args.EPut(GetExpression(e));
  } else if (sense>0) {
    auto args = ew.OPut2(nl::OR);
    auto args_lhs = args.OPut2(nl::EQ);
    args_lhs.VPut( GetNewVarIndex(GetVariable(e)) );
    args_lhs.NPut(1.0);
    auto args_rhs = args.OPut1(nl::NOT);
    args_rhs.EPut(GetExpression(e));
  } else {
    auto args = ew.OPut2(nl::IFF);
    auto args_lhs = args.OPut2(nl::EQ);
    args_lhs.VPut( GetNewVarIndex(GetVariable(e)) );
    args_lhs.NPut(1.0);
    args.EPut(GetExpression(e));
  }
}

template <class MPExpr, class ArgWriter>
void MP2NLModelAPI::FdArgs(
    const MPExpr& e, ArgWriter ew_arg) {
  for (int i=0; i<GetNumArguments(e); ++i)
    ew_arg.EPut(GetArgExpression(e, i));
  for (int i=0; i<GetNumParameters(e); ++i)
    ew_arg.NPut(GetParameter(e, i));
}

template <class MPExpr, class ArgWriter>
void MP2NLModelAPI::FdArgs_Params1st(
    const MPExpr& e, ArgWriter ew_arg) {
  for (int i=0; i<GetNumParameters(e); ++i)
    ew_arg.NPut(GetParameter(e, i));
  for (int i=0; i<GetNumArguments(e); ++i)
    ew_arg.EPut(GetArgExpression(e, i));
}

template <class MPExpr, class ArgWriter>
void MP2NLModelAPI::FdLogA(
    const MPExpr& e, ArgWriter& aw) {
  auto A = GetParameter(e, 0);
  if (std::fabs(10.0-A) < 1e-15) {
    auto aw1 = aw.OPut1(nl::LOG10);
    aw1.EPut(GetArgExpression(e, 0));
  } else if (std::fabs(exp(1.0)-A) < 1e-15) {
    auto aw1 = aw.OPut1(nl::LOG);
    aw1.EPut(GetArgExpression(e, 0));
  } else {
    MP_ASSERT_ALWAYS(std::fabs(1.0-A)>1e-15,
                     "logarithm with base 1");
    auto aw1 = aw.OPut2(nl::MUL);
    aw1.NPut(1.0/std::log(A));
    auto aw2 = aw1.OPut1(nl::LOG);
    aw2.EPut(GetArgExpression(e, 0));
  }
}

template <class MPExpr, class ArgWriter>
void MP2NLModelAPI::FdLogicArgs(
    const MPExpr& e, ArgWriter ew_arg) {
  // std::printf("   Feed '%s':  %d args\n",
  //             GetItemName(e), GetNumArguments(e));
  for (int i=0; i<GetNumArguments(e); ++i) {
    FeedLogicalExpression(
        GetArgExpression(e, i), ew_arg);
  }
  for (int i=0; i<GetNumParameters(e); ++i)
    ew_arg.NPut(GetParameter(e, i));
}

template <class MPExpr, class ArgWriter>
void MP2NLModelAPI::FdIfArgs(
    const MPExpr& e, ArgWriter ew_arg) {
  for (int i=0; i<GetNumArguments(e); ++i) {
    if (!i)                            // 1st arg.
      FeedLogicalExpression(           // @todo others too?
          GetArgExpression(e, i), ew_arg);
    else
      ew_arg.EPut(GetArgExpression(e, i));
  }
}

template <class ArgWriter>
void MP2NLModelAPI::FeedLogicalExpression(
    MP2NL_Expr mp2nle, ArgWriter& aw) {
  if (mp2nle.IsVariable()) {
    auto arg2 = aw.OPut2(nl::NE);
    arg2.EPut(mp2nle);
    arg2.NPut(0.0);
  } else {
    bool f_logical = false;
    if (mp2nle.IsExpression()) {
      f_logical
          = expr_info_.at(mp2nle.GetExprIndex()).IsLogical();
    }
    if (!f_logical) {      // Wrap as (expr != 0)
      auto argw = aw.OPut2(nl::NE);
      argw.EPut(mp2nle);
      argw.NPut(0.0);
    } else
      aw.EPut(mp2nle);
  }
}

template <class ColSizeWriter>
void MP2NLModelAPI::FeedColumnSizes(ColSizeWriter& csw) {
  if (WantColumnSizes())
    for (int i=0; i < var_lbs_.size()-1; ++i)        // use old ordering
      csw.Write(mark_data_.col_sizes_orig_[ GetOldVarIndex( i ) ]);
}


/** Initial primal guesses.
   *
   *  Implementation: write all meaningfuls entries (incl. zeros.)
   *      if (ini_guess.size()) {
   *        auto ig = igw.MakeVectorWriter(ini_guess.size());
   *        for (size_t i=0; i<ini_guess.size(); ++i)
   *          ig.Write(ini_guess[i].index_, ini_guess[i].value_);
   *      }
   */
template <class IGWriter>
void MP2NLModelAPI::FeedInitialGuesses(IGWriter& igw) {
  auto x0 = GetNLSolver().GetCallbacks()->GetInitialGuesses();
  if (x0.size()) {
    auto ig = igw.MakeVectorWriter(x0.size());
    for (size_t i=0; i<x0.size(); ++i) {
      ig.Write(GetNewVarIndex(x0[i].first), x0[i].second);
    }
  }
}

/** Initial dual guesses. */
template <class IDGWriter>
void MP2NLModelAPI::FeedInitialDualGuesses(IDGWriter& igw) {
  auto y0 = GetNLSolver().GetCallbacks()->GetInitialDualGuesses();
  if (y0.size()) {
    auto ig = igw.MakeVectorWriter(y0.size());
    for (size_t i=0; i<y0.size(); ++i) {
      ig.Write(GetNewAlgConIndex( i ), y0[i]);
    }
  }
}


template <class SuffixWriterFactory>
void MP2NLModelAPI::FeedSuffixes(SuffixWriterFactory& swf) {
  PrepareSOSSuffixes();
  FeedOriginalSuffixes(swf);
  // Custom suffixes: SOS
  Feed1Suffix(suf_sosno_, swf);
  Feed1Suffix(suf_ref_, swf);
}

void MP2NLModelAPI::PrepareSOSSuffixes() {
  suf_sosno_ = MP2NLModelSuffix
      {"sosno", suf::VAR, "", { std::vector<double>(var_lbs_.size(), 0.0) }};
  suf_ref_ = MP2NLModelSuffix
      {"ref", suf::VAR | suf::FLOAT, "",
                              { std::vector<double>(var_lbs_.size(), 0.0) }};
  int sosno = 0;
  for (const auto& sos: sos_info_) {
    ++sosno;
    if (StaticItemTypeID::ID_SOS1Constraint == sos.GetStaticTypeID())
      PrepareSOSConstraint(*(const SOS1Constraint*)sos.GetPItem(), sosno);
    else {
      assert(StaticItemTypeID::ID_SOS2Constraint == sos.GetStaticTypeID());
      PrepareSOSConstraint(*(const SOS2Constraint*)sos.GetPItem(), -sosno);
    }
  }
}

template <int SOSType>
void MP2NLModelAPI::PrepareSOSConstraint(
    const SOS_1or2_Constraint<SOSType>& sos, int sosno) {
  assert(sosno);
  assert(sosno<0 == (2==SOSType));
  for (size_t i=0; i<sos.get_vars().size(); ++i) {
    auto v = sos.get_vars()[i];
    assert(!suf_sosno_.values_[0][v]);   // single .sosno for all SOS
    suf_sosno_.values_[0][v] = sosno;
    suf_ref_.values_[0][v] = sos.get_weights()[i];
  }
}


template <class SuffixWriterFactory>
void MP2NLModelAPI::FeedOriginalSuffixes(SuffixWriterFactory& swf) {
  auto p_qc = p_nlsi_->GetCallbacks();
  auto suffixnames = p_qc->GetSuffixNames();
  for (const auto& sufname: suffixnames) {
    if ("sosno" != sufname && "ref" != sufname) {
      auto modelsuffix = p_qc->GetModelSuffix(sufname);
      Feed1Suffix(modelsuffix, swf);
    }
  }
}

template <class SuffixWriterFactory>
void MP2NLModelAPI::Feed1Suffix(
    const MP2NLModelSuffix& modelsuffix, SuffixWriterFactory& swf) {
  auto write1suf = [this](auto& sw, int kind, const auto& vals) {
    for (size_t i=0; i<vals.size(); ++i) {
      if (auto val = vals[i]) {
        int i0=i;
        if (suf::VAR==kind)
          i0 = this->GetNewVarIndex(i);
        else if (suf::CON==kind && i<alg_con_info_.size())
          i0 = this->GetNewAlgConIndex(i);
        sw.Write(i0, val);
      }
    }
  };
  for (int kind=0; kind<(int)modelsuffix.values_.size(); ++kind) {
    const auto& vals = modelsuffix.values_[kind];
    if (vals.size()) {                 // even all-0 suffixes
      auto nnz = std::count_if(vals.begin(), vals.end(),
                               [](auto n){ return bool(n); });
      if (modelsuffix.flags_ & suf::FLOAT) {
        auto sw = swf.StartDblSuffix(modelsuffix.name_.c_str(),
                                     kind | suf::FLOAT, nnz);
        write1suf(sw, kind, vals);
      } else {
        auto sw = swf.StartIntSuffix(modelsuffix.name_.c_str(),
                                     kind, nnz);
        write1suf(sw, kind, vals);
      }
    }
  }
}

template <class ColNameWriter>
void MP2NLModelAPI::FeedRowAndObjNames(ColNameWriter& wrt) {
	auto has_name0 = [](const auto& infos) {
		if (infos.size()) {
				const auto nm = infos.front().GetDispatcher().GetName(
					infos.front().GetPItem());
				return (nm && std::strlen(nm));     // remove padding?
		}
		return false;
  };
  if ((has_name0(alg_con_info_)
       || has_name0(log_con_info_) || has_name0(obj_info_)) && wrt) {
		auto write_names = [&wrt](const auto& infos) {
      for (size_t i=0; i<infos.size(); ++i) {
        const auto* nm
            = infos[i].GetDispatcher().GetName(infos[i].GetPItem());
        wrt << (nm ? nm : "..");
      }
    };
    for (size_t i=0; i<alg_con_info_.size(); ++i) {
      auto i0 = GetOldAlgConIndex(i);
      const auto* nm
					= alg_con_info_[i0].GetDispatcher().GetName(
						alg_con_info_[i0].GetPItem());
      wrt << (nm ? nm : "..");
    }
    write_names(log_con_info_);
    write_names(obj_info_);
  }
}

template <class ColNameWriter>
void MP2NLModelAPI::FeedColNames(ColNameWriter& wrt) {
  if (var_names_.size() && wrt) {
    assert(var_names_.size() == var_lbs_.size());
    for (size_t i=0; i<var_names_.size(); ++i)
      wrt << var_names_[ GetOldVarIndex(i) ];
  }
}


/// Implement NLSolverIntf
class MP2NLSolverImpl
    : public MP2NLSolverIntf,
      public SOLHandler {
public:
  /// Construct
  MP2NLSolverImpl(MP2NLModelAPI& mapi) : mapi_(mapi) { }

  /// Solve
  void Solve(const char* solver, const char* sopts) override {
    if_solve_attempted_ = true;
    nlsol_.SetFileStub(mapi_.GetFileStub());
    if_solve_ok_ = nlsol_.Solve(mapi_, *this, solver, sopts);
  }

  /// AMPL solve result code
  int GetSolveResult() const override {
    return (if_solve_ok_ || !if_solve_attempted_)
        ? sresult_ : 500;             // 500: failure
  }

  /// Solve result message or error message
  const char* GetSolveMessage() const override {
    if (if_solve_ok_)
      solve_message_final_ = solve_message_;
    else {
      solve_message_final_ = nlsol_.GetErrorMessage();
      if (solve_message_.size()) {
        solve_message_final_ +=
            "\nOriginal solve message: ";
        solve_message_final_ += solve_message_;
      }
    }
    return solve_message_final_.c_str();
  }

  /// Number of backspaces to print
  /// if printing the solve message right here,
  /// or skip so many symbols first.
  int GetSolveMessageNbs() const override { return nbs_; }

  /// Stub file used
  const char* GetFileStub() const override
  { return nlsol_.GetFileStub().c_str(); }

  /// Objno used
  int GetObjnoUsed() const override { return objno_used_; }

  /// Primal solution
  ArrayRef<double> GetX() const override { return primals_; }

  /// Dual solution
  ArrayRef<double> GetY() const override { return duals_; }

  /// Suffix names
  std::set<std::string> GetSuffixNames() override {
    std::set<std::string> result;
    for (const auto& suf: suf_map_)
      result.insert(suf.first);
    return result;
  }

  /// Get model suffix with given name
  const MP2NLModelSuffix& GetModelSuffix(
      const std::string& name) override {
    return suf_map_.at(name);
  }

  /// Num alg cons
  int GetNumAlgCons() const override {
    return n_alg_con_;
  }

  /////////////////////////////////////////////////////////////////////////
  //////////////////////// SOLHandler implementation //////////////////////
  /////////////////////////////////////////////////////////////////////////
public:
  /** The NLHeader used to write the NL file. */
  NLHeader Header() const { return mapi_.Header(); }

  /** Receive solve message.
   *  The message always ends with '\n'.
   *
   *  @param nbs: number of backspaces
   *  in the original solve message.
   *  So many characters should be skipped
   *  from the message if printed straightaway.
   *  AMPL solver drivers can supply the message
   *  with initial backspaces to indicate
   *  that so many characters should be skipped
   *  when printing. For example, if the driver prints
   *  MINOS 5.51:
   *  and exits, and the message starts with that again,
   *  this part should be skipped.
   */
  void OnSolveMessage(const char* s, int nbs) {
    solve_message_ = s;
    nbs_ = nbs;
  }

  /**
   * Can be ignored by external systems.
   * @return non-zero to stop solution input.
   */
  int OnAMPLOptions(const AMPLOptions& ) {
    suf_map_.clear();              // clear for new solution
    if_suf_data_registered_ = false;
    return 0;
  }

  /**
   * Dual values for algebraic constraints,
   * if provided in the solution.
   * Number of values <= NumAlgCons().
   * Implementation:
   *
   *   duals.reserve(rd.Size());
   *   while (rd.Size())
   *     duals.push_back(rd.ReadNext());
   */
  template <class VecReader>
  void OnDualSolution(VecReader& rd) {
    duals_.clear();
    if (int nac_sol = rd.Size()) {
      auto n_alg_cons = Header().num_algebraic_cons;
      if (nac_sol > n_alg_cons) {
        mapi_.GetEnv().AddWarning(
            "MP2NL_subsolver_solution_more_alg_cons",
            "The subsolver reported more duals than algebraic constraints");
      } else if (nac_sol < n_alg_cons) {
        mapi_.GetEnv().AddWarning(
            "MP2NL_subsolver_solution_fewer_alg_cons",
            "The subsolver reported fewer duals than algebraic constraints");
      }
      duals_.resize(nac_sol);
      int j=0;
      for ( ; rd.Size(); ++j ) {
        int j0 = mapi_.GetOldAlgConIndex(j);
        assert(j0>=0 && j0 < n_alg_cons);
        duals_[j0] = rd.ReadNext();
      }
    }
  }

  /**
   * Variable values, if provided.
   * Number of values <= NumVars().
   */
  template <class VecReader>
  void OnPrimalSolution(VecReader& rd) {
    primals_.clear();
    if (int nv_sol = rd.Size()) {
      auto n_vars = Header().num_vars;
      if (nv_sol > n_vars) {
        mapi_.GetEnv().AddWarning(
            "MP2NL_subsolver_solution_more_vars",
            "The subsolver reported more variables");
      } else if (nv_sol < n_vars) {
        mapi_.GetEnv().AddWarning(
            "MP2NL_subsolver_solution_fewer_vars",
            "The subsolver reported fewer variables");
      }
      primals_.resize(n_vars);
      int j=0;
      for ( ; rd.Size(); ++j ) {
        int j0 = mapi_.GetOldVarIndex(j);
        assert(j0>=0 && j0 < n_vars);
        primals_[j0] = rd.ReadNext();
      }
    }
  }

  /**
   * Receive notification of the objective index
   * used by the driver (solver option 'objno'-1).
   */
  void OnObjno(int objno) { objno_used_ = objno; }

  /**
   * Receive notification of the solve code.
   * Solve result codes docu:
   * https://mp.ampl.com/features-guide.html#solve-result-codes
   */
  void OnSolveCode(int sr) { sresult_ = sr; }

  /**
   * OnIntSuffix().
   *
   * For constraints, can include values for
   * logical constraints (after algebraic.)
   * Sparse representation - can be empty
   * (i.e., all values zero.)
   *
   * const auto& si = sr.SufInfo();
   * int kind = si.Kind();
   * int nmax = nitems_max[kind & 3];
   * const std::string& name = si.Name();
   * const std::string& table = si.Table();
   * while (sr.Size()) {
   *   std::pair<int, int> val = sr.ReadNext();
   *   if (val.first<0 || val.first>=nmax) {
   *     sr.SetError(NLW2_SOLRead_Bad_Suffix,
   *       "bad suffix element index");
   *     return;
   *   }
   *   suf[val.first] = val.second;
   * }
   * if (NLW2_SOLRead_OK == sr.ReadResult())    // Can check
   *   RegisterSuffix(kind, name, table, suf);
   */
  template <class SuffixReader>
  void OnIntSuffix(SuffixReader& sr) {
    OnSuffix(sr);
  }

  /**
   * Same as OnIntSuffix(), but
   * sr.ReadNext() returns pair<int, double>
   */
  template <class SuffixReader>
  void OnDblSuffix(SuffixReader& sr) {
    OnSuffix(sr);
  }


protected:
  template <class SuffixReader>
  void OnSuffix(SuffixReader& sr) {
    if (!if_suf_data_registered_) {
      if_suf_data_registered_ = true;
      RegisterSuffixData();
    }
    const auto& si = sr.SufInfo();
    int kind = si.Kind();
    int nmax = nitems_[kind & 3];
    const std::string& name = si.Name();
    const std::string& table = si.Table();
    auto& modelsuf = suf_map_[name];
    modelsuf.name_ = name;
    // set to FLOAT if at least one
    modelsuf.flags_ |= (kind & suf::FLOAT);
    if (modelsuf.table_.size() < table.size())
      modelsuf.table_ = table;
    auto& suf = modelsuf.values_.at(kind & 3);
    suf.clear();
    suf.resize(nitems_[kind & 3]);
    auto n_alg_cons = Header().num_algebraic_cons;
    while (sr.Size()) {
      auto sparse_entry = sr.ReadNext();
      if (sparse_entry.first<0 || sparse_entry.first>=nmax) {
        sr.SetError(NLW2_SOLRead_Bad_Suffix,
                    fmt::format(
          "bad suffix({}, kind={}) element index: {}, should be in [{}, {})",
                        name, kind, sparse_entry.first, 0, nmax));
        return;
      }
      int i0 = sparse_entry.first;
      if (0 == (kind & 3))                   // variable suffix
        i0 = mapi_.GetOldVarIndex(i0);
      else
        if (1 == (kind & 3) && i0 < n_alg_cons)
          i0 = mapi_.GetOldAlgConIndex(i0);
      suf[i0] = sparse_entry.second;
    }
  }

  /// Register some data for suffix reporting
  void RegisterSuffixData() {
    const auto& hdr = Header();

    nitems_[0] = hdr.num_vars;
    nitems_[1] = hdr.num_algebraic_cons + hdr.num_logical_cons;
    nitems_[2] = hdr.num_objs;
    nitems_[3] = 1;            // N problems

    n_alg_con_ = hdr.num_algebraic_cons;
  }


private:
  MP2NLModelAPI& mapi_;

  mp::NLUtils utils_;

  mp::NLSolver nlsol_ { &utils_ };

  bool if_suf_data_registered_ {};
  std::array<int, 4> nitems_ {};
  int n_alg_con_ {};               // For Backend to split alg + log cons

  std::unordered_map<std::string, MP2NLModelSuffix> suf_map_;

  /// Solution
  bool if_solve_attempted_ {};
  bool if_solve_ok_ {};     // return NLSolver's error message if not
  std::string solve_message_;
  mutable std::string solve_message_final_;
  int nbs_ {};
  std::vector<double> duals_,
      primals_;
  int objno_used_ {-1},
      sresult_ {-1};
};

template <class Lambda>
void MP2NLModelAPI::DispatchStaticItem(
    const ItemInfo& info, Lambda lambda) {
  switch (info.GetStaticTypeID()) {
  case StaticItemTypeID::ID_None:
    MP_RAISE("MP2NL static item type not specified");
  case StaticItemTypeID::ID_LinearObjective:
    lambda(*(const LinearObjective*)info.GetPItem());
    break;
  default:
    MP_RAISE("MP2NL: unknown static item type");
  }
        /*
      ID_NLObjective,

      ID_LinConRange,
      ID_LinConLE,
      ID_LinConEQ,
      ID_LinConGE,
      ID_NLConstraint,
      ID_NLAssignLE,
      ID_NLAssignEQ,
      ID_NLAssignGE,
      ID_NLLogical,
      ID_NLEquivalence,
      ID_NLImpl,
      ID_NLRimpl,
      ID_IndicatorConstraintLinLE,
      ID_IndicatorConstraintLinEQ,
      ID_IndicatorConstraintLinGE,
      ID_SOS1Constraint,
      ID_SOS2Constraint,
      ID_NLComplementarity  */
}

void MP2NLModelAPI::CreateInterfaces() {
  p_nls_ = std::make_unique<MP2NLSolverImpl>(*this);
  this->MP2NLCommonInfo::p_nlsi_ = p_nls_.get();
}

} // namespace mp
