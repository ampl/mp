#ifndef CONVERTER_MODEL_H
#define CONVERTER_MODEL_H

#include <memory>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cfloat>

#include "mp/utils-vec.h"

#include "mp/flat/obj_std.h"
#include "mp/flat/constr_std.h"
#include "mp/flat/constr_keeper.h"
#include "mp/flat/converter_model_base.h"
#include "mp/flat/model_info.h"

namespace mp {

class BasicConstraintKeeper;

struct DefaultFlatModelParams {
  using Var = int;
  static constexpr Var VoidVar() { return -1; }
};

/// Class FlatModel stores vars, objs, custom constraints
/// to be used internally in a FlatConverter
template <class ModelParams = DefaultFlatModelParams>
class FlatModel
    : public ConstraintManager, public BasicFlatModel {
public:
  using Params = ModelParams;
  using Var = typename Params::Var;
  static constexpr Var VoidVar() { return Params::VoidVar(); }

  using VarArray = std::vector<Var>;
  using VarTypeVec = std::vector<var::Type>;
  using VarNameVec = std::vector<const char*>;

  using ConstraintManager::GetFileAppender;

  /// Destruct
  virtual ~FlatModel() { }

  /// Converter info
  virtual const ConverterInfo* GetConverterInfo() const = 0;

  ///////////////////////////// VARIABLES ////////////////////////////////
  /// Add variable, return its index
  Var AddVar__basic(double lb=MinusInf(), double ub=Inf(),
             var::Type type=var::CONTINUOUS) {
    assert(check_vars());
    var_lb_.push_back(lb);
    var_ub_.push_back(ub);
    var_lb_best_.push_back(lb);
    var_ub_best_.push_back(ub);
    var_type_.push_back(type);
    ExportVars(var_type_.size()-1, {lb}, {ub}, {type});
    return var_type_.size()-1;
  }

  /// Add several variables.
  /// This is only done once when we add original variables.
  void AddVars__basic(const VarBndVec& lbs, const VarBndVec& ubs,
               const VarTypeVec& types) {
    assert(check_vars());
    var_lb_.insert(var_lb_.end(), lbs.begin(), lbs.end());
    var_ub_.insert(var_ub_.end(), ubs.begin(), ubs.end());
    var_lb_best_.insert(var_lb_best_.end(), lbs.begin(), lbs.end());
    var_ub_best_.insert(var_ub_best_.end(), ubs.begin(), ubs.end());
    var_type_.insert(var_type_.end(), types.begin(), types.end());
    assert(check_vars());
    assert(!num_vars_orig_);
    num_vars_orig_ = var_lb_.size();
    ExportVars(var_type_.size()-lbs.size(), lbs, ubs, types);
  }

protected:
  /// Export new variables
  template <class BndVec=std::array<double, 1>,
            class TypeVec=std::array<var::Type, 1> >
  void ExportVars(int i_start, const BndVec& lbs, const BndVec& ubs,
                  const TypeVec types,
                  const char* comment
                  = "Initial model information. "
                    "Can be updated later with new bounds, names, etc.")
  const {
    for (int i=0;
         GetFileAppender().IsOpen() && i<(int)lbs.size(); ++i) {
      fmt::MemoryWriter wrt;
      if (0==i && i_start==0) {
        {
          MiniJSONWriter jw(wrt);
          jw["COMMENT"] = comment;
        }
        wrt.write("\n");                      // with EOL
      }
      {
        MiniJSONWriter jw(wrt);
        int i_actual = i+i_start;
        jw["VAR_index"] = i_actual;
        {
					auto name = GetVarNamer().at(i_actual);
          jw["name"] = name;
          fmt::MemoryWriter pr;
          WriteVar(pr, name, lbs[i], ubs[i], types[i]);
          jw["printed"] = pr.c_str();
        }
        jw["bounds"]
            << (lbs[i] < -DBL_MAX ? -DBL_MAX : lbs[i])
            << (ubs[i] > DBL_MAX ? DBL_MAX : ubs[i]);
        jw["type"] = (int)types[i];
        jw["is_from_nl"] = (int)is_var_original(i_actual);
      }
      wrt.write("\n");                      // with EOL
      GetFileAppender().Append(wrt);
    }
  }

public:
  /// Set var names vector
  void AddVarNames(const std::vector<std::string>& names) {
    var_names_storage_ = names;
  }

  /// To be used only after names presolve
  const char* var_name(int i) const {
    return i<(int)var_names_storage_.size()
        ? var_names_storage_[i].c_str() : nullptr;
  }

  /// N vars
  int num_vars() const
  { assert(check_vars()); return (int)var_lb_.size(); }

  /// Not auxiliary?
  bool is_var_original(int i) const
  { return i<num_vars_orig_; }

  /// Actually, lb_best
  double lb(Var v) const {
    assert(0<=v && v<num_vars());
    return var_lb_best_[v];
  }

  /// Actually, ub_best
  double ub(Var v) const {
    assert(0<=v && v<num_vars());
    return var_ub_best_[v];
  }

  /// Hard LB.
  /// We should not inline the init expression when this
  /// bound on the result variable is significant #237 #266.
  double lb_hard(Var v) const {
    assert(0<=v && v<num_vars());
    return var_lb_[v];
  }

  /// Hard UB.
  /// We should not inline the init expression when this
  /// bound on the result variable is significant.
  double ub_hard(Var v) const {
    assert(0<=v && v<num_vars());
    return var_ub_[v];
  }

  var::Type var_type(Var v) const {
    assert(0<=v && v<num_vars());
    return var_type_[v];
  }

  const std::vector<double>& var_lb_vec() const
  { return var_lb_; }
  const std::vector<double>& var_ub_vec() const
  { return var_ub_; }
  const std::vector<double>& var_lb_best_vec() const
  { return var_lb_best_; }
  const std::vector<double>& var_ub_best_vec() const
  { return var_ub_best_; }
  const std::vector<var::Type>& var_type_vec() const
  { return var_type_; }

  /// Best-known bounds
  template <class VarArray>
  double lb_array(const VarArray& va) const {
    double result = Inf();
    for (auto v: va) {
      result = std::min( result, lb(v) );
    }
    return result;
  }

  template <class VarArray>
  double lb_max_array(const VarArray& va) const {
    double result = MinusInf();
    for (auto v: va) {
      result = std::max( result, lb(v) );
    }
    return result;
  }

  template <class VarArray>
  double ub_array(const VarArray& va) const {
    double result = MinusInf();
    for (auto v: va) {
      result = std::max( result, ub(v) );
    }
    return result;
  }

  template <class VarArray>
  double ub_min_array(const VarArray& va) const {
    double result = Inf();
    for (auto v: va) {
      result = std::min( result, ub(v) );
    }
    return result;
  }

  /// Uses best-known bounds
  bool is_fixed(int v) const {
    return lb(v)==ub(v);
  }

  double fixed_value(int v) const {
    assert(is_fixed(v));
    return 0.5*((long double)lb(v) + ub(v));
  }

  /// Uses hard bounds
  bool is_free(int v) const {
    return
        lb_hard(v) <= MinusInf()
           && ub_hard(v) >= Inf();
  }

  bool is_integer_var(int v) const {
    return var::Type::INTEGER==var_type(v);
  }

  /// Returns true also when fixed
  bool is_binary_var(int v) const {
    return
        (0.0==lb(v) && 1.0==ub(v)    // best bounds
            && is_integer_var(v))
        || (is_fixed(v)
            && (0.0==fixed_value(v) || 1.0==fixed_value(v)));
  }

  template <class VarArray=std::initializer_list<int> >
  var::Type common_type(const VarArray& va) const {
    auto type = var::Type::INTEGER;
    for (auto v: va) {
      if (!is_integer_var(v) &&
          (!is_fixed(v) || !is_integer_value(fixed_value(v)))) {
        type = var::Type::CONTINUOUS;
        break;
      }
    }
    return type;
  }

  /// Improve lb, incl. best-known
  void set_lb(int v, double l) {
    assert(0<=v && v<num_vars());
    if (var_lb_[v]<l)
      var_lb_[v] = l;
    set_best_lb(v, l);
  }

  /// Improve ub, incl. best-known
  void set_ub(int v, double u) {
    assert(0<=v && v<num_vars());
    if (var_ub_[v]>u)
      var_ub_[v] = u;
    set_best_ub(v, u);
  }

  /// Improve best-known lb
  void set_best_lb(int v, double l) {
    assert(0<=v && v<num_vars());
    if (var_lb_best_[v]<l)
      var_lb_best_[v] = l;
  }

  /// Improve best-known ub
  void set_best_ub(int v, double u) {
    assert(0<=v && v<num_vars());
    if (var_ub_best_[v]>u)
      var_ub_best_[v] = u;
  }

  /// To be called first
  void MarkAllResultVarsAsVars() {
    var_result_.clear();
    var_result_.resize(num_vars(), true);
  }

  /// Mark as an explicit result variable
  /// @todo AutoExpand fills 'false' for new elements...
  void MarkAsResultVar(int v) {
    AutoExpand(var_result_, v, true);
  }

  /// Mark as a proper expression
  void MarkAsExpression(int v) {
    AutoExpand(var_result_, v, false);
  }

  /// Variable has expr marking?
  bool VarHasMarking(int v) const {
    assert(v>=0);
    return v < (int)var_result_.size();
  }

  /// Is the variable proper - an explicit result var of a flat con,
  /// or just a primary variable?
  /// (Otherwise, it's marked as implicit
  ///   - the init expr will be an expression)
  /// @note this is also used to prevent AMPL defined variables
  ///   from inlining in ConsiderInliningAlgExpr().
  bool IsProperVar(int v) const {
    return !VarHasMarking(v) || var_result_[v];
  }

  /// Get var proper flags
  const std::vector<bool>& GetVarProperFlags() const
  { return var_result_; }

  /// Mark as an explicit defined variable
  /// @todo AutoExpand fills 'false' for new elements...
  void MarkAsExplicitDV(int v) {
    AutoExpand(var_dv_explicit_, v, true);
  }

  bool IsExplicitDV(int v) const {
    return (int)var_dv_explicit_.size()>v
           && var_dv_explicit_[v];
  }

  /// Mark var as eliminated.
  void MarkVarAsEliminated(int v) {
    AutoExpand(var_elim_, v, true);
  }

  /// Get the flags of eliminated variables
  /// (currenlty they are fixed)
  const std::vector<bool>& GetVarElimFlags() const
  { return var_elim_; }


  /// Add function definitions
  void AddFuncDefs(std::vector<FuncDef> funcs)
  { funcs_ = std::move(funcs); }

public:
  /// N functions
  int GetNumFuncs() const {
    return (int) funcs_.size();
  }

  /// Get function \a i
  const FuncDef& GetFunction(int i) const {
    return funcs_.at(i);
  }

  ///////////////////////////// OBJECTIVES ////////////////////////////
public:
  /// List of objectives
  using ObjList = std::vector<QuadraticObjective>;
  /// Get list of objectives, const
  const ObjList& get_objectives() const { return objs_; }
  /// Get list of objectives
  ObjList& get_objectives() { return objs_; }
  /// N obj
  int num_objs() const { return (int)objs_.size(); }
  /// Skip pushing objectives?
  bool if_skip_pushing_objs() const { return if_skip_push_objs_; }
  /// Set skip pushing objs
  void set_skip_pushing_objs(bool v=true) { if_skip_push_objs_=v; }
  /// Submit best-known bounds? (cvt:pre:boundsbest: 0 no, 1 all,
  /// 2 all but smooth-nonlinear variables, 3 only those)
  int if_submit_best_known_bounds() const { return if_submit_best_bounds_; }
  /// Submit best-known bounds?
  int& if_submit_best_known_bounds() { return if_submit_best_bounds_; }
  /// Submit the best-known bounds of variable \a v?
  /// Value 4 is decided by the submitted model's class, known
  /// in PrepareVariables(); before that it answers as 0.
  bool if_submit_best_bounds_var(int v) const {
    switch (if_submit_best_bounds_) {
    case 1: return true;
    case 2: return !is_var_nl(v);
    case 3: return is_var_nl(v);
    case 4: return model_nl_known_ && !model_nl_;
    default: return false;
    }
  }
  /// Is the submitted model nonlinear (nonlinear or quadratic
  /// constraints or objectives)? Valid after PrepareVariables()
  bool is_submitted_model_nonlinear() const
  { assert(model_nl_known_); return model_nl_; }
  /// Has a nonlinear or quadratic objective?
  bool HasNonlinearObjective() const {
    for (int i=0; i<num_objs(); ++i)
      if (get_obj(i).HasExpr() || get_obj(i).GetQPTerms().size())
        return true;
    return false;
  }
  /// Mark \a v as the result of a smooth nonlinear expression
  void MarkVarNL(int v) { SetVarFlag(var_nl_, v); }
  /// Is \a v the result of a smooth nonlinear expression?
  bool is_var_nl(int v) const { return GetVarFlag(var_nl_, v); }
  /// Make fixed vars continuous?
  int if_make_fixed_vars_continuous() const { return if_make_fixed_vars_continuous_; }
  /// Make fixed vars continuous?
  int& if_make_fixed_vars_continuous() { return if_make_fixed_vars_continuous_; }
  /// Submit implied-integer auxiliary variables as integer?
  int if_submit_aux_int() const { return if_submit_aux_int_; }
  /// Submit implied-integer auxiliary variables as integer?
  int& if_submit_aux_int() { return if_submit_aux_int_; }
  /// Mark \a v as an auxiliary variable whose integrality
  /// is implied by its defining constraint (cvt:pre:auxint)
  void MarkVarImpliedInt(int v) { SetVarFlag(var_implied_int_, v); }
  /// Mark \a v as a variable the solver needs declared
  /// with its type, e.g., an indicator binary
  void MarkVarTypeRequired(int v) { SetVarFlag(var_type_required_, v); }
  /// Is \a v marked implied-integer?
  bool is_var_implied_int(int v) const
  { return GetVarFlag(var_implied_int_, v); }
  /// Is \a v marked type-required?
  bool is_var_type_required(int v) const
  { return GetVarFlag(var_type_required_, v); }
  /// Submitted type (var_type() before PrepareVariables())
  var::Type var_type_subm(int v) const {
    assert(0<=v && v<num_vars());
    return var_type_subm_.size() ? var_type_subm_[v] : var_type_[v];
  }
  /// Submitted types (var_type_vec() before PrepareVariables())
  const VarTypeVec& var_type_subm_vec() const
  { return var_type_subm_.size() ? var_type_subm_ : var_type_; }
  /// Number of implied-integer variables submitted as continuous
  int num_aux_int_relaxed() const { return n_aux_int_relaxed_; }
  /// Get obj [i]
  const QuadraticObjective& get_obj(int i) const
  { return get_objectives().at(i); }
  /// Has a QP objective?
  bool HasQPObjective() const {
    for (const auto& obj: get_objectives())
      if (!obj.GetQPTerms().empty())
        return true;
    return false;
  }
  /// Add an objective
  void AddObjective(QuadraticObjective&& obj) {
    get_objectives().push_back(std::move(obj));
    ExportObjective(num_objs()-1, get_objectives().back());
  }

protected:
  void ExportObjective(
      int i_obj, const QuadraticObjective& obj) const {
    if (GetFileAppender().IsOpen()) {
      fmt::MemoryWriter wrt;
      {
        MiniJSONWriter jw(wrt);
        jw["OBJECTIVE_index"] = i_obj;
        if (obj.name() && *obj.name()) {
          jw["name"] = obj.name();
          fmt::MemoryWriter pr;
          WriteModelItem(pr, obj, GetVarNamer());
          jw["printed"] = pr.c_str();
        }
        jw["sense"] = (int)obj.obj_sense();
        WriteJSON(jw["qp_terms"], obj.GetQPTerms());
        WriteJSON(jw["lin_terms"], obj.GetLinTerms());
      }
      wrt.write("\n");                     // EOL
      GetFileAppender().Append(wrt);
    }
  }

  ///////////////////////////// FLAT CONSTRAINTS ////////////////////////////
protected:

public:

  /////////////////////////////// UTILITIES //////////////////////////////////
public:
  static constexpr double Inf()
  { return INFINITY; }

  static constexpr double MinusInf()
  { return -INFINITY; }

  template <class Num>
  static bool is_integer_value(Num n)
  { return std::floor(n)==n; }

  /// Provide variable lower bounds
  const VarBndVec& GetVarLBs() const { return var_lb_; }
  /// Provide variable upper bounds
  const VarBndVec& GetVarUBs() const { return var_ub_; }

  /// Get dimension of SDP var \a sdpi
  int GetSDPVarDim(int sdpi) const
  { return SDP_var_dims_.at(sdpi); }

  /// Add an SDP variable of dimension \a n
  /// @return its index
  int AddSDPVar(int dim) {
    SDP_var_dims_.push_back(dim);
    return (int) SDP_var_dims_.size()-1;
  }

  /// Describes to which SDP var a usual var belongs
  struct SDPVarElementInfo   {
    int sdp_var_index_ {-1};
    int sdp_var_row_ {-1};
    int sdp_var_col_ {-1};
  };

  /// Get SDP var info for scalar variable \a i
  const SDPVarElementInfo& GetSDPVarInfo(int i) const {
    assert(i>=0);
    if (i>=(int)var_SDP_var_info_.size())
      return sdp_var_elem_dummy_;
    return var_SDP_var_info_[i];
  }

  /// Set SDP var info for scalar variable \a i
  void SetSDPVarInfo(int i, SDPVarElementInfo vei) {
    if (var_SDP_var_info_.size() < var_lb_.size())
      var_SDP_var_info_.resize(var_lb_.size());
    assert(i>=0 && i<(int)var_lb_.size());
    var_SDP_var_info_[i] = vei;
  }

  /// Model state
  struct ModelState {
    /// If we emulate MO
    bool ifMOEmulator_ {false};
    /// If !=0, current MOE objective
    const QuadraticObjective* pObj_{nullptr};
  };

  /// Create default model state
  static ModelState MakeDefaultModelState()
  { return {}; }

  /////////////// EXPORT THE INSTANCE TO A ModelAPI //////////////
  ///
  /// Pushing the whole instance to the mapi.
  template <class ModelAPI>
  void PushModelTo(ModelAPI& mapi) {
    ModelState model_state {if_skip_pushing_objs()};
    CreateFlatModelInfo(mapi, model_state);     // If no previous
                                   // stats output, create in any case
    mapi.PassFlatModelInfo(GetModelInfo());

    mapi.InitProblemModificationPhase(GetModelInfo());
    mapi.AddFuncDefs(funcs_);
    PushVariablesTo(mapi);
    PushCustomConstraintsTo(mapi);
    if (!if_skip_pushing_objs())
      PushObjectivesTo(mapi);
    mapi.FinishProblemModificationPhase();
  }

protected:  
  void CreateFlatModelInfo(
      const BasicFlatModelAPI& mapi,
      const ModelState& model_state = {}) {
    if (!pfmi_)       // Create after construction of the Converter
      pfmi_ = mp::CreateFlatModelInfo(GetConverterInfo());
    FillVarStats(GetModelInfoWrt());
    FillObjStats(GetModelInfoWrt(), model_state);
    FillConstraintCounters(mapi, *GetModelInfoWrt());
  }

  void FillVarStats(FlatModelInfo* pfmi) const {
    int n=0;
    FlatModelInfo::VarInfo vi{0,0,0,0,0,0};
    const auto& vlb = var_lb_subm_.size() ?  // when computed
                          var_lb_subm_ : var_lb_;
    const auto& vub = var_ub_subm_.size() ?
                          var_ub_subm_ : var_ub_;
    for (auto i=var_lb_.size(); i--; ) {
      bool is_int = var::Type::CONTINUOUS != var_type_subm(i);
      if (vlb[i] < vub[i] && is_int)
        ++n;
      int index_add = 3*(!is_var_original(i));
      ++vi[index_add];
      if (is_int) {
        bool is_bin = !vlb[i] && 1.0==vub[i];
        ++vi[index_add + 1 + is_bin];
      }
    }
    pfmi->SetNumUnfixedIntVars(n);
    pfmi->SetVarInfo(vi);
  }

public:
  void FillObjStats(
      FlatModelInfo* pfmi, const ModelState& model_state = {}) const {
    FlatModelInfo::ObjInfo oi{0,0,0};
    auto IncrObjCount = [&oi](const QuadraticObjective& obj) {
      if (obj.HasExpr()) {
        ++ oi[2];
      } else if (obj.GetQPTerms().size())
        ++ oi[1];
      else
        ++ oi[0];
    };
    if (model_state.ifMOEmulator_) {
      if (auto* pObj = model_state.pObj_) {
        IncrObjCount(*pObj);
      }  // else, leave all 0
    } else   // Count all normally
    if (int n_objs = num_objs()) {
      for (int i = 0; i < n_objs; ++i) {
        const auto& obj = get_obj(i);
        IncrObjCount(obj);
      }
    }
    pfmi->SetObjInfo(oi);
  }

public:
  void PrepareVariables() const {
    assert(check_vars());
    model_nl_ = HasFinalNonlinearConstraints() || HasNonlinearObjective();
    model_nl_known_ = true;
    PrepareVarTypes();
    var_lb_subm_ = var_lb_;
    var_ub_subm_ = var_ub_;
    for (auto i=var_lb_subm_.size(); i--; )
      if (if_submit_best_bounds_var(i)) {
        var_lb_subm_[i] = var_lb_best_[i];
        var_ub_subm_[i] = var_ub_best_[i];
      }
    // Fix 'eliminated' variables - no proper deletion
    for (auto i=std::min(var_elim_.size(), var_lb_subm_.size());
         i--; ) {
      if (var_elim_[i]) {
        if (var_lb_subm_[i] > -1e20)
          var_ub_subm_[i] = var_lb_subm_[i];
        else if (var_ub_subm_[i] < 1e20)
          var_lb_subm_[i] = var_ub_subm_[i];
        else
          var_lb_subm_[i] = var_ub_subm_[i] = 0.0;
      }
    }
    // Fix variables with deduced equal (or lb>ub) bounds.
    for (auto i=var_lb_subm_.size(); i--; ) {
      if (var_lb_best_[i] >= var_ub_best_[i]) {
        var_lb_subm_[i] = var_lb_best_[i];
        var_ub_subm_[i] = var_ub_best_[i];
      }
    }
    // Make fixed vars continuous
    if (if_make_fixed_vars_continuous())
      for (auto i=var_lb_subm_.size(); i--; ) {
        if (var_lb_subm_[i] >= var_ub_subm_[i])
          var_type_subm_[i] = var::CONTINUOUS;   // avoid MIP classification
      }
  }

  /// Submitted types: implied-integer auxiliaries become continuous
  /// unless type-required or cvt:pre:auxint=1
  void PrepareVarTypes() const {
    var_type_subm_ = var_type_;
    n_aux_int_relaxed_ = 0;
    if (if_submit_aux_int())
      return;
    for (auto i=var_type_subm_.size(); i--; )
      if (var::CONTINUOUS != var_type_subm_[i]
          && is_var_implied_int(i) && !is_var_type_required(i)) {
        var_type_subm_[i] = var::CONTINUOUS;
        ++n_aux_int_relaxed_;
      }
  }

protected:
  template <class Backend>
  void PushVariablesTo(Backend& backend) const {
    assert(var_lb_subm_.size() == var_lb_.size()
           && var_ub_subm_.size() == var_ub_.size()
           && var_lb_.size() == var_ub_.size());
    if (var_names_storage_.size() > 0) {
      // Convert names to c-str if needed
      var_names_.reserve(var_names_storage_.size());
      for (const std::string& s : var_names_storage_)
        var_names_.push_back(s.c_str());
      backend.AddVariables({ var_lb_subm_, var_ub_subm_,
                             var_type_subm_, var_names_ });
    } else {
      backend.AddVariables({ var_lb_subm_, var_ub_subm_, var_type_subm_ });
    }
    ExportVars(0, var_lb_subm_, var_ub_subm_, var_type_subm_,
               "Updated model information.");
  }

  template <class Backend>
  void PushObjectivesTo(Backend& backend) const {
    if (int n_objs = num_objs()) {
      for (int i = 0; i < n_objs; ++i) {
        const auto& obj = get_obj(i);
        SetObjectiveTo(backend, i, obj);
        ExportObjective(i, obj);      // can change e.g., for SOCP
      }
    }
  }

  template <class Backend>
  void PushCustomConstraintsTo(Backend& backend) const {
    this->AddUnbridgedConstraintsToBackend(
        backend, GetVarNamer());
    this->LogConstraintGroups(backend);
  }


public:
  /// Model info
  const FlatModelInfo* GetModelInfo() const { return pfmi_.get(); }
  /// Model info, writable
  FlatModelInfo* GetModelInfoWrt() const { return pfmi_.get(); }

  /// Set given objective
  template <class Backend>
  void SetObjectiveTo(
      Backend& backend, int i, const QuadraticObjective& obj) const {
    if (obj.HasExpr()) {
      assert(obj.GetQPTerms().empty());     // not mixing quadratics and expr
      backend.SetNLObjective(i, obj);
    } else if (obj.GetQPTerms().size())
      backend.SetQuadraticObjective(i, obj);
    else
      backend.SetLinearObjective(i, obj);
		ExportObjective(i, obj);
  }

  /// Var namer
  ItemNamer& GetVarNamer() const { return var_namer_; }

private:
  /// Variables' "hard" bounds
  /// - solver should know them to have correct solution
  VarBndVec var_lb_, var_ub_;
  /// Variables' best-known bounds, see cvt:pre:boundsbest.
  /// By default, we don't submit such bounds e.g. for aux result vars,
  /// so solvers can eliminate them.
  VarBndVec var_lb_best_, var_ub_best_;
  /// Which bounds to submit
  int if_submit_best_bounds_ {0};
  int if_make_fixed_vars_continuous_ {1};
  /// Variables' submitted bounds - what goes to the solver
  mutable VarBndVec var_lb_subm_, var_ub_subm_;
  /// Variables' types, as known to MP
  mutable VarTypeVec var_type_;
  /// Variables' submitted types - what goes to the solver,
  /// see cvt:pre:auxint
  mutable VarTypeVec var_type_subm_;
  /// Number of implied-integer variables submitted as continuous
  mutable int n_aux_int_relaxed_ {0};
  /// Auxiliary variables whose integrality is implied
  /// by their defining constraint
  std::vector<bool> var_implied_int_;
  /// Variables which the solver needs declared with their type
  std::vector<bool> var_type_required_;
  /// Results of smooth nonlinear expressions (cvt:pre:boundsbest 2/3)
  std::vector<bool> var_nl_;
  /// Submitted model nonlinear? (cvt:pre:boundsbest 4)
  mutable bool model_nl_ {false}, model_nl_known_ {false};
  /// Submit implied-integer auxiliaries as integer?
  int if_submit_aux_int_ {0};
  /// Set a per-variable flag, growing the vector on demand
  static void SetVarFlag(std::vector<bool>& flags, int v) {
    assert(v>=0);
    if ((int)flags.size()<=v)
      flags.resize(v+1);
    flags[v] = true;
  }
  /// Get a per-variable flag
  static bool GetVarFlag(const std::vector<bool>& flags, int v)
  { return v<(int)flags.size() && flags[v]; }
  /// Whether the variable, being the result variable of a functional constraint,
  /// needs to stay a variable (vs being eliminated because the constraint
  /// is becoming an expression.)
  /// Normal variables are marked too.
  /// @note Only used in the expression interface.
  std::vector<bool> var_result_;
  /// In contrast to var_result_, this is used
  /// to mark if AMPL defined variable stays explicit.
  std::vector<bool> var_dv_explicit_;
  /// Eliminated variables.
  /// Currently they are just fixed.
  std::vector<bool> var_elim_;
  ///  Variables' names
  mutable VarNameVec var_names_;
  std::vector<std::string> var_names_storage_;
  mutable ItemNamer var_namer_ {var_names_storage_, "_svar"};
  /// Number of original NL variables
  int num_vars_orig_ {0};

  /// SDP variable dimensions
  std::vector<int> SDP_var_dims_;

  /// SDP var info for usual variables
  std::vector<SDPVarElementInfo> var_SDP_var_info_;
  /// dummy
  SDPVarElementInfo sdp_var_elem_dummy_ {};

  /// FuncDef defs
  std::vector<FuncDef> funcs_;
  /// Objectives
  ObjList objs_;
  /// Whether to skip pushing objectives
  /// (can be set by the MOManager.)
  bool if_skip_push_objs_ {false};

  /// Flat model info
  std::unique_ptr<FlatModelInfo> pfmi_;

public:
  /// Check var arrays
  bool check_vars() const {
    return var_lb_.size() == var_ub_.size() &&
           var_lb_best_.size() == var_ub_.size() &&
           var_ub_best_.size() == var_ub_.size() &&
           var_type_.size() == var_ub_.size();
  }

  void RelaxIntegrality() {
    std::fill(var_type_.begin(), var_type_.end(), var::CONTINUOUS);
  }
};

} // namespace mp

#endif // CONVERTER_MODEL_H
