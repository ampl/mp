#ifndef CONVERTER_MODEL_H
#define CONVERTER_MODEL_H

#include <memory>
#include <vector>
#include <algorithm>
#include <cmath>

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

/// Class BasicFlatModel stores vars, objs, custom constraints
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

  ///////////////////////////// VARIABLES ////////////////////////////////
  /// Add variable, return its index
  Var AddVar__basic(double lb=MinusInf(), double ub=Inf(),
             var::Type type=var::CONTINUOUS) {
    assert(check_vars());
    var_lb_.push_back(lb);
    var_ub_.push_back(ub);
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
                  const TypeVec types) {
    for (int i=0;
         GetFileAppender().IsOpen() && i<(int)lbs.size(); ++i) {
      fmt::MemoryWriter wrt;
      {
        MiniJSONWriter jw(wrt);
        jw["var_index"] = i+i_start;
        jw["bounds"] << lbs[i] << ubs[i];
        jw["type"] = (int)types[i];
        jw["is_from_nl"] = (int)is_var_original(i);
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

  double lb(Var v) const {
    assert(0<=v && v<num_vars());
    return var_lb_[v];
  }

  double ub(Var v) const {
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
  const std::vector<var::Type>& var_type_vec() const
  { return var_type_; }

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

  bool is_fixed(int v) const {
    return lb(v)==ub(v);
  }

  double fixed_value(int v) const {
    assert(is_fixed(v));
    return lb(v);
  }

  bool is_integer_var(int v) const {
    return var::Type::INTEGER==var_type(v);
  }

  /// Returns true also when fixed
  bool is_binary_var(int v) const {
    return
        (0.0==lb(v) && 1.0==ub(v) && is_integer_var(v))
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

  void set_lb(int v, double l) {
    assert(0<=v && v<num_vars());
    var_lb_[v] = l;
  }

  void set_ub(int v, double u) {
    assert(0<=v && v<num_vars());
    var_ub_[v] = u;
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
  void ExportObjective(int i_obj, const QuadraticObjective& obj) {
    if (GetFileAppender().IsOpen()) {
      fmt::MemoryWriter wrt;
      {
        MiniJSONWriter jw(wrt);
        jw["obj_index"] = i_obj;
        jw["sense"] = (int)obj.obj_sense();
        wrt.write("\"qp_terms\": ");
        WriteJSON(wrt, obj.GetQPTerms());
        wrt.write("\"lin_terms\": ");
        WriteJSON(wrt, obj.GetLinTerms());
        wrt.write(" {}\n", '}');                      // with EOL
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
  { return std::floor(n)==std::ceil(n); }

  /// Provide variable lower bounds
  const VarBndVec& GetVarLBs() const { return var_lb_; }
  /// Provide variable upper bounds
  const VarBndVec& GetVarUBs() const { return var_ub_; }


  /////////////// EXPORT THE INSTANCE TO A ModelAPI //////////////
  ///
  /// Pushing the whole instance to the mapi.
  template <class ModelAPI>
  void PushModelTo(ModelAPI& mapi) const {
    CreateFlatModelInfo(mapi);
    mapi.PassFlatModelInfo(GetModelInfo());

    mapi.InitProblemModificationPhase(GetModelInfo());
    PushVariablesTo(mapi);
    PushObjectivesTo(mapi);
    PushCustomConstraintsTo(mapi);
    mapi.FinishProblemModificationPhase();
  }

protected:  
  void CreateFlatModelInfo(const BasicFlatModelAPI& mapi) const {
    FillVarStats(GetModelInfoWrt());
    FillConstraintCounters(mapi, *GetModelInfoWrt());
  }

  void FillVarStats(FlatModelInfo* pfmi) const {
    int n=0;
    for (auto i=var_lb_.size(); i--; ) {
      if (var_lb_[i] < var_ub_[i]
          && var::Type::CONTINUOUS != var_type(i))
        ++n;
    }
    pfmi->SetNumUnfixedIntVars(n);
  }

  template <class Backend>
  void PushVariablesTo(Backend& backend) const {
    if (var_names_storage_.size() > 0) {
      // Convert names to c-str if needed
      for (const std::string& s : var_names_storage_)
        var_names_.push_back(s.c_str());
      backend.AddVariables({ var_lb_, var_ub_, var_type_, var_names_ });
    } else
      backend.AddVariables({ var_lb_, var_ub_, var_type_ });
  }

  template <class Backend>
  void PushObjectivesTo(Backend& backend) const {
    if (int n_objs = num_objs()) {
      for (int i = 0; i < n_objs; ++i) {
        const auto& obj = get_obj(i);
        if (obj.GetQPTerms().size())
          backend.SetQuadraticObjective(i, obj);
        else
          backend.SetLinearObjective(i, obj);
      }
    }
  }

  template <class Backend>
  void PushCustomConstraintsTo(Backend& backend) const {
    this->AddUnbridgedConstraintsToBackend(backend);
  }

  FlatModelInfo* GetModelInfoWrt() const { return pfmi_.get(); }


public:
  /// Model info
  const FlatModelInfo* GetModelInfo() const  { return pfmi_.get(); }


private:
  /// Variables' bounds
  VarBndVec var_lb_, var_ub_;
  /// Variables' types
  VarTypeVec var_type_;
  ///  Variables' names
  mutable VarNameVec var_names_;
  std::vector<std::string> var_names_storage_;
  /// Number of original NL variables
  int num_vars_orig_ {0};
  /// Objectives
  ObjList objs_;

  /// Flat model info
  std::unique_ptr<FlatModelInfo>
      pfmi_ {mp::CreateFlatModelInfo()};

public:
  /// Check var arrays
  bool check_vars() const {
    return var_lb_.size() == var_ub_.size() &&
        var_type_.size() == var_ub_.size();
  }

  void RelaxIntegrality() {
    std::fill(var_type_.begin(), var_type_.end(), var::CONTINUOUS);
  }
};

} // namespace mp

#endif // CONVERTER_MODEL_H
