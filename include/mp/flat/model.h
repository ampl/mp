#ifndef MODEL_H
#define MODEL_H

#include <memory>
#include <vector>
#include <algorithm>
#include <cmath>

#include "mp/flat/obj_std.h"
#include "mp/flat/constraints_std.h"
#include "mp/flat/constraint_keeper.h"

namespace mp {

class BasicConstraintKeeper;

struct DefaultFlatModelParams {
  using Var = int;
  static constexpr Var VoidVar() { return -1; }
};

/// Class BasicFlatModel stores vars, objs, custom constraints
/// to be used internally in a FlatConverter
template <class ModelParams = DefaultFlatModelParams>
class BasicFlatModel : public ConstraintManager {
public:
  using Params = ModelParams;
  using Var = typename Params::Var;
  static constexpr Var VoidVar() { return Params::VoidVar(); }

  using VarArray = std::vector<Var>;
  using VarBndVec = std::vector<double>;
  using VarTypeVec = std::vector<var::Type>;

  ///////////////////////////// VARIABLES ////////////////////////////////
  /// Add variable, return its index
  Var AddVar__basic(double lb=MinusInf(), double ub=Inf(),
             var::Type type=var::CONTINUOUS) {
    assert(check_vars());
    var_lb_.push_back(lb);
    var_ub_.push_back(ub);
    var_type_.push_back(type);
    return var_type_.size()-1;
  }

  /// Add several variables
  void AddVars__basic(const VarBndVec& lbs, const VarBndVec& ubs,
               const VarTypeVec& types) {
    assert(check_vars());
    var_lb_.insert(var_lb_.end(), lbs.begin(), lbs.end());
    var_ub_.insert(var_ub_.end(), ubs.begin(), ubs.end());
    var_type_.insert(var_type_.end(), types.begin(), types.end());
    assert(check_vars());
  }

  int num_vars() const
  { assert(check_vars()); return (int)var_lb_.size(); }

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
    return 0.0<=lb(v) && 1.0>=ub(v) && is_integer_var(v);
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
protected:
  using ObjList = std::vector<QuadraticObjective>; // TODO just an item
  const ObjList& get_objectives() const { return objs_; }
  ObjList& get_objectives() { return objs_; }
  int num_objs() const { return (int)objs_.size(); }
  const QuadraticObjective& get_obj(int i) const
  { return get_objectives().at(i); }

public:
  void AddObjective(QuadraticObjective&& obj)
  { get_objectives().push_back(std::move(obj)); }


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
  static bool is_integer_value(Num n) { return std::floor(n)==std::ceil(n); }


  //////////////////////////////// EXPORT INSTANCE TO A BACKEND ///////////////////////////////
  ///
  /// Pushing the whole instance to a backend or converter.
  /// A responsible backend should handle all essential items
  template <class Backend>
  void PushModelTo(Backend& backend) const {
    backend.InitProblemModificationPhase();
    PushVariablesTo(backend);
    PushObjectivesTo(backend);
    PushCustomConstraintsTo(backend);
    backend.FinishProblemModificationPhase();
  }

protected:
  template <class Backend>
  void PushVariablesTo(Backend& backend) const {
    backend.AddVariables( { var_lb_, var_ub_, var_type_ } );
  }

  template <class Backend>
  void PushObjectivesTo(Backend& backend) const {
    if (int n_objs = num_objs()) {
      for (int i = 0; i < n_objs; ++i) {
        const auto& obj = get_obj(i);
        if (obj.GetQPTerms().size())  // TODO make objectives just items
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

private:
  /// Variables' bounds
  VarBndVec var_lb_, var_ub_;
  /// Variables' types
  VarTypeVec var_type_;

  /// Objectives
  /// TODO storing QuadraticObjective now, make it just an item
  ObjList objs_;


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

#endif // MODEL_H
