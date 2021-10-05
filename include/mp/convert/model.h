#ifndef MODEL_H
#define MODEL_H

#include <algorithm>
#include <limits>
#include <memory>
#include <vector>

#include "mp/problem.h"
#include "mp/convert/quad_expr.h"
#include "mp/convert/constraint_keeper.h"

namespace mp {

class BasicConstraintKeeper;

/// Storing extra info in BasicModel's items
struct FlatConverterModelExraItemInfo : public DefaultExtraItemInfo {
  struct AlgConExtraInfo {
    QuadTerms qt_;
    AlgConExtraInfo() { }
    template <class QT>
    AlgConExtraInfo(QT&& qt) noexcept : qt_(std::forward<QT>(qt)) { }
  };
  struct ObjExtraInfo {
    double obj_const_term_ = 0.0;
    QuadTerms qt_;
    ObjExtraInfo() { }
    template <class QT>
    ObjExtraInfo(double c, QT&& qt) noexcept :
      obj_const_term_(c), qt_(std::forward<QT>(qt)) { }
  };
};

struct DefaultFlatConverterModelParams : public BasicProblemParams<> {
  using ExtraItemInfo = FlatConverterModelExraItemInfo;
};

/// class Model extends Problem to store custom constraints
template <class ModelParams = DefaultFlatConverterModelParams >
class BasicModel : public BasicProblem<ModelParams> {
  using BaseClass = BasicProblem<ModelParams>;
public:
  using Params = ModelParams;
protected:
  using PConstraintKeeper = std::unique_ptr<BasicConstraintKeeper>;

private:
  std::vector<PConstraintKeeper> custom_constr_;


public:

  /** Returns the number of custom constraints. */
  int num_custom_cons() const { return static_cast<int>(custom_constr_.size()); }

  /** Returns custom constraint i */
  BasicConstraintKeeper* custom_con(int i) const {
    internal::CheckIndex(i, num_custom_cons());
    return custom_constr_[i].get();
  }

  /// Add custom constraint. Takes ownership
  void AddConstraint(BasicConstraintKeeper* pbc) {
    PConstraintKeeper pc;
    pc.reset(pbc);
    AddConstraint(std::move(pc));
  }
  void AddConstraint(PConstraintKeeper pc) {
    custom_constr_.push_back(std::move(pc));
  }


  /////////////////////////////// UTILITIES //////////////////////////////////
  static constexpr double PlusInfinity() { return std::numeric_limits<double>::infinity(); }
  static constexpr double MinusInfinity() { return -std::numeric_limits<double>::infinity(); }
  template <class Num>
  static bool is_integer_value(Num n) { return std::floor(n)==std::ceil(n); }

  template <class VarArray>
  double lb_array(const VarArray& va) const {
    double result = PlusInfinity();
    for (auto v: va) {
      result = std::min( result, this->var(v).lb() );
    }
    return result;
  }
  template <class VarArray>
  double lb_max_array(const VarArray& va) const {
    double result = MinusInfinity();
    for (auto v: va) {
      result = std::max( result, this->var(v).lb() );
    }
    return result;
  }
  template <class VarArray>
  double ub_array(const VarArray& va) const {
    double result = MinusInfinity();
    for (auto v: va) {
      result = std::max( result, this->var(v).ub() );
    }
    return result;
  }
  template <class VarArray>
  double ub_min_array(const VarArray& va) const {
    double result = PlusInfinity();
    for (auto v: va) {
      result = std::min( result, this->var(v).ub() );
    }
    return result;
  }

  bool is_fixed(int v) const {
    auto vv = this->var(v);
    return vv.lb()==vv.ub();
  }

  double fixed_value(int v) const {
    if (!is_fixed(v))
      throw std::logic_error("Variable is not fixed");
    return this->var(v).lb();
  }

  void narrow_var_bounds(int v, double lb, double ub) {
    auto vv = this->var(v);
    vv.set_lb(std::max(vv.lb(), lb));
    vv.set_ub(std::min(vv.ub(), ub));
    if (vv.lb()>vv.ub())
      throw std::logic_error("infeasibility: empty variable domain");
  }

  bool is_integer_var(int v) const {
    auto vv = this->var(v);
    return var::Type::INTEGER==vv.type();
  }

  bool is_binary_var(int v) const {
    auto vv = this->var(v);
    return 0.0==vv.lb() && 1.0==vv.ub() && var::Type::INTEGER==vv.type();
  }

  template <class VarArray=std::initializer_list<int> >
  var::Type common_type(const VarArray& va) const {
    auto type = var::Type::INTEGER;
    for (auto v: va) {
      auto vv = this->var(v);
      if (var::Type::INTEGER!=vv.type() && (!is_fixed(v) || !is_integer_value(fixed_value(v)))) {
        type = var::Type::CONTINUOUS;
        break;
      }
    }
    return type;
  }


  //////////////////////////////// EXPORT INSTANCE TO A BACKEND ///////////////////////////////
  ///
  /// Pushing the whole instance to a backend or converter.
  /// A responsible backend should handle all essential items
  template <class Backend>
  void PushModelTo(Backend& backend) const {
    this->InitProblemModificationPhase(backend);
    this->PushVariablesTo(backend);
    this->PushObjectivesTo(backend);
    this->PushAlgebraicConstraintsTo(backend);
    PushCustomConstraintsTo(backend);
    this->FinishProblemModificationPhase(backend);
  }

protected:
  template <class Backend>
  void PushCustomConstraintsTo(Backend& backend) const {
    if (int n_ccons = num_custom_cons()) {
      for (int i = 0; i < n_ccons; ++i) {
        const auto* pConstraint = custom_con(i);
        if (!pConstraint->IsRemoved())
          pConstraint->AddToBackend(backend);
      }
    }
  }


};

} // namespace mp

#endif // MODEL_H
