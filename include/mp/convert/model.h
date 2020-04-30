#ifndef MODEL_H
#define MODEL_H

#include <algorithm>
#include <limits>
#include <memory>
#include <vector>

#include "mp/problem.h"

namespace mp {

/// class Model extends Problem to store custom constraints
template <class Allocator>
class BasicModel : public BasicProblem<Allocator> {
  using BaseClass = BasicProblem<Allocator>;
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

  double fix(int v) const {
    if (!is_fixed(v))
      throw std::logic_error("Variable is not fixed");
    return this->var(v).lb();
  }

  template <class VarArray>
  var::Type common_type(const VarArray& va) const {
    auto type = var::Type::INTEGER;
    for (auto v: va) {
      auto vv = this->var(v);
      if (var::Type::INTEGER!=vv.type() && (!is_fixed(v) || !is_integer_value(fix(v)))) {
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
    this->PushStandardMPItemsTo(backend);
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
