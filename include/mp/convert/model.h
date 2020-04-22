#ifndef MODEL_H
#define MODEL_H

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
  void AddConstraint(PConstraintKeeper&& pc) {
    custom_constr_.push_back(std::move(pc));
  }

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
