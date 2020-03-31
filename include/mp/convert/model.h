#ifndef MODEL_H
#define MODEL_H

#include <memory>

#include "mp/problem.h"
#include "mp/convert/constraint.h"

namespace mp {

/// class Model extends Problem to store custom constraints
template <class Converter, class Backend>
class CustomModel : public Problem {
  using TConstraint = Constraint<Converter, Backend>;
  using PConstraint = std::unique_ptr<TConstraint>;

  std::vector<PConstraint> custom_constr_;
};

} // namespace mp

#endif // MODEL_H
