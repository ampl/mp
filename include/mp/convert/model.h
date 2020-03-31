#ifndef MODEL_H
#define MODEL_H

#include <memory>

#include "mp/problem.h"
#include "mp/convert/constraint.h"

namespace mp {

/// class Model extends Problem to store custom constraints
template <class Allocator>
class BasicModel : public BasicProblem<Allocator> {
  using PConstraint = std::unique_ptr<BasicConstraint>;

  std::vector<PConstraint> custom_constr_;
};

} // namespace mp

#endif // MODEL_H
