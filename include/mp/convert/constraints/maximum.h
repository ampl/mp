#ifndef MAXIMUM_H
#define MAXIMUM_H

#include <vector>

#include "mp/convert/constraint.h"
#include "mp/convert/expr2constraint.h"

namespace mp {

template <class Converter, class Backend>
class MaximumConstraint : public BasicConstraint {
  std::vector<EExpr> args_;
  int result_var_;
};

} // namespace mp

#endif // MAXIMUM_H
