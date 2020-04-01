#ifndef MAXIMUM_H
#define MAXIMUM_H

#include <vector>

#include "mp/convert/constraint.h"
#include "mp/convert/expr2constraint.h"

namespace mp {

class MaximumConstraint : public BasicConstraint {
  std::vector<EExpr> args_;
  int result_var_;
public:
  const std::vector<EExpr>& GetArguments() const { return args_; }
  int GetResultVar() const { return result_var_; }
  MaximumConstraint(std::vector<EExpr>&& a, int r) : args_(a), result_var_(r) { }
};

} // namespace mp

#endif // MAXIMUM_H
