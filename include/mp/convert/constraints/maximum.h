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
public:
  const std::vector<EExpr>& GetArguments() const { return args_; }
  int GetResultVar() const { return result_var_; }
  MaximumConstraint(std::vector<EExpr>&& a, int r) : args_(a), result_var_(r) { }
  void ConvertWith(BasicConstraintConverter& cvt) override {
    throw std::runtime_error("Cannot redefine Maximum yet");
  }
  void AddToBackend(BasicConstraintAdder& be) const override {
    try {
      static_cast<Backend&>(be).AddConstraint(*this);
    } catch (const std::exception& exc) {
      throw std::logic_error(std::string("MaximumConstraint: ") + exc.what());
    }
  }
};

} // namespace mp

#endif // MAXIMUM_H
