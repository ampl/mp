#ifndef PREPRO_ARGS_H
#define PREPRO_ARGS_H

#include <limits>
#include <cmath>

#include "mp/problem.h"
#include "mp/convert/std_constr.h"

namespace mp {

template <class Num>
bool is_integer(Num n) { return std::floor(n)==std::ceil(n); }

template <class FuncConstraint>
struct BasicPreprocessInfo {
  double lb_=-std::numeric_limits<double>::max(),
    ub_=std::numeric_limits<double>::max();
  var::Type type_=var::CONTINUOUS;
  int result_var_ = -1;

  BasicPreprocessInfo() { }
  BasicPreprocessInfo(double l, double u, var::Type t) : lb_(l), ub_(u), type_(t) { }
  bool is_constant() const { return lb_==ub_; }
  bool is_result_var_known() const { return result_var_>=0; }
  void set_result_var(int r) { result_var_ = r; }
  int get_result_var() const {
    assert(result_var_>=0);
    return result_var_;
  }
};

/// Typical preprocess info
using PreprocessInfo = BasicPreprocessInfo<int>;

/// Default arguments prepro
template <class Model, class Constraint, class PreproInfo>
void PreprocessConstraint(
    Model& , Constraint&, PreproInfo& ) {
  // ...
}

template <class Model>
void ComputeBoundsAndType(Model& model, AffineExpr& ae, PreprocessInfo& result) {
  result.lb_ = result.ub_ = ae.constant_term();
  result.type_ = var::INTEGER;
  for (const auto& term: ae) {
    auto v = model.var(term.var_index());
    if (term.coef() >= 0.0) {
      result.lb_ += term.coef() * v.lb();
      result.ub_ += term.coef() * v.ub();
    } else {
      result.lb_ += term.coef() * v.ub();
      result.ub_ += term.coef() * v.lb();
    }
    if (var::INTEGER!=v.type() || !is_integer(term.coef()))
      result.type_=var::CONTINUOUS;
  }
}

/// Preprocess minimum
template <class Model>
void PreprocessConstraint(
    Model& m, MinimumConstraint& c, BasicPreprocessInfo<MinimumConstraint>& prepro) {
  prepro.lb_ = m.lb_array(c.GetArguments());
  prepro.ub_ = m.ub_min_array(c.GetArguments());
  prepro.type_ = m.common_type(c.GetArguments());
}

template <class Model>
void PreprocessConstraint(
    Model& m, MaximumConstraint& c, BasicPreprocessInfo<MaximumConstraint>& prepro) {
  prepro.lb_ = m.lb_max_array(c.GetArguments());
  prepro.ub_ = m.ub_array(c.GetArguments());
  prepro.type_ = m.common_type(c.GetArguments());
}

} // namespace mp

#endif // PREPRO_ARGS_H
