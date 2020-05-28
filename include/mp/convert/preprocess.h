#ifndef PREPROCESS_H
#define PREPROCESS_H

#include <limits>
#include <cmath>

#include "mp/problem.h"
#include "mp/convert/std_constr.h"

namespace mp {

template <class Num>
bool is_integer(Num n) { return std::floor(n)==std::ceil(n); }

template <class FuncConstraint>
struct PreprocessInfo {
  double lb_=-std::numeric_limits<double>::max(),
    ub_=std::numeric_limits<double>::max();
  var::Type type_=var::CONTINUOUS;
  int result_var_ = -1;

  PreprocessInfo() { }
  PreprocessInfo(double l, double u, var::Type t) : lb_(l), ub_(u), type_(t) { }
  double lb() const { return lb_; }
  double ub() const { return ub_; }
  void narrow_result_bounds(double l, double u) {
    lb_ = std::max(lb_, l);
    ub_ = std::min(ub_, u);
  }
  var::Type get_result_type() const { return type_; }
  void set_result_type(var::Type t) { type_=t; }
  bool is_constant() const { return lb_==ub_; }
  bool is_result_var_known() const { return result_var_>=0; }
  void set_result_var(int r) { result_var_ = r; }
  int get_result_var() const {
    assert(result_var_>=0);
    return result_var_;
  }
};


/// Typical preprocess info
using PreprocessInfoStd = PreprocessInfo<int>;

template <class Model>
PreprocessInfoStd ComputeBoundsAndType(Model& model, AffineExpr& ae) {
  PreprocessInfoStd result;
  ComputeBoundsAndType(model, ae, result);
  return result;
}

template <class Model>
void ComputeBoundsAndType(Model& model, AffineExpr& ae, PreprocessInfoStd& result) {
  result.lb_ = result.ub_ = ae.constant_term();           // TODO reuse bounds if supplied
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


} // namespace mp

#endif // PREPROCESS_H
