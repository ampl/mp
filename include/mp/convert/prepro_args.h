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
struct PreprocessInfo {
  double lb_=-std::numeric_limits<double>::max(),
    ub_=std::numeric_limits<double>::max();
  var::Type type_=var::CONTINUOUS;
  int result_var_ = -1;

  PreprocessInfo() { }
  PreprocessInfo(double l, double u, var::Type t) : lb_(l), ub_(u), type_(t) { }
  void narrow_result_bounds(double l, double u) {
    lb_ = std::max(lb_, l);
    ub_ = std::min(ub_, u);
  }
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

/// Default arguments prepro
/// All parameters are 'in-out'
template <class Converter, class Constraint, class PreproInfo>
void PreprocessConstraint(
    Converter& , Constraint&, PreproInfo& ) {
  // ... do nothing by default
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

/// Preprocess minimum's arguments
template <class Converter>
void PreprocessConstraint(
    Converter& cvt, MinimumConstraint& c, PreprocessInfo<MinimumConstraint>& prepro) {
  auto& m = cvt.GetModel();
  auto& args = c.GetArguments();
  prepro.narrow_result_bounds( m.lb_array(args),
                        m.ub_min_array(args) );
  prepro.set_result_type( m.common_type(args) );
}

/// Preprocess maximum's arguments
template <class Converter>
void PreprocessConstraint(
    Converter& cvt, MaximumConstraint& c, PreprocessInfo<MaximumConstraint>& prepro) {
  auto& m = cvt.GetModel();
  auto& args = c.GetArguments();
  prepro.narrow_result_bounds( m.lb_max_array(args),
                        m.ub_array(args) );
  prepro.set_result_type( m.common_type(args) );
}

/// Preprocess EQ's arguments
template <class Converter>
void PreprocessConstraint(
    Converter& cvt, EQConstraint& c, PreprocessInfo<EQConstraint>& prepro) {
  auto& m = cvt.GetModel();
  auto& args = c.GetArguments();
  prepro.narrow_result_bounds(0.0, 1.0);
  prepro.set_result_type( var::INTEGER );
  if (m.is_fixed(args[0]) && m.is_fixed(args[1])) {
    auto res = (double)int(m.fixed_value(args[0])==m.fixed_value(args[1]));
    prepro.narrow_result_bounds(res, res);
    return;
  }
  if (m.is_fixed(args[0])) {                 // Constant on the right
    std::swap(args[0], args[1]);
  }
  if (m.is_fixed(args[1])) {                 // See if this is binary var==const
    if (m.is_binary_var(args[0])) {
      if (1.0==std::fabs(m.fixed_value(args[1])))
        prepro.set_result_var( args[0] );
      else if (0.0==m.fixed_value(args[1]))
        prepro.set_result_var( cvt.MakeComplementVar(args[0]) );
      else
        prepro.narrow_result_bounds(0.0, 0.0);    // not 0/1 value, result false
      return;
    }
  }
}

/// Preprocess NE's arguments
template <class Converter>
void PreprocessConstraint(
    Converter& cvt, LEConstraint& c, PreprocessInfo<LEConstraint>& prepro) {
  prepro.narrow_result_bounds(0.0, 1.0);
  prepro.set_result_type( var::INTEGER );
  auto& m = cvt.GetModel();
  auto& args = c.GetArguments();
  // TODO special cases
}

/// Preprocess Disjunction's arguments
template <class Converter>
void PreprocessConstraint(
    Converter& cvt, DisjunctionConstraint& c, PreprocessInfo<DisjunctionConstraint>& prepro) {
  prepro.narrow_result_bounds(0.0, 1.0);
  prepro.set_result_type( var::INTEGER );
  auto& m = cvt.GetModel();
  auto& args = c.GetArguments();
  // TODO special cases
}

/// Preprocess Not's arguments
template <class Converter>
void PreprocessConstraint(
    Converter& cvt, NotConstraint& c, PreprocessInfo<NotConstraint>& prepro) {
  prepro.narrow_result_bounds(0.0, 1.0);
  prepro.set_result_type( var::INTEGER );
  auto& m = cvt.GetModel();
  auto& args = c.GetArguments();
  // TODO special cases
}

} // namespace mp

#endif // PREPRO_ARGS_H
