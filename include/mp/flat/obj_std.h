#ifndef STD_OBJ_H
#define STD_OBJ_H

#include <string>
#include <vector>

#include "mp/common.h"
#include "mp/flat/expr_quadratic.h"
#include "mp/flat/constr_eval.h"

namespace mp {

/// Linear objective incl. sense and name
class LinearObjective {
  obj::Type sense_;
  LinTerms lt_;
  std::string name_;
public:
  /// Construct
  template <class CoefVec=std::initializer_list<double>,
            class VarVec=std::initializer_list<int> >
  LinearObjective(obj::Type s, CoefVec&& c, VarVec&& v,
                  std::string nm = {}) noexcept :
    sense_(s),
    lt_(std::forward<CoefVec>(c), std::forward<VarVec>(v)),
    name_(std::move(nm)){ }
  /// Get sense
  obj::Type obj_sense() const { return sense_; }
  /// Get lin terms
  const LinTerms& GetLinTerms() const { return lt_; }
  /// Get N terms
  int num_terms() const { assert(check()); return lt_.size(); }
  /// Validate
  bool check() const { return lt_.check(); }
  /// Coefs vector
  const std::vector<double>& coefs() const { return lt_.coefs(); }
  /// Var vector
  const std::vector<int>& vars() const { return lt_.vars(); }
  /// Name
  const char* name() const { return name_.c_str(); }
  /// Set name
  void set_name(std::string nm) { name_ = std::move(nm); }

  /// Testing API
  bool operator==(const LinearObjective& lc) const {
    return sense_==lc.sense_ && lt_==lc.lt_;
  }
};

/// Quadragtic objective
class QuadraticObjective : public LinearObjective {
  QuadTerms qt_;
public:
  /// Construct
  QuadraticObjective(LinearObjective&& lc, QuadTerms&& qt) :
    LinearObjective(std::move(lc)), qt_(std::move(qt)) { sort_qp_terms(); }

  /// Get QP terms
  const QuadTerms& GetQPTerms() const { return qt_; }

  /// Sort QP terms
  void sort_qp_terms() {
    qt_.sort_terms();
  }

  /// Testing API
  bool operator==(const QuadraticObjective& qc) const {
    return LinearObjective::operator==(qc) && qt_==qc.qt_;
  }
};

/// Compute value of an objective.
template <class VarVec>
double ComputeValue(
    const QuadraticObjective& obj, const VarVec& x) {
  return
      obj.GetLinTerms().ComputeValue(x)
      + obj.GetQPTerms().ComputeValue(x);
}

} // namespace mp

#endif // STD_OBJ_H
