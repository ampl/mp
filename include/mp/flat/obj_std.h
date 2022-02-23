#ifndef STD_OBJ_H
#define STD_OBJ_H

#include <vector>

#include "mp/common.h"
#include "mp/flat/expr_quadratic.h"

namespace mp {

class LinearObjective {
  obj::Type sense_;
  std::vector<double> coefs_;
  std::vector<int> vars_;
public:
  template <class CoefVec=std::initializer_list<double>,
            class VarVec=std::initializer_list<int> >
  LinearObjective(obj::Type s, CoefVec&& c, VarVec&& v) noexcept :
    sense_(s),
    coefs_(std::forward<CoefVec>(c)), vars_(std::forward<VarVec>(v)) { }
  obj::Type obj_sense() const { return sense_; }
  int num_terms() const { assert(check()); return (int)vars_.size(); }
  bool check() const { return coefs_.size()==vars_.size(); }
  const std::vector<double>& coefs() const { return coefs_; }
  const std::vector<int>& vars() const { return vars_; }

  /// Testing API
  bool operator==(const LinearObjective& lc) const {
    return sense_==lc.sense_ && coefs_==lc.coefs_ && vars_==lc.vars_;
  }
};

class QuadraticObjective : public LinearObjective {
  QuadTerms qt_;
public:
  QuadraticObjective(LinearObjective&& lc, QuadTerms&& qt) noexcept :
    LinearObjective(std::move(lc)), qt_(std::move(qt)) { sort_qp_terms(); }

  const QuadTerms& GetQPTerms() const { return qt_; }

  void sort_qp_terms() {
    qt_.sort_terms();
  }

  /// Testing API
  bool operator==(const QuadraticObjective& qc) const {
    return LinearObjective::operator==(qc) && qt_==qc.qt_;
  }
};

} // namespace mp

#endif // STD_OBJ_H
