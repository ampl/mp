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
  obj::Type sense_, sense_true_;
  LinTerms lt_;
  std::string name_;
public:
  /// Construct
  template <class CoefVec=std::initializer_list<double>,
            class VarVec=std::initializer_list<int> >
  LinearObjective(obj::Type s, CoefVec&& c, VarVec&& v,
                  std::string nm = {}) noexcept :
      sense_(s), sense_true_(s),
    lt_(std::forward<CoefVec>(c), std::forward<VarVec>(v)),
    name_(std::move(nm)){ }
  /// Get original sense
  obj::Type obj_sense() const { return sense_; }
  /// Get true sense (in multiobj mode)
  obj::Type obj_sense_true() const { return sense_true_; }
  /// Set original sense
  void set_sense(obj::Type s) { sense_ = s; }
  /// Set true sense
  void set_sense_true(obj::Type s) { sense_true_ = s; }
  /// Get lin terms, const
  const LinTerms& GetLinTerms() const { return lt_; }
  /// Get lin terms
  LinTerms& GetLinTerms() { return lt_; }
  /// Get N terms
  int num_terms() const { assert(check()); return (int)lt_.size(); }
  /// Validate
  bool check() const { return lt_.check(); }
  /// Coefs vector
  ArrayRef<double> coefs() const { return lt_.coefs(); }
  /// Var vector
  ArrayRef<int> vars() const { return lt_.vars(); }
  /// Name
  const char* name() const { return name_.c_str(); }
  /// Set name
  void set_name(std::string nm) { name_ = std::move(nm); }

  /// Testing API
  bool operator==(const LinearObjective& lc) const {
    return sense_==lc.sense_ && lt_==lc.lt_;
  }
};


/// NL objective.
/// We have no proper transformation mechanism for objectives,
/// so mixing them all.
class NLObjective : public LinearObjective {
  int expr_ {-1};
public:
  /// Construct
  NLObjective(LinearObjective&& lc, int expr=-1) :
      LinearObjective(std::move(lc)), expr_(expr) { }

  /// Has expression term?
  bool HasExpr() const { return expr_>=0; }

  /// Expression index.
  /// @note ModelAPI should call self.HasExpression()
  ///   and self.GetExpression()
  ///   to obtain the expression term.
  int ExprIndex() const { assert(HasExpr()); return expr_; }

  /// Set expression index
  void SetExprIndex(int ei) { expr_ = ei; }
};


/// Quadratic objective.
/// @note Should have no expression
///   but we check HasExpr() to see if it is an NLObjective really
class QuadraticObjective : public NLObjective {
  QuadTerms qt_;
public:
  /// Construct
  QuadraticObjective(LinearObjective&& lc, QuadTerms&& qt) :
    NLObjective(std::move(lc)), qt_(std::move(qt)) { sort_qp_terms(); }

  /// Get QP terms, const
  const QuadTerms& GetQPTerms() const { return qt_; }
  /// Get QP terms
  QuadTerms& GetQPTerms() { return qt_; }

  /// Sort QP terms
  void sort_qp_terms() {
    qt_.sort_terms();
  }

  /// Testing API
  bool operator==(const QuadraticObjective& qc) const {
    return LinearObjective::operator==(qc) && qt_==qc.qt_;
  }
};

/// Objective argument visitor
inline void VisitArguments(
    const QuadraticObjective& obj, std::function<void (int)> argv) {
  VisitArgumentsOnce(obj.GetLinTerms(), obj.GetQPTerms(), argv);
  // not for Expr - it is done as reformulation
}



/// Write objective
void WriteModelItem(fmt::MemoryWriter& wrt, const QuadraticObjective& obj,
                    ItemNamer& vnam);

/// Compute value of an objective.
template <class VarVec>
double ComputeValue(
    const QuadraticObjective& obj, const VarVec& x) {
  double linp = obj.GetLinTerms().ComputeValue(x);
  double qp = obj.GetQPTerms().ComputeValue(x);
  double nlp = (obj.HasExpr() ? x[obj.ExprIndex()] : 0.0);
  return linp + qp + nlp;
}

} // namespace mp

#endif // STD_OBJ_H
