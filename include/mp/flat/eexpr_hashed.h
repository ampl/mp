#ifndef EEXPR_HASHED_H
#define EEXPR_HASHED_H

#include <unordered_map>

#include "mp/flat/eexpr.h"
#include "mp/utils-hash-stream.h"

namespace mp {

/// linear or quadratic var/coef list
template <class VarTuple>
class EExpr_Hashed_CoefVarList {
public:
};

/// Hashed EExpr for temporary storage
class EExpr_Hashed {
protected:
  /// QP key
  using QPKey = std::array<int, 2>;

public:
  /// Get constant
  double GetConstant() const { return const_term_; }
  /// Set constant
  void SetConstant(double ct) { const_term_ = ct; }
  /// Add linear term
  void AddTerm(double c, int v) {
    auto& val = map_l_[v];
    val += c;
    if (!val)                 // don't need 0 coefs
      map_l_.erase(v);
  }
  /// Add QP term
  void AddTerm(double c, int v1, int v2) {
    QPKey qpk {(v1<=v2) ? QPKey{v1, v2} : QPKey{v2, v1}};
    auto& val = map_q_[qpk];
    val += c;
    if (!val)                 // don't need 0 coefs
      map_q_.erase(qpk);
  }
  /// Add LinTerms
  void Add(const LinTerms& lt) {
    for (int i=0; i<lt.size(); ++i)
      AddTerm(lt.coef(i), lt.var(i));
  }
  /// Add QPTerms
  void Add(const QuadTerms& qt) {
    for (int i=0; i<qt.size(); ++i)
      AddTerm(qt.coef(i), qt.var1(i), qt.var2(i));
  }
  /// Add an EExpr
  void Add(const EExpr& ee) {
    SetConstant(GetConstant() + ee.constant_term());
    Add(ee.GetLinTerms());
    Add(ee.GetQPTerms());
  }
  /// Produce EExpr.
  /// @note no sorting but unique.
  EExpr ToEExpr() const {
    EExpr result;
    result.constant_term(GetConstant());
    result.GetLinTerms().reserve(map_l_.size());
    for (const auto v: map_l_)
      result.GetLinTerms().add_term(v.second, v.first);
    result.GetQPTerms().reserve(map_q_.size());
    for (const auto v: map_q_)
      result.GetQPTerms().add_term(v.second, v.first[0], v.first[1]);
    return result;

  }

private:
  double const_term_ {};
  std::unordered_map<int, double> map_l_;
  std::unordered_map<std::array<int, 2>, double> map_q_;
};

}  // namespace mp

#endif // EEXPR_HASHED_H
