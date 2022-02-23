#include <map>

#include "mp/flat/expr_quadratic.h"

namespace mp {

void LinTerms::sort_terms(bool force_sort)  {
  std::map<int, double> var_coef_map;
  for (size_t i=0; i<size(); ++i)
    if (0.0!=std::fabs(coefs_[i]))
      var_coef_map[vars_[i]] += coefs_[i];
  if (force_sort ||                    // force sorting for tests
      var_coef_map.size() < size()) {
    coefs_.clear();
    vars_.clear();
    for (const auto& vc: var_coef_map) {
      if (0.0!=std::fabs(vc.second)) {         // Need tolerance?
        coefs_.push_back(vc.second);
        vars_.push_back(vc.first);
      }
    }
  }
}


void QuadTerms::sort_terms()  {
  auto sort_pair = [](int a, int b) {
    return a<b ? std::pair<int, int>(a, b) : std::pair<int, int>(b, a);
  };
  std::map<std::pair<int, int>, double> var_coef_map;
  for (size_t i=0; i<size(); ++i)
    if (0.0!=std::fabs(coefs_[i]))
      var_coef_map[sort_pair(vars1_[i], vars2_[i])] += coefs_[i];
  if (true) {
    coefs_.clear();
    vars1_.clear();
    vars2_.clear();
    for (const auto& vc: var_coef_map) {
      if (0.0!=std::fabs(vc.second))         // Need tolerance?
        add_term(vc.second, vc.first.first, vc.first.second);
    }
  }
}

/// TODO modular printing of algebraic constraints

//void print(std::ostream& os) const {
//  os << lb_ << " <= ";
//  for (int i=0; i<nnz(); ++i) {
//    os << coefs_[i] << "*[" << vars_[i] << ']';
//    if (i<nnz()-1)
//      os << " + ";
//  }
//  os << " <= " << ub_;
//}

  //  void print(std::ostream& os) const {
  //    os << lb() << " <= ";
  //    for (int i=0; i<nnz(); ++i) {
  //      os << coefs()[i] << "*[" << vars()[i] << ']';
  //      if (i<nnz()-1)
  //        os << " + ";
  //    }
  //    for (int i=0; i<qt_.num_terms(); ++i) {
  //      os << " + "
  //         << qt_.coef(i) << "*[" << qt_.var1(i) << "]*[" << qt_.var2(i) << "]";
  //    }
  //    os << " <= " << ub();
  //  }


} // namespace mp
