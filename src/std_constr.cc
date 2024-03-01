#include <map>

#include "mp/format.h"
#include "mp/util-json-write.hpp"

#include "mp/flat/expr_quadratic.h"
#include "mp/flat/model_info.hpp"

namespace mp {

void LinTerms::sort_terms(bool force_sort) {
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
  for (int i=0; i<size(); ++i)
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

template <>
void WriteJSON(JSONW jw, const QuadTerms& qt) {
  jw["coefs"] = qt.coefs();
  jw["vars1"] = qt.vars1();
  jw["vars2"] = qt.vars2();
}

template <>
void WriteJSON(JSONW jw, const LinTerms& qt) {
  jw["coefs"] = qt.coefs();
  jw["vars"] = qt.vars();
}

template <>
void WriteJSON(JSONW jw, const QuadAndLinTerms& qlt) {
  WriteJSON(jw["qp_terms"], qlt.GetQPTerms());
  WriteJSON(jw["lin_terms"], qlt.GetLinTerms());
}

} // namespace mp
