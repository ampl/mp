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

template <class Writer, class Vec>
void WriteJSONVec(Writer& wrt, const Vec& vec) {
  wrt.write("{} ", '[');
  for (size_t i=0; i<vec.size(); ++i) {
    if (i)
      wrt.write(", ");
    wrt.write("{}", vec[i]);
  }
  wrt.write(" {}", ']');
}

template <>
void WriteJSON(fmt::MemoryWriter& wrt, const QuadTerms& qt) {
  wrt.write("{} ", '{');
  wrt.write("\"coefs\": ");
  WriteJSONVec(wrt, qt.coefs());
  wrt.write(", \"vars1\": ");
  WriteJSONVec(wrt, qt.vars1());
  wrt.write(", \"vars2\": ");
  WriteJSONVec(wrt, qt.vars2());
  wrt.write(" {}", '}');
}

template <>
void WriteJSON(fmt::MemoryWriter& wrt, const LinTerms& qt) {
  wrt.write("{} ", '{');
  wrt.write("\"coefs\": ");
  WriteJSONVec(wrt, qt.coefs());
  wrt.write(", \"vars\": ");
  WriteJSONVec(wrt, qt.vars());
  wrt.write(" {}", '}');
}

} // namespace mp
