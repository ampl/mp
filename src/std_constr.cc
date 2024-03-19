#include <map>
#include <cfloat>
#include <cassert>

#include "mp/format.h"
#include "mp/util-json-write.hpp"
#include "mp/common.h"

#include "mp/flat/expr_quadratic.h"
#include "mp/flat/obj_std.h"
#include "mp/flat/model_info.hpp"

namespace mp {

/// @todo Keep consistent with the \a ConstraintGroups enum.
static const char* const congroup_names[]
= {
 "Default",
 "All",
 "Algebraic",
 "Linear",
 "Quadratic",
 "Conic",
 "General",
 "Piecewiselinear",
 "SOS",
 "SOS1",
 "SOS2",
 "Logical"
};

const char* ConGroupName(int cg) {
  assert(0<=cg
         && cg<int(sizeof (congroup_names)/sizeof(congroup_names[0])));
  return congroup_names[cg];
}

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

template <class Writer>
void WriteVar(Writer& pr, const char* name,
              double lb, double ub, var::Type ty) {
  assert(*name);
  pr << "var " << name;
  if (!lb && 1.0==ub && var::INTEGER==ty)
    pr << " binary";
  else if (lb==ub)
    pr << " = " << lb;
  else {
    if (lb > -DBL_MAX)
    pr << " >=" << lb;
    if (ub < DBL_MAX)
    pr << " <=" << ub;
    if (var::INTEGER == ty)
    pr << " integer";
  }
}

void WriteModelItem(fmt::MemoryWriter& wrt, const LinTerms& lt,
                    const std::vector<std::string>& vnam) {
  for (int i=0; i<(int)lt.size(); ++i) {
    if (i) {
      wrt << (lt.coef(i)>=0.0 ? " + " : " - ");
    }
    wrt << std::fabs(lt.coef(i)) << '*' << vnam.at(lt.var(i));
  }
}

void WriteModelItem(fmt::MemoryWriter& wrt, const QuadTerms& qt,
                    const std::vector<std::string>& vnam) {
  for (int i=0; i<(int)qt.size(); ++i) {
    if (i) {
      wrt << (qt.coef(i)>=0.0 ? " + " : " - ");
    }
    wrt << std::fabs(qt.coef(i))
        << '*' << vnam.at(qt.var1(i))
        << '*' << vnam.at(qt.var2(i));
  }
}

void WriteModelItem(fmt::MemoryWriter& wrt, const QuadAndLinTerms& qlt,
                    const std::vector<std::string>& vnam) {
  WriteModelItem(wrt, qlt.GetLinTerms(), vnam);
  if (qlt.GetQPTerms().size()) {
    if (qlt.GetLinTerms().size())
      wrt << " + ";
    wrt << '(';
    WriteModelItem(wrt, qlt.GetQPTerms(), vnam);
    wrt << ')';
  }
}

void WriteModelItem(fmt::MemoryWriter& wrt, const QuadraticObjective& obj,
                    const std::vector<std::string>& vnam) {
  wrt << (obj.obj_sense() ? "maximize " : "minimize ");
  assert(obj.name() && *obj.name());
  wrt << obj.name() << ": ";
  WriteModelItem(wrt, obj.GetLinTerms(), vnam);
  if (obj.GetQPTerms().size()) {
    if (obj.GetLinTerms().size())
      wrt << " + ";
    wrt << '(';
    WriteModelItem(wrt, obj.GetQPTerms(), vnam);
    wrt << ')';
  }
}

// Generate
template
void WriteVar(fmt::MemoryWriter& pr, const char* name,
              double lb, double ub, var::Type ty);

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
