#include <map>
#include <cfloat>
#include <cmath>
#include <cassert>

#include "mp/format.h"
#include "mp/common.h"
#include "mp/utils-vec.h"
#include "mp/utils-string.h"

#include "mp/flat/eexpr.h"
#include "mp/flat/obj_std.h"
#include "mp/flat/constr_keeper.h"
#include "mp/flat/model_info.hpp"

#include "mp/flat/bucketaccum.hpp"

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
        "Nonlinear",
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


//////////////////////////// SORTING /////////////////////////////

template <class Vec>
void LinTerms::fold_into(Vec& vec) {
  vec.resize(size());
  for (size_t i=0; i<size(); ++i)
    vec.push_back({ var(i), coef(i) });
}

template <class Vec>
void LinTerms::unfold_from(const Vec& vec) {
  clear();
  reserve(vec.size());
  for (const auto& v: vec)
    add_term(v.second, v.first);
}

bool LinTerms::is_sorted() const {
  if (size()) {           // empty expr ==> "sorted" ???
    if (!coefs_.back())       // last coef == 0
      return false;
    for (auto i = size()-1; (i--)>0; ) {
      if (vars_[i] >= vars_[i+1]
          || !coefs_[i])      // coef == 0
        return false;
    }
  }
  return true;            // Check emptyness elsewhere? @todo
}

/// Sort, leave only unique keys with non-0 values
/// @param vec: some_vector< std::pair<Key, Value> >
template <class Vec>
void SortUnifyNon0(Vec& vec) {
  assert(vec.size() >= 1);
  if (vec.size() < 1)
    return;
  // Sort by (Key, abs(Value)) for numerics
  auto Cmp = [](const auto& a, const auto& b) {
    return a.first<b.first ? true
                             : a.first==b.first
               ? std::fabs(a.second) < std::fabs(b.second)
        : false;
  };
  std::sort(vec.begin(), vec.end(), Cmp);
  // Merge same keys, leaving non-0 values
  auto i2=vec.begin(), i1=i2;
  while (++i2!=vec.end()) {
    if (i1->first == i2->first)
      i1->second += i2->second;
    else {
      if (i1->second)
        ++i1;
      *i1 = *i2;
    }
  }
  vec.resize(i1-vec.begin()
             +bool(i1->second));   // last target element non-0
}

void LinTerms::sort_terms(bool force_sort) {
  if (1==size()) {
    if (!coef(0))
      clear();
  } else {
    if (1<size() && !is_sorted()) {
      SmallVec< std::pair<int, double>, 256 > fold;
      fold_into(fold);
      SortUnifyNon0(fold);
      assert(fold.size() <= size());
      if (force_sort || fold.size() < size())
        unfold_from(fold);
    }
  }
  assert(!force_sort || is_sorted());
}


template <class Vec>
void QuadTerms::fold_into(Vec& vec) {
  vec.resize(size());
  for (size_t i=0; i<size(); ++i) {
    auto key = std::pair {var1(i), var2(i)};
    if (key.first > key.second)       // index pair ordered
      std::swap(key.first, key.second);
    vec.push_back({ key, coef(i) });
  }
}

template <class Vec>
void QuadTerms::unfold_from(const Vec& vec) {
  clear();
  reserve(vec.size());
  for (const auto& v: vec)
    add_term(v.second, v.first.first, v.first.second);
}

bool QuadTerms::is_sorted() const {
  if (size()) {           // empty expr ==> "sorted" ???
    if (!coefs_.back())       // last coef == 0
      return false;
    if (vars1_.back() > vars2_.back())   // v1>v2
      return false;
    for (auto i = size()-1; (i--)>0; ) {
      if (vars1_[i] > vars2_[i]          // v1>v2
          || vars1_[i] > vars1_[i+1]
          || (vars1_[i]==vars1_[i+1]
              && vars2_[i]>=vars2_[i+1]) // also when ==
          || !coefs_[i])                 // coef == 0
        return false;
    }
  }
  return true;            // Check emptyness elsewhere? @todo
}

void QuadTerms::sort_terms()  {
  if (1==size()) {
    if (!coef(0))
      clear();
    else
      if (var1(0) > var2(0))              // order index pair
        std::swap(vars1_[0], vars2_[0]);
  } else {
    if (1<size() && !is_sorted()) {
      SmallVec< std::pair<std::pair<int, int>, double>, 256 >
          fold;
      fold_into(fold);
      SortUnifyNon0(fold);
      assert(fold.size() <= size());
      unfold_from(fold);
    }
  }
  assert(is_sorted());
}


template <class Terms>
Terms MergeSorted(const Terms& t1, const Terms& t2) {
  assert(t1.is_sorted());
  assert(t2.is_sorted());
  Terms result;
  result.reserve(t1.size() + t2.size());

  size_t i1=0, i2=0;
  while (true) {
    if (i1<t1.size() && i2<t2.size()) {
      auto t1i = t1.IndexValue(i1);
      auto t2i = t2.IndexValue(i2);
      if (t1i.first < t2i.first) {
        result.add_index_value(t1i);
        ++i1;
      } else if (t1i.first > t2i.first) {
        result.add_index_value(t2i);
        ++i2;
      } else {
        t1i.second += t2i.second;
        if (t1i.second)
          result.add_index_value(t1i);
        ++i1;
        ++i2;
      }
    } else if (i1<t1.size()) {
      for ( ; i1<t1.size(); ++i1) {
        auto t1i = t1.IndexValue(i1);
        result.add_index_value(t1i);
      }
      break;
    } else if (i2<t2.size()) {
      for ( ; i2<t2.size(); ++i2) {
        auto t2i = t2.IndexValue(i2);
        result.add_index_value(t2i);
      }
      break;
    } else
      break;
  }

  result.shrink_to_fit();
  return result;
}

/// Merge 2 sorted LinTerms
LinTerms Merge(const LinTerms& t1, const LinTerms& t2) {
  return MergeSorted(t1, t2);
}

/// Merge 2 sorted QuadTerms
QuadTerms Merge(const QuadTerms& t1, const QuadTerms& t2) {
  return MergeSorted(t1, t2);
}


QuadTerms MultiplyOut(const LinTerms& ae1, const LinTerms& ae2) {
  QuadTerms result;
  result.reserve(ae1.size()*ae2.size());

#ifndef OLD_MULTOUT
  // Dave style
  size_t i=0, j=0;
  while (i<ae1.size() && j<ae2.size()) {
    auto vi = ae1.var(i), vj = ae2.var(j);
    if (vi < vj) {
      for (auto j1=j; j1<ae2.size(); ++j1)   // follow from j
        result.add_term(
            ae1.coef(i)*ae2.coef(j1), vi, ae2.var(j1));
      ++i;
    } else if (vj < vi) {
      for (auto i1=i; i1<ae1.size(); ++i1)   // follow from i
        result.add_term(
            ae1.coef(i1)*ae2.coef(j), vj, ae1.var(i1));
      ++j;
    } else {                                 // vi==vj
      result.add_term(
          ae1.coef(i)*ae2.coef(j), vi, vi);
      size_t i1=i, j1=j;
      ++i1;
      ++j1;
      while (i1<ae1.size() && j1<ae2.size()) {
        auto vi1 = ae1.var(i1), vj1 = ae2.var(j1);
        if (vi1 < vj1) {
          result.add_term(
              ae1.coef(i1)*ae2.coef(j), vi, vi1);
          ++i1;
        } else if (vj1 < vi1) {
          result.add_term(
              ae1.coef(i)*ae2.coef(j1), vi, vj1);
          ++j1;
        } else {                             // vi1==vj1
          if (auto v = ae1.coef(i)*ae2.coef(j1) + ae1.coef(i1)*ae2.coef(j))
            result.add_term(v, vi, vi1);
          ++i1;
          ++j1;
        }
      }
      ++i;
      ++j;
    }
  }
#else
  for (auto i1 = ae1.size(); i1--; ) {
    for (auto i2 = ae2.size(); i2--; ) {
      result.add_term(ae1.coef(i1) * ae2.coef(i2),
                      ae1.var(i1), ae2.var(i2) );
    }
  }
  result.sort_terms();      // eliminate 0's and duplicates
#endif

  result.shrink_to_fit();
  assert(result.is_sorted());
  return result;
}


EExpr MultiplyOut(const EExpr& el, const EExpr& er) {
  assert((el.is_affine() && er.is_affine()) ||
         (el.is_constant() || er.is_constant()));
  EExpr result;
  if (er.constant_term()) {
    result.GetLinTerms().add(el.GetLinTerms());
    result.GetLinTerms() *= er.constant_term();
    result.GetQPTerms().add(el.GetQPTerms());
    result.GetQPTerms() *= er.constant_term();
    result.constant_term(
        er.constant_term() * el.constant_term());
  }
  if (el.constant_term()) {
    if (er.GetLinTerms().size()) {
      if (result.GetLinTerms().size()) {
        auto ae2 = er.GetLinTerms();    // @todo in-place if el!=er
        ae2 *= el.constant_term();
        result.GetLinTerms().sort_terms();  // @todo only in Release
        ae2.sort_terms();
        result.GetLinTerms()
            = Merge(result.GetLinTerms(), ae2);
      } else {
        result.GetLinTerms().add(er.GetLinTerms());
        result.GetLinTerms() *= el.constant_term();
      }
    }      // at most 1 factor has QP:
    result.GetQPTerms().add(er.GetQPTerms());
    result.GetQPTerms() *= el.constant_term();
  }
  if (el.is_affine() && er.is_affine()) {
    assert(result.GetQPTerms().empty());
    result.GetQPTerms()
        = MultiplyOut(el.GetLinTerms(), er.GetLinTerms());
  }
  return result;
}


// Instantiate BucketAccumulator
template
class BucketAccumulator<EExpr>;

// For some reason we need explicit:
template
    class BucketAccum1Type<LinTerms>;
template
    class BucketAccum1Type<QuadTerms>;

///////////////////////// END SORTING /////////////////////////////


const char*
BasicConstraintKeeper::GetShortTypeName() const {
  if (type_name_short_.empty()) {
    std::string acc_opt = GetAcceptanceOptionNames();
    assert(acc_opt.size());
    auto word_end = std::min(acc_opt.find(' '),
                             acc_opt.size());
    auto colon_pos = acc_opt.find(':');
    if (colon_pos>word_end)
      colon_pos = 0;
    type_name_short_ = acc_opt.substr(
        colon_pos, word_end-colon_pos);
    for (auto& c: type_name_short_)
      if (':'==c)
        c = '_';                // Markdown
    assert(type_name_short_.size());
  }
  return type_name_short_.c_str();
}


/// acceptance when constraint only
static const mp::OptionValueInfo values_con_acceptance[] = {
    { "0", "Not accepted natively, automatic redefinition will be attempted", 0},
    { "1", "Accepted but automatic redefinition will be used where possible", 1},
    { "2", "Accepted natively and preferred", 2}
};

/// acceptance when expression only
static const mp::OptionValueInfo values_expr_acceptance[] = {
    { "0", "Not accepted natively, automatic redefinition will be attempted", 0},
    { "3", "Accepted but automatic redefinition will be used where possible", 3},
    { "4", "Accepted natively and preferred", 4}
};

/// acceptance when both constraint and expression are possible
/// (ANY SOLVER DOING THIS? Ilog CP?)
static const mp::OptionValueInfo values_universal_acceptance[] = {
    { "0", "Not accepted natively, automatic redefinition will be attempted", 0},
    { "1", "Accepted as constraint but automatic redefinition will be used where possible", 1},
    { "2", "Accepted as constraint natively and preferred", 2},
    { "3", "Accepted as expression but automatic redefinition will be used where possible", 3},
    { "4", "Accepted as expression natively and preferred", 4}
};



void BasicConstraintKeeper::DoAddAcceptanceOptions(
    BasicFlatConverter& ,
    const BasicFlatModelAPI& ma,
    Env& env) {
  auto cal = GetModelAPIAcceptance(ma);
  auto eal = GetModelAPIAcceptanceEXPR(ma);
  auto eial = GetModelAPIAcceptance_EXPR_INTF(ma);
  const bool conacc = (ConstraintAcceptanceLevel::NotAccepted != cal);
  const bool expracc = (ExpressionAcceptanceLevel::NotAccepted != eal);
  const bool expr_intf_acc = (ExpressionAcceptanceLevel::NotAccepted != eial);
  acc_level_item_ = -1;           // user: unset
  acc_level_default_ = 0;         // default: not accepted
  if (conacc)
    acc_level_default_
        = std::underlying_type_t<ConstraintAcceptanceLevel>(cal);
  // we prefer expressions, if ModelAPI accepts expression interface
  if (expracc && expr_intf_acc)
    acc_level_default_      // Won't be taken however, if acc:_expr==0
        = std::underlying_type_t<ExpressionAcceptanceLevel>(eal) + 2;
  if (conacc && expracc) {
    env.AddStoredOption(GetAcceptanceOptionNames(),
                        fmt::format(
                            "Solver acceptance level for '{}' as either constraint or expression, "
                            "default {}:\n\n.. value-table::",
                            GetConstraintName(), acc_level_default_).c_str(),
                        acc_level_item_, values_universal_acceptance);
  } else
    if (conacc) {
      env.AddStoredOption(GetAcceptanceOptionNames(),
                          fmt::format(
                              "Solver acceptance level for '{}' as flat constraint, "
                              "default {}:\n\n.. value-table::",
                              GetConstraintName(), acc_level_default_).c_str(),
                          acc_level_item_, values_con_acceptance);
    } else
      if (expracc) {
        env.AddStoredOption(GetAcceptanceOptionNames(),
                            fmt::format(
                                "Solver acceptance level for '{}' as expression, "
                                "default {}:\n\n.. value-table::",
                                GetConstraintName(), acc_level_default_).c_str(),
                            acc_level_item_, values_expr_acceptance);
      } else {
        env.AddStoredOption(GetAcceptanceOptionNames(),
                            "HIDDEN",
                            acc_level_item_, 0, 4);
      }
}

void BasicConstraintKeeper::DoPopulateConstraintList(
    BasicFlatConverter& cvt,
    const BasicFlatModelAPI& ma,
    Env& env) {
  auto cancvt = IfConverterConverts(cvt);
  auto cal = GetModelAPIAcceptance(ma);
  auto eal = GetModelAPIAcceptanceEXPR(ma);
  // Description table
  env.SetConstraintListHeader(
      "List of flat constraints and corresponding expressions.\n"
      "For each constraint/expression, the following are given:\n"
      "\n"
      "  - name,\n"
      "  - convertibility into simpler forms,\n"
      "  - solver acceptance natively as flat constraint,\n"
      "  - solver acceptance natively as expression,\n"
      "  - driver option(s) to modify acceptance\n"
      "    (effective if both convertible and accepted).");
  std::string con_descr = (cancvt) ? "Convertible" : "NonConvertible";
  con_descr += "; ";
  const char * const acc_lev_nam[] = {
       "NotAccepted", "NativeAcceptedButNotRecommended", "NativeRecommended"
  };
  con_descr += acc_lev_nam[std::underlying_type_t<ConstraintAcceptanceLevel>(cal)];
  con_descr += "; ";
  con_descr += acc_lev_nam[std::underlying_type_t<ExpressionAcceptanceLevel>(eal)];
  con_descr += "; ";
  con_descr += GetAcceptanceOptionNames();
  env.AddConstraintDescr(GetConstraintName(), con_descr);
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
                    ItemNamer& vnam) {
  for (int i=0; i<(int)lt.size(); ++i) {
    auto coef = lt.coef(i);
    bool ifpos = coef>=0.0;
    if (i) {
      wrt << (ifpos ? " + " : " - ");
    } else {
      if (!ifpos)
        wrt << "-";
    }
    auto abscoef = std::fabs(coef);
    if (1.0 != abscoef)
      wrt << abscoef << '*';
    wrt << vnam.at(lt.var(i));
  }
}

void WriteModelItem(fmt::MemoryWriter& wrt, const QuadTerms& qt,
                    ItemNamer& vnam) {
  for (int i=0; i<(int)qt.size(); ++i) {
    auto coef = qt.coef(i);
    bool ifpos = coef>=0.0;
    if (i) {
      wrt << (qt.coef(i)>=0.0 ? " + " : " - ");
    } else
      if (!ifpos)
        wrt << "-";
    auto abscoef = std::fabs(coef);
    if (1.0 != abscoef)
      wrt << abscoef << '*';
    if (qt.var1(i)==qt.var2(i))
      wrt << vnam.at(qt.var1(i)) << "^2";
    else
      wrt << vnam.at(qt.var1(i))
          << '*' << vnam.at(qt.var2(i));
  }
}

void WriteModelItem(fmt::MemoryWriter& wrt, const QuadAndLinTerms& qlt,
                    ItemNamer& vnam) {
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
                    ItemNamer& vnam) {
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

void VisitArguments(const LinTerms& lt, std::function<void (int)> argv) {
  for (auto v: lt.vars())
    argv(v);
}

void VisitArguments(const QuadTerms& lt, std::function<void (int)> argv) {
  for (auto v: lt.vars1())
    argv(v);
  for (auto v: lt.vars2())
    argv(v);
}

void VisitArguments(const QuadAndLinTerms& qlt, std::function<void (int)> argv) {
  VisitArguments(qlt.GetLinTerms(), argv);
  VisitArguments(qlt.GetQPTerms(), argv);
}

std::unique_ptr<FlatModelInfo> CreateFlatModelInfo() {
  return std::unique_ptr<FlatModelInfo>{new FlatModelInfoImpl()};
}

void PrintModelInfo(const FlatModelInfo& fmi,
    const char* header, bool aux_vars) {
  auto vi = fmi.GetVarInfo();
  int nv = vi[0] + aux_vars*vi[3];
  int nvi = vi[1] + aux_vars*vi[4];
  int nvb = vi[2] + aux_vars*vi[5];
  fmt::print(fmt::format(
      "{} has {} variables ({} integer, {} binary);\n",
      header, nv, nvi, nvb));
  auto oi = fmi.GetObjInfo();
  if (true || 1!=oi[0] || oi[1] || oi[2]) {
    if (!oi[0] && !oi[1] && !oi[2])
      fmt::print("No objectives;\n");
    else {
      fmt::print("Objectives: ");
      if (oi[0])
        fmt::print(fmt::format("{} linear; ", oi[0]));
      if (oi[1])
        fmt::print(fmt::format("{} quadratic; ", oi[1]));
      if (oi[2])
        fmt::print(fmt::format("{} nonlinear; ", oi[2]));
      fmt::print("\n");
    }
  }
  const auto& coninfo = fmi.GetConstraintTypes();
  int n_lin = 0;
  int n_quad = 0;
  int n_nl = 0;
  int n_cones = 0;
  int n_condlin = 0, n_condquad = 0;
  int n_indlin = 0, n_indquad = 0;
  int n_sos1 = 0, n_sos2 = 0;
  std::map<std::string, int> expr_alg, expr_logic, cones;
  auto NewName = [](const char* old) { return old+1; };  // skip _
  for (const auto& val: coninfo) {
    if (val.second.n_) {
      if (begins_with(val.first, "_lin"))
        n_lin += val.second.n_;
      else if (begins_with(val.first, "_quad"))
        n_quad += val.second.n_;
      else if (ends_with(val.first, "cone")) {
        n_cones += val.second.n_;
        cones[NewName(val.second.name_)] = val.second.n_;
      }
      else if (begins_with(val.first, "_condlin"))
        n_condlin += val.second.n_;
      else if (begins_with(val.first, "_condquad"))
        n_condquad += val.second.n_;
      else if (begins_with(val.first, "_nl")) {
        if (val.first == "_nlcon")
          n_nl += val.second.n_;
      }  // else, NL assignment or logical - skip
      else if (begins_with(val.first, "_sos1"))
        n_sos1 += val.second.n_;
      else if (begins_with(val.first, "_sos2"))
        n_sos2 += val.second.n_;
      else if (begins_with(val.first, "_indlin"))
        n_indlin += val.second.n_;
      else if (begins_with(val.first, "_indquad"))
        n_indquad += val.second.n_;
      else if (begins_with(val.first, "_uenc"))
      { }
      // Else, it's expressions
      else if (val.second.is_logical_)
        expr_logic[NewName(val.second.name_)] = val.second.n_;
      else
        expr_alg[NewName(val.second.name_)] = val.second.n_;
    }
  }
  auto PrnType = [](const char* descr, int n) {
    if (n)
      fmt::print(fmt::format(" {} {};", n, descr));
  };
  if (n_lin + n_quad + n_nl + n_cones + n_sos1 + n_sos2) {
    fmt::print("Constraints: ");
    PrnType("linear", n_lin);
    PrnType("quadratic", n_quad);
    PrnType("nonlinear", n_nl);
    if (n_cones) {
      PrnType("conic", n_cones);
      for (const auto& cone: cones)
        PrnType(cone.first.c_str(), cone.second);
      fmt::print(");");
    }
    PrnType("SOS1", n_sos1);
    PrnType("SOS2", n_sos2);
    fmt::print("\n");
  }
  if (expr_alg.size()) {
    fmt::print("Algebraic expressions: ");
    for (const auto& expr: expr_alg)
      PrnType(expr.first.c_str(), expr.second);
    fmt::print("\n");
  }
  if (expr_logic.size()
      || n_condlin || n_condquad
      || n_indlin || n_indquad) {
    fmt::print("Logical expressions: ");
    PrnType("indicator(s)", n_indlin);
    PrnType("quadratic indicator(s)", n_indquad);
    PrnType("conditional (in)equalitie(s)", n_condlin);
    PrnType("conditional quadratic (in)equalitie(s)", n_condquad);
    for (const auto& expr: expr_logic)
      PrnType(expr.first.c_str(), expr.second);
    fmt::print("\n");
  }
  fmt::print("\n");
}

void ReportModelInfoSuffixes(const FlatModelInfo& fmi,
    std::string suf_prefix, SuffixGetterSetter sgs) {
  auto PutIntSuf = [suf_prefix, sgs](const char* name_extra, int val) {
    sgs.ssi_( {suf_prefix + name_extra, suf::PROBLEM}, {&val, 1} );
  };
  {
    auto vi = fmi.GetVarInfo();
    PutIntSuf("var_orig", vi[0]);
    PutIntSuf("var_orig_int", vi[1]);
    PutIntSuf("var_orig_bin", vi[2]);
    PutIntSuf("var_aux", vi[3]);
    PutIntSuf("var_aux_int", vi[4]);
    PutIntSuf("var_aux_bin", vi[5]);
  }
  {
    auto oi = fmi.GetObjInfo();
    PutIntSuf("obj_lin", oi[0]);
    PutIntSuf("obj_quad", oi[1]);
    PutIntSuf("obj_nonlin", oi[2]);
  }
  const auto& coninfo = fmi.GetConstraintTypes();
  for (const auto& val: coninfo) {
    PutIntSuf(val.second.name_, val.second.n_);
  }
}

} // namespace mp
