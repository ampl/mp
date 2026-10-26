#include <cmath>
#include <limits>
#include <stdexcept>

#include "cpsatmodelapi.h"

namespace mp {

namespace sat = operations_research::sat;

namespace {

/// Build a CP-SAT linear expression from terms coef(i)*var(i).
/// The coefficients should be integer.
template <class CoefFn, class VarFn, class MakeVar>
sat::LinearExpr BuildLinExpr(size_t n, CoefFn coef, VarFn var,
                             MakeVar make_var) {
  std::vector<sat::IntVar> vs;
  std::vector<int64_t> cs;
  vs.reserve(n);
  cs.reserve(n);
  for (size_t i = 0; i < n; ++i) {
    vs.push_back(make_var(var(i)));
    cs.push_back(ToInt64(coef(i), "coefficients"));
  }
  return sat::LinearExpr::WeightedSum(vs, cs);
}

}  // namespace


void CpsatModelAPI::InitProblemModificationPhase(
    const FlatModelInfo*) { }

void CpsatModelAPI::FinishProblemModificationPhase() { }


//////////////////////////// helpers ////////////////////////////

CpsatModelAPI::BoolVar CpsatModelAPI::Lit(int i) const {
  const auto* d = data();
  if (d->var_lb.at(i) < 0 || d->var_ub.at(i) > 1)
    throw std::runtime_error(fmt::format(
        "CP-SAT: variable {} is used as a logical value "
        "but its domain is [{}, {}], not within [0, 1]",
        i, d->var_lb.at(i), d->var_ub.at(i)));
  return d->vars[i].ToBoolVar();
}

std::vector<CpsatModelAPI::IntVar> CpsatModelAPI::Vars(
    const std::vector<int>& idx) const {
  std::vector<IntVar> r;
  r.reserve(idx.size());
  for (auto i: idx)
    r.push_back(Var(i));
  return r;
}

CpsatModelAPI::Domain CpsatModelAPI::MakeDomain(double lb, double ub) {
  const auto kMin = std::numeric_limits<int64_t>::min();
  const auto kMax = std::numeric_limits<int64_t>::max();
  const int64_t l = (lb <= -1e18) ? kMin : (int64_t)std::ceil(lb - 1e-9);
  const int64_t u = (ub >= 1e18) ? kMax : (int64_t)std::floor(ub + 1e-9);
  return Domain(l, u);
}

bool CpsatModelAPI::IsFixedTo(int var, int64_t value) const {
  return data()->var_lb.at(var) == value && data()->var_ub.at(var) == value;
}

void CpsatModelAPI::AddAnd(int r, const std::vector<BoolVar>& lits) {
  auto& b = data()->builder;
  const auto rl = Lit(r);
  if (lits.empty()) {                       // empty conjunction is true
    b.AddBoolAnd(std::vector<BoolVar>{rl});
    return;
  }
  std::vector<BoolVar> neg;
  neg.reserve(lits.size());
  for (const auto& l: lits)
    neg.push_back(l.Not());
  b.AddBoolAnd(lits).OnlyEnforceIf(rl);     // r => all
  b.AddBoolOr(neg).OnlyEnforceIf(rl.Not()); // !r => some false
}


//////////////////////////// variables, objective ////////////////////////////

void CpsatModelAPI::AddVariables(const VarArrayDef& v) {
  auto* d = data();
  d->vars.reserve(d->vars.size() + v.size());
  for (int i = 0; i < (int)v.size(); ++i) {
    // MP represents constants as fixed (continuous) variables;
    // those are fine if the value is integer
    const bool fixed_int = v.plb()[i] == v.pub()[i] &&
        IsIntegerValue(v.plb()[i]);
    if (v.ptype()[i] != var::Type::INTEGER && !fixed_int)
      throw std::runtime_error(fmt::format(
          "CP-SAT supports integer variables only; "
          "variable {} is continuous", d->vars.size()));
    const auto lb = LowerBoundToInt64(v.plb()[i]);
    const auto ub = UpperBoundToInt64(v.pub()[i]);
    d->var_lb.push_back(lb);
    d->var_ub.push_back(ub);
    d->vars.push_back(d->builder.NewIntVar(Domain(lb, ub)));
  }
}

void CpsatModelAPI::SetLinearObjective(int iobj, const LinearObjective& lo) {
  if (iobj > 0)
    throw std::runtime_error("CP-SAT: multiple objectives not supported");
  auto& b = data()->builder;
  const auto n = (size_t)lo.num_terms();
  const auto coefs = lo.coefs();
  const auto vars = lo.vars();
  bool all_int = true;
  for (size_t j = 0; j < n; ++j)
    if (!IsIntegerValue(coefs[j]))
      all_int = false;
  const bool maximize = (lo.obj_sense() == obj::Type::MAX);
  if (all_int) {
    auto e = BuildLinExpr(n,
        [&](size_t j) { return coefs[j]; },
        [&](size_t j) { return vars[j]; },
        [this](int v) { return Var(v); });
    if (maximize) b.Maximize(e); else b.Minimize(e);
  } else {                                   // floating-point objective
    sat::DoubleLinearExpr e;
    for (size_t j = 0; j < n; ++j)
      e.AddTerm(Var(vars[j]), coefs[j]);
    if (maximize) b.Maximize(e); else b.Minimize(e);
  }
}

void CpsatModelAPI::SetQuadraticObjective(int, const QuadraticObjective&) {
  throw std::runtime_error("CP-SAT: quadratic objectives not supported");
}


//////////////////////////// linear constraints ////////////////////////////

#define CPSAT_ADD_LINCON(ConType)                                          \
  void CpsatModelAPI::AddConstraint(const ConType& lc) {                   \
    const auto e = BuildLinExpr(lc.size(),                                 \
        [&](size_t i) { return lc.coef(i); },                              \
        [&](size_t i) { return lc.var(i); },                               \
        [this](int v) { return Var(v); });                                 \
    data()->builder.AddLinearConstraint(e, MakeDomain(lc.lb(), lc.ub()));  \
  }
CPSAT_ADD_LINCON(LinConRange)
CPSAT_ADD_LINCON(LinConLE)
CPSAT_ADD_LINCON(LinConEQ)
CPSAT_ADD_LINCON(LinConGE)
#undef CPSAT_ADD_LINCON


//////////////////////////// arithmetic ////////////////////////////

void CpsatModelAPI::AddConstraint(const MaxConstraint& mc) {
  data()->builder.AddMaxEquality(
      Var(mc.GetResultVar()), Vars(mc.GetArguments()));
}

void CpsatModelAPI::AddConstraint(const MinConstraint& mc) {
  data()->builder.AddMinEquality(
      Var(mc.GetResultVar()), Vars(mc.GetArguments()));
}

void CpsatModelAPI::AddConstraint(const AbsConstraint& ac) {
  data()->builder.AddAbsEquality(
      Var(ac.GetResultVar()), Var(ac.GetArguments().at(0)));
}

void CpsatModelAPI::AddConstraint(const MulConstraint& mc) {
  const auto& args = mc.GetArguments();
  data()->builder.AddMultiplicationEquality(
      Var(mc.GetResultVar()), Var(args[0]), Var(args[1]));
}

void CpsatModelAPI::AddConstraint(const PowConstExpConstraint& pc) {
  auto& b = data()->builder;
  const double e = pc.GetParameters().at(0);
  if (!IsIntegerValue(e) || e < 0.0)
    throw std::runtime_error(fmt::format(
        "CP-SAT: only powers with integer exponents >= 0 are "
        "supported, got exponent {}", e));
  const auto r = Var(pc.GetResultVar());
  const auto x = Var(pc.GetArguments().at(0));
  const auto k = (size_t)std::llround(e);
  if (0 == k)                                   // x^0 = 1
    b.AddEquality(r, 1);
  else if (1 == k)
    b.AddEquality(r, x);
  else
    b.AddMultiplicationEquality(r, std::vector<LinearExpr>(k, x));
}


//////////////////////////// logic ////////////////////////////

void CpsatModelAPI::AddConstraint(const AndConstraint& cc) {
  std::vector<BoolVar> lits;
  for (auto a: cc.GetArguments())
    lits.push_back(Lit(a));
  AddAnd(cc.GetResultVar(), lits);
}

void CpsatModelAPI::AddConstraint(const OrConstraint& cc) {
  auto& b = data()->builder;
  const auto rl = Lit(cc.GetResultVar());
  std::vector<BoolVar> lits, neg;
  for (auto a: cc.GetArguments()) {
    lits.push_back(Lit(a));
    neg.push_back(lits.back().Not());
  }
  if (lits.empty()) {                       // empty disjunction is false
    b.AddBoolAnd(std::vector<BoolVar>{rl.Not()});
    return;
  }
  b.AddBoolOr(lits).OnlyEnforceIf(rl);        // r => some
  b.AddBoolAnd(neg).OnlyEnforceIf(rl.Not());  // !r => all false
}

template <class Con>
void CpsatModelAPI::AddIndicator(const Con& ic) {
  const auto b = Lit(ic.get_binary_var());
  const auto& c = ic.get_constraint();
  const auto e = BuildLinExpr(c.size(),
      [&](size_t i) { return c.coef(i); },
      [&](size_t i) { return c.var(i); },
      [this](int v) { return Var(v); });
  data()->builder.AddLinearConstraint(e, MakeDomain(c.lb(), c.ub()))
      .OnlyEnforceIf(ic.get_binary_value() ? b : b.Not());
}

void CpsatModelAPI::AddConstraint(const IndicatorConstraintLinLE& ic) {
  AddIndicator(ic);
}
void CpsatModelAPI::AddConstraint(const IndicatorConstraintLinEQ& ic) {
  AddIndicator(ic);
}
void CpsatModelAPI::AddConstraint(const IndicatorConstraintLinGE& ic) {
  AddIndicator(ic);
}


//////////////////////////// constraint programming ////////////////////////////

void CpsatModelAPI::AddConstraint(const AllDiffConstraint& ac) {
  auto& b = data()->builder;
  const auto& args = ac.GetArguments();
  const int r = ac.GetResultVar();
  const auto xs = Vars(args);
  if (xs.size() < 2) {                      // trivially true
    b.AddEquality(Var(r), 1);
    return;
  }
  if (IsFixedTo(r, 1)) {                    // the common case: just alldiff
    b.AddAllDifferent(xs);
    return;
  }
  // Reified alldiff: r <==> all pairs different
  std::vector<BoolVar> diff;
  for (size_t i = 0; i < xs.size(); ++i) {
    for (size_t j = i + 1; j < xs.size(); ++j) {
      const auto d = b.NewBoolVar();
      b.AddNotEqual(xs[i], xs[j]).OnlyEnforceIf(d);
      b.AddEquality(xs[i], xs[j]).OnlyEnforceIf(d.Not());
      diff.push_back(d);
    }
  }
  AddAnd(r, diff);
}

void CpsatModelAPI::AddConstraint(const CountConstraint& cc) {
  // The arguments are logical (0/1) values
  data()->builder.AddEquality(
      Var(cc.GetResultVar()), LinearExpr::Sum(Vars(cc.GetArguments())));
}

void CpsatModelAPI::AddConstraint(const NumberofConstConstraint& nc) {
  auto& b = data()->builder;

  constexpr double feasTol = 1e-6;
  const double kd = nc.GetParameters().at(0);
  if (!IsInteger(kd, feasTol)) {                     // never equal
    b.AddEquality(Var(nc.GetResultVar()), 0);
    return;
  }
  const int64_t k = (int64_t)std::round(kd);
  std::vector<BoolVar> flags;
  for (auto a: nc.GetArguments()) {
    const auto x = Var(a);
    const auto f = b.NewBoolVar();                // f <==> (x == k)
    b.AddEquality(x, k).OnlyEnforceIf(f);
    b.AddNotEqual(x, k).OnlyEnforceIf(f.Not());
    flags.push_back(f);
  }
  b.AddEquality(Var(nc.GetResultVar()), LinearExpr::Sum(flags));
}

void CpsatModelAPI::AddConstraint(const NumberofVarConstraint& nc) {
  auto& b = data()->builder;
  const auto& args = nc.GetArguments();
  const auto xx = Var(args.at(0));
  std::vector<BoolVar> flags;
  for (size_t i = 1; i < args.size(); ++i) {
    const auto x = Var(args[i]);
    const auto f = b.NewBoolVar();                // f <==> (x == xx)
    b.AddEquality(x, xx).OnlyEnforceIf(f);
    b.AddNotEqual(x, xx).OnlyEnforceIf(f.Not());
    flags.push_back(f);
  }
  b.AddEquality(Var(nc.GetResultVar()), LinearExpr::Sum(flags));
}

void CpsatModelAPI::AddConstraint(const IfThenConstraint& ic) {
  auto& b = data()->builder;
  const auto& a = ic.GetArguments();
  const auto cond = Lit(a.at(0));
  const auto r = Var(ic.GetResultVar());
  b.AddEquality(r, Var(a.at(1))).OnlyEnforceIf(cond);
  b.AddEquality(r, Var(a.at(2))).OnlyEnforceIf(cond.Not());
}

void CpsatModelAPI::AddConstraint(const ImplicationConstraint& ic) {
  // Logic-valued if-then-else: same encoding on 0/1 values
  auto& b = data()->builder;
  const auto& a = ic.GetArguments();
  const auto cond = Lit(a.at(0));
  const auto r = Var(ic.GetResultVar());
  b.AddEquality(r, Var(a.at(1))).OnlyEnforceIf(cond);
  b.AddEquality(r, Var(a.at(2))).OnlyEnforceIf(cond.Not());
}

} // namespace mp
