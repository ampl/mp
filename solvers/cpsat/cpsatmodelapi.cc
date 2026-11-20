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
  if (std::isnan(lb) || std::isnan(ub))
    throw std::runtime_error("CP-SAT: a constraint bound is not a number");
  // The expression is integer-valued, so rounding a fractional bound
  // inward keeps the row's solutions; a near-integer bound is that integer.
  // Bounds beyond int64 are clamped: Domain(l, u) is empty when l > u.
  return Domain(ClampToInt64(RoundBoundInward(lb, true)),
                ClampToInt64(RoundBoundInward(ub, false)));
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
    const bool fixed = v.plb()[i] == v.pub()[i];
    const bool fixed_int = fixed && IsIntegerValue(v.plb()[i]);
    // Contradictory bounds (lb>ub) make the model infeasible, whatever
    // the type: do not complain about the type
    const bool contradictory = v.plb()[i] > v.pub()[i];
    d->var_const.push_back(NAN);
    if (fixed && !fixed_int && v.ptype()[i] != var::Type::INTEGER) {
      d->var_const.back() = v.plb()[i];  // see Var(); a placeholder at 0
      d->var_lb.push_back(0);
      d->var_ub.push_back(0);
      d->vars.push_back(d->builder.NewIntVar(Domain(0, 0)));
      continue;
    }
    if (v.ptype()[i] != var::Type::INTEGER && !fixed_int && !contradictory)
      throw std::runtime_error(fmt::format(
          "CP-SAT supports integer variables only; "
          "variable {} is continuous", d->vars.size()));
    auto lb = LowerBoundToInt64(v.plb()[i]);
    auto ub = UpperBoundToInt64(v.pub()[i]);
    if (lb > ub) {
      // An empty domain (also: an integer variable fixed at a fractional
      // value) is rejected by CP-SAT as an invalid model, not reported
      // as infeasible. Keep a valid domain and state the infeasibility.
      d->builder.AddLinearConstraint(sat::LinearExpr(int64_t{0}), Domain(1, 1));
      lb = ub;
    }
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
  // Variables fixed at non-integer values (constants) go to the offset
  std::vector<double> cs;
  std::vector<int> vs;
  double offset = 0.0;
  for (size_t j = 0; j < n; ++j) {
    const double c = data()->var_const.at(vars[j]);
    if (std::isnan(c)) {
      cs.push_back(coefs[j]);
      vs.push_back(vars[j]);
    } else
      offset += coefs[j] * c;
  }
  bool all_int = IsIntegerValue(offset);
  for (double c: cs)
    if (!IsIntegerValue(c))
      all_int = false;
  const bool maximize = (lo.obj_sense() == obj::Type::MAX);
  if (all_int) {
    auto e = BuildLinExpr(cs.size(),
        [&](size_t j) { return cs[j]; },
        [&](size_t j) { return vs[j]; },
        [this](int v) { return Var(v); });
    e += sat::LinearExpr(ToInt64(offset, "objective constant"));
    if (maximize) b.Maximize(e); else b.Minimize(e);
  } else {                                   // floating-point objective
    sat::DoubleLinearExpr e;
    for (size_t j = 0; j < cs.size(); ++j)
      e.AddTerm(Var(vs[j]), cs[j]);
    e += offset;
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

CpsatModelAPI::IntVar CpsatModelAPI::NonzeroDivisor(int y) {
  auto* d = data();
  const auto ly = d->var_lb.at(y), uy = d->var_ub.at(y);
  if (ly > 0 || uy < 0)
    return Var(y);
  if (0 == ly && 0 == uy)
    throw std::runtime_error("CP-SAT: division by 0");
  const auto y1 = d->builder.NewIntVar(
        Domain(ly, -1).UnionWith(Domain(1, uy)));
  d->builder.AddEquality(y1, Var(y));
  return y1;
}

void CpsatModelAPI::AddConstraint(const IntDivConstraint& dc) {
  const auto& a = dc.GetArguments();
  data()->builder.AddDivisionEquality(      // CP-SAT truncates, as AMPL
        Var(dc.GetResultVar()), Var(a[0]), NonzeroDivisor(a[1]));
}

void CpsatModelAPI::AddConstraint(const ModConstraint& mc) {
  // r = x - y*q, q = x div y; |q| <= |x| and |y*q| <= |x|
  auto* d = data();
  auto& b = d->builder;
  const int x = mc.GetArguments()[0], y = mc.GetArguments()[1];
  const int64_t m = std::max(std::abs(d->var_lb.at(x)),
                             std::abs(d->var_ub.at(x)));
  const auto q = b.NewIntVar(Domain(-m, m));
  b.AddDivisionEquality(q, Var(x), NonzeroDivisor(y));
  const auto p = b.NewIntVar(Domain(-m, m));
  b.AddMultiplicationEquality(p, Var(y), q);
  b.AddEquality(Var(mc.GetResultVar()), LinearExpr(Var(x)) - p);
}

void CpsatModelAPI::AddConstraint(const PowConstExpConstraint& pc) {
  auto& b = data()->builder;
  const double e = pc.GetParameters().at(0);
  if (!IsIntegerValue(e) || e < 0.0)
    throw std::runtime_error(fmt::format(
        "CP-SAT: only powers with integer exponents >= 0 are "
        "supported, got exponent {:.17g}", e));
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
