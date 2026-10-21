#ifndef DISJ2IND_H
#define DISJ2IND_H

/**
 * Absorbing static disjunctions into indicator constraints
 * (option cvt:pre:disj2ind).
 *
 * AMPL sends an implication  b==1 ==> x<=5  as the disjunction
 * b!=1 || x<=5, flattened into  Or(Not(b), r)  with  r = (x<=5).
 * Converted as is, this costs the binary r, the row Not(b)+r>=1,
 * and the indicator r==1 ==> x<=5.  When r is used only here,
 * the disjunction is equivalent to the single indicator
 * b==1 ==> x<=5, with the same projected LP relaxation.
 *
 * The absorption is the disjunction's own conversion
 * (OrConverter_MIP::Convert()). It is triggered by whichever item
 * of the chain is converted first: the Or itself, or a member
 * (conditional comparison, And, Not), via
 * ConstraintKeeper::ConvertConstraintStd() -> TryConvertParentFirst().
 * Bridging the Or then leaves the members unused.
 * A chain whose comparison awaits the unary encoding decision
 * (CondLinConEQ x==k) is delayed as a whole and
 * decided in GoWithoutEqEnc(), or by the final reassessment.
 */

#include <algorithm>
#include <tuple>
#include <type_traits>
#include <unordered_set>
#include <vector>

#include "mp/flat/constr_std.h"
#include "mp/flat/redef/MIP/cond_literal.h"
#include "mp/flat/redef/MIP/cond_ineq.h"
#include "mp/flat/redef/MIP/cond_eq.h"

namespace mp {

/// Bits of option cvt:pre:disj2ind
enum Disj2IndFlags {
  D2I_LINEAR = 1,       ///< linear conditional comparisons
  D2I_QUADRATIC = 2,    ///< quadratic conditional comparisons
  D2I_TWO_CHAINS = 4,   ///< both disjuncts absorbable: absorb one
  D2I_UENC = 8,         ///< x==k subject to the unary encoding decision
  D2I_AND_POS = 16,     ///< And of comparisons in positive context
  D2I_DISTRIBUTE = 32   ///< static l1 || .. || (y1 && ..): clauses li || yj
};

/// Linear conditional comparisons which can be absorbed
using Disj2IndLinLeafTypes = std::tuple<
  CondLinConLT, CondLinConLE, CondLinConEQ, CondLinConGE, CondLinConGT>;

/// Quadratic conditional comparisons which can be absorbed
using Disj2IndQuadLeafTypes = std::tuple<
  CondQuadConLT, CondQuadConLE, CondQuadConEQ,
  CondQuadConGE, CondQuadConGT>;

/// All absorbable conditional comparisons
using Disj2IndLeafTypes = decltype(std::tuple_cat(
  std::declval<Disj2IndLinLeafTypes>(),
  std::declval<Disj2IndQuadLeafTypes>()));

/// Is \a Con one of the types in tuple \a Tuple?
template <class Con, class Tuple>
struct TupleHasType;

/// Specialization for std::tuple
template <class Con, class... Ts>
struct TupleHasType<Con, std::tuple<Ts...> >
  : std::disjunction<std::is_same<Con, Ts>...> { };

/// Is \a Con an absorbable conditional comparison type?
template <class Con>
constexpr bool IsDisj2IndLeaf()
{ return TupleHasType<Con, Disj2IndLeafTypes>::value; }

/// Can \a Con be a chain member below the disjunction?
template <class Con>
constexpr bool IsDisj2IndMember() {
  return IsDisj2IndLeaf<Con>()
      || std::is_same_v<Con, AndConstraint>
      || std::is_same_v<Con, NotConstraint>;
}

/// Is \a Con an item of a chain, including the disjunction?
template <class Con>
constexpr bool IsDisj2IndItem()
{ return IsDisj2IndMember<Con>() || std::is_same_v<Con, OrConstraint>; }


/// A conditional comparison to be absorbed into an indicator
struct Disj2IndLeaf {
  int res_ {-1};        ///< its result variable, -1 if none
  bool neg_ {false};    ///< under Not: imply its negation
  bool literal_ {false};  ///< And member kept as a literal: imply res_==1
};

/// Analysis of a disjunction l || C
/// for absorption into indicators (!l) ==> C
struct Disj2IndChain {
  /// Absorption state
  enum class State {
    None,        ///< not absorbable
    Ready,       ///< absorbable now
    Blocked      ///< absorbable, but a comparison is delayed
  };

  State state_ {State::None};          ///< absorption state
  CondLiteral lit_ {-1, 0};            ///< literal triggering C: !l
  std::vector<Disj2IndLeaf> leaves_;   ///< the comparisons of C

  /// Has absorbable comparisons?
  bool HasLeaves() const { return !leaves_.empty(); }

  /// Contains a comparison with result variable \a res?
  bool Contains(int res) const {
    return std::any_of(leaves_.begin(), leaves_.end(),
                       [res](const Disj2IndLeaf& lf) { return lf.res_==res; });
  }
};


/// Add the indicator (lit) ==> \a cc, or (lit) ==> !cc if \a f_neg,
/// for any conditional comparison type
template <class MC, class AlgCon>
void AddLeafAsIndicator(MC& mc, const ConditionalConstraint<AlgCon>& cc,
                        CondLiteral lit, bool f_neg) {
  if constexpr (0 == AlgCon::kind()) {
    if (f_neg)
      AddCondNeAsIndicators(mc, cc, lit);
    else
      AddCondEqAsIndicator(mc, cc, lit);
  } else
    AddCondIneqAsIndicator(mc, cc, lit, f_neg);
}


/// Mix-in for the MIP converter \a Impl:
/// absorbing static disjunctions into indicators.
///
/// A disjunction Or(l, x) is absorbable when
/// - its result is fixed to 1 and it has two arguments;
/// - x is a chain used only here: a conditional comparison C,
///   or Not(C) for an equality C (a != consequent),
///   or And(...) of these; each node used once and not fixed;
/// - each comparison is still to be converted by this backend,
///   is in the appropriate context, and its type is enabled.
/// The literal is l, or b for l = Not(b).
/// If both arguments qualify, D2I_TWO_CHAINS selects one.
///
/// Parents are recorded per variable (allocated on demand)
/// and validated when used.
template <class Impl>
class Disj2IndMixin {
public:
  /// Register option cvt:pre:disj2ind
  void InitDisj2IndOptions() {
    MC().GetEnv().AddOption("cvt:pre:disj2ind disj2ind",
        "Absorb static disjunctions l || C, where C is a conditional "
        "comparison used only there (or an AND of such), into "
        "indicator constraints !l ==> C. This is how implications "
        "b==1 ==> C are received from AMPL. "
        "Bitwise OR of the following values:\n"
        "\n"
        "|  1 - Linear comparisons\n"
        "|  2 - Quadratic comparisons\n"
        "|  4 - Both disjuncts absorbable: absorb one into an "
        "indicator on the other's result (if-else shape, "
        "complementarity)\n"
        "|  8 - Comparisons x==k awaiting the unary encoding decision "
        "for x (absorbed only if no encoding is chosen)\n"
        "| 16 - Any AND of such comparisons (used only there) in positive "
        "context, e.g., a disjunct (C1 && C2) of a multi-term "
        "disjunction: indicators on the AND's result variable\n"
        "| 32 - Static disjunctions l1 || .. || lm || (y1 && .. && yk) "
        "with the AND used only there: clauses l1 || .. || lm || yj, "
        "without the AND's result variable.\n"
        "\n"
        "Requires cvt:pre:all=1. Default 11.",
        disj2ind_flags_, 0, 63);
  }

  /// Hook: record \a con as the parent of its arguments
  /// (Or, And, Not only)
  template <class Con>
  void NoteConstraintAdded(const Con& con, int ) {
    if constexpr (std::is_same_v<Con, OrConstraint>
                  || std::is_same_v<Con, AndConstraint>
                  || std::is_same_v<Con, NotConstraint>)
      if (IfDisj2Ind())
        NoteParentOf(con.GetArguments(), con.GetResultVar());
  }

  /// Hook: delay a chain item while its chain is blocked
  template <class Con>
  bool IfDelayCvtForParent(const Con& con, int ) {
    if constexpr (IsDisj2IndItem<Con>())
      return IfDisj2Ind() && DelayIfChainBlocked(con.GetResultVar());
    return false;
  }

  /// Hook: convert an accepted disjunction or conjunction
  /// if absorbable
  template <class Con>
  bool IfNeedsCvtForParent(const Con& con, int ) {
    if constexpr (std::is_same_v<Con, OrConstraint>)
      return IfDisj2Ind()
          && (Disj2IndChain::State::Ready
              == AnalyzeChain(con.GetResultVar()).state_
              || FindDistributableAnd(con) >= 0);
    if constexpr (std::is_same_v<Con, AndConstraint>)
      return IfDisj2Ind()
          && Disj2IndChain::State::Ready
             == AnalyzeConjunction(con.GetResultVar()).state_;
    return false;
  }

  /// Hook: before converting a chain member,
  /// convert (absorb) its disjunction if ready
  template <class Con>
  bool TryConvertParentFirst(const Con& con, int ) {
    if constexpr (IsDisj2IndMember<Con>()) {
      if (!IfDisj2Ind())
        return false;
      if (ConvertTopOrIfReady(con.GetResultVar()))
        return true;
      if constexpr (std::is_same_v<Con, AndConstraint>)
        return ConvertDistributingOr(con.GetResultVar());
      else
        return ConvertParentAndIfReady(con.GetResultVar());
    }
    return false;
  }

  /// Absorb the comparisons of conjunction \a conj into indicators
  /// on its result, if possible (D2I_AND_POS).
  /// Called by AndConverter_MIP::Convert().
  /// @return whether absorbed
  bool AbsorbConjunction(const AndConstraint& conj) {
    if (!IfDisj2Ind())
      return false;
    return EmitIfReady(AnalyzeConjunction(conj.GetResultVar()));
  }

  /// Absorb disjunction \a disj into indicators, if possible.
  /// Called by OrConverter_MIP::Convert().
  /// @return whether absorbed
  bool AbsorbDisjunction(const OrConstraint& disj) {
    if (!IfDisj2Ind())
      return false;
    return EmitIfReady(AnalyzeChain(disj.GetResultVar()))
        || DistributeOverAnd(disj);
  }


protected:
  /// The final converter
  Impl& MC() { return static_cast<Impl&>(*this); }

  /// Enabled option bits (0 also if cvt:pre:all=0)
  int Disj2IndFlags() const {
    return static_cast<const Impl&>(*this).CanPreprocess(disj2ind_flags_);
  }

  /// Is the absorption on?
  bool IfDisj2Ind() const { return 0 != Disj2IndFlags(); }


  ////////////////////////// Parents //////////////////////////

  /// Record \a parent as the parent of arguments with
  /// a defining constraint
  template <class Args>
  void NoteParentOf(const Args& args, int parent) {
    for (auto v: args)
      if (MC().HasActiveOrInactiveInitExpression(v))
        SetParent(v, parent);
  }

  /// Set the parent result variable of variable \a v
  void SetParent(int v, int parent) {
    if (v >= (int)parent_res_.size())
      parent_res_.resize(std::max(v+1, MC().num_vars()), -1);
    parent_res_[v] = parent;
  }

  /// Validated parent of variable \a v: the result variable
  /// of an active Or/And/Not with \a v as argument, or -1
  int GetParent(int v) {
    if (v < 0 || v >= (int)parent_res_.size())
      return -1;
    const int p = parent_res_[v];
    return p >= 0 && IsActiveParentOf(p, v) ? p : -1;
  }

  /// Does variable \a p have an active Or/And/Not definition
  /// with argument \a v?
  bool IsActiveParentOf(int p, int v) {
    return ArgsContain<OrConstraint>(p, v)
        || ArgsContain<AndConstraint>(p, v)
        || ArgsContain<NotConstraint>(p, v);
  }

  /// Does variable \a p have an active definition of type
  /// \a Con with argument \a v?
  template <class Con>
  bool ArgsContain(int p, int v) {
    const auto* pc = MC().template GetActiveInitExpressionOfType<Con>(p);
    if (!pc)
      return false;
    const auto& args = pc->GetArguments();
    return std::find(args.begin(), args.end(), v) != args.end();
  }

  /// The disjunction on top of \a v's chain, or -1:
  /// \a v itself if it is an Or, otherwise up to 3 parents up
  /// (a comparison under Not under And under Or).
  int FindTopOr(int v) {
    for (int depth=0; depth<=3 && v>=0; ++depth, v=GetParent(v))
      if (MC().template GetActiveInitExpressionOfType<OrConstraint>(v))
        return v;
    return -1;
  }


  ////////////////////////// Chain analysis //////////////////////////

  /// If the chain above \a v is blocked, note \a v as delayed
  /// for its chain (not for its own sake, see IsLeafDelayed()).
  /// @return whether blocked
  bool DelayIfChainBlocked(int v) {
    const bool f_blocked = Disj2IndChain::State::Blocked
        == AnalyzeChainAbove(v).state_;
    if (f_blocked)
      chain_delayed_.insert(v);
    return f_blocked;
  }

  /// Analyze the chain on top of variable \a v
  Disj2IndChain AnalyzeChainAbove(int v) {
    const int top = FindTopOr(v);
    return top >= 0 ? AnalyzeChain(top) : Disj2IndChain{};
  }

  /// Analyze disjunction \a or_res for absorption
  Disj2IndChain AnalyzeChain(int or_res) {
    const auto* pOr = GetStaticOr2(or_res);
    if (!pOr)
      return {};
    const auto& args = pOr->GetArguments();
    auto chain = SelectChain(ChainAt(args[0], args[1]),
                             ChainAt(args[1], args[0]));
    if (chain.HasLeaves())
      chain.state_ = IsBlocked(chain)
          ? Disj2IndChain::State::Blocked : Disj2IndChain::State::Ready;
    return chain;
  }

  /// The active Or defining \a v if it is fixed to 1
  /// and has two arguments, otherwise nullptr
  const OrConstraint* GetStaticOr2(int v) {
    const auto* pOr =
        MC().template GetActiveInitExpressionOfType<OrConstraint>(v);
    const bool f_static =
        pOr && 2 == pOr->GetArguments().size()
        && MC().is_fixed(v) && 1.0 == MC().fixed_value(v);
    return f_static ? pOr : nullptr;
  }

  /// Chain for consequent \a x and literal from disjunct \a l.
  /// Without D2I_TWO_CHAINS, \a l must not be a comparison's
  /// result: neither as a chain nor, after its own conversion,
  /// as a literal.
  Disj2IndChain ChainAt(int x, int l) {
    Disj2IndChain chain;
    if (!(Disj2IndFlags() & D2I_TWO_CHAINS) && IsComparisonResult(l))
      return chain;
    chain.leaves_ = ConsequentLeaves(x);
    chain.lit_ = MakeLiteral(l);
    return chain;
  }

  /// Choose between the chains of both disjuncts
  Disj2IndChain SelectChain(Disj2IndChain&& c0, Disj2IndChain&& c1) {
    if (!c0.HasLeaves())
      return std::move(c1);
    if (!c1.HasLeaves())
      return std::move(c0);
    if (!(Disj2IndFlags() & D2I_TWO_CHAINS))
      return {};
    return c1.Contains(trigger_) ? std::move(c1) : std::move(c0);
  }

  /// Literal triggering the consequent: !l, i.e.,
  /// b==1 for l = Not(b), otherwise l==0
  CondLiteral MakeLiteral(int l) {
    const auto* pNot =
        MC().template GetActiveInitExpressionOfType<NotConstraint>(l);
    if (pNot)
      return { pNot->GetArguments()[0], 1 };
    return { l, 0 };
  }

  /// Is the chain blocked by a delayed comparison
  /// (other than the trigger)?
  bool IsBlocked(const Disj2IndChain& chain) {
    for (const auto& leaf: chain.leaves_)
      if (!leaf.literal_ && leaf.res_ != trigger_
          && IsLeafDelayed(leaf.res_))
        return true;
    return false;
  }

  /// Absorbable comparisons of consequent \a x, or empty
  std::vector<Disj2IndLeaf> ConsequentLeaves(int x) {
    if (!IsSingleUseFree(x))
      return {};
    const auto leaf = MemberLeaf(x);
    if (leaf.res_ >= 0)
      return { leaf };
    return AndLeaves(x);
  }

  /// \a x is an absorbable comparison, or Not of an absorbable
  /// equality. Otherwise returns a leaf with res_ == -1.
  Disj2IndLeaf MemberLeaf(int x) {
    if (IsAbsorbableLeaf(x, false))
      return { x, false };
    const int y = NotArgIfEquality(x);
    if (y >= 0 && IsSingleUseFree(y) && IsAbsorbableLeaf(y, true))
      return { y, true };
    return {};
  }

  /// If \a x = Not(y) with y a conditional equality, return y,
  /// otherwise -1
  int NotArgIfEquality(int x) {
    const auto* pNot =
        MC().template GetActiveInitExpressionOfType<NotConstraint>(x);
    if (!pNot)
      return -1;
    const int y = pNot->GetArguments()[0];
    const bool f_eq =
        MC().template GetActiveInitExpressionOfType<CondLinConEQ>(y)
        || MC().template GetActiveInitExpressionOfType<CondQuadConEQ>(y);
    return f_eq ? y : -1;
  }

  /// If \a x = And(...) of absorbable members, their comparisons
  /// (all or nothing)
  /// If \a x = And(...), its members: absorbable comparisons,
  /// and other arguments as literals (implied by clauses).
  /// Empty if no member is an absorbable comparison.
  std::vector<Disj2IndLeaf> AndLeaves(int x) {
    const auto* pAnd =
        MC().template GetActiveInitExpressionOfType<AndConstraint>(x);
    if (!pAnd)
      return {};
    std::vector<Disj2IndLeaf> leaves;
    bool f_any_comparison = false;
    for (auto y: pAnd->GetArguments()) {
      auto leaf = IsSingleUseFree(y) ? MemberLeaf(y) : Disj2IndLeaf{};
      f_any_comparison = f_any_comparison || leaf.res_ >= 0;
      if (leaf.res_ < 0)
        leaf = { y, false, true };
      leaves.push_back(leaf);
    }
    return f_any_comparison ? leaves : std::vector<Disj2IndLeaf>{};
  }

  /// Is \a x used once and not fixed?
  /// (A fixed result would carry meaning beyond this chain.)
  bool IsSingleUseFree(int x) {
    return 1 == MC().VarUsage(x) && !MC().is_fixed(x);
  }


  /// Analyze And \a and_res for absorption into indicators
  /// (and_res == 1) ==> each comparison: it must be in positive
  /// context and not fixed, with all arguments absorbable.
  Disj2IndChain AnalyzeConjunction(int and_res) {
    Disj2IndChain chain;
    if (!(Disj2IndFlags() & D2I_AND_POS) || !IsPositiveFreeAnd(and_res))
      return chain;
    chain.leaves_ = AndLeaves(and_res);
    chain.lit_ = { and_res, 1 };
    if (chain.HasLeaves())
      chain.state_ = IsBlocked(chain)
          ? Disj2IndChain::State::Blocked : Disj2IndChain::State::Ready;
    return chain;
  }

  /// Is \a v the result of an active And in positive context,
  /// not fixed?
  bool IsPositiveFreeAnd(int v) {
    const auto* pAnd =
        MC().template GetActiveInitExpressionOfType<AndConstraint>(v);
    return pAnd && pAnd->GetContext().IsPositive() && !MC().is_fixed(v);
  }

  /// The And above comparison result \a v (directly, or via Not
  /// for an equality), or -1
  int FindParentAnd(int v) {
    for (int depth=0; depth<2 && v>=0; ++depth) {
      v = GetParent(v);
      if (v >= 0
          && MC().template GetActiveInitExpressionOfType<AndConstraint>(v))
        return v;
    }
    return -1;
  }

  /// If the And above \a v can absorb its comparisons,
  /// convert it, with \a v as the trigger.
  /// @return whether converted
  bool ConvertParentAndIfReady(int v) {
    const int and_res = FindParentAnd(v);
    if (and_res < 0)
      return false;
    TriggerScope trigger_scope(trigger_, v);
    if (Disj2IndChain::State::Ready != AnalyzeConjunction(and_res).state_)
      return false;
    MC().GetConstraintKeeper((AndConstraint*)nullptr).ConvertConstraintStd(
          MC().GetInitExpression(and_res).GetIndex());
    return true;
  }

  ////////////////////////// Distribution over And //////////////////////////

  /// Static disjunction \a disj with a distributable And argument
  /// (D2I_DISTRIBUTE): add a clause per And argument instead.
  /// @return whether distributed
  bool DistributeOverAnd(const OrConstraint& disj) {
    const int x = FindDistributableAnd(disj);
    if (x < 0)
      return false;
    std::vector<int> others;
    for (auto v: disj.GetArguments())
      if (v != x)
        others.push_back(v);
    for (auto y: MC().template GetActiveInitExpressionOfType<AndConstraint>(x)
         ->GetArguments())
      AddClause(others, y);
    return true;
  }

  /// The first And argument of static \a disj used only there,
  /// positive and not fixed, or -1
  int FindDistributableAnd(const OrConstraint& disj) {
    if (!(Disj2IndFlags() & D2I_DISTRIBUTE)
        || !IsStatic(disj.GetResultVar()))
      return -1;
    for (auto v: disj.GetArguments())
      if (IsSingleUseFree(v) && IsPositiveFreeAnd(v))
        return v;
    return -1;
  }

  /// Is \a v fixed to 1?
  bool IsStatic(int v) {
    return MC().is_fixed(v) && 1.0 == MC().fixed_value(v);
  }

  /// Add the static clause  others || y
  void AddClause(std::vector<int> others, int y) {
    others.push_back(y);
    MC().FixAsTrue( MC().AssignResultVar2Args(
                      OrConstraint{ std::move(others) } ) );
  }

  /// If the static Or above And \a and_res distributes over it,
  /// convert the Or first.
  /// @return whether converted
  bool ConvertDistributingOr(int and_res) {
    const int or_res = GetParent(and_res);
    const auto* pOr = or_res >= 0
        ? MC().template GetActiveInitExpressionOfType<OrConstraint>(or_res)
        : nullptr;
    if (!pOr || FindDistributableAnd(*pOr) != and_res)
      return false;
    MC().GetConstraintKeeper((OrConstraint*)nullptr).ConvertConstraintStd(
          MC().GetInitExpression(or_res).GetIndex());
    return true;
  }


  /// Emit the indicators of \a chain if it is ready.
  /// @return whether emitted
  bool EmitIfReady(const Disj2IndChain& chain) {
    if (Disj2IndChain::State::Ready != chain.state_)
      return false;
    for (const auto& leaf: chain.leaves_)
      EmitLeaf(leaf, chain.lit_);
    return true;
  }


  ////////////////////////// Comparisons //////////////////////////

  /// Call \a f(cc, i) for the active conditional comparison
  /// defining \a y, if any.
  /// @return whether there is one
  template <class F>
  bool VisitLeaf(int y, F&& f) {
    return VisitLeafOfTypes(y, f, (Disj2IndLeafTypes*)nullptr);
  }

  /// VisitLeaf() over the types of a tuple
  template <class F, class... Cs>
  bool VisitLeafOfTypes(int y, F& f, std::tuple<Cs...>* ) {
    return ( VisitLeafOfType<Cs>(y, f) || ... );
  }

  /// VisitLeaf() for type \a Con
  template <class Con, class F>
  bool VisitLeafOfType(int y, F& f) {
    const auto* pc = MC().template GetActiveInitExpressionOfType<Con>(y);
    if (pc)
      f(*pc, MC().GetInitExpression(y).GetIndex());
    return nullptr != pc;
  }

  /// Is \a l defined by a conditional comparison,
  /// whether still active or not?
  bool IsComparisonResult(int l) {
    return IsComparisonResultOfTypes(l, (Disj2IndLeafTypes*)nullptr);
  }

  /// IsComparisonResult() over the types of a tuple
  template <class... Cs>
  bool IsComparisonResultOfTypes(int l, std::tuple<Cs...>* ) {
    return ( (nullptr !=
              MC().template GetInitExpressionOfType<Cs>(l)) || ... );
  }

  /// Is \a y an absorbable comparison, in positive context
  /// (negative if \a f_neg)?
  bool IsAbsorbableLeaf(int y, bool f_neg) {
    bool f_ok = false;
    VisitLeaf(y, [this, f_neg, &f_ok](const auto& cc, int i) {
      f_ok = this->IsAbsorbable(cc, i, f_neg);
    });
    return f_ok;
  }

  /// Is comparison \a cc (index \a i) absorbable?
  template <class Con>
  bool IsAbsorbable(const Con& cc, int i, bool f_neg) {
    const auto ctx = cc.GetContext();
    return IfTypeEnabled<Con>()
        && !cc.GetConstraint().empty()
        && (f_neg ? ctx.IsNegative() : ctx.IsPositive())
        && MC().GetConstraintKeeper((Con*)nullptr).IfConvertsItems()
        && ((Disj2IndFlags() & D2I_UENC) || !IsLeafDelayed(cc, i));
  }

  /// Is the comparison type enabled by the option?
  template <class Con>
  bool IfTypeEnabled() {
    return Disj2IndFlags() &
        (TupleHasType<Con, Disj2IndQuadLeafTypes>::value
         ? D2I_QUADRATIC : D2I_LINEAR);
  }

  /// Is comparison \a y's conversion delayed?
  bool IsLeafDelayed(int y) {
    bool f_delayed = false;
    VisitLeaf(y, [this, &f_delayed](const auto& cc, int i) {
      f_delayed = this->IsLeafDelayed(cc, i);
    });
    return f_delayed;
  }

  /// Is comparison \a cc's conversion delayed by its own
  /// converter (not by its chain), now or before?
  /// (Before: e.g., x==k delayed in the main loop
  /// awaits the unary encoding decision in the map phase.)
  template <class Con>
  bool IsLeafDelayed(const Con& cc, int i) {
    return MC().IfDelayCvt_impl(cc, i)
        || (MC().GetConstraintKeeper((Con*)nullptr).IsDelayed(i)
            && !chain_delayed_.count(cc.GetResultVar()));
  }


  ////////////////////////// Conversion //////////////////////////

  /// If the chain above \a v is ready, convert its disjunction,
  /// with \a v as the trigger.
  /// @return whether converted
  bool ConvertTopOrIfReady(int v) {
    const int top = FindTopOr(v);
    if (top < 0)
      return false;
    TriggerScope trigger_scope(trigger_, v);
    if (Disj2IndChain::State::Ready != AnalyzeChain(top).state_)
      return false;
    MC().GetConstraintKeeper((OrConstraint*)nullptr).ConvertConstraintStd(
          MC().GetInitExpression(top).GetIndex());
    return true;
  }

  /// Add the indicator for comparison \a leaf triggered by \a lit
  void EmitLeaf(const Disj2IndLeaf& leaf, CondLiteral lit) {
    if (leaf.literal_)
      return AddImpliedLiteral(lit, leaf.res_);
    VisitLeaf(leaf.res_, [this, &leaf, lit](const auto& cc, int ) {
      AddLeafAsIndicator(this->MC(), cc, lit, leaf.neg_);
    });
  }

  /// Add (lit) ==> y==1 as the clause row
  /// [lit does not trigger] + y >= 1
  void AddImpliedLiteral(CondLiteral lit, int y) {
    auto row = NonTriggerExpr(lit);
    row.add_term(1.0, y);
    MC().AddConstraint( LinConGE( std::move(row.GetBody()),
                                  1.0 - row.constant_term() ) );
  }

  /// Sets the trigger variable for the scope
  class TriggerScope {
  public:
    /// Set \a trigger to \a v until destruction
    TriggerScope(int& trigger, int v) : trigger_(trigger), saved_(trigger)
    { trigger_ = v; }
    /// Restore
    ~TriggerScope() { trigger_ = saved_; }
  private:
    int& trigger_;
    const int saved_;
  };


private:
  int disj2ind_flags_ = D2I_LINEAR | D2I_QUADRATIC | D2I_UENC;
  std::vector<int> parent_res_;   // parent result var per variable
  int trigger_ = -1;              // chain item whose turn it is
  std::unordered_set<int> chain_delayed_;  // delayed for their chains
};

} // namespace mp

#endif // DISJ2IND_H
