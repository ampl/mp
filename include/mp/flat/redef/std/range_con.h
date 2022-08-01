#ifndef RANGE_CON_H
#define RANGE_CON_H

/*
 * Convert range constraints
 */

#include "mp/flat/redef/redef_base.h"
#include "mp/valcvt-link.h"
#include "mp/flat/constr_std.h"

namespace mp {

namespace pre {

/// Presolve link between RangeCon and either
/// Con(LE/GE) or ConEQ+Slack.
/// Parameters of the base class: self, N value nodes (3),
/// N indexes (3), N source nodes / indexes in the list (1).
template <class ModelConverter, class RangeCon>
class RangeCon2Slack :
    public BasicStaticIndivEntryLink<
      RangeCon2Slack<ModelConverter, RangeCon>, 3, 3, 1> {
public:
  /// Base class
  using Base = BasicStaticIndivEntryLink<
    RangeCon2Slack<ModelConverter, RangeCon>, 3, 3, 1>;

  /// Constructor.
  /// @param cvt: the model converter
  /// @param ndl: list of 3 value nodes
  /// (orig ranges, new constraints, new slacks)
  RangeCon2Slack(ModelConverter& cvt,
                 const typename Base::NodeList& ndl) :
    Base(cvt.GetValuePresolver(), ndl), cvt_(cvt) { }

  /// Type name
  const char* GetTypeName() const override { return "Range2Slk<...>"; }


  /// Names for indexes in a Base::LinkEntry
  enum LinkEntryIndexes {
    /// Source range constraint index
    CON_SRC = 0,
    /// Target equality constraint index
    CON_TARGET = 1,
    /// Target slack index
    VAR_SLK = 2
  };

  /// Define pre- / postsolve methods for individual link entries.
  /// Typedef LinkEntry is created in Base as std::array<int, 3>
  /// (3 is the base class template parameter)
  /// and means the following indexes: {range_con, lin_con, slack_index}

  /// PresolveGeneric double
  void PresolveGenericDblEntry(const typename Base::LinkEntry& be) {
    auto src_val = GetDbl(be, CON_SRC);
    SetDbl(be, CON_TARGET, src_val);
    SetDbl(be, VAR_SLK, src_val);
  }

  /// PostsolveGeneric double
  void PostsolveGenericDblEntry(const typename Base::LinkEntry& be) {
    SetDbl(be, CON_SRC, std::max(
             GetDbl(be, CON_TARGET), GetDbl(be, VAR_SLK)));
  }

  /// PresolveGeneric int
  void PresolveGenericIntEntry(const typename Base::LinkEntry& be) {
    auto src_val = GetInt(be, CON_SRC);
    SetInt(be, CON_TARGET, src_val);
    SetInt(be, VAR_SLK, src_val);
  }

  /// PostsolveGeneric int
  void PostsolveGenericIntEntry(const typename Base::LinkEntry& be) {
    SetInt(be, CON_SRC, std::max(
             GetInt(be, CON_TARGET), GetInt(be, VAR_SLK)));
  }

  /// Presolve solution (primal + dual).
  /// Duals: just copy.
  /// Primals: compute slack.
  /// Slacks are probably only needed for Gurobi LP warmstart.
  /// Who does not need them, can use unpresolved primal solution.
  void PresolveSolutionEntry(const typename Base::LinkEntry& be) {
    SetDbl(be, CON_TARGET, GetDbl(be, CON_SRC));
    const auto& orig_cons =
        GetMC().template GetConstraint<RangeCon>(be[CON_SRC]);
    SetDbl(be, VAR_SLK, orig_cons.ComputeLowerSlack(GetNode(VAR_SLK)));
  }

  /// Postsolve solution (primal + dual)
  void PostsolveSolutionEntry(const typename Base::LinkEntry& be) {
    SetDbl(be, CON_SRC, GetDbl(be, CON_TARGET));
  }

  /// Presolve basis
  ///
  /// From a range constraint's basis status,
  /// transfer it to the slack.
  /// Set the new constraint's status to 'equ'.
  void PresolveBasisEntry(const typename Base::LinkEntry& be) {
    SetInt(be, VAR_SLK, GetInt(be, CON_SRC));
    SetInt(be, CON_TARGET, (int)BasicStatus::equ);
  }
  /// Postsolve basis
  ///
  /// The reverse (forget solver's constraint status)
  void PostsolveBasisEntry(const typename Base::LinkEntry& be) {
    SetInt(be, CON_SRC, GetInt(be, VAR_SLK));
  }

  /// Presolve IIS
  void PresolveIISEntry(const typename Base::LinkEntry& ) {
    /// Should not need
  }

  /// Postsolve IIS
  ///
  /// Take slack's if set, otherwise the constraint's
  void PostsolveIISEntry(const typename Base::LinkEntry& be) {
    if (auto slk_iis = GetInt(be, VAR_SLK))
      SetInt(be, CON_SRC, slk_iis);
    else
      SetInt(be, CON_SRC, GetInt(be, CON_TARGET));
  }

  /// Mark Lazy/user cut: copy flag
  void PresolveLazyUserCutFlagsEntry(const typename Base::LinkEntry& be) {
    SetInt(be, CON_TARGET, GetInt(be, CON_SRC));
  }

  void PostsolveLazyUserCutFlagsEntry(const typename Base::LinkEntry& ) {
    /// Should not need
  }

protected:
  /// Get the model converter
  ModelConverter& GetMC() { return cvt_; }

  using Base::GetInt;
  using Base::SetInt;
  using Base::GetDbl;
  using Base::SetDbl;

  using Base::GetNode;

private:
  ModelConverter& cvt_;
};


/// Typedef RangeLinCon2Slack
template <class MC>
using RangeLinCon2Slack = RangeCon2Slack<MC, LinConRange>;

/// Typedef RangeQuadCon2Slack
template <class MC>
using RangeQuadCon2Slack = RangeCon2Slack<MC, QuadConRange>;

} // namespace pre


/// Converts proper range constraints to body+slack=ub,
/// otherwise to body </=/> rhs.
/// @param ModelConverter: the main converter class
/// @param AlgConBody: The \a Body argument of the
/// \a AlgebraicConstraint
template <class ModelConverter, class AlgConBody>
class RangeConstraintConverter :
    public BasicItemConverter<ModelConverter> {
public:
  /// Base class
  using Base = BasicItemConverter<ModelConverter>;
  /// Constructor
  RangeConstraintConverter(ModelConverter& mc) : Base(mc) { }
  /// Converted item type
  using ItemType = AlgebraicConstraint<AlgConBody, AlgConRange>;
  /// Resulting constraint type: LE
  using AlgConLE = AlgebraicConstraint< AlgConBody, AlgConRhs<-1> >;
  /// Resulting constraint type: equality
  using AlgConEQ = AlgebraicConstraint< AlgConBody, AlgConRhs<0> >;
  /// Resulting constraint type: GE
  using AlgConGE = AlgebraicConstraint< AlgConBody, AlgConRhs<1> >;

  /// Conversion
  ///
  /// Responsible for adding presolve links
  /// @param item: the item to be converted
  /// @param i: item index, used to create a presolve link
  void Convert(const ItemType& item, int i) {
    auto rr = Relate(item.lb(), item.ub());
    if (rr[0] && rr[1] && rr[2])
      ConvertRange(item, i);
    else
      ConvertWithRhs(item, rr);
  }

protected:
  using Base::GetMC;
  using RangeRelations = std::array<bool, 3>;
  RangeRelations Relate(double lb, double ub) {
    return {lb!=ub, lb>GetMC().MinusInfty(), ub<GetMC().Infty()};
  }
  void ConvertRange(const ItemType& item, int i) {
    GetMC().TurnOffAutoLinking();     // for range only
    auto slk = GetMC().AddVar(0.0, item.ub()-item.lb());
    auto body = item.GetBody();
    body.add_term(1.0, slk);
    AlgConEQ lceq { std::move(body), item.ub() };
    int i1 = GetMC().AddConstraint(std::move(lceq));
    GetSlackLink().AddEntry({i, i1, slk});
  }
  void ConvertWithRhs(const ItemType& item, RangeRelations rr) {
    if (rr[1] && !rr[2]) {
      GetMC().AddConstraint(
              AlgConGE( item.GetBody(), item.lb() ) );
    } else if (!rr[1] && rr[2]) {
      GetMC().AddConstraint(
              AlgConLE( item.GetBody(), item.ub() ) );
    } else if (rr[1] && rr[2]) {
      assert(item.lb()>=item.ub());
      GetMC().AddConstraint(
            AlgConEQ( item.GetBody(),
                      (item.lb()+item.ub()) / 2.0 ) );
    } // else, both are inf, forget
  }

  using SlackLink = pre::RangeLinCon2Slack<ModelConverter>;
  SlackLink& GetSlackLink() { return link_rng2slk_; }

private:
  SlackLink link_rng2slk_ {
    GetMC(), {
          &GET_CONSTRAINT_VALUE_NODE(ItemType),
          &GET_CONSTRAINT_VALUE_NODE(AlgConEQ),
          &this->GetMC().GetVarValueNode()
    } };
};


/// Typedef RangeLinearConstraintConverter
template <class MC>
using RangeLinearConstraintConverter =
    RangeConstraintConverter<MC, LinTerms>;

/// Typedef RangeQuadraticConstraintConverter
template <class MC>
using RangeQuadraticConstraintConverter =
    RangeConstraintConverter<MC, QuadAndLinTerms>;

} // namespace mp

#endif // RANGE_CON_H
