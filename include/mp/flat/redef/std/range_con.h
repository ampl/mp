#ifndef RANGE_CON_H
#define RANGE_CON_H

/*
 * Convert range constraints
 */

#include "mp/flat/redef_base.h"
#include "mp/presolve-bridge.h"
#include "mp/flat/constraints_std.h"

namespace mp {

namespace pre {

/// Presolve bridge between RangeLinCon and LinCon+Slack.
/// TODO dependency inversion #164
/// (FlatConverter needs just a BasicBridge*? Entry type?)
template <class ModelConverter>
class RangeLinearConstraint2Slack :
    public BasicStaticIndivEntryBridge<
      RangeLinearConstraint2Slack<ModelConverter>, 3, 3> {
public:
  /// Base class
  using Base = BasicStaticIndivEntryBridge<
    RangeLinearConstraint2Slack<ModelConverter>, 3, 3>;
  /// Constructor
  RangeLinearConstraint2Slack(ModelConverter& cvt,
                              const typename Base::NodeList& ndl) :
    Base(cvt.GetPresolver(), ndl), cvt_(cvt) { }

  /// Define pre- / postsolve methods for individual bridge entries
  /// Typedef BridgeEntry is created in Base as std::array<int, 3>
  /// (3 is the base class template parameter)
  /// and means the following indexes: {range_con, lin_con, slack_index}

  /// Presolve solution (primal + dual)
  /// Duals: just copy
  /// Primals: compute slack.
  /// Slacks are probably only needed for Gurobi LP warmstart
  /// Who does not need them, can use unpresolved primal solution
  void PresolveSolutionEntry(const typename Base::BridgeEntry& be) {
    SetDbl(be, 1, GetDbl(be, 0));
    const auto& orig_cons =
        GetMC().template GetConstraint<RangeLinCon>(be[0]);
    SetDbl(be, 2, orig_cons.ComputeLowerSlack(GetNode(2)));
  }
  /// Postsolve solution (primal + dual)
  void PostsolveSolutionEntry(const typename Base::BridgeEntry& be) {
    SetDbl(be, 0, GetDbl(be, 1));
  }

  /// Presolve basis
  ///
  /// From a range constraint's basis status,
  /// transfer it to the slack.
  /// Set the new constraint's status to 'equ'
  void PresolveBasisEntry(const typename Base::BridgeEntry& be) {
    SetInt(be, 2, GetInt(be, 0));
    SetInt(be, 1, (int)BasicStatus::equ);
  }
  /// Postsolve basis
  ///
  /// The reverse (forget solver's constraint status)
  void PostsolveBasisEntry(const typename Base::BridgeEntry& be) {
    SetInt(be, 0, GetInt(be, 2));
  }

  /// Presolve IIS
  void PresolveIISEntry(const typename Base::BridgeEntry& ) {
    /// Should not need
  }
  /// Postsolve IIS
  ///
  /// Take slack's if set, otherwise the constraint's
  void PostsolveIISEntry(const typename Base::BridgeEntry& be) {
    if (auto slk_iis = GetInt(be, 2))
      SetInt(be, 0, slk_iis);
    else
      SetInt(be, 0, GetInt(be, 1));
  }

  /// Mark Lazy/user cut: copy flag
  void PresolveLazyUserCutFlagsEntry(const typename Base::BridgeEntry& be) {
    SetInt(be, 1, GetInt(be, 0));
  }
  void PostsolveLazyUserCutFlagsEntry(const typename Base::BridgeEntry& ) {
    /// Should not need
  }

protected:
  ModelConverter& GetMC() { return cvt_; }

  using Base::GetInt;
  using Base::SetInt;
  using Base::GetDbl;
  using Base::SetDbl;

  using Base::GetNode;

private:
  ModelConverter& cvt_;
};

} // namespace pre


/// Converts proper range linear constraints to c'x-slack=ub,
/// otherwise to c'x ? rhs.
template <class ModelConverter>
class RangeConstraintConverter :
    public BasicItemConverter<ModelConverter> {
public:
  /// Base class
  using Base = BasicItemConverter<ModelConverter>;
  /// Constructor
  RangeConstraintConverter(ModelConverter& mc) : Base(mc) { }
  /// Converted item type
  using ItemType = RangeLinCon;

  /// Conversion
  ///
  /// Responsible for adding presolve bridges
  /// @param item: the item to be converted
  /// @param i: item index, used to create a presolve bridge
  void Convert(const ItemType& item, int i) {
    auto rr = Relate(item.lb(), item.ub());
    if (rr[0] && rr[1] && rr[2])
      ConvertRange(item, i);
    else
      ConvertWithRhs(item, i, rr);
  }

protected:
  using Base::GetMC;
  using RangeRelations = std::array<bool, 3>;
  RangeRelations Relate(double lb, double ub) {
    return {lb<ub, lb>GetMC().MinusInfty(), ub<GetMC().Infty()};
  }
  void ConvertRange(const ItemType& item, int i) {
    auto slk = GetMC().AddVar(0.0, item.ub()-item.lb());
    LinConEQ lceq { {item.coefs(), item.vars()}, item.ub() };
    lceq.add_term(1.0, slk);
    int i1 = GetMC().AddConstraint(std::move(lceq));
    GetSlackBridge().AddEntry({i, i1, slk});
  }
  void ConvertWithRhs(const ItemType& item, int i, RangeRelations rr) {
    pre::NodeRange nr;              // target node+index
    if (rr[1] && !rr[2]) {
      nr = GetMC().AddConstraint(
              LinConGE( {item.coefs(), item.vars()}, item.lb() ) );
    } else if (!rr[1] && rr[2]) {
      nr = GetMC().AddConstraint(
              LinConLE( {item.coefs(), item.vars()}, item.ub() ) );
    } else if (rr[1] && rr[2]) {
      assert(item.lb()>=item.ub()); // TODO have an option for eps tolerance
      nr = GetMC().AddConstraint(
            LinConEQ( {item.coefs(), item.vars()},
                      (item.lb()+item.ub()) / 2.0 ) );
    } // else, both are inf, forget
    GetMC().GetCopyBridge().AddEntry(
          { GET_CONSTRAINT_VALUE_NODE(ItemType).Select(i), nr });
  }

  using SlackBridge = pre::RangeLinearConstraint2Slack<ModelConverter>;
  SlackBridge& GetSlackBridge() { return bridge_rng2slk_; }

private:
  SlackBridge bridge_rng2slk_ {
    GetMC(), {
          &GET_CONSTRAINT_VALUE_NODE(ItemType), // just some constraints for now
          &GET_CONSTRAINT_VALUE_NODE(LinConEQ),
          &this->GetMC().GetVarValueNode()
    } };
};

} // namespace mp

#endif // RANGE_CON_H
