#ifndef RANGE_CON_H
#define RANGE_CON_H

/*
 * Convert range constraints
 */

#include "mp/flat/item_cvt/item_cvt_base.h"
#include "mp/flat/bridges/range2slack.h"
#include "mp/flat/std_constr.h"
#include "mp/flat/constraint_keeper.h"

namespace mp {

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
