#ifndef RANGE_CON_H
#define RANGE_CON_H

/*
 * Convert range constraints
 */

#include "mp/flat/convert/base.h"
#include "mp/flat/bridges/range2slack.h"

namespace mp {

template <class ModelConverter>
class RangeConstraintConverter :
    public BasicConverter<ModelConverter> {
public:
  /// Base class
  using Base = BasicConverter<ModelConverter>;
  /// Constructor
  RangeConstraintConverter(ModelConverter& mc) : Base(mc) { }

private:
  pre::RangeLinearConstraint2Slack bridge_rng2slk_ {
    this->GetPresolver(), {
          &GET_CONSTRAINT_VALUE_NODE(LinearConstraint), // just some constraints for now
          &GET_CONSTRAINT_VALUE_NODE(LinearDefiningConstraint),
          &this->GetMC().GetVarValueNode()
    } };
};

} // namespace mp

#endif // RANGE_CON_H
