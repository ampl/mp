#ifndef RANGE2SLACK_H
#define RANGE2SLACK_H

#include "mp/presolve_bridge.h"

namespace mp {

namespace pre {

class RangeLinearConstraint2Slack :
    public BasicStaticIndivEntryBridge<RangeLinearConstraint2Slack, 3, 3> {
public:
  /// Base class
  using Base = BasicStaticIndivEntryBridge<RangeLinearConstraint2Slack, 3, 3>;
  /// Constructor
  RangeLinearConstraint2Slack(Presolver& pre, const Base::NodeList& ndl) :
    Base(pre, ndl) { }

  /// Define pre- / postsolve methods for individual bridge entries

  /// Presolve solution (primal + dual)
  void PresolveSolutionEntry(const BridgeEntry&) {}
  /// Postsolve solution (primal + dual)
  void PostsolveSolutionEntry(const BridgeEntry&) {}

  /// Presolve basis
  void PresolveBasisEntry(const BridgeEntry&) {}
  /// Postsolve basis
  void PostsolveBasisEntry(const BridgeEntry&) {}

};

} // namespace pre

} // namespace mp

#endif // RANGE2SLACK_H
