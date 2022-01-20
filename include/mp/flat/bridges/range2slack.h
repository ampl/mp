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
  /// Typedef BridgeEntry is created in Base as std::array<int, 3>
  /// (3 is the base class template parameter)
  /// and means the following indexes: {range_con, lin_con, slack_index}

  /// Presolve solution (primal + dual)
  /// Duals: just copy
  /// Primals: compute slack TODO
  void PresolveSolutionEntry(const BridgeEntry& be) {
    SetDbl(be, 1, GetDbl(be, 0));
  }
  /// Postsolve solution (primal + dual)
  void PostsolveSolutionEntry(const BridgeEntry& be) {
    SetDbl(be, 0, GetDbl(be, 1));
  }

  /// Presolve basis
  ///
  /// From a range constraint's basis status,
  /// transfer it to the slack.
  /// Set the new constraint's status to 'equ'
  void PresolveBasisEntry(const BridgeEntry& be) {
    SetInt(be, 2, GetInt(be, 0));
    SetInt(be, 1, (int)BasicStatus::equ);
  }
  /// Postsolve basis
  ///
  /// The reverse (forget solver's constraint status)
  void PostsolveBasisEntry(const BridgeEntry& be) {
    SetInt(be, 0, GetInt(be, 2));
  }

  /// Presolve IIS
  void PresolveIISEntry(const BridgeEntry& ) {
    /// Should not need
  }
  /// Postsolve IIS
  ///
  /// Take slack's if set, otherwise the constraint's
  void PostsolveIISEntry(const BridgeEntry& be) {
    if (auto slk_iis = GetInt(be, 2))
      SetInt(be, 0, slk_iis);
    else
      SetInt(be, 0, GetInt(be, 1));
  }

};

} // namespace pre

} // namespace mp

#endif // RANGE2SLACK_H
