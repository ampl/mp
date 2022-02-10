#ifndef RANGE2SLACK_H
#define RANGE2SLACK_H

#include "mp/flat/std_constr.h"
#include "mp/presolve_bridge.h"

namespace mp {

namespace pre {

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

} // namespace mp

#endif // RANGE2SLACK_H
