#ifndef PRESOLVE_BRIDGE_H
#define PRESOLVE_BRIDGE_H

#include <vector>
#include <array>

#include "mp/common.h"
#include "presolve-node.h"


namespace mp {

namespace pre {

/// Index range in a bridge
using BridgeIndexRange = IndexRange;

/// Declare Presolver
class Presolver;

/// Macro for a list of pre- / postsolve method definitions
/// in a bridge.
/// Requires PRESOLVE_KIND defined to declare / define corr. methods
#define LIST_PRESOLVE_METHODS \
  PRESOLVE_KIND(Solution) \
  PRESOLVE_KIND(Basis) \
  PRESOLVE_KIND(IIS) \
  PRESOLVE_KIND(LazyUserCutFlags)
// ...

/// Bridge interface
///
/// A bridge is an array of value converters between nodes.
/// All converters are of the same type
class BasicBridge {
public:
  /// Constructor
  BasicBridge(Presolver& pre) : presolver_(pre) { }
  /// Virtual destructor
  virtual ~BasicBridge() = default;

  /// The below pre- / postsolves work on a range of bridge entries
  /// Postsolves should usually loop the range backwards

#define PRESOLVE_KIND(name) \
  virtual void Presolve ## name (BridgeIndexRange ) = 0; \
  virtual void Postsolve ## name(BridgeIndexRange ) = 0;

  LIST_PRESOLVE_METHODS

protected:
  /// Add a range of bridge entries to the Presolver's list
  /// Every derived bridge should call either of the next two
  /// whenever adding a bridge entry
  void RegisterBridgeIndex(int i)
  { RegisterBridgeIndexRange( {i, i+1} ); }
  void RegisterBridgeIndexRange(BridgeIndexRange );

private:
  Presolver& presolver_;
};


/// Bridge range: range of conversion specifiers of certain type
/// The bridge is specified as well
struct BridgeRange {
  /// Try & extend the range
  /// @return true iff extension worked,
  /// otherwise the caller would have to add the new range
  /// separately
  bool ExtendRange(BridgeRange br) {
    if (&b_ == &br.b_) {                  // same bridge
      if (ir_.end == br.ir_.beg) {        // and consecutive ranges
        ir_.end = br.ir_.end;
        return true;
      }
    }
    return false;
  }

  BasicBridge& b_;
  BridgeIndexRange ir_;
};


/// A specific bridge: each entry just copies a range of values
/// Useful to transfer values between NL and internal model,
/// internal model and solver, and in conversions preserving
/// all values
class CopyBridge : public BasicBridge {
public:
  /// Constructor
  CopyBridge(Presolver& pre) : BasicBridge(pre) { }

  /// Single bridge entry
  using BridgeEntry = std::pair<NodeRange, NodeRange>;

  /// Add entry
  /// Instead of a new entry, tries to extend the last one
  /// if exists
  void AddEntry(BridgeEntry be) {
    if (entries_.empty() ||
        !entries_.back().first.ExtendableBy(be.first) ||
        !entries_.back().second.ExtendableBy(be.second)) {
      entries_.push_back(be);             // Add new entry
      RegisterBridgeIndex(entries_.size()-1);
    } else {                              // Extend last entry
      entries_.back().first.ExtendBy(be.first);
      entries_.back().second.ExtendBy(be.second);
    }
  }

#undef PRESOLVE_KIND
#define PRESOLVE_KIND(name) \
  void Presolve ## name (BridgeIndexRange ir) override \
  { CopySrcDest(ir); } \
  void Postsolve ## name(BridgeIndexRange ir) override \
  { CopyDestSrc(ir); }

  LIST_PRESOLVE_METHODS

protected:
  /// Copy src -> dest for index range ir
  void CopySrcDest(BridgeIndexRange ir) {
    for (int i=ir.beg; i!=ir.end; ++i) {
      const auto& br = entries_[i];
      Copy(br.first, br.second);
    }
  }
  /// Copy src <- dest for index range ir. Loop backwards
  void CopyDestSrc(BridgeIndexRange ir) {
    for (int i=ir.end; (i--)!=ir.beg; ) {
      const auto& br = entries_[i];
      Copy(br.second, br.first);
    }
  }

private:
  std::vector<BridgeEntry> entries_;
};


/// A stub for bridges which process each entry individually,
/// in contrast to Presolve...(BridgeIndexRange ...)
/// Some of such bridges could be optimized to store
/// a range of transformations in each entry if that helps
///
/// Usage: a derived class should define types
/// BridgeEntry and NodeList, and some methods
/// Presolve...(const BridgeEntry& ) and Postsolve...(const BridgeEntry&).
/// TODO: Those methods which are not defined, just copy values
/// (which might be correct in _some_ cases).
/// Need a "default copy" method.
///
/// Using CRTP: https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
template <class Impl, class BridgeEntry>
class BasicIndivEntryBridge : public BasicBridge {
public:
  /// Constructor
  BasicIndivEntryBridge(Presolver& pre) :
    BasicBridge(pre) { }

  /// Add entry
  /// Instead of a new entry, tries to extend the last one
  /// if exists
  void AddEntry(BridgeEntry be) {
    entries_.push_back(be);             // Add new entry
    RegisterBridgeIndex(entries_.size()-1);
  }

#undef PRESOLVE_KIND
#define PRESOLVE_KIND(name) \
  void Presolve ## name (BridgeIndexRange ir) override { \
    for (int i=ir.beg; i!=ir.end; ++i) \
      MPD( Presolve ## name ## Entry(entries_.at(i)) ); } \
  void Postsolve ## name(BridgeIndexRange ir) override { \
    for (int i=ir.end; i--!=ir.beg; ) \
      MPD( Postsolve ## name ## Entry(entries_.at(i)) ); }

  LIST_PRESOLVE_METHODS


private:
  std::vector<BridgeEntry> entries_;
};


/// A static indiv entry bridge has a fixed number of ValueNodes
/// and indexes into them.
/// Generally NNodes==NIndexes
template <class Impl, int NNodes, int NIndexes>
class BasicStaticIndivEntryBridge :
    public BasicIndivEntryBridge<Impl, std::array<int, NIndexes> > {
public:
  /// Base class
  using Base = BasicIndivEntryBridge<Impl, std::array<int, NIndexes> >;
  /// Typedef NodeList
  using NodeList = std::array<ValueNode*, NNodes>;
  /// Typedef: BridgeEntry is just an array if node indexes
  using BridgeEntry = std::array<int, NIndexes>;

  /// Constructor
  BasicStaticIndivEntryBridge(Presolver& pre, const NodeList& ndl) :
    Base(pre), ndl_(ndl) { }

protected:
  /// Access whole node at specific index, const
  const ValueNode& GetNode(size_t i) const { return *ndl_.at(i); }
  /// Access whole node at specific index
  ValueNode& GetNode(size_t i) { return *ndl_.at(i); }

  /// Access int value at the node \a pos at the index
  /// stored in \a be[pos]
  ///
  /// @param be: a BridgeEntry
  /// @param pos: node number from 0..NNodes-1
  int GetInt(const BridgeEntry& be, size_t pos) const
  { return GetNode(pos).GetInt(be[pos]); }
  /// SetInt
  void SetInt(const BridgeEntry& be, size_t pos, int i)
  { GetNode(pos).SetInt(be[pos], i); }

  /// Access double value at the node \a pos at the index
  /// stored in \a be[pos]
  ///
  /// @param be: a BridgeEntry
  /// @param pos: node number from 0..NNodes-1
  double GetDbl(const BridgeEntry& be, size_t pos) const
  { return GetNode(pos).GetDbl(be[pos]); }
  /// SetDbl
  void SetDbl(const BridgeEntry& be, size_t pos, double i)
  { GetNode(pos).SetDbl(be[pos], i); }

private:
  NodeList ndl_;
};

} // namespace pre

} // namespace mp

#endif // PRESOLVE_BRIDGE_H
