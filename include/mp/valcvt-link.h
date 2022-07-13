#ifndef VALUE_PRESOLVE_LINK_H
#define VALUE_PRESOLVE_LINK_H

#include <array>
#include <vector>
#include <deque>
#include <algorithm>
#include <mp/arrayref.h>

#include "mp/common.h"
#include "valcvt-node.h"


namespace mp {

namespace pre {

/// Index range in a link
using LinkIndexRange = IndexRange;

/// Declare ValuePresolver
class ValuePresolver;

/// Macro for a list of pre- / postsolve method definitions
/// in a link.
/// Requires PRESOLVE_KIND defined to declare / define corr. methods.
/// Generic(Dbl/Int): maximizes suffix value among non-0.
/// Example 1 (presolve): expression exp(y) is used in various places
/// and marked with different values of .funcpieces.
/// The largest is chosen.
/// Example 2 (postsolve): IIS membership value for a converted
/// high-level constraint can be reported as the maximum of those
/// for its low-level representation.
#define LIST_PRESOLVE_METHODS \
  PRESOLVE_KIND(GenericDbl, double) \
  PRESOLVE_KIND(GenericInt, int) \
  PRESOLVE_KIND(Solution, double) \
  PRESOLVE_KIND(Basis, int) \
  PRESOLVE_KIND(IIS, int) \
  PRESOLVE_KIND(LazyUserCutFlags, int)
// ...


/// ValuePresolveLink interface
///
/// A link is a node of the structural conversion graph.
/// It contains an array of concrete nodes,
/// each describes a transformation from
/// some source vars+cons+objs values into/from some target ones.
/// All concrete nodes are of the same type.
class BasicLink {
public:
  /// Constructor
  BasicLink(ValuePresolver& pre) : value_presolver_(pre) { }
  /// Virtual destructor
  virtual ~BasicLink() = default;

  /// The below pre- / postsolve methods
  /// work on a range of link entries.
  /// Postsolves should usually loop the range backwards
#define PRESOLVE_KIND(name, ValType) \
  virtual void Presolve ## name (LinkIndexRange ) = 0; \
  virtual void Postsolve ## name(LinkIndexRange ) = 0;

  LIST_PRESOLVE_METHODS

  /// Get source/target nodes.
  /// BasicLink only defines an interface for this.

  /// Typedef variable container
  using VarList = ArrayRef<int>;

  /// Typedef constraint container
  using ConList = ArrayRef<NodeRange>;

  /// Typedef objective container
  using ObjList = ArrayRef<int>;

  /// Get source variables
  virtual VarList GetSrcVars() const { return {}; }
  /// Get source constraint
  virtual ConList GetSrcCon() const { return {}; }
  /// Get source objective
  virtual ObjList GetSrcObj() const { return {}; }

  /// Get target variables
  virtual VarList GetTargetVars() const { return {}; }
  /// Get target constraint
  virtual ConList GetTargetCon() const { return {}; }
  /// Get target objective
  virtual ObjList GetTargetObj() const { return {}; }

  /// Add source/target nodes.
  /// BasicLink only defines an interface for this.

  /// Add source variables
  virtual void AddSrcVars(ArrayRef<int> ) { }
  /// Add source constraint
  virtual void AddSrcCon(NodeRange ) { }
  /// Add source objective
  virtual void AddSrcObj(ArrayRef<int> ) { }

  /// Add target variables
  virtual void AddTargetVars(ArrayRef<int> ) { }
  /// Add target constraint
  virtual void AddTargetCon(NodeRange ) { }
  /// Add target objective
  virtual void AddTargetObj(ArrayRef<int> ) { }


protected:
  /// Add a range of link entries to the Presolver's list.
  /// Every derived link should call either of the next two methods
  /// whenever adding a link entry.
  /// Version 1: add single link
  void RegisterLinkIndex(int i)
  { RegisterLinkIndexRange( {i, i+1} ); }
  /// Version 2: add a range of links
  void RegisterLinkIndexRange(LinkIndexRange );

private:
  ValuePresolver& value_presolver_;
};


/// Link range: range of conversion specifiers of certain type.
/// The link is specified as well
struct LinkRange {
  /// Try & extend the range
  /// @return true iff extension worked,
  /// otherwise the caller would have to add the new range
  /// separately
  bool ExtendRange(LinkRange br) {
    if (&b_ == &br.b_) {                  // same link?
      if (ir_.end == br.ir_.beg) {        // and consecutive ranges?
        ir_.end = br.ir_.end;
        return true;
      }
    }
    return false;
  }

  /// Reference to BasicLink
  BasicLink& b_;
  /// Link index range
  LinkIndexRange ir_;
};


/// A specific link: each entry just copies a range of values.
/// Useful to transfer values between NL and internal model,
/// internal model and solver, and in conversions preserving
/// all values
class CopyLink : public BasicLink {
public:
  /// Constructor
  CopyLink(ValuePresolver& pre) : BasicLink(pre) { }

  /// Single link entry,
  /// stores src + dest ranges
  using LinkEntry = std::pair<NodeRange, NodeRange>;

  /// Collection of entries
  using CollectionOfEntries = std::deque<LinkEntry>;

  /// Add entry.
  /// Instead of a new entry, tries to extend the last one
  /// if exists
  void AddEntry(LinkEntry be) {
    if (entries_.empty() ||
        !entries_.back().first.ExtendableBy(be.first) ||
        !entries_.back().second.ExtendableBy(be.second)) {
      entries_.push_back(be);             // Add new entry
      RegisterLinkIndex(entries_.size()-1);
    } else {                              // Extend last entry
      entries_.back().first.ExtendBy(be.first);
      entries_.back().second.ExtendBy(be.second);
    }
  }

  /// Retrieve entries
  const CollectionOfEntries& GetEntries() const { return entries_; }

  /// Copy everything, MaxAmongNon0 should not apply
#undef PRESOLVE_KIND
#define PRESOLVE_KIND(name, ValType) \
  void Presolve ## name (LinkIndexRange ir) override \
  { CopySrcDest <ValType>(ir); } \
  void Postsolve ## name(LinkIndexRange ir) override \
  { CopyDestSrc <ValType>(ir); }

  LIST_PRESOLVE_METHODS

protected:
  /// Copy src -> dest for index range ir
  template <class T>
  void CopySrcDest(LinkIndexRange ir) {
    for (int i=ir.beg; i!=ir.end; ++i) {
      const auto& br = entries_[i];
      Copy<T>(br.first, br.second);
    }
  }
  /// Copy src <- dest for index range ir. Loop backwards
  template <class T>
  void CopyDestSrc(LinkIndexRange ir) {
    for (int i=ir.end; (i--)!=ir.beg; ) {
      const auto& br = entries_[i];
      Copy<T>(br.second, br.first);
    }
  }

private:
  CollectionOfEntries entries_;
};


/// A specific link: each entry just copies from 1 source value
/// into a range of values.
/// Useful to transfer values into expression subtree,
/// in conversions, e.g., form 1 var / constraint / objective
/// into several new ones
class One2ManyLink : public BasicLink {
public:
  /// Constructor
  One2ManyLink(ValuePresolver& pre) : BasicLink(pre) { }

  /// Single link entry,
  /// stores src + dest ranges.
  /// HOWEVER the source range keeps just 1 node.
  using LinkEntry = std::pair<NodeRange, NodeRange>;

  /// Collection of entries
  using CollectionOfEntries = std::deque<LinkEntry>;

  /// Add entry.
  /// Instead of a new entry, tries to extend the last one
  /// if exists
  void AddEntry(LinkEntry be) {
    assert(be.first.IsSingleIndex());
    if (entries_.empty() ||
        entries_.back().first!=be.first ||
        !entries_.back().second.ExtendableBy(be.second)) {
      entries_.push_back(be);             // Add new entry
      RegisterLinkIndex(entries_.size()-1);
    } else {                              // Extend last entry
      entries_.back().second.ExtendBy(be.second);
    }
  }

  /// Retrieve entries
  const CollectionOfEntries& GetEntries() const { return entries_; }

  /// All pre- / postsolves just take max from non-0
#undef PRESOLVE_KIND
#define PRESOLVE_KIND(name, ValType) \
  void Presolve ## name (LinkIndexRange ir) override \
  { DistributeFromSrc2Dest <ValType>(ir); } \
  void Postsolve ## name(LinkIndexRange ir) override \
  { CollectFromDest2Src <ValType>(ir); }

  LIST_PRESOLVE_METHODS

protected:
  /// Distribute values of type T from nr1 to nr2
  template <class T>
  void Distr(NodeRange nr1, NodeRange nr2) {
    assert(nr1.IsSingleIndex());
    auto val = nr1.GetValueNode()->
        GetVal<T>(nr1.GetSingleIndex());
    for (auto i=nr2.GetIndexRange().beg;
         i!=nr2.GetIndexRange().end; ++i)
      nr2.GetValueNode()->SetVal(i, val);
  }

  /// Collect values of type T from nr2 to nr1
  template <class T>
  void Collect(NodeRange nr1, NodeRange nr2) {
    assert(nr1.IsSingleIndex());
    auto nr1_idx = nr1.GetSingleIndex();
    auto& vec2 = nr2.GetValueNode()->GetValVec<T>();
    for (auto i=nr2.GetIndexRange().beg;
         i!=nr2.GetIndexRange().end; ++i)
      nr1.GetValueNode()->SetVal(nr1_idx, vec2.at(i));
  }

  /// Src -> dest for the entries index range ir
  template <class T>
  void DistributeFromSrc2Dest(LinkIndexRange ir) {
    for (int i=ir.beg; i!=ir.end; ++i) {
      const auto& br = entries_[i];
      Distr<T>(br.first, br.second);
    }
  }

  /// Collect src <- dest for index range ir. Loop backwards
  template <class T>
  void CollectFromDest2Src(LinkIndexRange ir) {
    for (int i=ir.end; (i--)!=ir.beg; ) {
      const auto& br = entries_[i];
      Collect<T>(br.first, br.second);
    }
  }


private:
  CollectionOfEntries entries_;
};


/// A stub for links which process each entry individually,
/// in contrast to Presolve...(LinkIndexRange ...).
/// Some of such links could be optimized to store
/// a range of transformations in each entry if that helps
/// (but this could be an overoptimization?)
///
/// Usage: a derived class should define some methods
/// Presolve...(const LinkEntry& ) and Postsolve...(const LinkEntry&).
///
/// TODO: we could default as follows:
/// Those methods which are not defined, just copy values
/// (which might be correct in _some_ cases).
/// Then, need a "default copy" method.
///
/// Using CRTP:
/// https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
template <class Impl, class LinkEntry>
class BasicIndivEntryLink : public BasicLink {
public:
  /// Constructor
  BasicIndivEntryLink(ValuePresolver& pre) :
    BasicLink(pre) { }

  /// Add entry
  /// Instead of a new entry, tries to extend the last one
  /// if exists
  void AddEntry(LinkEntry be) {
    entries_.push_back(be);             // Add new entry
    RegisterLinkIndex(entries_.size()-1);
  }

  /// Pre- / postsolve loops over link entries
  /// and calls the derived class' method for each.
#undef PRESOLVE_KIND
#define PRESOLVE_KIND(name, ValType) \
  void Presolve ## name (LinkIndexRange ir) override { \
    for (int i=ir.beg; i!=ir.end; ++i) \
      MPD( Presolve ## name ## Entry(entries_.at(i)) ); } \
  void Postsolve ## name(LinkIndexRange ir) override { \
    for (int i=ir.end; i--!=ir.beg; ) \
      MPD( Postsolve ## name ## Entry(entries_.at(i)) ); }

  LIST_PRESOLVE_METHODS


private:
  std::deque<LinkEntry> entries_;
};


/// A static indiv entry link has a fixed number (\a NNodes)
/// of ValueNodes and indexes (\a NIndexes) into them.
/// Generally NNodes==NIndexes
template <class Impl, int NNodes, int NIndexes>
class BasicStaticIndivEntryLink :
    public BasicIndivEntryLink<Impl, std::array<int, NIndexes> > {
public:
  /// Base class
  using Base = BasicIndivEntryLink<Impl, std::array<int, NIndexes> >;
  /// Typedef NodeList.
  /// Remember, a node is a node of the transformation graph
  /// and can store arrays of values for a list of model items
  using NodeList = std::array<ValueNode*, NNodes>;
  /// Typedef: LinkEntry is just an array if node indexes
  using LinkEntry = std::array<int, NIndexes>;

  /// Constructor
  BasicStaticIndivEntryLink(ValuePresolver& pre, const NodeList& ndl) :
    Base(pre), ndl_(ndl) { }

protected:
  /// Access whole node at specific index, const
  const ValueNode& GetNode(size_t i) const { return *ndl_.at(i); }
  /// Access whole node at specific index
  ValueNode& GetNode(size_t i) { return *ndl_.at(i); }

  /// Access int value at the node \a pos at the index
  /// stored in \a be[pos]
  ///
  /// @param be: a LinkEntry
  /// @param pos: node number from 0..NNodes-1
  int GetInt(const LinkEntry& be, size_t pos) const
  { return GetNode(pos).GetInt(be[pos]); }
  /// SetInt
  void SetInt(const LinkEntry& be, size_t pos, int i)
  { GetNode(pos).SetInt(be[pos], i); }

  /// Access double value at the node \a pos at the index
  /// stored in \a be[pos]
  ///
  /// @param be: a LinkEntry
  /// @param pos: node number from 0..NNodes-1
  double GetDbl(const LinkEntry& be, size_t pos) const
  { return GetNode(pos).GetDbl(be[pos]); }
  /// SetDbl
  void SetDbl(const LinkEntry& be, size_t pos, double i)
  { GetNode(pos).SetDbl(be[pos], i); }

private:
  NodeList ndl_;
};

} // namespace pre

} // namespace mp

#endif // VALUE_PRESOLVE_LINK_H
