#ifndef VALUE_PRESOLVE_H
#define VALUE_PRESOLVE_H

/**
  @file valcvt.h.
  Implementation of value presolver
  */

#include <deque>
#include <unordered_set>

#include "valcvt-node.h"
#include "valcvt-link.h"


namespace mp {

namespace pre {


/// Array of link ranges
///
/// This defines a chain of transformations
class LinkRangeList : private std::deque<LinkRange> {
public:
  /// Base class typedef
  using Base = std::deque<LinkRange>;

  /// Add another range to the list,
  /// tries to extend the last range if possible
  void Add(LinkRange br) {
    if (empty() || !back().ExtendRange(br))
      push_back(br);
  }

  /// using Base::empty()
  using Base::empty;
  /// iterator access: begin
  using Base::begin;
  /// iterator access: end
  using Base::end;
  /// iterator access: rbegin
  using Base::rbegin;
  /// iterator access: rend
  using Base::rend;

protected:
  /// using Base::back()
  using Base::back;
};


/// Value presolver implementation.
/// A ValuePresolver converts solutions and suffixes
/// between the original model and the presolved one.
class ValuePresolver : public BasicValuePresolver {
public:
  /// Source nodes of the conversion graph, const
  const ModelValuesTerminal& GetSourceNodes() const { return src_; }
  /// Source nodes of the conversion graph
  ModelValuesTerminal& GetSourceNodes() { return src_; }

  /// Target nodes of the conversion graph, const
  const ModelValuesTerminal& GetTargetNodes() const { return dest_; }
  /// Target nodes of the conversion graph
  ModelValuesTerminal& GetTargetNodes() { return dest_; }

  /// Add (register) a link range
  /// This is normally called automatically from a link
  /// when an entry is added
  void Add(LinkRange br) { brl_.Add(br); }

  /// Presolve solution (primal + dual)
  ModelValuesDbl PresolveSolution(const ModelValuesDbl& mvd) override
  { return RunPresolve(&BasicLink::PresolveSolution, mvd); }
  /// Postsolve solution (primal + dual)
  ModelValuesDbl PostsolveSolution(const ModelValuesDbl& mvd) override
  { return RunPostsolve(&BasicLink::PostsolveSolution, mvd); }

  /// Presolve basis (vars + cons)
  ModelValuesInt PresolveBasis(const ModelValuesInt& mvi) override
  { return RunPresolve(&BasicLink::PresolveBasis, mvi); }
  /// Postsolve solution (vars + cons)
  ModelValuesInt PostsolveBasis(const ModelValuesInt& mvi) override
  { return RunPostsolve(&BasicLink::PostsolveBasis, mvi); }

  /// Presolve IIS (vars + cons)
  ModelValuesInt PresolveIIS(const ModelValuesInt& mvi) override
  { return RunPresolve(&BasicLink::PresolveIIS, mvi); }
  /// Postsolve IIS (vars + cons)
  ModelValuesInt PostsolveIIS(const ModelValuesInt& mvi) override
  { return RunPostsolve(&BasicLink::PostsolveIIS, mvi); }

  /// Presolve LazyUserCutFlags (vars + cons)
  ModelValuesInt PresolveLazyUserCutFlags(const ModelValuesInt& mvi) override
  { return RunPresolve(&BasicLink::PresolveLazyUserCutFlags, mvi); }

  /// Register a ValueNode*
  void Register(ValueNode* pvn) override
  { assert(val_nodes_.insert(pvn).second); }

  /// Deregister a ValueNode*
  void Deregister(ValueNode* pvn) override { assert(val_nodes_.erase(pvn)); }


protected:
  /// Helper type: virtual member function pointer
  using LinkFn = void (BasicLink::*)(LinkIndexRange);

  /// Generic value presolve loop: forward
  template <class ModelValues>
  ModelValues RunPresolve(LinkFn fn, const ModelValues& mv) const {
    CleanUpValueNodes();
    src_ = mv;
    for (const auto& br: brl_)
      (br.b_.*fn)(br.ir_);
    return dest_;
  }

  /// Generic value postsolve loop: backward
  template <class ModelValues>
  ModelValues RunPostsolve(LinkFn fn, const ModelValues& mv) const {
    CleanUpValueNodes();
    dest_ = mv;
    for (auto rit=brl_.rbegin(); rit!=brl_.rend(); ++rit)
      (rit->b_.*fn)(rit->ir_);
    return src_;
  }

  /// Clean up value nodes for new propagation
  void CleanUpValueNodes() const {
    for (auto pvn: val_nodes_)
      pvn->CleanUp();
  }

private:
  /// val_nodes_ should be before src_ / dest_
  std::unordered_set<ValueNode*> val_nodes_;
  /// val_nodes_ should be before src_ / dest_
  mutable ModelValuesTerminal
    src_{*this, "source_item_nodes__"},
    dest_{*this, "target_item_nodes__"};
  LinkRangeList brl_;
};


/// Implement link range registration in global Presolver
inline void
BasicLink::RegisterLinkIndexRange(LinkIndexRange bir)
{ value_presolver_.Add( { *this, bir } ); }

} // namespace pre

} // namespace mp

#endif // VALUE_PRESOLVE_H
