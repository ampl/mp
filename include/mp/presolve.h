#ifndef PRESOLVE_H
#define PRESOLVE_H

/**
  @file presolve.h
  Implementation of value presolver
  */

#include <deque>

#include "presolve-node.h"
#include "presolve-bridge.h"


namespace mp {

namespace pre {


/// Array of bridge ranges
///
/// This defines a chain of transformations
class BridgeRangeList : private std::deque<BridgeRange> {
public:
  /// Base class typedef
  using Base = std::deque<BridgeRange>;

  /// Add another range to the list,
  /// tries to extend the last range if possible
  void Add(BridgeRange br) {
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


/// Value presolver implementation
class Presolver : public BasicPresolver {
public:
  /// Source nodes for model items
  const ModelValuesTerminal& GetSourceNodes() const { return src_; }
  ModelValuesTerminal& GetSourceNodes() { return src_; }

  /// Target nodes for model items
  const ModelValuesTerminal& GetTargetNodes() const { return dest_; }
  ModelValuesTerminal& GetTargetNodes() { return dest_; }

  /// Add (register) a bridge range
  /// This is normally called automatically from a bridge
  /// when an entry is added
  void Add(BridgeRange br) { brl_.Add(br); }

  /// Presolve solution (primal + dual)
  ModelValuesDbl PresolveSolution(const ModelValuesDbl& mvd) override
  { return RunPresolve(&BasicBridge::PresolveSolution, mvd); }
  /// Postsolve solution (primal + dual)
  ModelValuesDbl PostsolveSolution(const ModelValuesDbl& mvd) override
  { return RunPostsolve(&BasicBridge::PostsolveSolution, mvd); }

  /// Presolve basis (vars + cons)
  ModelValuesInt PresolveBasis(const ModelValuesInt& mvi) override
  { return RunPresolve(&BasicBridge::PresolveBasis, mvi); }
  /// Postsolve solution (vars + cons)
  ModelValuesInt PostsolveBasis(const ModelValuesInt& mvi) override
  { return RunPostsolve(&BasicBridge::PostsolveBasis, mvi); }

  /// Presolve IIS (vars + cons)
  ModelValuesInt PresolveIIS(const ModelValuesInt& mvi) override
  { return RunPresolve(&BasicBridge::PresolveIIS, mvi); }
  /// Postsolve IIS (vars + cons)
  ModelValuesInt PostsolveIIS(const ModelValuesInt& mvi) override
  { return RunPostsolve(&BasicBridge::PostsolveIIS, mvi); }

  /// Presolve LazyUserCutFlags (vars + cons)
  ModelValuesInt PresolveLazyUserCutFlags(const ModelValuesInt& mvi) override
  { return RunPresolve(&BasicBridge::PresolveLazyUserCutFlags, mvi); }

protected:
  /// Helper type: virtual member function pointer
  using BridgeFn = void (BasicBridge::*)(BridgeIndexRange);

  /// Generic value presolve loop: forward
  template <class ModelValues>
  ModelValues RunPresolve(BridgeFn fn, const ModelValues& mv) const {
    src_ = mv;
    for (const auto& br: brl_)
      (br.b_.*fn)(br.ir_);
    return dest_;
  }

  /// Generic value postsolve loop: backward
  template <class ModelValues>
  ModelValues RunPostsolve(BridgeFn fn, const ModelValues& mv) const {
    dest_ = mv;
    for (auto rit=brl_.rbegin(); rit!=brl_.rend(); ++rit)
      (rit->b_.*fn)(rit->ir_);
    return src_;
  }

private:
  mutable ModelValuesTerminal src_{std::string("src")},
    dest_{std::string("dest")};
  BridgeRangeList brl_;
};


/// Implement bridge range registration in global Presolver
inline void
BasicBridge::RegisterBridgeIndexRange(BridgeIndexRange bir)
{ presolver_.Add( { *this, bir } ); }

} // namespace pre

} // namespace mp

#endif // PRESOLVE_H
