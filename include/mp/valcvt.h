#ifndef VALUE_PRESOLVE_H
#define VALUE_PRESOLVE_H

/**
  @file valcvt.h.
  Implementation of value presolver
  */

#include <deque>

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
  /// Source nodes for model items
  const ModelValuesTerminal& GetSourceNodes() const { return src_; }
  ModelValuesTerminal& GetSourceNodes() { return src_; }

  /// Target nodes for model items
  const ModelValuesTerminal& GetTargetNodes() const { return dest_; }
  ModelValuesTerminal& GetTargetNodes() { return dest_; }

  /// Add (register) a link range
  /// This is normally called automatically from a link
  /// when an entry is added
  void Add(LinkRange br) { brl_.Add(br); }

  /// Presolve solution (primal + dual)
  ModelValuesDbl PresolveSolution(const ModelValuesDbl& mvd) override
  { return RunPresolve(&BasicValcvtLink::PresolveSolution, mvd); }
  /// Postsolve solution (primal + dual)
  ModelValuesDbl PostsolveSolution(const ModelValuesDbl& mvd) override
  { return RunPostsolve(&BasicValcvtLink::PostsolveSolution, mvd); }

  /// Presolve basis (vars + cons)
  ModelValuesInt PresolveBasis(const ModelValuesInt& mvi) override
  { return RunPresolve(&BasicValcvtLink::PresolveBasis, mvi); }
  /// Postsolve solution (vars + cons)
  ModelValuesInt PostsolveBasis(const ModelValuesInt& mvi) override
  { return RunPostsolve(&BasicValcvtLink::PostsolveBasis, mvi); }

  /// Presolve IIS (vars + cons)
  ModelValuesInt PresolveIIS(const ModelValuesInt& mvi) override
  { return RunPresolve(&BasicValcvtLink::PresolveIIS, mvi); }
  /// Postsolve IIS (vars + cons)
  ModelValuesInt PostsolveIIS(const ModelValuesInt& mvi) override
  { return RunPostsolve(&BasicValcvtLink::PostsolveIIS, mvi); }

  /// Presolve LazyUserCutFlags (vars + cons)
  ModelValuesInt PresolveLazyUserCutFlags(const ModelValuesInt& mvi) override
  { return RunPresolve(&BasicValcvtLink::PresolveLazyUserCutFlags, mvi); }

protected:
  /// Helper type: virtual member function pointer
  using LinkFn = void (BasicValcvtLink::*)(ValcvtLinkIndexRange);

  /// Generic value presolve loop: forward
  template <class ModelValues>
  ModelValues RunPresolve(LinkFn fn, const ModelValues& mv) const {
    src_ = mv;
    for (const auto& br: brl_)
      (br.b_.*fn)(br.ir_);
    return dest_;
  }

  /// Generic value postsolve loop: backward
  template <class ModelValues>
  ModelValues RunPostsolve(LinkFn fn, const ModelValues& mv) const {
    dest_ = mv;
    for (auto rit=brl_.rbegin(); rit!=brl_.rend(); ++rit)
      (rit->b_.*fn)(rit->ir_);
    return src_;
  }

private:
  mutable ModelValuesTerminal src_{std::string("src")},
    dest_{std::string("dest")};
  LinkRangeList brl_;
};


/// Implement link range registration in global Presolver
inline void
BasicValcvtLink::RegisterLinkIndexRange(ValcvtLinkIndexRange bir)
{ vsalue_presolver_.Add( { *this, bir } ); }

} // namespace pre

} // namespace mp

#endif // VALUE_PRESOLVE_H
