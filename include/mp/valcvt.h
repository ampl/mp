#ifndef VALUE_PRESOLVE_H
#define VALUE_PRESOLVE_H

/**
  Implementation of value presolver
  */

#include <deque>
#include <unordered_set>
#include <functional>

#include "mp/utils-string.h"
#include "mp/util-json-write.hpp"

#include "valcvt-node.h"
#include "valcvt-link.h"
#include "mp/flat/converter_model_base.h"


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
    if (empty() || !back().TryExtendBy(br))
      push_back(br);
  }

  /// using Base::empty()
  using Base::empty;
  /// using Base::size()
  using Base::size;
  /// using operator[]
  using Base::operator[];
  /// iterator access: begin
  using Base::begin;
  /// iterator access: end
  using Base::end;
  /// iterator access: rbegin
  using Base::rbegin;
  /// iterator access: rend
  using Base::rend;
  /// using Base::back()
  using Base::back;
};


/// Value presolver intermediate implementation.
/// A ValuePresolver converts solutions and suffixes
/// between the original model and the presolved one.
class ValuePresolverImpl : public BasicValuePresolver {
public:
  /// Constructor
  ValuePresolverImpl(Env& env, BasicLogger& bts) :
    BasicValuePresolver(env), bts_(bts) { }

  /// Source nodes of the conversion graph, const
  const ModelValuesTerminal& GetSourceNodes() const { return src_; }
  /// Source nodes of the conversion graph
  ModelValuesTerminal& GetSourceNodes() { return src_; }

  /// Target nodes of the conversion graph, const
  const ModelValuesTerminal& GetTargetNodes() const { return dest_; }
  /// Target nodes of the conversion graph
  ModelValuesTerminal& GetTargetNodes() { return dest_; }

  /// Add (register) a link range.
  /// This is normally called automatically from a link
  /// when an entry is added.
  /// Exports previous entries if exporter provided.
  void Add(LinkRange br) {
    if (brl_.empty() || !brl_.back().TryExtendBy(br)) {
      ExportRemainingEntries();       // new entry, export existing
      brl_.Add(br);
    }
  }

  /// Want Export?
  bool GetExport() const { return bts_.IsOpen(); }

  /// Finish exporting entries, if exporting=ON.
  /// This should be called after model conversions are finished,
  /// because the exports are lazy.
  void FinishExportingLinkEntries() { ExportRemainingEntries(); }

  /// Check that the whole entries list has been exported
  bool AllEntriesExported() const { return brl_.size()==(size_t)i_exported_; }


  /// Pre- / postsolve loops over link entries
  /// and calls the link's method for each.
#undef PRESOLVE_KIND
#define PRESOLVE_KIND(name, ValType) \
  MVOverEl<ValType> \
    Presolve ## name ( \
      const MVOverEl<ValType> & mv) override { \
    return RunPresolve(&BasicLink::Presolve ## name, mv); \
  } \
  MVOverEl<ValType> \
    Postsolve ## name ( \
      const MVOverEl<ValType> & mv) override { \
    return RunPostsolve(&BasicLink::Postsolve ## name, mv); \
  }

  LIST_PRESOLVE_METHODS


  /// Register a ValueNode*
  void Register(ValueNode* pvn) override {
    auto res = val_nodes_.insert(pvn).second;
    assert(res);
  }

  /// Deregister a ValueNode*
  void Deregister(ValueNode* pvn) override {
    auto res = val_nodes_.erase(pvn);
    assert(res);
  }


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
    for (auto rit=brl_.rbegin(); rit!=brl_.rend(); ++rit) {
      (rit->b_.*fn)(rit->ir_);
    }
    return src_;
  }


public:
  /// Clean up value nodes for new propagation.
  /// Numeric arrays only.
  void CleanUpValueNodes() const {
    for (auto pvn: val_nodes_)
      pvn->CleanUpAndRealloc();
  }
  /// Clean up name nodes for new propagation.
  void CleanUpNameNodes() const {
    for (auto pvn: val_nodes_)
      pvn->CleanUpAndRealloc_Names();
  }


protected:
  /// Export unexported entries, if exporter provided
  void ExportRemainingEntries() {
    if (GetExport()) {
      for ( ; i_exported_<(int)brl_.size(); ++i_exported_) {
        const auto& ln_rng = brl_[i_exported_];
        for (auto i_entry=ln_rng.ir_.beg_; i_entry!=ln_rng.ir_.end_;
             ++i_entry) {
          ExportLinkEntry(ln_rng.b_, i_entry);
        }
      }
    }
  }

  /// Export certain link entry
  void ExportLinkEntry(const BasicLink& ln, int i_entry) {
    if (GetExport()) {
      ln.ExportEntryItems(entry_items_, i_entry);
      fmt::MemoryWriter wrt;
      {
        MiniJSONWriter jw(wrt);
        jw["link_index"] = std::make_tuple(i_exported_, i_entry);
        jw["link_type"] = ln.GetTypeName();
        WriteNodes(jw["source_nodes"], entry_items_.src_items_);
        WriteNodes(jw["dest_nodes"], entry_items_.dest_items_);
      }
      /// Export record
      wrt.write("\n");                     // EOL
      bts_.Append(wrt);
    }
  }

  /// Write a vector of nodes
  template <class JSON>
  void WriteNodes(JSON jw, const std::vector<NodeRange>& nodes) {
    for (size_t i=0; i<nodes.size(); ++i) {
      const auto& nd = nodes[i];
      auto jelement = ++jw;       // We write an array
      if (nd.IsSingleIndex())
        jelement[nd.GetValueNode()->GetName()] = (int)nd;
      else
        jelement[nd.GetValueNode()->GetName()]
            << nd.GetIndexRange().beg_
            << (nd.GetIndexRange().end_-1);  // print the exact last value
    }
  }

private:
  /// val_nodes_ should be before src_ / dest_
  std::unordered_set<ValueNode*> val_nodes_;

  /// val_nodes_ should be before src_ / dest_
  mutable ModelValuesTerminal
    src_{*this, "src"},
    dest_{*this, "dest"};

  /// The link ranges
  LinkRangeList brl_;

  /// Exporter functor
  BasicLogger& bts_;

  /// 1-after-last exported entry
  int i_exported_=0;

  /// Link entry items temporary storage for exporting
  BasicLink::EntryItems entry_items_;
};


/// How to call a solution checker
using SolCheckerCall = bool(
    ArrayRef<double> x,
    const ValueMapDbl& y,
    ArrayRef<double> obj,
    void* p_extra);
/// Function wrapper
using SolCheckerType = std::function<SolCheckerCall>;


/// Final ValuePresolver.
/// It specializes PresolveSolution() to update fixed variables
class ValuePresolver : public ValuePresolverImpl {
public:
  ValuePresolver(BasicFlatModel& m, Env& env,
                 BasicLogger& bts, SolCheckerType sc={})
    : ValuePresolverImpl(env, bts), model_(m), solchk_(sc) { }

  /// Override PresolveSolution().
  /// Move warm start values into the bounds.
  MVOverEl<double> PresolveSolution (
      const MVOverEl<double>& mv) override {
    auto result = ValuePresolverImpl::PresolveSolution(mv);
    auto& x = result.GetVarValues()();
    const auto& lbs = model_.GetVarLBs();
    const auto& ubs = model_.GetVarUBs();
    assert(x.size() == lbs.size());
    assert(lbs.size() == ubs.size());
    for (auto i=x.size(); i--; ) {
      if (x[i]<lbs[i])
        x[i]=lbs[i];
      else
        if (x[i]>ubs[i])
          x[i]=ubs[i];
    }
    return result;
  }

  /// Override PostsolveSolution().
  /// Check solution if checker provided.
  /// mv's ExtraData() is passed to the checker.
  MVOverEl<double> PostsolveSolution (
      const MVOverEl<double>& mv) override {
    if (solchk_) {
      const auto& mx = mv.GetVarValues();
      const auto& mo = mv.GetObjValues();
      if (mx.IsSingleKey()
          && mx().size()) {          // solution available
        ArrayRef<double> objs;
        if (mo.IsSingleKey())
          objs = mo();
        solchk_(mx(), mv.GetConValues(), objs, mv.ExtraData());
      }
    }
    return ValuePresolverImpl::PostsolveSolution(mv);
  }


private:
  BasicFlatModel& model_;
  SolCheckerType solchk_;
};


/// Implement link range registration in global Presolver
inline void
BasicLink::RegisterLinkIndexRange(LinkIndexRange bir)
{ value_presolver_.Add( { *this, bir } ); }

} // namespace pre

} // namespace mp

#endif // VALUE_PRESOLVE_H
