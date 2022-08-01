#ifndef VALUE_PRESOLVE_H
#define VALUE_PRESOLVE_H

/**
  @file valcvt.h.
  Implementation of value presolver
  */

#include <deque>
#include <unordered_set>
#include <functional>

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


/// Value presolver implementation.
/// A ValuePresolver converts solutions and suffixes
/// between the original model and the presolved one.
class ValuePresolver : public BasicValuePresolver {
public:
  /// Exporter functor type
  using ExporterFn = std::function< void (const char*) >;

  /// Constructor
  ValuePresolver(Env& env, ExporterFn bts={}) :
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

  /// Switch exporting on/off.
  void SetExport(bool onoff) { f_export_=onoff; }

  /// Want Export?
  bool GetExport() const { return f_export_; }

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

  /// Clean up value nodes for new propagation
  void CleanUpValueNodes() const {
    for (auto pvn: val_nodes_)
      pvn->CleanUpAndRealloc();
  }

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
      wrt.write("{} ", '{');
      wrt.write("\"link_index\": [{}, {}], ", i_exported_, i_entry);
      wrt.write("\"link_type\": {}, ", ln.GetTypeName());
      wrt.write("\"source_nodes\": [");
        WriteNodes(wrt, entry_items_.src_items_);
      wrt.write(" ], ");                      // end source nodes
      wrt.write("\"dest_nodes\": [ ");
        WriteNodes(wrt, entry_items_.dest_items_);
      wrt.write(" ]");                        // end dest nodes
      wrt.write(" {}\n", '}');                      // with EOL
      /// Export record
      assert(bts_);
      bts_(wrt.c_str());
    }
  }

  /// Write a vector of nodes
  void WriteNodes(fmt::MemoryWriter& wrt, const std::vector<NodeRange>& nodes) {
    for (size_t i=0; i<nodes.size(); ++i) {
      const auto& nd = nodes[i];
      if (nd.IsSingleIndex())
        wrt.write("{} \"{}\": {} {}",
                  '{', nd.GetValueNode()->GetName(), (int)nd, '}');
      else
        wrt.write("{} \"{}\": [{}, {}] {}",
                  '{', nd.GetValueNode()->GetName(),
                  nd.GetIndexRange().beg_,
                  nd.GetIndexRange().end_-1,     // print the exact last value
                  '}');
      if (i != nodes.size()-1)
        wrt.write(",");
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

  /// Export on/off
  bool f_export_ { false };

  /// Exporter functor
  ExporterFn const bts_{};

  /// 1-after-last exported entry
  int i_exported_=0;

  /// Link entry items temporary storage for exporting
  BasicLink::EntryItems entry_items_;
};


/// Implement link range registration in global Presolver
inline void
BasicLink::RegisterLinkIndexRange(LinkIndexRange bir)
{ value_presolver_.Add( { *this, bir } ); }

} // namespace pre

} // namespace mp

#endif // VALUE_PRESOLVE_H
