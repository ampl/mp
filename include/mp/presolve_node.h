#ifndef PRESOLVE_NODE_H
#define PRESOLVE_NODE_H

#include "presolve_base.h"


namespace mp {

namespace pre {

/// Index range for a single node
using NodeIndexRange = IndexRange;


/// Value node, contains arrays of int's and double's
/// corresponding to variables, or a constraint type, or objectives
class ValueNode {
public:
  /// Assign from ArrayRef<int>
  ValueNode& operator=(ArrayRef<int> ai)
  { vi_ = ai.move_or_copy(); return *this; }
  /// Assign from ArrayRef<double>
  ValueNode& operator=(ArrayRef<double> ai)
  { vd_ = ai.move_or_copy(); return *this; }

  /// Retrieve whole ArrayRef<int>
  operator ArrayRef<int> () const { return vi_; }
  /// Retrieve whole ArrayRef<double>
  operator ArrayRef<double> () const { return vd_; }

  /// Retrieve int[i]
  int GetInt(size_t i) const { assert(i<vi_.size()); return vi_[i]; }
  /// Set int[i]
  void SetInt(size_t i, int v)
  { if (i>=vi_.size()) vi_.resize(i+1); vi_[i]=v; }

  /// Retrieve double[i]
  double GetDbl(size_t i) const { assert(i<vd_.size()); return vd_[i]; }
  /// Set double[i]
  void SetDbl(size_t i, double v)
  { if (i>=vd_.size()) vd_.resize(i+1); vd_[i]=v; }

  /// Copy to another node
  void Copy(NodeIndexRange ir, ValueNode& dest, int index1) {
    if (vi_.size()>=ir.end) {
      if (dest.vi_.size() < index1 + ir.size())
        dest.vi_.resize(index1 + ir.size());
      std::copy(vi_.begin()+ir.beg, vi_.begin()+ir.end,
                dest.vi_.begin()+index1);
    }
    if (vd_.size()>=ir.end) {
      if (dest.vd_.size() < index1 + ir.size())
        dest.vd_.resize(index1 + ir.size());
      std::copy(vd_.begin()+ir.beg, vd_.begin()+ir.end,
                dest.vd_.begin()+index1);
    }
  }

private:
  std::vector<int> vi_;
  std::vector<double> vd_;
};


/// Typedef map of nodes
using NodeMap = ValueMap< ValueNode >;
/// Typedef ModelValues stored in NodeMaps
using ModelValuesTerminal = ModelValues< NodeMap >;


/// Node range: range of node entries
/// The node is specified as well
struct NodeRange {
  ValueNode& vn_;
  NodeIndexRange ir_;
};


} // namespace pre

} // namespace mp


#endif // PRESOLVE_NODE_H
