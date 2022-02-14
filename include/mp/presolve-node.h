#ifndef PRESOLVE_NODE_H
#define PRESOLVE_NODE_H

#include "presolve-base.h"


namespace mp {

namespace pre {

/// Index range for a single node
using NodeIndexRange = IndexRange;

/// Declare ValueNode
class ValueNode;

/// Node range: range of node entries
/// The node is specified as well
class NodeRange {
public:
  /// Construct
  NodeRange() = default;
  /// Construct
  NodeRange(const NodeRange& nr) noexcept : pvn_(nr.pvn_), ir_(nr.ir_) { }

  /// Assign
  NodeRange& operator=(const NodeRange& nr)
  { pvn_=nr.pvn_; ir_=nr.ir_; return *this; }
  /// Get pvn
  ValueNode* GetValueNode() const { assert(pvn_); return pvn_; }
  /// Get index range
  NodeIndexRange GetIndexRange() const { assert(ir_.check()); return ir_; }
  /// Get single index if it is
  operator int() const { return ir_; }

  /// Check extendability
  bool ExtendableBy(NodeRange nr) const
  { return pvn_==nr.pvn_ && ir_.end==nr.ir_.beg; }
  /// Extend
  void ExtendBy(NodeRange nr)
  { assert(ExtendableBy(nr)); ir_.end = nr.ir_.end; }

protected:
  /// Declare ValueNode
  friend class ValueNode;
  /// Assign members, only accessible to ValueNode
  void Assign(ValueNode* pvn, NodeIndexRange ir) { pvn_=pvn; ir_=ir; }

private:
  ValueNode* pvn_ = nullptr;
  NodeIndexRange ir_;
};


/// Value node, contains arrays of int's and double's
/// corresponding to variables, or a constraint type, or objectives
class ValueNode {
public:
  /// Default constructor
  ValueNode() = default;
  /// Constructor
  ValueNode(std::string nm) : name_(nm) { }
  /// Declared size (what is being used by bridges)
  size_t size() const { return sz_; }
  /// Create entry (range) pointer: add n elements
  NodeRange Add(int n=1) {
    NodeRange nr;
    nr.Assign(this, {(int)sz_, (int)sz_+n});
    sz_ += n;
    return nr;
  }
  /// Create entry (range) pointer: select n elements at certain pos
  NodeRange Select(int pos, int n=1) {
    NodeRange nr;
    nr.Assign(this, {pos, pos+n});
    if ((int)sz_<pos+n)
      sz_ = pos+n;
    return nr;
  }

  /////////////////////// Access value vectors ///////////////////////

  /// Assign from ArrayRef<int>. Always copy the values.
  /// TODO Check that values fit the declared size?
  /// Problematic since e.g. Gurobi adds more variables in feasrelax
  ValueNode& operator=(std::vector<int> ai)
  {
    // assert(ai.size() <= size());
    vi_ = std::move(ai);
    return *this;
  }
  /// Assign from vector. Always copy
  ValueNode& operator=(std::vector<double> ad)
  {
    // assert(ad.size() <= size());
    vd_ = std::move(ad);
    return *this;
  }

  /// Retrieve whole ArrayRef<int>
  operator ArrayRef<int> () const { return vi_; }
  /// Retrieve whole ArrayRef<double>
  operator ArrayRef<double> () const { return vd_; }

  /// Retrieve whole vector<int>&
  operator const std::vector<int>& () const { return vi_; }
  /// Retrieve whole vector<double>&
  operator const std::vector<double>& () const { return vd_; }

  /////////////////////// Access individual values ///////////////////////

  /// Retrieve int[i]
  int GetInt(size_t i) const { assert(i<vi_.size()); return vi_[i]; }
  /// Set int[i]
  void SetInt(size_t i, int v) {
    assert(i<size());
    if (vi_.size()<=i)  // can happen after CopySrcDest / CopyDestSrc
      vi_.resize(size());
    vi_[i]=v;
  }

  /// Retrieve double[i]
  double GetDbl(size_t i) const { assert(i<vd_.size()); return vd_[i]; }
  /// Set double[i]
  void SetDbl(size_t i, double v) {
    assert(i<size());
    if (vd_.size()<=i)  // can happen after CopySrcDest / CopyDestSrc
      vd_.resize(size());
    vd_[i]=v;
  }

  /// Copy node to another node
  friend void Copy(NodeRange ir1, NodeRange ir2);

  /// SetName
  void SetName(std::string s) { name_ = std::move(s); }

private:
  std::vector<int> vi_;
  std::vector<double> vd_;
  size_t sz_;
  std::string name_ = "default_value_node";
};

/// Spec for ValueNode
inline void
SetValueNodeName(ValueNode& vn, std::string nm) { vn.SetName(nm); }


/// Copy int or double range only
/// @return always true currently
template <class Vec> inline
bool CopyRange(Vec& src, Vec& dest, NodeIndexRange ir, int i1) {
  if ((int)src.size() < ir.end) {
    src.resize(ir.end);
    // assert(src.empty());      // attempt to process a smaller vector
    /// We cannot do this as long there exist "rootless" constraints
    /// No value is copied into them
  }
  if ((int)dest.size() < i1 + ir.size())
    dest.resize(i1 + ir.size());
  std::copy(src.begin()+ir.beg, src.begin()+ir.end,
            dest.begin()+i1);
  return true;
}

/// Copy node to another node
inline
void Copy(NodeRange ir1, NodeRange ir2) {
  assert(ir1.GetIndexRange().size() == ir2.GetIndexRange().size());
  auto fi = CopyRange(ir1.GetValueNode()->vi_, ir2.GetValueNode()->vi_,
                      ir1.GetIndexRange(), ir2.GetIndexRange().beg);
  auto fd = CopyRange(ir1.GetValueNode()->vd_, ir2.GetValueNode()->vd_,
                      ir1.GetIndexRange(), ir2.GetIndexRange().beg);
  assert(fi || fd);
}


/// Typedef map of nodes
using NodeMap = ValueMap< ValueNode >;
/// Typedef ModelValues stored in NodeMaps
using ModelValuesTerminal = ModelValues< NodeMap >;


} // namespace pre

} // namespace mp


#endif // PRESOLVE_NODE_H
