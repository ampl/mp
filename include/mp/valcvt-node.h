#ifndef VALUE_PRESOLVE_NODE_H
#define VALUE_PRESOLVE_NODE_H

#include "valcvt-base.h"


namespace mp {

namespace pre {

/// Index range for a single node
using NodeIndexRange = IndexRange;

/// Declare ValueNode
class ValueNode;

/// Node range: range of node entries in a specific value node.
/// The node is specified as well
class NodeRange {
public:
  /// Get pvn
  ValueNode* GetValueNode() const { assert(pvn_); return pvn_; }

  /// Get index range
  NodeIndexRange GetIndexRange() const { assert(ir_.check()); return ir_; }

  /// check if the range represents just 1 index
  bool IsSingleIndex() const { return ir_.IfSingleIndex(); }

  /// Get the single index if range is just that
  int GetSingleIndex() const { return ir_; }

  /// Get single index if it is
  operator int() const { return ir_; }

  /// Check extendability
  bool ExtendableBy(NodeRange nr) const
  { return pvn_==nr.pvn_ && ir_.end==nr.ir_.beg; }

  /// Extend the range by another one
  void ExtendBy(NodeRange nr)
  { assert(ExtendableBy(nr)); ir_.end = nr.ir_.end; }


protected:
  /// Friend ValueNode
  friend class ValueNode;

  /// Assign members, only accessible to ValueNode
  void Assign(ValueNode* pvn, NodeIndexRange ir) { pvn_=pvn; ir_=ir; }


private:
  ValueNode* pvn_ = nullptr;
  NodeIndexRange ir_;
};


/// Value node, a node of the conversion graph.
/// Stores arrays of int's and double's
/// corresponding to variables, or a constraint type, or objectives.
/// The data is stored temporarily during a conversion run.
class ValueNode {
public:
  /// Default constructor
  ValueNode() = default;

  /// Constructor
  ValueNode(std::string nm) : name_(nm) { }

  /// Declared size (what is being used by links)
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
  /// Problematic since, e.g., Gurobi adds more variables in feasrelax
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

  /// Retrieve vec<T>&
  template <class T>
  std::vector<T>& GetValVec() { assert(false); return {}; }

  /// Retrieve vec<double>&
  template <>
  std::vector<double>& GetValVec<double>() { return vd_; }

  /// Retrieve vec<int>&
  template <>
  std::vector<int>& GetValVec<int>() { return vi_; }

  /// Retrieve whole const vector<int>&
  operator const std::vector<int>& () const { return vi_; }

  /// Retrieve whole vector<double>&
  operator const std::vector<double>& () const { return vd_; }


  /////////////////////// Access individual values ///////////////////////

  /// Retrieve T
  template <class T>
  T GetVal(size_t i) const { return {}; }

  /// Retrieve double
  template <>
  double GetVal<double>(size_t i) const { return GetDbl(i); }

  /// Retrieve int
  template <>
  int GetVal<int>(size_t i) const { return GetInt(i); }



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

  /// SetName
  void SetName(std::string s) { name_ = std::move(s); }

private:
  std::vector<int> vi_;
  std::vector<double> vd_;
  size_t sz_;
  std::string name_ = "default_value_node";
};


/// Specialize for ValueNode
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
    /// (SOS?)
    /// No value is copied into them
  }
  if ((int)dest.size() < i1 + ir.size())
    dest.resize(i1 + ir.size());
  std::copy(src.begin()+ir.beg, src.begin()+ir.end,
            dest.begin()+i1);
  return true;
}

/// Copy range of node entries to another node
template <class T>
inline
void Copy(NodeRange ir1, NodeRange ir2) {
  assert(ir1.GetIndexRange().size() == ir2.GetIndexRange().size());
  auto fd = CopyRange(ir1.GetValueNode()->GetValVec<T>(),
                      ir2.GetValueNode()->GetValVec<T>(),
                      ir1.GetIndexRange(), ir2.GetIndexRange().beg);
  assert(fd);
}


/// Typedef map of nodes
using NodeMap = ValueMap< ValueNode >;

/// Terminal nodes of a conversion graph
using ModelValuesTerminal = ModelValues< NodeMap >;


} // namespace pre

} // namespace mp


#endif // VALUE_PRESOLVE_NODE_H
