#ifndef VALUE_PRESOLVE_NODE_H
#define VALUE_PRESOLVE_NODE_H

#include <cmath>

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
  /// Constructor
  NodeRange(ValueNode* pvn=nullptr, NodeIndexRange ir={}) :
    pvn_(pvn), ir_(ir) { }

  /// Valid?
  bool IsValid() const { return pvn_ && ir_.IsValid(); }

  /// invalidate
  void Invalidate() { pvn_=nullptr; ir_.Invalidate(); }

  /// Get pvn
  ValueNode* GetValueNode() const { assert(pvn_); return pvn_; }

  /// Get index range
  NodeIndexRange GetIndexRange() const { assert(ir_.IsValid()); return ir_; }

  /// check if the range represents just 1 index
  bool IsSingleIndex() const { return ir_.IsSingleIndex(); }

  /// Get the single index if range is just that
  int GetSingleIndex() const { return ir_; }

  /// Get single index if it is
  operator int() const { return ir_; }

  /// Check extendability
  bool ExtendableBy(NodeRange nr) const
  { return pvn_==nr.pvn_ && ir_.end_==nr.ir_.beg_; }

  /// Try & extend the range by another one
  bool TryExtendBy(NodeRange nr) {
    if (!ExtendableBy(nr))
      return false;
    ExtendBy(nr);
    return true;
  }

  /// Extend the range by another one
  void ExtendBy(NodeRange nr)
  { assert(ExtendableBy(nr)); ir_.end_ = nr.ir_.end_; }


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
  /// Destructor
  ~ValueNode() { DeregisterMe(); }

  /// Constructor.
  /// Need ValuePresolver to register itself.
  ValueNode(BasicValuePresolver& pre, std::string nm={}) :
    pre_(pre), name_(nm) { RegisterMe(); }

  /// Move constructor
  ValueNode(ValueNode&& vn) : pre_(vn.pre_) {
    vi_ = std::move(vn.vi_);
    vd_ = std::move(vn.vd_);
    sz_ = std::move(vn.sz_);
    name_ = std::move(vn.name_);
    RegisterMe();
  }

  /// Copy constructor
  ValueNode(const ValueNode& vn) : pre_(vn.pre_) {
    vi_ = (vn.vi_);
    vd_ = (vn.vd_);
    sz_ = (vn.sz_);
    name_ = (vn.name_);
    RegisterMe();
  }

  /// Declared size (what is being used by links)
  size_t Size() const { return sz_; }

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
  /// Some (target) node assignments can be longer vectors:
  /// e.g., Gurobi  adds variables for FeasRelax.
  ValueNode& operator=(std::vector<int> ai)
  {
    // assert(ai.size() <= size());
    vi_ = std::move(ai);
    vi_.resize(Size());       // cut off / complement values
    return *this;
  }

  /// Assign from vector. Always copy
  ValueNode& operator=(std::vector<double> ad)
  {
    // assert(ad.size() <= size());
    vd_ = std::move(ad);
    vd_.resize(Size());       // cut off / complement values
    return *this;
  }

  /// Retrieve whole ArrayRef<int>
  operator ArrayRef<int> () const { return vi_; }

  /// Retrieve whole ArrayRef<double>
  operator ArrayRef<double> () const { return vd_; }

  /// Retrieve vec<T>& - dummy version
  template <class T>
  std::vector<T>& GetValVec()
  { assert(false); static std::vector<T> dum; return dum; }

  /// Retrieve whole const vector<int>&
  operator const std::vector<int>& () const { return vi_; }

  /// Retrieve whole vector<double>&
  operator const std::vector<double>& () const { return vd_; }


  /////////////////////// Access individual values ///////////////////////

  /// Retrieve T, dummy version
  template <class T>
  T GetVal(size_t ) const { return {}; }

  /// Set T, dummy version
  template <class T>
  void SetVal(size_t i, T ) { assert(i<Size()); }

  /// Retrieve int[i]
  int GetInt(size_t i) const { assert(i<vi_.size()); return vi_[i]; }

  /// Set int[i].
  /// If existing value non-0, only allow larger value.
  void SetInt(size_t i, int v) { SetNum(vi_, i, v); }

  /// Retrieve double[i]
  double GetDbl(size_t i) const { assert(i<vd_.size()); return vd_[i]; }

  /// Set double[i].
  /// If existing value non-0, only allow larger value.
  void SetDbl(size_t i, double v) { SetNum(vd_, i, v); }

  /// SetName
  void SetName(std::string s) { name_ = std::move(s); }

  /// Clean up and realloc with current size, fill by 0's
  void CleanUpAndRealloc() {
    vi_.clear(); vd_.clear();
    vi_.resize(Size()); vd_.resize(Size());
  }

protected:
  /// Set int[i] or dbl[i].
  /// If existing value non-0, only allow larger value,
  /// BUT only if that is non-0 too.
  /// This way, the order of propagations is irrelevant.
  template <class Vec>
  void SetNum(Vec& vec, size_t i, typename Vec::value_type v) {
    assert(i<Size());   // index into the originally declared suffix size
    if (vec.size()<=i)  // can happen after CopySrcDest / CopyDestSrc
      vec.resize(Size());
    /// As we initialize vectors (by 0's), don't need std::fpclassify().
    /// Moreover on Windows: https://github.com/microsoft/STL/issues/519
    if ( vec[i] ) {
      if (v>vec[i] && v)
        vec[i]=v;
      // TODO warning for unequal values
    } else
      vec[i]=v;
  }

  /// Register with the ValuePresolver
  void RegisterMe() {
    pre_.Register(this);
  }

  /// Deregister with the ValuePresolver
  void DeregisterMe() {
    pre_.Deregister(this);
  }

private:
  // Move & copy constructors should copy all members!
  BasicValuePresolver& pre_;
  std::vector<int> vi_;
  std::vector<double> vd_;
  size_t sz_=0;
  std::string name_ = "default_value_node";
};


template <>
std::vector<double>& ValueNode::GetValVec<double>() { return vd_; }

template <>
std::vector<int>& ValueNode::GetValVec<int>() { return vi_; }

template <>
double ValueNode::GetVal<double>(size_t i) const { return GetDbl(i); }

template <>
int ValueNode::GetVal<int>(size_t i) const { return GetInt(i); }

template <>
void ValueNode::SetVal<double>(size_t i, double v) { SetDbl(i, v); }

template <>
void ValueNode::SetVal<int>(size_t i, int v) { SetInt(i, v); }


/// Specialize SetValueNodeName() for ValueNode
inline void
SetValueNodeName(ValueNode& vn, std::string nm) { vn.SetName(nm); }


/// Default CreateArray() for ValueNode
template <>
ValueNode CreateArray(BasicValuePresolver& vp) { return ValueNode{vp}; }


/// Copy int or double range only
/// @return always true currently
template <class Vec>
inline
bool CopyRange(Vec& src, Vec& dest, NodeIndexRange ir, int i1) {
  assert(ir.end_ <= (int)src.size());
  assert(i1 + ir.Size() <= (int)dest.size());
  std::copy(src.begin()+ir.beg_, src.begin()+ir.end_,
            dest.begin()+i1);
  return true;
}

/// Copy range of node entries to another node
template <class T>
inline
void Copy(NodeRange ir1, NodeRange ir2) {
  assert(ir1.GetIndexRange().Size() == ir2.GetIndexRange().Size());
  auto fd = CopyRange(ir1.GetValueNode()->GetValVec<T>(),
                      ir2.GetValueNode()->GetValVec<T>(),
                      ir1.GetIndexRange(), ir2.GetIndexRange().beg_);
  assert(fd);
}


/// Typedef BasicValuePresolverRef
using BasicValuePresolverRef = BasicValuePresolver&;

/// Typedef map of nodes
using NodeMap = ValueMap< ValueNode, BasicValuePresolverRef >;

/// Terminal nodes of a conversion graph
using ModelValuesTerminal = ModelValues< NodeMap >;


} // namespace pre

} // namespace mp


#endif // VALUE_PRESOLVE_NODE_H
