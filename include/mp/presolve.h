#ifndef PRESOLVE_H
#define PRESOLVE_H

/**
  @file presolve.h
  Implementation of value presolver
  */

#include <deque>

#include "presolve_base.h"


namespace mp {

namespace pre {

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

private:
  std::vector<int> vi_;
  std::vector<double> vd_;
};


/// Typedef map of nodes
using NodeMap = ValueMap< ValueNode >;
/// Typedef ModelValues stored in NodeMaps
using ModelValuesTerminal = ModelValues< NodeMap >;


/// index range for a specific bridge
struct BridgeIndexRange {
  /// Construct, possibly from a single index
  BridgeIndexRange(int b, int e=-1) : beg(b), end(e<0 ? b+1 : e) { }

  int beg=0, end=0;
};

/// Bridge interface
///
/// A bridge is an array of value converters between nodes
/// As such, a bridge is an array of conversion specifiers
/// of the same type
class BasicBridge {
public:
  /// Virtual destructor
  virtual ~BasicBridge() = default;

  /// The below pre- / postsolves work on a range of bridges

  /// Presolve solution (primal + dual)
  virtual void PresolveSolution(BridgeIndexRange) = 0;
  /// Postsolve solution (primal + dual)
  virtual void PostsolveSolution(BridgeIndexRange) = 0;

  /// Presolve basis (primal + dual)
  virtual void PresolveBasis(BridgeIndexRange) = 0;
  /// Postsolve solution (primal + dual)
  virtual void PostsolveBasis(BridgeIndexRange) = 0;
};


/// Bridge range: range of conversion specifiers of certain type
/// The bridge is specified as well
struct BridgeRange {
  /// Try & extend the range
  /// @return true iff extension worked,
  /// otherwise the caller would have to add the new range
  /// separately
  bool ExtendRange(BridgeRange br) {
    if (&b_ == &br.b_) {
      if (ir_.end == br.ir_.beg) {
        ir_.end = br.ir_.end;
        return true;
      }
    }
    return false;
  }

  BasicBridge& b_;
  BridgeIndexRange ir_;
};


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

protected:
  /// using Base::back()
  using Base::back;
};

/// Value presolver implementation
class Presolver : public BasicPresolver {
public:
  /// Presolve solution (primal + dual)
  ModelValuesDbl PresolveSolution(const ModelValuesDbl& mvd) override
  { return RunPresolve(&BasicBridge::PresolveSolution, mvd); }
  /// Postsolve solution (primal + dual)
  ModelValuesDbl PostsolveSolution(const ModelValuesDbl& ) override;

  /// Presolve basis (primal + dual)
  ModelValuesInt PresolveBasis(const ModelValuesInt& ) override;
  /// Postsolve solution (primal + dual)
  ModelValuesInt PostsolveBasis(const ModelValuesInt& ) override;

protected:
  /// Helper type: virtual member function pointer
  using BridgeFn = void (BasicBridge::*)(BridgeIndexRange);
  /// Actual value presolve cycle: forward
  template <class ModelValues>
  ModelValues RunPresolve(BridgeFn fn, const ModelValues& mv) const {
    src_ = mv;
    for (const auto& br: brl_)
      (br.b_.*fn)(br.ir_);
    return dest_;
  }

private:
  mutable ModelValuesTerminal src_, dest_;
  BridgeRangeList brl_;
};


} // namespace pre

} // namespace mp

#endif // PRESOLVE_H
