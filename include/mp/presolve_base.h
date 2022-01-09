#ifndef PRESOLVE_BASE_H
#define PRESOLVE_BASE_H

/**
  @file presolve_base.h
  Interface for value presolver: transforms solutions etc.
  between original and presolved model
  */

#include <map>
#include <cassert>

#include "mp/arrayref.h"

namespace mp {

namespace pre {


/// Value map, contains a map of arrays of int's and/or double's
/// corresponding to variables, or a constraint type, or objectives
template <class Array>
class ValueMap {
public:
  using MapType = std::map< int, Array >;

  /// Default constructor
  ValueMap() = default;

  /// Construct from a single SomeArray
  template <class SomeArray>
  ValueMap(SomeArray r) { map_[0] = r; }

  /// Construct from the low-level MapType
  ValueMap(MapType m) : map_{m} { }

  /// Construct from another map
  template <class Array2>
  ValueMap(const ValueMap<Array2>& m) {
    for (const auto& a2: m.GetMap())
      map_[a2.first] = a2.second;
  }

  /// Check if we have only the single array
  bool IfSingleArray() const
  { return 1==map_.size() && 0==map_.begin()->first; }

  /// Retrieve the single array, const
  const Array& operator()() const
  { assert(IfSingleArray()); return map_.at(0); }
  /// Retrieve the single array (create if need)
  Array& operator()()
  { assert(map_.empty() || IfSingleArray()); return map_[0]; }

  /// Retrieve the array with index i, const
  const Array& operator()(int i) const { return map_.at(i); }
  /// Retrieve the array with index i (create if need)
  Array& operator()(int i)
  { return map_[i]; }

  const MapType& GetMap() const { return map_; }
private:
  MapType map_;
};


/// Group of values for variables, constraints, and objectives
template <class VMap>
class ModelValues {
public:
  /// Construct from values for vars, cons, objs
  /// (last 2 can be omitted)
  ModelValues(VMap v, VMap c = {}, VMap o = {}) :
    vars_(v), cons_(c), objs_(o) { }
  /// Construct from ModelValues<AnotherVMap>
  template <class MV2>
  ModelValues(const MV2& vm) :
    vars_(vm.GetVarValues()),
    cons_(vm.GetConValues()),
    objs_(vm.GetObjValues()) { }

  /// Retrieve vars map, const
  const VMap& GetVarValues() const { return vars_; }
  /// Retrieve vars map
  VMap& GetVarValues() { return vars_; }

  /// Retrieve cons map, const
  const VMap& GetConValues() const { return cons_; }
  /// Retrieve const map
  VMap& GetConValues() { return cons_; }

  /// Retrieve objs map, const
  const VMap& GetObjValues() const { return objs_; }
  /// Retrieve objs map
  VMap& GetObjValues() { return objs_; }

private:
  VMap vars_, cons_, objs_;
};


/// Specialize ModelValues<> for int
using ModelValuesInt = ModelValues< ValueMap< mp::ArrayRef<int> > >;
/// Specialize ModelValues<> for double
using ModelValuesDbl = ModelValues< ValueMap< mp::ArrayRef<double> > >;


/// Presolver interface
class BasicPresolver {
public:
  /// Virtual destructor
  virtual ~BasicPresolver() = default;

  /// Presolve solution (primal + dual)
  virtual ModelValuesDbl PresolveSolution(const ModelValuesDbl& ) = 0;
  /// Postsolve solution (primal + dual)
  virtual ModelValuesDbl PostsolveSolution(const ModelValuesDbl& ) = 0;

  /// Presolve basis (primal + dual)
  virtual ModelValuesInt PresolveBasis(const ModelValuesInt& ) = 0;
  /// Postsolve solution (primal + dual)
  virtual ModelValuesInt PostsolveBasis(const ModelValuesInt& ) = 0;
};

} // namespace pre

} // namespace mp

#endif // PRESOLVE_BASE_H
