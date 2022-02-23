#ifndef PRESOLVE_BASE_H
#define PRESOLVE_BASE_H

/**
  Interface for value presolver: transforms solutions etc.
  between original and presolved model
  */

#include <map>
#include <string>
#include <cassert>

#include "mp/arrayref.h"

namespace mp {

namespace pre {

/// Debugging template
template <class Any>
void SetValueNodeName(Any& , std::string ) { }


/// Value map, contains a map of arrays of int's and/or double's
/// corresponding to variables, or a constraint type, or objectives
template <class Array>
class ValueMap {
public:
  /// Typedef to map item type number to value array
  /// For example, distinguish values for linear constraints etc
  using MapType = std::map< int, Array >;

  /// Default constructor
  ValueMap() = default;

  /// Construct from a single SomeArray
  template <class SomeArray>
  ValueMap(SomeArray r) {
    map_[0] = std::move(r);
    SetValueNodeName(map_[0], name_ + std::to_string(0));
  }

  /// Construct from the low-level MapType
  ValueMap(MapType m) : map_{std::move(m)} { }

  /// Construct from a low-level map<AnotherArray>
  template <class Array2>
  ValueMap(const std::map<int, Array2>& llm) {
    for (const auto& a2: llm)
      map_[a2.first] = a2.second;
  }

  /// Construct from another map
  template <class Array2>
  ValueMap(const ValueMap<Array2>& m) {
    for (const auto& a2: m.GetMap())
      map_[a2.first] = a2.second;
  }

  /// Assign from another map
  template <class Array2>
  ValueMap& operator=(const ValueMap<Array2>& m) {
    for (const auto& a2: m.GetMap())
      map_[a2.first] = a2.second;
    return *this;
  }

  /// Check if we have only the single key
  bool IfSingleKey() const
  { return 1==map_.size() && 0==map_.begin()->first; }

  /// Make single key, or check one, and return
  Array& MakeSingleKey() { return operator()(); }

  /// Retrieve the single array, const
  const Array& operator()() const
  { assert(IfSingleKey()); return map_.at(0); }
  /// Retrieve the single array (create if need)
  Array& operator()() {
    if (map_.empty())
      SetValueNodeName(map_[0], name_ + std::to_string(0));
    else
      assert(IfSingleKey());
    return map_[0];
  }

  /// Retrieve the array with index i, const
  const Array& operator()(int i) const { return map_.at(i); }
  /// Retrieve the array with index i (create if need)
  Array& operator()(int i) {
    if (map_.end() == map_.find(i))
      SetValueNodeName(map_[i], name_ + std::to_string(i));
    return map_[i];
  }

  const MapType& GetMap() const { return map_; }

  void SetName(std::string s) { name_ = std::move(s); }

private:
  std::string name_ { "default_value_map" };
  MapType map_;
};


/// Convenience typedef
using ValueMapInt = ValueMap< std::vector<int> >;
/// Convenience typedef
using ValueMapDbl = ValueMap< std::vector<double> >;

/// Group of values for variables, constraints, and objectives
template <class VMap>
class ModelValues {
public:
  /// Constructor
  ModelValues(std::string nm) : name_{nm} {
    vars_.SetName(nm+"_vars");
    cons_.SetName(nm+"_cons");
    objs_.SetName(nm+"_objs");
  }
  /// Default move-copy
  ModelValues(ModelValues&& ) = default;
  /// Default move-assign
  ModelValues& operator=(ModelValues&& ) = default;
  /// Construct from values for vars, cons, objs
  /// (last 2 can be omitted)
  ModelValues(VMap v, VMap c = {}, VMap o = {}) :
    vars_(std::move(v)), cons_(std::move(c)), objs_(std::move(o)) { }

  /// Construct from ModelValues<AnotherVMap>
  template <class VM2>
  ModelValues(const ModelValues<VM2>& vm) :
    vars_(vm.GetVarValues()),
    cons_(vm.GetConValues()),
    objs_(vm.GetObjValues()) { }

  /// Assign from ModelValues<AnotherVMap>
  template <class VM2>
  ModelValues& operator=(const ModelValues<VM2>& vm) {
    vars_ = vm.GetVarValues();
    cons_ = vm.GetConValues();
    objs_ = vm.GetObjValues();
    return *this;
  }

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
  std::string name_;
  VMap vars_, cons_, objs_;
};


/// Specialize ModelValues<> for int
using ModelValuesInt = ModelValues< ValueMapInt >;
/// Specialize ModelValues<> for double
using ModelValuesDbl = ModelValues< ValueMapDbl >;


/// Presolver interface
/// Currently only addresses value pre- / postsolve (solutions, basis, ...)
class BasicPresolver {
public:
  /// Virtual destructor
  virtual ~BasicPresolver() = default;

  /// Presolve solution (primal + dual)
  virtual ModelValuesDbl PresolveSolution(const ModelValuesDbl& ) = 0;
  /// Postsolve solution (primal + dual)
  virtual ModelValuesDbl PostsolveSolution(const ModelValuesDbl& ) = 0;

  /// Presolve basis (vars + cons)
  virtual ModelValuesInt PresolveBasis(const ModelValuesInt& ) = 0;
  /// Postsolve solution (vars + cons)
  virtual ModelValuesInt PostsolveBasis(const ModelValuesInt& ) = 0;

  /// Presolve IIS (vars + cons)
  virtual ModelValuesInt PresolveIIS(const ModelValuesInt& ) = 0;
  /// Postsolve IIS (vars + cons)
  virtual ModelValuesInt PostsolveIIS(const ModelValuesInt& ) = 0;

  /// Presolve LazyUserCutFlags
  virtual ModelValuesInt PresolveLazyUserCutFlags(const ModelValuesInt& ) = 0;
};


/// index range for some bridge or node
struct IndexRange {
  /// Construct, possibly from a single index
  IndexRange(int b=0, int e=-1) : beg(b), end(e<0 ? b+1 : e) { }

  /// Validate
  bool check() const { return end>beg; }
  /// Size()
  int size() const { assert(check()); return end-beg; }
  /// Represents just 1 index
  bool IfSingleIndex() const { return beg==end-1; }
  /// Return single index if it is
  operator int() const { assert(IfSingleIndex()); return beg; }

  int beg=0, end=0;
};


} // namespace pre

} // namespace mp

#endif // PRESOLVE_BASE_H
