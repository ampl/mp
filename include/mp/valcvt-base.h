#ifndef VALUE_PRESOLVE_BASE_H
#define VALUE_PRESOLVE_BASE_H

/**
  Value presolver: namespace mp::pre
  */

#include <map>
#include <string>
#include <cassert>

#include "mp/arrayref.h"

namespace mp {


/// Namespace mp::pre: value presolve methods.
/// They transform solutions, etc.,
/// between the original and presolved models
namespace pre {

/// Default SetValueNodeName().
/// For std::vector<> it does nothing.
template <class Any>
void SetValueNodeName(Any& , std::string ) { }

/// Default CreateArray().
/// For std::vector<> it does not use the Param.
template <class Array, class Param>
Array CreateArray(Param ) { return Array{}; }


/// Value(Node) map.
///
/// Contains a map of either:
/// - concrete arrays of int's and/or double's (ValueMapInt/Dbl), or
/// - ValueNode's (conversion graph nodes) which manage such arrays.
///
/// The data in a single map corresponds to either variables, constraints,
/// or objectives,
/// and map keys correspond to different item subcategories (e.g., linear vs
/// quadratic constraints).
/// When IsSingleKey(), the map only has a single array (at key 0),
/// accessible via ().
///
/// @param Array: the type of array stored (std::vector<int> / <double>),
///         or ValueNode.
/// @param Param: type of parameter to call CreateArray<Array>(Param)
///         when we need to create a mapped value.
template <class Array, class Param=int>
class ValueMap {
public:
  /// Typedef ParamType
  using ParamType =Param;

  /// Typedef to map item type number to value array
  /// For example, distinguish values for linear constraints etc
  using MapType = std::map< int, Array >;

  /// Constructor
  ValueMap(Param p={}) : prm_(p) { }

  /// Construct from a single SomeArray
  template <class SomeArray>
  ValueMap(SomeArray r) :
      map_{ {0, std::move(r) } } {
    SetValueNodeName(map_.at(0), name_ + "[0]");
  }

  /// Construct from the low-level MapType
  ValueMap(MapType m) : map_{std::move(m)} { }

  /// Construct from a low-level map<AnotherArray>
  template <class Array2>
  ValueMap(const std::map<int, Array2>& llm) {
    for (const auto& a2: llm) {
      map_.insert({ a2.first, CreateArray<Array>(prm_) }).first->second =
          a2.second;
    }
  }

  /// Construct from another map
  template <class Array2, class Param2>
  ValueMap(const ValueMap<Array2, Param2>& m) {
    for (const auto& a2: m.GetMap())
      map_.insert({ a2.first, CreateArray<Array>(prm_) }).first->second =
          a2.second;
  }

  /// Assign from another map
  template <class Array2>
  ValueMap& operator=(const ValueMap<Array2>& m) {
    for (const auto& a2: m.GetMap())
      map_.insert({
                    a2.first, CreateArray<Array, Param>(prm_)
                  }).first->second =
          a2.second;
    return *this;
  }

  /// Check if we have only the single key
  bool IsSingleKey() const
  { return 1==map_.size() && 0==map_.begin()->first; }

  /// Make single key, or check one, and return
  Array& MakeSingleKey() { return operator()(); }

  /// Retrieve the single array, const.
  /// WARNING:
  /// When returning the single array from a temporary map,
  /// use MoveOut()
  const Array& operator()() const
  { assert(IsSingleKey()); return map_.at(0); }

  /// Retrieve the single array (creating if need)
  Array& operator()() {
    if (map_.empty())
      SetValueNodeName(
            *map_.insert({ 0,
                           CreateArray<Array, Param>(prm_) }).first,
            name_ + "[0]");
    else
      assert(IsSingleKey());
    return map_.at(0);
  }

  /// Retrieve the array with key \a i, const.
  /// WARNING:
  /// When returning array from a temporary map, use MoveOut(i)
  const Array& operator()(int i) const { return map_.at(i); }

  /// Retrieve the array with key \a i (creating if need)
  Array& operator()(int i) {
    if (map_.end() == map_.find(i)) {
      Array arr = CreateArray<Array, Param>(prm_);
      SetValueNodeName(
            *map_.insert({ i, std::move(arr) }).first,
            name_ + std::to_string(i));
    }
    return map_.at(i);
  }

  /// Move out the single array
  Array MoveOut()
  { assert(IsSingleKey()); return std::move(map_.at(0)); }

  /// Move out array \a i
  Array MoveOut(int i) { return std::move(map_.at(i)); }

  /// Retrieve the whole map
  const MapType& GetMap() const { return map_; }

  /// Set map name
  void SetName(std::string s) { name_ = std::move(s); }


private:
  Param prm_;
  std::string name_ { "default_value_map" };
  MapType map_;
};


/// Convenience typedef
using ValueMapInt = ValueMap< std::vector<int> >;
/// Convenience typedef
using ValueMapDbl = ValueMap< std::vector<double> >;


/// Group of values or value nodes
/// for variables, constraints, and objectives
template <class VMap>
class ModelValues {
public:
  /// The extra parameter type, taken from VMap
  using ParamType = typename VMap::ParamType;

  /// Constructor
  ModelValues(ParamType prm, std::string nm) :
      name_{nm}, vars_(prm), cons_(prm), objs_(prm) {
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


/// Specialize ModelValues<> for concrete int data
using ModelValuesInt = ModelValues< ValueMapInt >;
/// Specialize ModelValues<> for concrete double data
using ModelValuesDbl = ModelValues< ValueMapDbl >;


class ValueNode;

/// ValuePresolver interface.
/// Addresses value pre- / postsolve (solutions, basis, etc)
class BasicValuePresolver {
public:
  /// Virtual destructor
  virtual ~BasicValuePresolver() = default;

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

  /// Register a ValueNode*
  virtual void Register(ValueNode* ) = 0;

  /// Deregister a ValueNode*
  virtual void Deregister(ValueNode* ) = 0;
};


/// index range for some link or node
struct IndexRange {
  /// Construct, possibly from a single index
  IndexRange(int b=0, int e=-1) : beg_(b), end_(e<0 ? b+1 : e) { }

  /// Validate
  bool IsValid() const { return end_>beg_; }

  /// Invalidate
  void Invalidate() { beg_=end_=0; }

  /// Size()
  int Size() const { assert(IsValid()); return end_-beg_; }

  /// Represents just 1 index
  bool IsSingleIndex() const { return beg_==end_-1; }

  /// Return single index if it is
  operator int() const { assert(IsSingleIndex()); return beg_; }

  /// begin index
  int beg_=0;

  /// end (= last+1) index
  int end_=0;
};

} // namespace pre


/// Some backends need to keep reference to the value presolver
/// for pre- / postsolving of suffix values etc.
class BasicValuePresolverKeeper {
protected:
  void SetValuePresolver(pre::BasicValuePresolver* pPre) {
    pPresolver_ = pPre;
  }

  const pre::BasicValuePresolver& GetValuePresolver() const
  { assert(pPresolver_); return *pPresolver_; }
  pre::BasicValuePresolver& GetValuePresolver()
  { assert(pPresolver_); return *pPresolver_; }


private:
  pre::BasicValuePresolver* pPresolver_ = nullptr;
};

} // namespace mp

#endif // VALUE_PRESOLVE_BASE_H
