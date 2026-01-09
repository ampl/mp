#ifndef VALUE_PRESOLVE_BASE_H
#define VALUE_PRESOLVE_BASE_H

/**
  Value presolver: namespace mp::pre
  */

#include <map>
#include <string>
#include <algorithm>
#include <cassert>

#include "mp/arrayref.h"
#include "mp/env.h"

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
/// @param Array: the type of array stored (std::vector(int) / (double)),
///         or ValueNode.
/// @param Param: type of parameter to call CreateArray &lt; Array &gt; (Param)
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
  ValueMap(const SomeArray& r) {
    Array a(r.begin(), r.end());
    map_[0] = std::move(a);
    SetValueNodeName(map_.at(0), name_ + "()");
  }

  /// Construct from the low-level MapType
  ValueMap(MapType m) : map_{std::move(m)} { }

  /// Construct from a low-level map<AnotherArray>
  template <class Array2>
  ValueMap(const std::map<int, Array2>& llm) {
    for (const auto& a2: llm) {
      map_.insert(
            { a2.first, CreateArray<Array>(prm_) }
            ).first->second =
          a2.second;
    }
  }

  /// Construct from another map
  template <class Array2, class Param2>
  ValueMap(const ValueMap<Array2, Param2>& m) {
    for (const auto& a2: m.GetMap())
      map_.insert(
            { a2.first, CreateArray<Array>(prm_) }
            ).first->second =
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

  /// Empty() if no values or all values empty()
  bool Empty() const {
    return map_.end() == std::find_if_not(
          map_.begin(), map_.end(),
          [](const typename MapType::value_type& mv) -> bool
    { return mv.second.empty(); });
  }

  /// Check if we have only the single key
  bool IsSingleKey() const
  { return 1==map_.size() && 0==map_.begin()->first; }

  /// Make single key, or check one, and return
  Array& MakeSingleKey() { return operator()(); }

  /// Retrieve the single array, const.
  const Array& operator()() const&
  { assert(IsSingleKey()); return map_.at(0); }

  /// Move out the single array, rvalue.
  Array operator()() &&  { return MoveOut(); }

  /// Retrieve the single array (creating if need)
  Array& operator()() & {
    if (map_.empty())
      SetValueNodeName(
            map_.insert({ 0,
                           CreateArray<Array, Param>(prm_) }).
            first->second,
            name_ + "()");
    else
      assert(IsSingleKey());
    return map_.at(0);
  }

  /// Retrieve the array with key \a i, const.
  const Array& operator()(int i) const& { return map_.at(i); }

  /// Move out the array with key \a i, rvalue.
  Array operator()(int i) && { return MoveOut(i); }

  /// Retrieve the array with key \a i (creating if need)
  Array& operator()(int i) & {
    if (map_.end() == map_.find(i)) {
      Array arr = CreateArray<Array, Param>(prm_);
      SetValueNodeName(
            map_.insert({ i, std::move(arr) }).first->second,
            name_ + '(' + std::to_string(i) + ')');
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
  std::string name_ { "VMapName__unset" };
  MapType map_;
};


/// ValueMap typedef template over the element type
template <class El>
using VMapOverElement = ValueMap< std::vector<El> >;


/// Specialize ValueMap storing int's
using ValueMapInt = VMapOverElement< int >;
/// Specialize ValueMap storing double's
using ValueMapDbl = VMapOverElement< double >;
/// Specialize ValueMap storing strings
using ValueMapStr = VMapOverElement< std::string >;


/// Group of values or value nodes
/// for variables, constraints, and objectives
template <class VMap>
class ModelValues {
public:
  /// The extra parameter type, taken from VMap
  using ParamType = typename VMap::ParamType;

  /// Constructor
  ModelValues(ParamType prm, const std::string& nm) :
      name_{nm}, vars_(prm), cons_(prm), objs_(prm) {
    vars_.SetName(nm+"_vars");
    cons_.SetName(nm+"_cons");
    objs_.SetName(nm+"_objs");
  }

  /// Default copy construct
  ModelValues(const ModelValues& ) = default;

  /// Default move-copy
  ModelValues(ModelValues&& ) = default;

  /// Default move-assign
  ModelValues& operator=(ModelValues&& ) = default;

  /// Construct from values for vars, cons, objs
  /// (last 2 maps can be even legally omitted)
  ModelValues(VMap v={}, VMap c = {}, VMap o = {},
              void* p_extra=nullptr) :
    vars_(std::move(v)), cons_(std::move(c)), objs_(std::move(o)),
      p_extra_(p_extra) { }

  /// Construct from ModelValues<AnotherVMap>
  template <class VM2>
  ModelValues(const ModelValues<VM2>& vm) :
      name_(vm.GetName()),
      vars_(vm.GetVarValues()),
      cons_(vm.GetConValues()),
      objs_(vm.GetObjValues()),
      p_extra_(vm.ExtraData())  { }

  /// Assign from ModelValues<AnotherVMap>
  template <class VM2>
  ModelValues& operator=(const ModelValues<VM2>& vm) {
    name_ = vm.GetName();
    vars_ = vm.GetVarValues();
    cons_ = vm.GetConValues();
    objs_ = vm.GetObjValues();
    p_extra_ = vm.ExtraData();
    return *this;
  }

  /// operator bool
  operator bool() const { return !Empty(); }

  /// Empty(). True when all VMaps are
  bool Empty() const
  { return vars_.Empty() && cons_.Empty() && objs_.Empty(); }

  /// Get name
  const std::string& GetName() const { return name_; }

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

  /// Extra data
  void* ExtraData() const { return p_extra_; }

private:
  std::string name_;
  VMap vars_, cons_, objs_;
  void* p_extra_ {nullptr};
};


/// ModelValues typedef template over vectors
/// of given element type
template <class El>
using MVOverEl = ModelValues< VMapOverElement<El> >;


/// Specialize ModelValues<> for concrete int data
using ModelValuesInt = MVOverEl< int >;
/// Specialize ModelValues<> for concrete double data
using ModelValuesDbl = MVOverEl< double >;
/// Specialize ModelValues<> for concrete string data
using ModelValuesStr = MVOverEl< std::string >;


/// Declare ValueNode
class ValueNode;


/// Typedef VCString.
/// This class stores node names
/// and counts references to a given name
/// and allows numbering derived items
/// 'Con01_1_', ...
class VCString {
public:
  /// Construct default
  VCString() = default;
  /// Construct from string
  VCString(std::string s) : s_(std::move(s)) { }
  /// Copy from VCString.
  /// Updates its counter,
  /// but inits own counter to 0.
  VCString(const VCString& vcs)
    : s_(vcs.MakeCountedName()) { }
  /// Assign VCString.
  /// Only if empty (then copy.)
  VCString& operator=(const VCString& vcs) {
    if (empty()) {
      s_ = vcs.MakeCountedName();
      assert(!n_);
    }
    return *this;
  }

  /// Return current name, unmodified.
  const std::string& MakeCurrentName() const
  { return s_; }

  /// Produce name s + '_<counter>_' if n_>0.
  /// Post-increment counter.
  std::string MakeCountedName() const {
    return
        n_++==0
        ? s_
        : s_ + '_' + std::to_string(n_) + '_';
  }
  /// Use MakeCountedName().
  operator std::string() const
  { return MakeCountedName(); }
  /// Produce another VCString with
  /// an added string or char.
  /// Does not change own counter,
  /// we assume the caller cares for name uniqueness.
  template <class T>
  VCString operator+(const T& arg) const
  { return {s_ + arg}; }

  /// empty?
  bool empty() const { return s_.empty(); }

private:
  std::string s_;
  mutable std::size_t n_ = 0;
};

/// Macro for a list of pre- / postsolve method definitions
/// in a ValuePresolver or a link.
/// Requires PRESOLVE_KIND defined to declare / define
/// corr. pre- and postsolve methods.
/// Generic(Dbl/Int): maximizes suffix value among non-0.
/// Example 1 (presolve): expression exp(y) is used in various places
/// and marked with different values of .funcpieces.
/// The largest is chosen.
/// Example 2 (postsolve): IIS membership value for a converted
/// high-level constraint can be reported as the maximum of those
/// for its low-level representation.
#define LIST_PRESOLVE_METHODS \
  PRESOLVE_KIND(GenericDbl, double) \
  PRESOLVE_KIND(GenericInt, int) \
  PRESOLVE_KIND(Solution, double) \
  PRESOLVE_KIND(Basis, int) \
  PRESOLVE_KIND(IIS, int) \
  PRESOLVE_KIND(LazyUserCutFlags, int) \
  PRESOLVE_KIND(Names, VCString)
// etc...


/// How to resolve the conflict
/// when several values are gathered in one
/// source or dest node
enum class ValueResolution {
  ValResAll,       // All values transferred
                   // and SetVal() decides what to do
                   // (can max out non-0 values.)
  ValResDefault = ValResAll,
  ValResLast       // just the last conection is transferred
};

/// ValuePresolver interface.
/// Addresses value pre- / postsolve (solutions, basis, etc)
class BasicValuePresolver : public EnvKeeper {
public:
  /// Virtual destructor
  virtual ~BasicValuePresolver() = default;

  /// Constructor
  BasicValuePresolver(Env& env) : EnvKeeper(env) { }

  /// Pre- / postsolve method declarations
#undef PRESOLVE_KIND
#define PRESOLVE_KIND(name, ValType) \
  virtual MVOverEl<ValType> \
    Presolve ## name ( \
      const MVOverEl<ValType> & mv, \
      ValueResolution = ValueResolution::ValResDefault) = 0; \
  virtual MVOverEl<ValType> \
    Postsolve ## name ( \
      const MVOverEl<ValType> & mv, \
      ValueResolution = ValueResolution::ValResDefault) = 0;

  LIST_PRESOLVE_METHODS


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

  /// operator==
  bool operator==(const IndexRange& ir) const {
    return beg_==ir.beg_ && end_==ir.end_;
  }

  /// operator!=
  bool operator!=(const IndexRange& ir) const
  { return !(*this==ir); }

  /// Represents just 1 index
  bool IsSingleIndex() const { return beg_==end_-1; }

  /// Return single index if it is
  int GetSingleIndex() const { assert(IsSingleIndex()); return beg_; }

  /// explicit operator int(): get single index
  explicit operator int() const { return GetSingleIndex(); }

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

public:
  const pre::BasicValuePresolver& GetValuePresolver() const
  { assert(pPresolver_); return *pPresolver_; }
  pre::BasicValuePresolver& GetValuePresolver()
  { assert(pPresolver_); return *pPresolver_; }


private:
  pre::BasicValuePresolver* pPresolver_ = nullptr;
};

} // namespace mp

#endif // VALUE_PRESOLVE_BASE_H
