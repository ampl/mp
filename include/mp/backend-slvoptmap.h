#ifndef BACKEND_SLVOPTMAP_H
#define BACKEND_SLVOPTMAP_H

#include <unordered_map>
#include <memory>
#include <typeinfo>
#include <typeindex>
#include <cassert>

#include "mp/error.h"

namespace mp {

/// Abstract map of solver option codes
class BasicSolverOptMap {
public:
  /// Destroy
  virtual ~BasicSolverOptMap() { }
  /// find option
  virtual const char* FindOption(const void* pkey) const = 0;
  /// add option reference
  virtual void AddOption(const void* pkey, const char* amplname) = 0;
};

/// Implementation
template <class Key>
class SolverOptMap : public BasicSolverOptMap {
  std::unordered_map<Key, const char*> opts_;
public:
  /// find option
  const char* FindOption(const void* pkey) const override {
    const auto& key = *static_cast<const Key*>(pkey);
    auto it = opts_.find(key);
    if (opts_.end() != it)
      return it->second;
    return nullptr;
  }
  /// add option
  void AddOption(const void* pkey, const char* name) override {
    const auto& key = *static_cast<const Key*>(pkey);
    auto it = opts_.find(key);
    MP_ASSERT_ALWAYS(opts_.end() == it,
                     "adding repeated solver option key, old name list:\n    "
                         + std::string(it->second) + "\nNew name list:\n    "
                         + std::string(name));
    opts_[key] = name;
  }
};

/// Solver option map manager
class SolverOptionMapManager {
std::unordered_map<std::type_index,
                     std::unique_ptr<BasicSolverOptMap> > opt_maps_;
public:
  /// Find option
  template <class Key, class Value>
  const char* FindOption(const Key& key) {
    auto& map = opt_maps_[std::type_index( typeid(Value) )];
    if (!map)
      return nullptr;
    return map->FindOption(&key);
  }
  /// Add option reference
  template <class Key, class Value>
  void AddOptionRef(const Key& key, const char* name) {
    auto& map = opt_maps_[std::type_index( typeid(Value) )];
    if (!map)
      map = std::make_unique< SolverOptMap<Key> >();
    map->AddOption(&key, name);
  }
};

}  // namespace mp

#endif // BACKEND_SLVOPTMAP_H
