#ifndef MODEL_INFO_HPP
#define MODEL_INFO_HPP

/**
 * Implementation of flat model info
 */

#include <unordered_map>
#include <functional>

#include "mp/flat/model_info.h"

namespace mp {

/// Implementation of flat model info
class FlatModelInfoImpl : public FlatModelInfo {
public:

  /// For hashing of type_info
  using TypeInfoRef = std::reference_wrapper<const std::type_info>;

  /// TypeInfoRefHasher
  struct TypeInfoRefHasher {
      std::size_t operator()(TypeInfoRef code) const
      {
          return code.get().hash_code();
      }
  };

  /// TypeInfoRef ==
  struct TypeInfoRefEqualTo {
      bool operator()(TypeInfoRef lhs, TypeInfoRef rhs) const
      {
          return lhs.get() == rhs.get();
      }
  };

  /// Hash map of ints by TypeInfoRef
  using TypeInfoRefIntMap =
    std::unordered_map<TypeInfoRef, int,
      TypeInfoRefHasher, TypeInfoRefEqualTo>;

  /// Hash map of ints by constraint groups
  using ConstrGroupIntMap = std::unordered_map<int, int>;

  /// Get number of constraints of certain group
  int GetNumberOfConstraintsOfGroup(int cg) const override {
    if (cg_map_.end() != cg_map_.find(cg))
      return cg_map_.at(cg);
    return 0;
  }

  /// Get number of constraints of single type
  int GetNumberOfConstraints(const std::type_info& it) const override {
    if (ti_map_.end() != ti_map_.find(it))
      return ti_map_.at(it);
    return 0;
  }

  /// Initialize constraint counting
  void InitConstraintCount() override
  { cg_map_.clear(); ti_map_.clear(); }

  /// Add number of constraints of single type
  void AddNumberOfConstraints(
      const std::type_info& ti, int igroup, int nc) override {
    cg_map_[igroup] += nc;
    ti_map_[ti] += nc;
  }


private:
  TypeInfoRefIntMap ti_map_;
  ConstrGroupIntMap cg_map_;
};


/// FlatModelInfo factory
std::unique_ptr<FlatModelInfo> CreateFlatModelInfo() {
  return std::unique_ptr<FlatModelInfo>{new FlatModelInfoImpl()};
}

} // namespace mp

#endif // MODEL_INFO_HPP
