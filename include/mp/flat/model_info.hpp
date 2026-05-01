#ifndef MODEL_INFO_HPP
#define MODEL_INFO_HPP

/**
 * Implementation of flat model info
 */

#include <unordered_map>
#include <map>
#include <functional>

#include "mp/flat/model_info.h"
#include "mp/flat/converter_info.h"

namespace mp {

/// Implementation of flat model info
class FlatModelInfoImpl : public FlatModelInfo {
public:
  FlatModelInfoImpl(const ConverterInfo* pci) : cvti_(*pci) { }

  /// Num unfixed int vars
  int NumUnfixedIntVars() const override { return nUnfxIntVars_; }

  /// Set N unfixed int vars
  void SetNumUnfixedIntVars(int n) override { nUnfxIntVars_ = n; }

  /// Get var info
  VarInfo GetVarInfo() const override { return var_info_; }

  /// Set var info
  void SetVarInfo(VarInfo vi) override { var_info_ = vi; }

  /// Get obj info
  ObjInfo GetObjInfo() const override { return obj_info_; }

  /// Set obj info
  void SetObjInfo(ObjInfo oi) override { obj_info_ = oi; }

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

  /// Obtain constraint types
  const ConstrTypeMapByName& GetConstraintTypes() const override
  { return coninfo_map_; }

  /// Initialize constraint counting
  void InitConstraintCount() override
  { cg_map_.clear(); ti_map_.clear(); coninfo_map_.clear(); }

  /// Add number of constraints of single type
  void AddNumberOfConstraints(
      const std::type_info& ti, const char* name,
      int igroup, bool is_logical, int nc, int nu, int na) override {
    cg_map_[igroup] += nc;
    ti_map_[ti] += nc;
    auto& ci = coninfo_map_[name];
    ci.name_ = name;
    ci.is_logical_ = is_logical;
    ci.n_ = nc;
    ci.n_used_ = nu;
    ci.n_total_ = na;
  }

  /// Value of option cvt:expr:refcountmax
  int RefCountMaxAlgebraic() const override
  { return cvti_.RefCountMaxAlgebraic(); }


private:
  const ConverterInfo& cvti_;

  TypeInfoRefIntMap ti_map_;
  ConstrGroupIntMap cg_map_;
  ConstrTypeMapByName coninfo_map_;

  int nUnfxIntVars_ = 0;
  VarInfo var_info_ {};
  ObjInfo obj_info_;
};

} // namespace mp

#endif // MODEL_INFO_HPP
