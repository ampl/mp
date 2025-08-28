#ifndef MODEL_INFO_H
#define MODEL_INFO_H

/**
 * Public interface for flat model info
 */

#include <array>
#include <typeinfo>
#include <memory>

#include "mp/suffix.h"

namespace mp {

/// Public interface for flat model info.
/// Information about the model
/// can be used to initialize storage.
class FlatModelInfo {
public:
  /// Virtual destructor
  virtual ~FlatModelInfo() { }

  /// Number of unfixed integer variables
  virtual int NumUnfixedIntVars() const =0;

  /// Set N unfixed int vars
  virtual void SetNumUnfixedIntVars(int n) =0;

  /// Full var info: for each of (original, auxiliary):
  /// all, int, binary
  using VarInfo = std::array<int, 6>;

  /// Get var info
  virtual VarInfo GetVarInfo() const = 0;

  /// Set var info
  virtual void SetVarInfo(VarInfo ) = 0;

  /// Obj info: N lin, QP, NL
  using ObjInfo = std::array<int, 3>;

  /// Get obj info
  virtual ObjInfo GetObjInfo() const = 0;

  /// Set obj info
  virtual void SetObjInfo(ObjInfo ) = 0;

  /// Get number of constraints of certain group
  virtual int GetNumberOfConstraintsOfGroup(int ng) const =0;

  /// Get number of constraints of single type
  virtual int GetNumberOfConstraints(const std::type_info& nt) const =0;

  /// Constraint type info
  struct ConTypeInfo {
    const char* name_ {nullptr};
    bool is_logical_ {0};
    int n_ {0};               // number of active
    int n_used_ {0};          // used (active + redefined)
    int n_total_ {0};         // total, including unused
    /// operator==, needed when checking if the model changed
    bool operator==(const ConTypeInfo& cti) const {
      assert(name_==cti.name_);
      return n_==cti.n_
             && n_used_==cti.n_used_
             && n_total_==cti.n_total_;
    }
  };

  /// Map by constraint name
  using ConstrTypeMapByName = std::map<std::string, ConTypeInfo>;

  /// Obtain constraint types
  virtual const ConstrTypeMapByName& GetConstraintTypes() const = 0;

  /// Initialize constraint counting
  virtual void InitConstraintCount() =0;

  /// Add number of constraints of single type
  /// to the counter.
  virtual void AddNumberOfConstraints(
      const std::type_info& ti, const char* name,
      int igroup, bool is_logical, int nc, int nu, int na) = 0;
};


/// FlatModelInfo factory
std::unique_ptr<FlatModelInfo> CreateFlatModelInfo();

/// Print on screen
void PrintModelInfo(const FlatModelInfo& fmi,
                    const char* header, bool aux_vars);

/// Report via suffixes
void ReportModelInfoSuffixes(const FlatModelInfo& fmi,
                             std::string suf_prefix,
                             SuffixGetterSetter sgs);

} // namespace mp

#endif // MODEL_INFO_H
