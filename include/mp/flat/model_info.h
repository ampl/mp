#ifndef MODEL_INFO_H
#define MODEL_INFO_H

/**
 * Public interface for flat model info
 */

#include <typeinfo>
#include <memory>

namespace mp {

/// Public interface for flat model info.
/// Can be used to initialize storage
class FlatModelInfo {
public:
  /// Virtual destructor
  virtual ~FlatModelInfo() { }

  /// Get number of constraints of certain group
  virtual int GetNumberOfConstraintsOfGroup(int ng) const =0;

  /// Get number of constraints of single type
  virtual int GetNumberOfConstraints(const std::type_info& nt) const =0;

  /// Initialize constraint counting
  virtual void InitConstraintCount() =0;

  /// Add number of constraints of single type
  virtual void AddNumberOfConstraints(
      const std::type_info& ti, int igroup, int nc) = 0;
};


/// FlatModelInfo factory
std::unique_ptr<const FlatModelInfo> CreateFlatModelInfo();

} // namespace mp

#endif // MODEL_INFO_H
