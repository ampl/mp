#ifndef CONVERTER_MODEL_BASE_H
#define CONVERTER_MODEL_BASE_H

#include <vector>

namespace mp {

/// Class BasicFlatModel.
/// Provides abstract interface to a few items in FlatModel
class BasicFlatModel {
public:
  /// Virtual Destruct
  virtual ~BasicFlatModel() { }

  /// Typedef VarBndVec
  using VarBndVec = std::vector<double>;

  /// Provide variable lower bounds
  virtual const VarBndVec& GetVarLBs() const = 0;
  /// Provide variable upper bounds
  virtual const VarBndVec& GetVarUBs() const = 0;
};

}  // namespace mp

#endif // CONVERTER_MODEL_BASE_H
