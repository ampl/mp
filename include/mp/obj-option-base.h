#ifndef OBJ_OPTION_BASE_H
#define OBJ_OPTION_BASE_H

#include <vector>

namespace mp {

/// Abstract obj option setter
class BasicObjOptionSetter {
public:
  /// Destroy
  virtual ~BasicObjOptionSetter() { }
  /// Get vector of multi-obj passes with options
  virtual std::vector<int> GetPassesWithOptions() = 0;
  /// Set solver options for given pass
  virtual void SetOptionsForMultiobjPass(int iPass) = 0;
};

}  // namespace mp

#endif // OBJ_OPTION_BASE_H
