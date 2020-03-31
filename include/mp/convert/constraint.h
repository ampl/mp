#ifndef CONSTRAINT_H
#define CONSTRAINT_H

namespace mp {

template <class Converter, class Backend>
class Constraint {
public:
  virtual ~Constraint() { }
  /// This normally dispatches conversion (decomposition) to the Converter
  virtual void ConvertWith(Converter& cvt) = 0;
  /// This adds the constraint to the backend without conversion
  virtual void AddToBackend(Backend& be) = 0;
};

} // namespace mp

#endif // CONSTRAINT_H
