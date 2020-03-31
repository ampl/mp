#ifndef CONSTRAINT_H
#define CONSTRAINT_H

namespace mp {

/// Converters handling custom constraints should derive from
class BasicConstraintConverter {
public:
};

/// Backends handling custom constraints should derive from
class BasicConstraintAdder {
public:
};

/// Custom constraints to derive from
class BasicConstraint {
public:
  virtual ~BasicConstraint() { }
  /// This normally dispatches conversion (decomposition) to the Converter
  virtual void ConvertWith(BasicConstraintConverter& cvt) = 0;
  /// This adds the constraint to the backend without conversion
  virtual void AddToBackend(BasicConstraintAdder& be) = 0;
};

} // namespace mp

#endif // CONSTRAINT_H
