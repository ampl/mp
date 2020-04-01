#ifndef CONSTRAINT_KEEPER_H
#define CONSTRAINT_KEEPER_H

namespace mp {

class BasicConstraint;

/// Converters handling custom constraints should derive from
class BasicConstraintConverter {
public:
};

/// Backends handling custom constraints should derive from
class BasicConstraintAdder {
public:
  void AddConstraint(const BasicConstraint& ) {
    throw std::logic_error("Adding this constraint has not been implemented");
  }
};

class BasicConstraintKeeper {
public:
  virtual ~BasicConstraintKeeper() { }
  /// This normally dispatches conversion (decomposition) to the Converter
  virtual void ConvertWith(BasicConstraintConverter& cvt) = 0;
  /// This adds the constraint to the backend without conversion
  virtual void AddToBackend(BasicConstraintAdder& be) const = 0;
};

template <class Converter, class Backend, class Constraint>
class ConstraintKeeper : public BasicConstraintKeeper {
  Constraint cons_;
public:
  template <class... Args>
  ConstraintKeeper(Args&&... args) : cons_(std::forward<Args>(args)...) { }
  void ConvertWith(BasicConstraintConverter& cvt) override {
    throw std::runtime_error(
          std::string("Not dispatching redefinition of ") + typeid(Constraint).name() + " yet");
  }
  void AddToBackend(BasicConstraintAdder& be) const override {
    try {
      static_cast<Backend&>(be).AddConstraint(cons_);
    } catch (const std::exception& exc) {
      throw std::logic_error(typeid(Constraint).name() + std::string(": ") + exc.what());
    }
  }
};

} // namespace mp

#endif // CONSTRAINT_KEEPER_H
