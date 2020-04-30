#ifndef CONSTRAINT_KEEPER_H
#define CONSTRAINT_KEEPER_H

namespace mp {

/// Converters handling custom constraints should derive from
class BasicConstraintConverter {
public:
  /// By default, we complain about someone trying to convert an unknown constraint
  template <class Constraint>
  void Convert(const Constraint& ) {
    throw std::logic_error(
          std::string("Not converting constraint ") + typeid(Constraint).name());
  }
  /// Derived converter classes have to tell C++ to use default handlers if they need them
  /// when they overload Convert(), due to C++ name hiding
#define USE_BASE_CONSTRAINT_CONVERTERS(BaseConverter) \
  using BaseConverter::Convert;
};

class BasicConstraint;

/// Backends handling custom constraints should derive from
class BasicConstraintAdder {
public:
  template <class Constraint>
  void AddConstraint(const Constraint& ) {
    throw std::logic_error(
          std::string("Not handling constraint ") + typeid(Constraint).name());
  }
  /// Level of acceptance of a constraint by a backend
  enum ConstraintAcceptanceLevel {
    NotAccepted,
    AcceptedButNotRecommended,
    Recommended
  };
  /// Derived backends have to tell C++ to use default handlers if they are needed
  /// when they overload AddConstraint(), due to C++ name hiding
#define USE_BASE_CONSTRAINT_HANDLERS(BaseBackend) \
  using BaseBackend::AddConstraint; \
  using BaseBackend::AcceptanceLevel;
  /// By default, we say constraint XYZ is not accepted
  static constexpr ConstraintAcceptanceLevel AcceptanceLevel(const BasicConstraint*) { return NotAccepted; }
  /// Then for a certain constraint it can be specified
#define ACCEPT_CONSTRAINT(ConstrType, level) \
  static constexpr ConstraintAcceptanceLevel AcceptanceLevel(const ConstrType*) { return level; }

};

class BasicConstraintKeeper {
public:
  virtual ~BasicConstraintKeeper() { }
  virtual bool IsRemoved() const = 0;
  virtual void Remove() = 0;
  /// This normally dispatches conversion (decomposition) to the Converter
  virtual void ConvertWith(BasicConstraintConverter& cvt) = 0;
  /// Checks backend's acceptance level for the constraint
  virtual BasicConstraintAdder::ConstraintAcceptanceLevel BackendAcceptance(
      const BasicConstraintAdder& ) const = 0;
  /// This adds the constraint to the backend without conversion
  virtual void AddToBackend(BasicConstraintAdder& be) const = 0;
};

template <class Converter, class Backend, class Constraint>
class ConstraintKeeper : public BasicConstraintKeeper {
  Constraint cons_;
  bool is_removed_ = false;
public:
  template <class... Args>
  ConstraintKeeper(Args&&... args) : cons_(std::move(args)...) { }
  bool IsRemoved() const override { return is_removed_; }
  void Remove() override { is_removed_=true; }
  void ConvertWith(BasicConstraintConverter& cvt) override {
    try {
      static_cast<Converter&>(cvt).Convert(cons_);
    } catch (const std::exception& exc) {
      throw std::logic_error(typeid(Converter).name() + std::string(": ") + exc.what());
    }
  }
  BasicConstraintAdder::ConstraintAcceptanceLevel BackendAcceptance(
      const BasicConstraintAdder& ba) const override {
    return static_cast<const Backend&>( ba ).AcceptanceLevel(&cons_);
  }
  void AddToBackend(BasicConstraintAdder& be) const override {
    try {
      static_cast<Backend&>(be).AddConstraint(cons_);
    } catch (const std::exception& exc) {
      throw std::logic_error(typeid(Backend).name() + std::string(": ") + exc.what());
    }
  }
};

/// Helper function constructing a constraint keeper
template <class Converter, class Constraint, class... Args>
BasicConstraintKeeper *makeConstraint(Args... args) {
  return
    new ConstraintKeeper<Converter, typename Converter::BackendType, Constraint>(
        std::forward<Args>(args)...);
}

} // namespace mp

#endif // CONSTRAINT_KEEPER_H
