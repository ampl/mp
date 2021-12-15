#ifndef CONSTRAINT_KEEPER_H
#define CONSTRAINT_KEEPER_H

#include <limits>

#include "mp/common.h"
#include "mp/flat/flat_model_api_basic.h"

namespace mp {

/// An array of constraints of a single type
class BasicConstraintKeeper;

/// A reference to a stored constraint of a known type
using ConstraintKeeperId = int;

/// Converters handling custom constraints should derive from
class BasicFlatConverter {
public:
  /// Default constraint prepro
  /// All parameters are 'in-out'
  template <class Constraint, class PreproInfo>
  void PreprocessConstraint( Constraint&, PreproInfo& ) {
    // ... do nothing by default
    // Should at least derive bounds & type for the result
  }
  /// TODO incapsulate parameters
  void PropagateResult(BasicConstraint& con, double lb, double ub, Context ctx) {
    internal::Unused(con, lb, ub, ctx);
    throw std::logic_error("This should not be called");
  }
  /// By default, we complain about someone trying to convert an unknown constraint
  template <class Constraint>
  void Convert(const Constraint& ) {
    throw std::logic_error(
          std::string("Not converting constraint ") + Constraint::GetConstraintName());
  }
  /// Derived converter classes have to tell C++ to use default handlers if they need them
  /// when they overload Convert(), due to C++ name hiding
#define USE_BASE_CONSTRAINT_CONVERTERS(BaseConverter) \
  using BaseConverter::PreprocessConstraint; \
  using BaseConverter::PropagateResult; \
  using BaseConverter::Convert


  /// For Common Subexpression Elimination, we can use maps
  const BasicConstraintKeeper* MapFind(const BasicConstraint& ) const { return nullptr; }
  /// Returns false when we do have a map and entry duplicated
  bool MapInsert(const BasicConstraintKeeper* ) { return true; }
#define USE_BASE_MAP_FINDERS(BaseConverter) \
  using BaseConverter::MapFind; \
  using BaseConverter::MapInsert;

  static constexpr double Infty() { return std::numeric_limits<double>::infinity(); }
  static constexpr double MinusInfty() { return -std::numeric_limits<double>::infinity(); }
};

/// Conversion failure helper
class ConstraintConversionFailure {
  const std::string msg_;
public:
  ConstraintConversionFailure(std::string msg) noexcept :
    msg_(std::move(msg)) { }
  const std::string& message() const { return msg_; }
};



/// A derived class stores an array of constraints of certain type
class BasicConstraintKeeper {
public:
  virtual ~BasicConstraintKeeper() { }
  /// Constraint keeper description
  virtual const std::string& GetDescription() const = 0;
  /// The item was bridged, don't pass it further
  virtual bool IsBridged() const = 0;
  /// Mark as bridged
  virtual void MarkAsBridged() = 0;
  /// Propagate expression result top-down
  virtual void PropagateResult(BasicFlatConverter& cvt,
                               double lb, double ub, Context ctx) = 0;
  /// Returns -1 if none
  virtual int GetResultVar() const = 0;
  /// This normally dispatches conversion (decomposition) to the Converter
  virtual void ConvertWith(BasicFlatConverter& cvt) = 0;
  /// Checks backend's acceptance level for the constraint
  virtual ConstraintAcceptanceLevel BackendAcceptance(
      const BasicFlatBackend& ) const = 0;
  /// This adds the constraint to the backend without conversion
  virtual void AddToBackend(BasicFlatBackend& be) const = 0;

  /// Pre- / postsolve: elementary interface
  /// 1: linear, 2: qcp
  virtual int ConstraintClass() const = 0;
};

/// Specialize ConstraintKeeper for a given constraint type
/// to store an array of such constraints
template <class Converter, class Backend, class Constraint>
class ConstraintKeeper : public BasicConstraintKeeper {
  Constraint cons_;
  bool is_removed_ = false;
  const std::string desc_ {
    std::string("ConstraintKeeper< ") +
        Converter::GetConverterName() + ", " +
        Backend::GetBackendName() + ", " +
        Constraint::GetConstraintName() + " >"};
public:
  template <class... Args>
  ConstraintKeeper(Args&&... args) noexcept
    : cons_(std::move(args)...) { }
  const std::string& GetDescription() const override {
    return desc_.c_str();
  }
  const Constraint& GetConstraint() const { return cons_; }
  Constraint& GetConstraint() { return cons_; }
  bool IsBridged() const override { return is_removed_; }
  void MarkAsBridged() override { is_removed_=true; }
  void PropagateResult(BasicFlatConverter& cvt,
                       double lb, double ub, Context ctx) override {
    try {
      static_cast<Converter&>(cvt).PropagateResult(cons_, lb, ub, ctx);
    } catch (const std::exception& exc) {
      throw std::logic_error(Converter::GetConverterName() +
                             std::string(": propagating result for constraint ") +
                             + Constraint::GetConstraintName() +
                             ":  " + exc.what());
    }
  }
  virtual int GetResultVar() const { return cons_.GetResultVar(); }
  void ConvertWith(BasicFlatConverter& cvt) override {
    try {
      static_cast<Converter&>(cvt).RunConversion(cons_);
    } catch (const std::exception& exc) {
      throw std::logic_error(Converter::GetConverterName() + std::string(": ")
                             + exc.what());
    }
  }
  ConstraintAcceptanceLevel BackendAcceptance(
      const BasicFlatBackend& ba) const override {
    return static_cast<const Backend&>( ba ).AcceptanceLevel(&cons_);
  }
  void AddToBackend(BasicFlatBackend& be) const override {
    try {
      static_cast<Backend&>(be).AddConstraint(cons_);
    } catch (const std::exception& exc) {
      throw std::logic_error(Backend::GetBackendName() + std::string(": ") +
                             exc.what());
    }
  }
  int ConstraintClass() const
  { return Converter::ConstraintClass((Constraint*)nullptr); }
};

/// Helper function constructing a constraint keeper
template <class Converter, class Constraint, class... Args>
ConstraintKeeper<Converter, typename Converter::BackendType, Constraint>
  *makeConstraintKeeper(Args... args) {
  return
    new ConstraintKeeper<Converter, typename Converter::BackendType, Constraint>(
        std::forward<Args>(args)...);
}

} // namespace mp

#endif // CONSTRAINT_KEEPER_H
