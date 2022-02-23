#ifndef CONSTRAINT_KEEPER_H
#define CONSTRAINT_KEEPER_H

#include <deque>
#include <map>
#include <cmath>

#include "mp/common.h"
#include "mp/flat/model_api_base.h"
#include "mp/presolve-node.h"

namespace mp {

/// Converters handling custom constraints should derive from
class BasicFlatConverter;

/// Interface for an array of constraints of certain type
class BasicConstraintKeeper {
public:
  /// Destructor
  virtual ~BasicConstraintKeeper() { }
  /// Constructor
  BasicConstraintKeeper(const char* nm) :
    value_node_(nm), constr_name_(nm) { }
  /// Constraint type
  using ConstraintType = BasicConstraint;
  /// Constraint keeper description
  virtual const std::string& GetDescription() const = 0;
  /// Propagate expression result of constraint \a i top-down
  virtual void PropagateResult(BasicFlatConverter& cvt,
                               int i,
                               double lb, double ub, Context ctx) = 0;
  /// Result variable of constraint \a i. Returns -1 if none
  virtual int GetResultVar(int i) const = 0;
  /// Convert all new items of this constraint.
  /// This normally dispatches conversion (decomposition) to the Converter
  /// @return whether any converted
  virtual bool ConvertAllNewWith(BasicFlatConverter& cvt) = 0;
  /// Checks backend's acceptance level for the constraint
  virtual ConstraintAcceptanceLevel BackendAcceptance(
      const BasicFlatModelAPI& ) const = 0;
  /// Checks backend's group number for the constraint
  virtual int BackendGroup(const BasicFlatModelAPI& ) const = 0;
  /// This adds all unbridged items to the backend (without conversion)
  virtual void AddUnbridgedToBackend(BasicFlatModelAPI& be) = 0;

  /// Value presolve node, const
  const pre::ValueNode& GetValueNode() const { return value_node_; }
  /// Value presolve node
  pre::ValueNode& GetValueNode() { return value_node_; }

private:
  pre::ValueNode value_node_;
  const char* constr_name_;
};


/// Full id of a constraint: CK + index
/// This helper class is parameterized by the Keeper
template <class ConstraintKeeper>
struct ConstraintLocationHelper {
  /// Default constructor: no valid Id
  ConstraintLocationHelper() = default;
  /// Normal constructor
  ConstraintLocationHelper(ConstraintKeeper* pck, int i) noexcept :
    pck_(pck), index_(i) { }

  /// Checks if we store a constraint's location
  operator bool() const { return HasId(); }
  /// Checks if we store a constraint's location
  bool HasId() const { return nullptr!=pck_; }
  /// High-level getter
  int GetResultVar() const { return GetCK()->GetResultVar(GetIndex()); }
  /// High-level getter
  const typename ConstraintKeeper::ConstraintType&
  GetConstraint() const { return GetCK()->GetConstraint(GetIndex()); }
  /// High-level getter
  typename ConstraintKeeper::ConstraintType&
  GetConstraint() { return GetCK()->GetConstraint(GetIndex()); }

  /// Get/set data
  ConstraintKeeper* GetCK() const { assert(HasId()); return pck_; }
  int GetIndex() const { return index_; }

  void SetCK(ConstraintKeeper* pck) { pck_ = pck; }
  void SetIndex(int i) { index_ = i; }

  ConstraintKeeper* pck_ = nullptr;
  int index_ = 0;        // constraint index
};

/// Without constraint type
using AbstractConstraintLocation =
  ConstraintLocationHelper<BasicConstraintKeeper>;

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

  /// Default conversion priority
  static constexpr double ConstraintCvtPriority(BasicConstraint*) { return 1.0; }

  /// By default, we complain about someone trying to convert an unknown constraint
  template <class Constraint>
  void Convert(const Constraint& ) {
    throw std::logic_error(
          std::string("Not converting constraint ") + Constraint::GetName());
  }

  /// Derived converter classes have to tell C++ to use
  /// default handlers if they need them
  /// when they overload Convert(), due to C++ name hiding
#define USE_BASE_CONSTRAINT_CONVERTERS(BaseConverter) \
  using BaseConverter::PreprocessConstraint; \
  using BaseConverter::PropagateResult; \
  using BaseConverter::Convert


  /// For Common Subexpression Elimination, we can use maps
  /// This stub returns empty Id
  AbstractConstraintLocation MapFind(const BasicConstraint& )
  { return { }; }

  /// Returns false when we do have a map and entry duplicated
  /// (should not happen).
  /// Can be conveniently overloaded for
  /// ConstraintLocation<Constraint> defined in Converters
  template <class CK>
  bool MapInsert(ConstraintLocationHelper<CK> ) { return true; }

  /// Similarly to Convert(),
  /// need to 'using' base class' map accessors in the Converter
#define USE_BASE_MAP_FINDERS(BaseConverter) \
  using BaseConverter::MapFind; \
  using BaseConverter::MapInsert; \
  template <class Constraint> \
  using ConstraintLocation = \
    ConstraintLocationHelper< ConstraintKeeper< Impl, Backend, Constraint > >;


  static constexpr double Infty() { return INFINITY; }
  static constexpr double MinusInfty() { return -INFINITY; }
  static constexpr double PracticallyInfty() { return 1e20; }  // TODO options
  static constexpr double PracticallyMinusInfty() { return -1e20; }
};

/// Conversion failure helper
class ConstraintConversionFailure {
  const char *key_, *msg_;
public:
  ConstraintConversionFailure(const char* key, const char* msg) noexcept :
    key_(key), msg_(msg) { }
  /// Failure type, used to display infos about failures
  const char* key() const { return key_; }
  /// Detailed message, should help improve model
  const char* message() const { return msg_; }
};


/// Specialize ConstraintKeeper for a given constraint type
/// to store an array of such constraints
template <class Converter, class Backend, class Constraint>
class ConstraintKeeper : public BasicConstraintKeeper {
public:
  /// Constraint type
  using ConstraintType = Constraint;
  /// Constrint Keeper description
  const std::string& GetDescription() const override
  { return desc_; }
  /// Assume Converter has the Backend
  Backend& GetBackend(BasicFlatConverter& cvt)
  { return static_cast<Converter&>(cvt).GetBackend(); }
  /// Constructor, adds this CK to the provided ConstraintManager
  /// Requires the CM to be already constructed
  ConstraintKeeper(Converter& cvt, const char* nm) :
    BasicConstraintKeeper(nm), cvt_(cvt) {
    GetConverter().AddConstraintKeeper(*this, ConversionPriority());
  }
  /// Add a pre-constructed constraint (or just arguments)
  /// @return index of the new constraint
  template <class... Args>
  int AddConstraint(Args&&... args) noexcept
  {
    cons_.emplace_back( std::move(args)... );
    return cons_.size()-1;
  }
  /// Get const constraint \a i
  const Constraint& GetConstraint(int i) const
  { assert(check_index(i)); return cons_[i].con_; }
  /// Get constraint \a i
  Constraint& GetConstraint(int i)
  { assert(check_index(i)); return cons_[i].con_; }
  /// Propagate expression result of constraint \a i top-down
  void PropagateResult(BasicFlatConverter& cvt,
                       int i,
                       double lb, double ub, Context ctx) override {
    try {
      static_cast<Converter&>(cvt).PropagateResult(
            GetConstraint(i), lb, ub, ctx);
    } catch (const std::exception& exc) {
      throw std::logic_error(Converter::GetConverterName() +
                             std::string(": propagating result for constraint ") +
                             std::to_string(i) + " of type '" +
                             Constraint::GetName() +
                             "':  " + exc.what());
    }
  }
  /// Result variable of constraint \a i. Returns -1 if none
  int GetResultVar(int i) const override
  { assert(check_index(i)); return cons_[i].con_.GetResultVar(); }
  /// Conversion priority. Uses that from Converter
  double ConversionPriority() const
  { return Converter::ConstraintCvtPriority((Constraint*)nullptr); }
  /// Convert all new items of this constraint.
  /// This normally dispatches conversion (decomposition) to the Converter
  /// @return whether any converted
  bool ConvertAllNewWith(BasicFlatConverter& cvt) override {
    assert(&cvt == &GetConverter());         // Using the same Converter
    try {
      return ConvertAllFrom(i_cvt_last_);
    } catch (const std::exception& exc) {
      throw std::logic_error(Converter::GetConverterName() + std::string(": ")
                             + exc.what());
    }
    return false;
  }
  ConstraintAcceptanceLevel BackendAcceptance(
      const BasicFlatModelAPI& ba) const override {
    return static_cast<const Backend&>( ba ).AcceptanceLevel((Constraint*)nullptr);
  }
  int BackendGroup(
      const BasicFlatModelAPI& ba) const override {
    return static_cast<const Backend&>( ba ).GroupNumber((Constraint*)nullptr);
  }
  void AddUnbridgedToBackend(BasicFlatModelAPI& be) override {
    try {
      AddAllUnbridged(be);
    } catch (const std::exception& exc) {
      throw std::logic_error(std::string("Adding constraint '") +
                             Constraint::GetName() + "' to " +
                             Backend::GetName() + std::string(": ") +
                             exc.what());
    }
  }

protected:
  /// Retrieve the Converter, const
  const Converter& GetConverter() const { return cvt_; }
  /// Retrieve the Converter
  Converter& GetConverter() { return cvt_; }

  /// Check constraint index
  bool check_index(int i) const { return i>=0 && i<(int)cons_.size(); }
  /// Container for a single constraint
  struct Container {
    Container(Constraint&& c) noexcept : con_(std::move(c)) { }

    bool IsBridged() const { return is_bridged_; }
    void MarkAsBridged() { is_bridged_=true; }

    Constraint con_;
    bool is_bridged_ = false;
  };
  bool ConvertAllFrom(int& i_last) {
    int i=i_last;
    ++i;
    const auto acceptanceLevel =
        BackendAcceptance(GetBackend(GetConverter()));
    if (NotAccepted == acceptanceLevel) {       // Convert all
      for (auto it=cons_.begin()+i; it!=cons_.end(); ++it, ++i)
        if (!it->IsBridged())
          ConvertConstraint(*it, i);
    }
    else if (AcceptedButNotRecommended == acceptanceLevel) {
      for (auto it=cons_.begin()+i; it!=cons_.end(); ++it, ++i) {
        if (!it->IsBridged()) {
          try {
            ConvertConstraint(*it, i);
          } catch (const ConstraintConversionFailure& ccf) {
            GetConverter().AddWarning( ccf.key(), ccf.message() );
          }
        }
      }
    }
    bool any_converted = i_last!=i-1;
    i_last = i-1;
    return any_converted;
  }
  /// Call Converter's RunConversion() and mark as "bridged"
  ///
  /// @param cnt the constraint container
  /// actually redundant as i is enough to find. But for speed
  /// @param i constraint index, needed for bridging
  void ConvertConstraint(Container& cnt, int i) {
    assert(!cnt.IsBridged());
    GetConverter().RunConversion(cnt.con_, i);
    cnt.MarkAsBridged();    // TODO should this be marked in Convert()?
  }
  void AddAllUnbridged(BasicFlatModelAPI& be) {
    int con_index=0;
    auto con_group = BackendGroup(be);
    for (const auto& cont: cons_) {
      if (!cont.IsBridged()) {
        static_cast<Backend&>(be).AddConstraint(cont.con_);
        GetConverter().GetCopyBridge().
            AddEntry({
                       GetValueNode().Select(con_index),
                       GetConverter().GetPresolver().GetTargetNodes().
                         GetConValues()(con_group).Add()
                     });
      }
      ++con_index;                      // increment index
    }
  }

private:
  Converter& cvt_;
  std::deque<Container> cons_;        // TODO see if vector is faster
  int i_cvt_last_ = -1;               // last converted constraint
  const std::string desc_ {
    std::string("ConstraintKeeper< ") +
        Converter::GetConverterName() + ", " +
        Backend::GetName() + ", " +
        Constraint::GetName() + " >"};
};

/// Macros to define / access constraint keepers
/// Assume ConstraintManager as public parent

/// Use this to obtain a certain keeper, const
#define GET_CONST_CONSTRAINT_KEEPER(Constraint) \
  MPCD( GetConstraintKeeper((Constraint*)nullptr) )
/// Use this to obtain a certain keeper
#define GET_CONSTRAINT_KEEPER(Constraint) \
  MPD( GetConstraintKeeper((Constraint*)nullptr) )

/// Define a constraint keeper
#define STORE_CONSTRAINT_TYPE(Constraint) \
  ConstraintKeeper<Impl, Backend, Constraint> \
    CONSTRAINT_KEEPER_VAR(Constraint) \
      {*static_cast<Impl*>(this), #Constraint}; \
  const ConstraintKeeper<Impl, Backend, Constraint>& \
  GetConstraintKeeper(Constraint* ) const { \
    return CONSTRAINT_KEEPER_VAR(Constraint); \
  } \
  ConstraintKeeper<Impl, Backend, Constraint>& \
  GetConstraintKeeper(Constraint* ) { \
    return CONSTRAINT_KEEPER_VAR(Constraint); \
  }

/// Internal use
#define CONSTRAINT_KEEPER_VAR(Constraint) \
  ck__ ## Constraint ## _


/// Manage ConstraintKeepers for different constraint types
class ConstraintManager {
  std::multimap<double, BasicConstraintKeeper&> con_keepers_;
public:
  /// Add a new CKeeper with given conversion priority (smaller = sooner)
  void AddConstraintKeeper(BasicConstraintKeeper& ck, double priority)
  { con_keepers_.insert( { priority, ck } ); }

  /// Convert all constraints (including any new appearing)
  void ConvertAllConstraints(BasicFlatConverter& cvt) {
    bool any_converted;
    do {
      any_converted = false;
      for (auto& ck: con_keepers_)
        any_converted = any_converted || ck.second.ConvertAllNewWith(cvt);
    } while (any_converted);
  }

  /// Add all unbridged constraints to Backend
  void AddUnbridgedConstraintsToBackend(BasicFlatModelAPI& be) const {
    for (const auto& ck: con_keepers_)
      ck.second.AddUnbridgedToBackend(be);
  }
};

} // namespace mp

#endif // CONSTRAINT_KEEPER_H
