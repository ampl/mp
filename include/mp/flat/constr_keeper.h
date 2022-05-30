#ifndef CONSTRAINT_KEEPER_H
#define CONSTRAINT_KEEPER_H

#include <deque>
#include <map>
#include <unordered_map>
#include <functional>
#include <cmath>

#include "mp/common.h"
#include "mp/flat/model_api_base.h"
#include "mp/flat/constr_hash.h"
#include "mp/flat/redef/redef_base.h"
#include "mp/valcvt-node.h"

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

  /// Backend's acceptance level for the constraint type
  virtual ConstraintAcceptanceLevel BackendAcceptance(
      const BasicFlatModelAPI& ) const = 0;

  /// Backend's group number for the constraint type
  virtual int BackendGroup(const BasicFlatModelAPI& ) const = 0;

  /// This adds all unbridged items to the backend (without conversion)
  virtual void AddUnbridgedToBackend(BasicFlatModelAPI& be) = 0;

  /// Value presolve node, const
  const pre::ValueNode& GetValueNode() const { return value_node_; }

  /// Value presolve node
  pre::ValueNode& GetValueNode() { return value_node_; }

  /// Constraint name
  const char* GetConstraintName() const { return constr_name_; }

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

  /// Get Keeper
  ConstraintKeeper* GetCK() const { assert(HasId()); return pck_; }
  /// Get index
  int GetIndex() const { return index_; }

  /// Set Keeper
  void SetCK(ConstraintKeeper* pck) { pck_ = pck; }
  /// Set index
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
  /// Default conversion priority
  static constexpr double ConstraintCvtPriority(BasicConstraint*) { return 1.0; }

  /// Derived converter classes have to tell C++ to use
  /// default handlers if they need them
  /// when they overload Convert(), due to C++ name hiding
#define USE_BASE_CONSTRAINT_CONVERTERS(BaseConverter) \
  using BaseConverter::PreprocessConstraint; \
  using BaseConverter::PropagateResult; \
  using BaseConverter::Convert


  /// For Common Subexpression Elimination, we can use maps
  /// This stub returns empty Id
  int MapFind(const BasicConstraint& ) { return -1; }

  /// Returns false when we do have a map and entry duplicated
  /// (should not happen).
  /// Can be conveniently overloaded
  template <class Con>
  bool MapInsert(const Con& , int ) { return true; }

  /// Similarly to Convert(),
  /// need to 'using' base class' map accessors in the Converter
#define USE_BASE_MAP_FINDERS(BaseConverter) \
  using BaseConverter::MapFind; \
  using BaseConverter::MapInsert; \
  template <class Constraint> \
  using ConstraintLocation = \
    ConstraintLocationHelper< \
      ConstraintKeeper< Impl, ModelAPI, Constraint > >;


  static constexpr double Infty() { return INFINITY; }
  static constexpr double MinusInfty() { return -INFINITY; }
  static constexpr double PracticallyInfty() { return 1e20; }  // TODO options
  static constexpr double PracticallyMinusInfty() { return -1e20; }
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
  { return static_cast<Converter&>(cvt).GetModelAPI(); }

  /// Constructor, adds this CK to the provided ConstraintManager
  /// Requires the CM to be already constructed
  ConstraintKeeper(Converter& cvt, const char* nm) :
    BasicConstraintKeeper(nm), cvt_(cvt) {
    GetConverter().AddConstraintKeeper(*this, ConversionPriority());
  }

  /// Add a pre-constructed constraint (or just arguments)
  /// @return index of the new constraint
  template <class... Args>
  int AddConstraint(Args&&... args)
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
      MP_RAISE(Converter::GetTypeName() +
                             std::string(": propagating result for constraint ") +
                             std::to_string(i) + " of type '" +
                             Constraint::GetTypeName() +
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
      MP_RAISE(Converter::GetTypeName() + std::string(": ")
                             + exc.what());
    }
    return false;
  }

  /// Acceptance level of this constraint type in the Backend
  ConstraintAcceptanceLevel BackendAcceptance(
      const BasicFlatModelAPI& ba) const override {
    return static_cast<const Backend&>( ba ).AcceptanceLevel((Constraint*)nullptr);
  }

  /// Group number of this constraint type in the Backend.
  /// This is needed for pre- / postsolve to group solution values
  int BackendGroup(
      const BasicFlatModelAPI& ba) const override {
    return static_cast<const Backend&>( ba ).GroupNumber((Constraint*)nullptr);
  }

  /// Add remaining constraints to Backend
  void AddUnbridgedToBackend(BasicFlatModelAPI& be) override {
    try {
      AddAllUnbridged(be);
    } catch (const std::exception& exc) {
      MP_RAISE(std::string("Adding constraint of type '") +
                             Constraint::GetTypeName() + "' to " +
                             Backend::GetTypeName() + std::string(": ") +
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
        GetConverter().GetCopyLink().
            AddEntry({
                       GetValueNode().Select(con_index),
                       GetConverter().GetValuePresolver().GetTargetNodes().
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
        Converter::GetTypeName() + ", " +
        Backend::GetTypeName() + ", " +
        Constraint::GetTypeName() + " >"};
};


////////////////////////////////////////////////////////////////////////////////////
/// Macros to define / access constraint keepers
/// Assume ConstraintManager as public parent

/// Use this to obtain a certain keeper, const
#define GET_CONST_CONSTRAINT_KEEPER(Constraint) \
  MPCD( GetConstraintKeeper((Constraint*)nullptr) )
/// Use this to obtain a certain keeper
#define GET_CONSTRAINT_KEEPER(Constraint) \
  MPD( GetConstraintKeeper((Constraint*)nullptr) )

/// Use this to obtain a certain map, const
#define GET_CONST_CONSTRAINT_MAP(Constraint) \
  MPCD( GetConstraintMap((Constraint*)nullptr) )
/// Use this to obtain a certain map
#define GET_CONSTRAINT_MAP(Constraint) \
  MPD( GetConstraintMap((Constraint*)nullptr) )

/// Define a constraint keeper
/// without a subexpression map
#define STORE_CONSTRAINT_TYPE__INTERNAL(Constraint) \
  ConstraintKeeper<Impl, ModelAPI, Constraint> \
    CONSTRAINT_KEEPER_VAR(Constraint) \
      {*static_cast<Impl*>(this), #Constraint}; \
  const ConstraintKeeper<Impl, ModelAPI, Constraint>& \
  GetConstraintKeeper(Constraint* ) const { \
    return CONSTRAINT_KEEPER_VAR(Constraint); \
  } \
  ConstraintKeeper<Impl, ModelAPI, Constraint>& \
  GetConstraintKeeper(Constraint* ) { \
    return CONSTRAINT_KEEPER_VAR(Constraint); \
  }

/// Define a constraint keeper
/// without a subexpression map.
/// Provide empty MapFind (returns -1) / MapInsert
#define STORE_CONSTRAINT_TYPE__NO_MAP(Constraint) \
  STORE_CONSTRAINT_TYPE__INTERNAL(Constraint) \
  int MapFind__Impl(const Constraint& ) { return -1; } \
  bool MapInsert__Impl(const Constraint&, int ) \
    { return true; }


/// Define a constraint keeper
/// with a subexpression map.
/// The Converter storing the Constraint
/// should define MapFind / MapInsert accessing
/// the GET_(CONST_)CONSTRAINT_MAP(Constraint)
#define STORE_CONSTRAINT_TYPE__WITH_MAP(Constraint) \
  STORE_CONSTRAINT_TYPE__INTERNAL(Constraint) \
  STORE_CONSTRAINT_MAP(Constraint)

/// Internal use. Name of the constraint container
#define CONSTRAINT_KEEPER_VAR(Constraint) \
  ck__ ## Constraint ## _

/// Create constraint map. Normally internal use
#define STORE_CONSTRAINT_MAP(Constraint) \
  ConstraintMap<Constraint> CONSTRAINT_MAP_VAR(Constraint); \
  const ConstraintMap<Constraint>& \
  GetConstraintMap(Constraint* ) const { \
    return CONSTRAINT_MAP_VAR(Constraint); \
  } \
  ConstraintMap<Constraint>& \
  GetConstraintMap(Constraint* ) { \
    return CONSTRAINT_MAP_VAR(Constraint); \
  }

/// Internal use. Name of the constraint map
#define CONSTRAINT_MAP_VAR(Constraint) \
  map__ ## Constraint ## _



/////////////////////////////////////////////////////////////////////////////
/// Subexpression map
///
/// Indexes constraint location
template <class Constraint>
using ConstraintMap = std::unordered_map<
    std::reference_wrapper< const Constraint >, int >;

/// Subexpression maps for an expression require
/// operator==(refwrap<expr>, refwrap<expr>).
///
/// The below one is for CustomFunctionalConstraint<>
template <class Args, class Params, class NumOrLogic, class Id>
inline
bool operator==(std::reference_wrapper<
                  const CustomFunctionalConstraint<Args, Params, NumOrLogic, Id> > c1,
                std::reference_wrapper<
                  const CustomFunctionalConstraint<Args, Params, NumOrLogic, Id> > c2) {
  return c1.get().GetArguments() == c2.get().GetArguments() &&
      c1.get().GetParameters() == c2.get().GetParameters();
}

/// operator==(refwrap<ConditionalConstraint<> >)
template <class Con>
inline
bool operator==(std::reference_wrapper<
                  const ConditionalConstraint<Con> > c1,
                std::reference_wrapper<
                  const ConditionalConstraint<Con> > c2) {
  return c1.get().GetConstraint() == c2.get().GetConstraint();
}

/// operator==(LFC)
inline
bool operator==(std::reference_wrapper<
                  const LinearFunctionalConstraint > c1,
                std::reference_wrapper<
                  const LinearFunctionalConstraint > c2) {
  return c1.get().GetAffineExpr() == c2.get().GetAffineExpr();
}

/// operator==(QFC)
inline
bool operator==(std::reference_wrapper<
                  const QuadraticFunctionalConstraint > c1,
                std::reference_wrapper<
                  const QuadraticFunctionalConstraint > c2) {
  return c1.get().GetQuadExpr() == c2.get().GetQuadExpr();
}


////////////////////////////////////////////////////////////////////////////////////
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
