#ifndef CONSTRAINT_KEEPER_H
#define CONSTRAINT_KEEPER_H

#include <deque>
#include <unordered_map>
#include <functional>
#include <cmath>

#include "mp/common.h"
#include "mp/format.h"
#include "mp/env.h"
#include "mp/flat/model_api_base.h"
#include "mp/flat/constr_hash.h"
#include "mp/flat/redef/redef_base.h"
#include "mp/valcvt-node.h"

namespace mp {

/// Converters handling custom constraints should derive from
class BasicFlatConverter;


static const mp::OptionValueInfo values_item_acceptance[] = {
  { "0", "Not accepted natively, automatic redefinition will be attempted", 0},
  { "1", "Accepted but automatic redefinition will be used where possible", 1},
  { "2", "Accepted natively and preferred", 2}
};


/// Violation summary for a class of vars/cons/objs
struct ViolSummary {
  /// Check if this violation should be counted
  void CheckViol(double val, double eps) {
    if (val > eps) {
      ++N_;
      if (epsMax_ < val)
        epsMax_ = val;
    }
  }
  /// Count violation
  void CountViol(double val) {
    ++N_;
    if (epsMax_ < val)
      epsMax_ = val;
  }
  int N_ {0};
  double epsMax_ {0.0};
};

/// Array of violation summaries.
/// For different kinds, e.g., original / aux vars.
template <int Nkinds>
using ViolSummArray = std::array<ViolSummary, Nkinds>;

/// Solution check data
struct SolCheck {
  /// Construct
  SolCheck(ArrayRef<double> x,
           const pre::ValueMapDbl& duals,
           ArrayRef<double> obj,
           double feastol, double inttol)
    : x_(x), y_(duals), obj_(obj),
      feastol_(feastol), inttol_(inttol) { }
  /// Any violations?
  bool HasAnyViols() const { return hasAnyViol_; }
  /// Summary
  const std::string& GetReport() const { return report_; }

  /// x
  ArrayRef<double>& x() { return x_; }
  /// x[i]
  double x(int i) const { return x_[i]; }
  /// Feasibility tolerance
  double GetFeasTol() const { return feastol_; }

  /// Var bnd violations
  ViolSummArray<2>& VarViolBnds() { return viol_var_bnds_; }
  /// Var int-ty violations
  ViolSummArray<2>& VarViolIntty() { return viol_var_int_; }

  /// Constraints: algebraic.
  /// Map by constraint type.
  /// Values: for original, intermediate,
  /// and solver-side constraints.
  std::map< std::string, ViolSummArray<3> >&
  ConViolAlg() { return viol_cons_alg_; }
  /// Constraints: logical.
  std::map< std::string, ViolSummArray<3> >&
  ConViolLog() { return viol_cons_log_; }

private:
  ArrayRef<double> x_;
  const pre::ValueMapDbl& y_;
  ArrayRef<double> obj_;
  double feastol_;
  double inttol_;

  bool hasAnyViol_ = false;
  std::string report_;

  /// Variable bounds: orig, aux
  ViolSummArray<2> viol_var_bnds_;
  /// Variable integrality: orig, aux
  ViolSummArray<2> viol_var_int_;
  /// Constraints: algebraic.
  std::map< std::string, ViolSummArray<3> > viol_cons_alg_;
  /// Constraints: logical.
  std::map< std::string, ViolSummArray<3> > viol_cons_log_;
  /// Objectives
  ViolSummary viol_obj_;
};


/// Interface for an array of constraints of certain type
class BasicConstraintKeeper {
public:
  /// Destructor
  virtual ~BasicConstraintKeeper() { }

  /// Constructor
  BasicConstraintKeeper(
      pre::BasicValuePresolver& pres,
      const char* nm, const char* optN) :
    value_node_(pres, nm),
    constr_name_(nm), solver_opt_nm_(optN) { }

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
  /// This normally dispatches conversion (decomposition)
  ///  to the Converter
  /// @return whether any converted
  virtual bool ConvertAllNewWith(BasicFlatConverter& cvt) = 0;

  /// Query (user-chosen, if sensible) constraint acceptance level
  virtual ConstraintAcceptanceLevel GetChosenAcceptanceLevel()
  const {
    assert(0<=acceptance_level_);       // has been initialized
    return ConstraintAcceptanceLevel(acceptance_level_);
  }

  /// Converter's ability to convert the constraint type
  virtual bool IfConverterConverts(
      BasicFlatConverter& cvt ) const = 0;

  /// ModelAPI's acceptance level for the constraint type.
  /// This should not be used directly, instead:
  /// GetChosenAcceptanceLevel()
  virtual ConstraintAcceptanceLevel GetModelAPIAcceptance(
      const BasicFlatModelAPI& ) const = 0;

  /// Constraint type_info
  virtual const std::type_info& GetTypeInfo() const =0;

  /// Backend's group number for the constraint type
  virtual int GetConstraintGroup(const BasicFlatModelAPI& ) const = 0;

  /// Report how many will be added to Backend
  virtual int GetNumberOfAddable() const = 0;

  /// This adds all unbridged items to the backend (without conversion)
  virtual void AddUnbridgedToBackend(BasicFlatModelAPI& be) = 0;

  /// Value presolve node, const
  const pre::ValueNode& GetValueNode() const { return value_node_; }

  /// Value presolve node
  pre::ValueNode& GetValueNode() { return value_node_; }

  /// Create ValueNode range pointer: add n elements
  pre::NodeRange AddValueNodeRange(int n=1)
  { return GetValueNode().Add(n); }

  /// Create ValueNode range pointer: select n elements at certain pos
  pre::NodeRange SelectValueNodeRange(int pos, int n=1)
  { return GetValueNode().Select(pos, n); }

  /// Constraint name
  const char* GetConstraintName() const { return constr_name_; }

  /// Acceptance option names
  virtual const char* GetAcceptanceOptionNames() const
  { return solver_opt_nm_; }

  /// See what options are available for this constraint:
  /// whether it is accepted natively by ModelAPI and/or can be
  /// converted by the Converter.
  /// If both, add constraint acceptance option.
  /// This should be called before using the class.
  virtual void ConsiderAcceptanceOptions(
      BasicFlatConverter& cvt,
      const BasicFlatModelAPI& ma,
      Env& env) {
    auto cancvt = IfConverterConverts(cvt);
    SetChosenAcceptanceLevel( GetModelAPIAcceptance(ma) );
    if (cancvt) {         // See if ModelAPI accepts too
      if (ConstraintAcceptanceLevel::Recommended ==
            GetChosenAcceptanceLevel()) {
        env.AddStoredOption(GetAcceptanceOptionNames(),
                            fmt::format(
                              "Solver acceptance level for '{}', "
                              "default 2:\n\n.. value-table::",
                              GetConstraintName()).c_str(),
                            GetAccLevRef(), values_item_acceptance);
      } else
        if (ConstraintAcceptanceLevel::AcceptedButNotRecommended ==
                       GetChosenAcceptanceLevel()) {
          env.AddStoredOption(GetAcceptanceOptionNames(),
                              fmt::format(
                                "Solver acceptance level for '{}', "
                                "default 1:\n\n.. value-table::",
                                GetConstraintName()).c_str(),
                              GetAccLevRef(), values_item_acceptance);
        }
    }
    env.SetConstraintListHeader(
          "List of flat constraints.\n"
          "For each constraint the following are given:\n"
          "\n"
          "  - name,\n"
          "  - convertibility into simpler forms,\n"
          "  - solver acceptance natively,\n"
          "  - driver option(s) to modify acceptance\n"
          "    (enabled if both convertible and accepted).");
    std::string con_descr = (cancvt) ? "Convertible" : "NonConvertible";
    con_descr += "; ";
    if (ConstraintAcceptanceLevel::Recommended ==
          GetChosenAcceptanceLevel())
      con_descr += "NativeRecommended";
    else if (ConstraintAcceptanceLevel::AcceptedButNotRecommended ==
          GetChosenAcceptanceLevel())
      con_descr += "NativeAcceptedButNotRecommended";
    else
      con_descr += "NotAccepted";
    con_descr += "; ";
    con_descr += GetAcceptanceOptionNames();
    env.AddConstraintDescr(GetConstraintName(), con_descr);
  }

  /// Set user preferred acceptance level
  virtual void SetChosenAcceptanceLevel(
      ConstraintAcceptanceLevel acc) { acceptance_level_ = acc;}

	/// Mark as deleted, use index only
	virtual void MarkAsDeleted(int i) = 0;

  /// Copy names from ValueNodes
  virtual void CopyNamesFromValueNodes() = 0;

  /// Compute violations
  virtual void ComputeViolations(SolCheck& ) = 0;


protected:
  int& GetAccLevRef() { return acceptance_level_; }


private:
  pre::ValueNode value_node_;
  const char* const constr_name_;
  const char* const solver_opt_nm_;
  int acceptance_level_ {-1};
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
  /// when they overload Convert() etc, due to C++ name hiding
#define USE_BASE_CONSTRAINT_CONVERTERS(BaseConverter) \
  using BaseConverter::PreprocessConstraint; \
  using BaseConverter::PropagateResult; \
  using BaseConverter::IfHasCvt_impl; \
  using BaseConverter::IfNeedsCvt_impl; \
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


  /// Value of Pi
  static constexpr double Pi() { return 3.14159265358979; }

  /// Infinity
  static constexpr double Infty() { return INFINITY; }
  /// -Infinity
  static constexpr double MinusInfty() { return -INFINITY; }
  /// Pract inf
	static constexpr double PracticallyInf() { return 1e20; }
  /// Pract -inf
	static constexpr double PracticallyMinusInf() { return -1e20; }
};


/// Specialize ConstraintKeeper for a given constraint type
/// to store an array of such constraints
template <class Converter, class Backend, class Constraint>
class ConstraintKeeper : public BasicConstraintKeeper {
public:
  /// Constructor, adds this CK to the provided ConstraintManager
  /// Requires the CM to be already constructed
  ConstraintKeeper(Converter& cvt, const char* nm, const char* optnm) :
      BasicConstraintKeeper(cvt.GetValuePresolver(), nm, optnm), cvt_(cvt)
  {
    GetConverter().AddConstraintKeeper(*this, ConversionPriority());
  }

  /// Constraint type
  using ConstraintType = Constraint;

  /// Constrint Keeper description
  const std::string& GetDescription() const override
  { return desc_; }

  /// Assume Converter has the Backend
  Backend& GetBackend(BasicFlatConverter& cvt)
  { return static_cast<Converter&>(cvt).GetModelAPI(); }

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
    MP_UNUSED(cvt);
    try {
      return ConvertAllFrom(i_cvt_last_);
    } catch (const std::exception& exc) {
      MP_RAISE(Converter::GetTypeName() + std::string(": ")
                             + exc.what());
    }
    return false;
  }

  /// Converter's ability to convert the constraint type
  bool IfConverterConverts(
      BasicFlatConverter& cvt ) const override {
    return static_cast<Converter&>(cvt).
        IfHasConversion((const Constraint*)nullptr);
  }

  /// Acceptance level of this constraint type in the ModelAPI
  ConstraintAcceptanceLevel GetModelAPIAcceptance(
      const BasicFlatModelAPI& ba) const override {
    return
        static_cast<const Backend&>( ba ).
        AcceptanceLevel((Constraint*)nullptr);
  }

  /// Constraint type_info
  const std::type_info& GetTypeInfo() const override
  { return typeid(ConstraintType); }

  /// Report how many will be added to Backend
  int GetNumberOfAddable() const override {
    return (int)cons_.size()-n_bridged_;
  }

  /// Group number of this constraint type in the Backend.
  /// This is needed for pre- / postsolve to group solution values
  int GetConstraintGroup(const BasicFlatModelAPI& ba) const override {
    return static_cast<const Backend&>( ba ).
        GroupNumber((Constraint*)nullptr);
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

    bool IsDeleted() const { return IsBridged(); }
    bool IsBridged() const { return is_bridged_; }
    void MarkAsBridged() { is_bridged_=true; }

    /// Depth in the redef tree
    int GetDepth() const { return depth_; }

    Constraint con_;
    bool is_bridged_ = false;
    int depth_ {0};
  };

	/// Convert all new constraints of this type
  bool ConvertAllFrom(int& i_last) {
    int i=i_last;
    const auto acceptanceLevel =
        GetChosenAcceptanceLevel();
    if (NotAccepted == acceptanceLevel) {
      for ( ; ++i!=(int)cons_.size(); )
        if (!cons_[i].IsBridged())
          ConvertConstraint(cons_[i], i);
    }
    else if (AcceptedButNotRecommended == acceptanceLevel) {
      for (; ++i != (int)cons_.size(); ) {
        if (!cons_[i].IsBridged()) {
          try {       // Try to convert all but allow failure
            ConvertConstraint(cons_[i], i);
          } catch (const ConstraintConversionGracefulFailure& ) {
            /// nothing
          } catch (const ConstraintConversionFailure& ccf) {
            GetConverter().AddWarning( ccf.key(), ccf.message() );
          }
        }
      }
    } else { // Recommended == acceptanceLevel &&
      for (; ++i != (int)cons_.size(); )
        if (!cons_[i].IsBridged() &&
            GetConverter().IfNeedsConversion(cons_[i].con_, i))
          ConvertConstraint(cons_[i], i);
    }
    bool any_converted = i_last!=i-1;
    i_last = i-1;
    return any_converted;
  }

	/// Call Converter's RunConversion() and mark as "bridged".
  ///
	/// @param cnt the constraint container -
	/// actually redundant, as \a i is enough to find it. But for speed.
  /// @param i constraint index, needed for bridging
  void ConvertConstraint(Container& cnt, int i) {
    assert(!cnt.IsBridged());
    GetConverter().RunConversion(cnt.con_, i);
		MarkAsDeleted(cnt, i);
  }

	/// Mark item as deleted
	void MarkAsDeleted(Container& cnt, int ) {
		cnt.MarkAsBridged();
		++n_bridged_;
	}


public:
	/// Mark as deleted, use index only
	void MarkAsDeleted(int i) override {
		MarkAsDeleted(cons_.at(i), i);
	}

  /// Copy names from ValueNodes
  void CopyNamesFromValueNodes() override {
    const auto& vn = GetValueNode().GetStrVec();
    assert(vn.size()==cons_.size());
    for (auto i=vn.size(); i--; )
      cons_[i].con_.SetName(vn[i].MakeCurrentName());
  }

	/// ForEachActive().
  /// Deletes every constraint where fn() returns true.
	template <class Fn>
	void ForEachActive(Fn fn) {
		for (int i=0; i<(int)cons_.size(); ++i)
			if (!cons_[i].IsBridged())
				if (fn(cons_[i].con_, i))
					MarkAsDeleted(cons_[i], i);
	}

  /// Compute violations for this constraint type.
  /// We do it for redefined ones too.
  void ComputeViolations(SolCheck& chk) {
    if (cons_.size()) {
      auto& conviolmap =
          cons_.front().con_.IsLogical() ?
            chk.ConViolAlg() :
            chk.ConViolLog();
      auto& conviolarray =
          conviolmap[cons_.front().con_.GetTypeName()];
      const auto& x = chk.x();
      for (int i=(int)cons_.size(); i--; ) {
        auto viol = cons_[i].con_.ComputeViolation(x);
        if (viol > chk.GetFeasTol()) {
          /// Solver-side?
          /// TODO also original NL constraints (index 0)
          int index = cons_[i].IsDeleted() ? 2 : 1;
          conviolarray[index].CountViol(viol);
        }
      }
    }
  }


protected:
	/// Add all non-converted items to ModelAPI
  void AddAllUnbridged(BasicFlatModelAPI& be) {
    int con_index=0;
    auto con_group = GetConstraintGroup(be);
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
  std::deque<Container> cons_;
  int i_cvt_last_ = -1;               // Last converted constraint.
  int n_bridged_ = 0;                 // Number of converted items,
                                      // they won't go to Backend
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
/// without a subexpression map.
/// @param optionNames: name (or a few, space-separated)
/// of the solver option(s) for acceptance of this
/// constraint
#define STORE_CONSTRAINT_TYPE__INTERNAL( \
    Constraint, optionNames) \
private: \
  ConstraintKeeper<Impl, ModelAPI, Constraint> \
    CONSTRAINT_KEEPER_VAR(Constraint) \
      {*static_cast<Impl*>(this), \
       #Constraint, optionNames}; \
public: \
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
#define STORE_CONSTRAINT_TYPE__NO_MAP( \
    Constraint, optionNames) \
  STORE_CONSTRAINT_TYPE__INTERNAL( \
    Constraint, optionNames) \
  int MapFind__Impl(const Constraint& ) { return -1; } \
  bool MapInsert__Impl(const Constraint&, int ) \
    { return true; }


/// Define a constraint keeper
/// with a subexpression map.
/// The Converter storing the Constraint
/// should define MapFind / MapInsert accessing
/// the GET_(CONST_)CONSTRAINT_MAP(Constraint)
#define STORE_CONSTRAINT_TYPE__WITH_MAP( \
    Constraint, optionNames) \
  STORE_CONSTRAINT_TYPE__INTERNAL( \
    Constraint, optionNames) \
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

  /// This should be called after adding all constraint keepers
  void ConsiderAcceptanceOptions(
      BasicFlatConverter& cvt,
      const BasicFlatModelAPI& ma,
      Env& env) {
    for (auto& ck: con_keepers_)
      ck.second.ConsiderAcceptanceOptions(cvt, ma, env);
  }

  /// Convert all constraints (including any new appearing)
  void ConvertAllConstraints(BasicFlatConverter& cvt) {
    bool any_converted;
    do {
      any_converted = false;
      for (auto& ck: con_keepers_)
        any_converted = any_converted || ck.second.ConvertAllNewWith(cvt);
    } while (any_converted);
  }

  /// Fill counters of unbridged constraints
  void FillConstraintCounters(
      const BasicFlatModelAPI& mapi, FlatModelInfo& fmi) const {
    fmi.InitConstraintCount();
    for (const auto& ck: con_keepers_) {
      fmi.AddNumberOfConstraints(
        ck.second.GetTypeInfo(),
            ck.second.GetConstraintGroup(mapi),
            ck.second.GetNumberOfAddable());
    }
  }

  /// Copy names from ValueNodes
  void CopyNamesFromValueNodes() {
    for (const auto& ck: con_keepers_)
      ck.second.CopyNamesFromValueNodes();
  }

  /// Add all unbridged constraints to Backend
  void AddUnbridgedConstraintsToBackend(
      BasicFlatModelAPI& be) const {
    for (const auto& ck: con_keepers_)
      ck.second.AddUnbridgedToBackend(be);
  }

  /// Compute violations
  void ComputeViolations(SolCheck& chk) {
    for (const auto& ck: con_keepers_)
      ck.second.ComputeViolations(chk);
  }
};

} // namespace mp

#endif // CONSTRAINT_KEEPER_H
