#ifndef ITEM_KEEPER_H
#define ITEM_KEEPER_H

#include <cmath>

#include "mp/flat/preprocess.h"
#include "mp/flat/model_api_base.h"
#include "mp/utils-file.h"

namespace mp {

/// Converters handling custom constraints should derive from
class BasicFlatConverter;


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

  /// Is a logical constraint?
  virtual bool IsLogical() const = 0;

  /// Get context of contraint \a i
  virtual Context GetContext(int i) const = 0;

  /// Add context of contraint \a i
  virtual void AddContext(int i, Context ctx) = 0;

  /// Set context of contraint \a i
  virtual void SetContext(int i, Context ctx) = 0;

  /// Propagate expression result of constraint \a i bottom-up
  virtual void PreprocessConstraint(int i, PreprocessInfoStd& preinfo) = 0;

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

  /// Mark whether we could keep result vars
  virtual void MarkExprResultVars(BasicFlatConverter& cvt) = 0;

  /// Then, mark flat constraint arguments as vars
  virtual void MarkArguments(BasicFlatConverter& cvt) = 0;

  /// Convert to use expressions
  virtual void ConvertAllWithExpressions(BasicFlatConverter& cvt) = 0;

  /// acc:_expr==1 ?
  virtual bool IfWantNLOutput() const = 0;

  /// Low-level combined final acceptance
  int GetFinalItemAcceptance() const {
    assert(acc_level_default_>=0 && acc_level_default_<=4);
    int al_user = acc_level_item_; // user option
    if (-2 == al_user)
      al_user = 0;                // not accepted
    int al = acc_level_default_;
    if (al_user>=0)               // al_user provided - prioritized
      al = al_user;
    else if (AccLevelCommon()>=0) // acc:_all provided - otherwise
      al = AccLevelCommon();
    return al;
  }

  /// Query (user-chosen) acceptance level.
  /// This is "combined" for constraint or expression
  ConstraintAcceptanceLevel GetChosenAcceptanceLevel() const {
    if (acceptance_level_<0) {      // not initialized
      int al = GetFinalItemAcceptance();
      std::array<int, 5> alv = {0, 1, 2, 1, 2};
      acceptance_level_ = alv.at(al);
    }
    return ConstraintAcceptanceLevel(acceptance_level_);
  }

  /// Query (user-chosen) expression acceptance level.
  ExpressionAcceptanceLevel GetChosenAcceptanceLevelEXPR() const {
    if (acc_level_expr_<0) {      // not initialized
      int al = GetFinalItemAcceptance();
      std::array<int, 5> alv = {0, 0, 0, 1, 2};
      acc_level_expr_ = alv.at(al);
    }
    return ExpressionAcceptanceLevel(acc_level_expr_);
  }

  /// Converter's ability to convert the constraint type
  virtual bool IfConverterConverts(
      BasicFlatConverter& cvt ) const = 0;

  /// ModelAPI's acceptance level for the constraint type.
  /// This should not normally be used directly, instead:
  /// GetChosenAcceptanceLevel().
  /// Only when initializing options, their user values
  /// are not available, then use this.
  virtual ConstraintAcceptanceLevel GetModelAPIAcceptance(
      const BasicFlatModelAPI& ) const = 0;

  /// ModelAPI's acceptance level for the expression type.
  /// This should not normally be used directly, instead:
  /// GetChosenAcceptanceLevelEXPR()
  virtual ExpressionAcceptanceLevel GetModelAPIAcceptanceEXPR(
      const BasicFlatModelAPI& ) const = 0;

  /// Acceptance level of the overall expression interface in the ModelAPI
  virtual ExpressionAcceptanceLevel GetModelAPIAcceptance_EXPR_INTF(
      const BasicFlatModelAPI& ba) const = 0;

  /// @return acc:_all
  virtual int AccLevelCommon() const = 0;

  /// Constraint type_info
  virtual const std::type_info& GetTypeInfo() const =0;

  /// Backend's group number for the constraint type
  virtual int GetConstraintGroup(const BasicFlatModelAPI& ) const = 0;

  /// Report how many will be added to Backend
  virtual int GetNumberOfAddable() const = 0;

  /// This adds all unbridged items to the backend (without conversion)
  virtual void AddUnbridgedToBackend(
      BasicFlatModelAPI& be, const std::vector<std::string>* vnames) = 0;

  /// Store the solver's native expression for constraint \a i.
  /// Have to abandon type safety - an alternative would be to
  /// parameterize BasicConstraintKeeper by ConverterType
  /// and ModelAPI incl. ExprType.
  virtual void StoreSolverExpression(
      BasicFlatModelAPI& be, int i, void* pexpr) = 0;

  /// This logs the constraint group
  virtual void LogConstraintGroup(BasicFlatModelAPI& be) = 0;

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

  /// Constraint type name, e.g., 'AbsConstraint'
  const char* GetConstraintName() const { return constr_name_; }

  /// Expression type, or, if appropriate, constraint type name,
  /// e.g., 'Abs'
  virtual const char* GetExprOrConstraintName() const =0;

  /// Acceptance option names
  virtual const char* GetAcceptanceOptionNames() const
  { return solver_opt_nm_; }

  /// Constraint type short name.
  /// Ideally should be in the constraint itself,
  /// but currently we derive it from acceptance options.
  virtual const char* GetShortTypeName() const;

  /// See what options are available for this constraint:
  /// whether it is accepted natively by ModelAPI,
  /// as flat constraint or expression.
  /// Add acceptance option(s) "acc:...".
  /// Populate constraint list for -c output.
  /// @note This should be called before using the class.
  void ConsiderAcceptanceOptions(
      BasicFlatConverter& cvt,
      const BasicFlatModelAPI& ma,
      Env& env) {
    DoAddAcceptanceOptions(cvt, ma, env);
    DoPopulateConstraintList(cvt, ma, env);  // for -c option
  }

  /// Mark as bridged. Use index only.
  virtual void MarkAsBridged(int i) = 0;

  /// Mark as unused. Use index only.
  virtual void MarkAsUnused(int i) = 0;

  /// Is constraint \a i unused?
  virtual bool IsUnused(int i) const = 0;

  /// Copy names from ValueNodes
  virtual void CopyNamesFromValueNodes() = 0;

  /// Compute result for constraint \a i
  /// (for functional constraints)
  virtual double ComputeValue(int i, const VarInfoRecomp& ) = 0;

  /// Compute violations
  virtual void ComputeViolations(SolCheck& ) = 0;

  /// Set logger
  void SetLogger(BasicLogger* lg) { exporter_=lg; }
  /// Get logger, if provided and open and ok.
  BasicLogger* GetLogger() const {
    return exporter_ && exporter_->IsOpen()
               ? exporter_ : nullptr;
  }

protected:
  void DoAddAcceptanceOptions(
      BasicFlatConverter& cvt,
      const BasicFlatModelAPI& ma,
      Env& env);
  /// For -c
  void DoPopulateConstraintList(
      BasicFlatConverter& cvt,
      const BasicFlatModelAPI& ma,
      Env& env);


private:
  pre::ValueNode value_node_;
  const char* const constr_name_;
  const char* const solver_opt_nm_;
  mutable std::string type_name_short_;
  mutable int acceptance_level_ {-1};     // combined, for either con or expr
  int acc_level_item_ {-1};               // solver option acc:... value
  int acc_level_default_ {-1};            // default value, if neither item_ nor acc:_all set.
  mutable int acc_level_expr_ {-1};       // expression only
  BasicLogger* exporter_{};
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

  /// Store native expression for result index \a i.
  void StoreSolverExpression(
      BasicFlatModelAPI& be, void* pexpr) const
  { GetCK()->StoreSolverExpression(be, GetIndex(), pexpr); }

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
  using BaseConverter::IfDelayCvt_impl; \
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


////////////////////////////////////////////////////////////////////////////////////
/// Manage ConstraintKeepers for different constraint types
class ConstraintManager {
public:
  /// Add a new CKeeper with given conversion priority (smaller = sooner)
  void AddConstraintKeeper(BasicConstraintKeeper& ck, double priority) {
    con_keepers_.insert( { priority, ck } );
    ck.SetLogger(&*graph_exporter_app_);
  }

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

  /// Mark which func cons can be expressions
  void MarkExprResultVars(BasicFlatConverter& cvt) {
    for (auto& ck: con_keepers_)
      ck.second.MarkExprResultVars(cvt);
  }

  /// Then, mark arguments of flat cons as vars
  void MarkArguments(BasicFlatConverter& cvt) {
    for (auto& ck: con_keepers_)
      ck.second.MarkArguments(cvt);
  }

  /// Convert to expression-based model
  void ConvertAllWithExpressions(BasicFlatConverter& cvt) {
    for (auto& ck: con_keepers_)
      ck.second.ConvertAllWithExpressions(cvt);
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
      BasicFlatModelAPI& be,
      const std::vector<std::string>* pvnam=nullptr) const {
    for (const auto& ck: con_keepers_)
      ck.second.AddUnbridgedToBackend(be, pvnam);
  }

  /// Log constraint groups
  void LogConstraintGroups(
      BasicFlatModelAPI& be) const {
    for (const auto& ck: con_keepers_)
      ck.second.LogConstraintGroup(be);
  }

  /// Compute violations
  void ComputeViolations(SolCheck& chk) {
    for (const auto& ck: con_keepers_)
      ck.second.ComputeViolations(chk);
  }

  /// Retrieve file logger
  BasicFileAppender& GetFileAppender() const
  { return *graph_exporter_app_; }

private:
  std::multimap<double, BasicConstraintKeeper&> con_keepers_;
  /// Conversion graph exporter file appender
  std::unique_ptr<BasicFileAppender>
      graph_exporter_app_{MakeFileAppender()};
};

}  // namespace mp

#endif // ITEM_KEEPER_H
