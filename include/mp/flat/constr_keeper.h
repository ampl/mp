#ifndef CONSTRAINT_KEEPER_H
#define CONSTRAINT_KEEPER_H

#include <deque>
#include <unordered_map>
#include <functional>
#include <cmath>
#include <climits>

#include "mp/common.h"
#include "mp/format.h"
#include "mp/env.h"
#include "mp/utils-math.h"
#include "mp/util-json-write.hpp"

#include "mp/flat/item_keeper.h"
#include "mp/flat/constr_hash.h"
#include "mp/flat/redef/redef_base.h"
#include "mp/flat/constr_eval.h"

namespace mp {

/// Specialize ConstraintKeeper for a given constraint type
/// to store an array of such constraints
template <class Converter, class Backend, class Constraint>
class ConstraintKeeper final
    : public BasicConstraintKeeper {
public:
  /// Constructor, adds this CK to the provided ConstraintManager
  /// Requires the CM to be already constructed
  ConstraintKeeper(Converter& cvt, const char* nm, const char* optnm) :
      BasicConstraintKeeper(cvt.GetValuePresolver(), nm, optnm), cvt_(cvt)
  {
    GetValueNode().SetName(GetShortTypeName());  // change value node name
    GetConverter().AddConstraintKeeper(*this, ConversionPriority());
  }

  /// Constraint type
  using ConstraintType = Constraint;

  /// The corresponding flat expression type
  using FlatExprType = ExprWrapper<Constraint>;

  /// Expression type, or, if appropriate, constraint type name,
  /// e.g., 'Abs'
  const char* GetExprOrConstraintName() const override
  { return Constraint::GetTypeName(); }

  /// Constrint Keeper description
  const std::string& GetDescription() const override
  { return desc_; }

  /// Is a logical constraint?
  bool IsLogical() const override
  { return Constraint::IsLogical(); }

  /// Assume Converter has the Backend
  Backend& GetBackend(BasicFlatConverter& cvt)
  { return static_cast<Converter&>(cvt).GetModelAPI(); }

  /// Add a pre-constructed constraint (or just arguments)
  /// @return index of the new constraint
  template <class... Args>
  int AddConstraint(int d, Args&&... args)
  {
    cons_.emplace_back( d, std::move(args)... );
    ++n_bridged_or_unused_;    // initially "unused"
    ExportConstraint(cons_.size()-1, cons_.back());
    // fmt::MemoryWriter wrt;
    // WriteCon2JSON(wrt, cons_.size()-1, cons_.back());
    // printf("%s\n", wrt.c_str());
    return cons_.size()-1;
  }

  /// Get const constraint \a i
  const Constraint& GetConstraint(int i) const
  { assert(check_index(i)); return cons_[i].GetCon(); }

  /// Get constraint \a i
  Constraint& GetConstraint(int i)
  { assert(check_index(i)); return cons_[i].GetCon(); }

  /// Is item \a i already bridged or abandoned?
  bool IsRedundant(int i) const override {
    return cons_[i].IsRedundant();
  }

  /// Get constraint depth in the reformulation tree
  int GetConstraintDepth(int i) const
  { assert(check_index(i)); return cons_[i].GetDepth(); }

  /// Get context of constraint \a i
  Context GetContext(int i) const override
  { assert(check_index(i)); return cons_[i].GetCon().GetContext(); }

  /// Add context of constraint \a i
  void AddContext(int i, Context ctx) override
  { assert(check_index(i)); cons_[i].GetCon().AddContext(ctx); }

  /// Set context of constraint \a i
  void SetContext(int i, Context ctx) override
  { assert(check_index(i)); cons_[i].GetCon().SetContext(ctx); }

  /// Propagate expression result of constraint \a i bottom-up
  void PreprocessConstraint(int i, PreprocessInfoStd& preinfo) override {
    if constexpr (std::is_base_of_v<FunctionalConstraint, Constraint>) {
      PreprocessInfo<Constraint> prepro;
      GetConverter().PreprocessConstraint(cons_[i].GetCon(), prepro);
      preinfo.narrow_result_bounds(prepro.lb(), prepro.ub());
      preinfo.set_result_type(prepro.get_result_type());
    }
  }

  /// Propagate expression result of constraint \a i top-down
  void PropagateResult(BasicFlatConverter& cvt,
                       int i,
                       double lb, double ub, Context ctx) override {
    try {
      // Too strong: instead, differentiate context
      // in which the redefinition happened #248.
      // MP_ASSERT_ALWAYS(ctx.IsSubsetOf(GetContext())
      //                  || IsBridgingToBeConsidered(i),
      auto ctx_redef = cons_[i].GetRedefContext();
      MP_ASSERT_ALWAYS(ctx_redef.IsNone()
                           || ctx.IsSubsetOf(ctx_redef),
                       "Part of expression redefinition\n"
                       "could be lost. Please contact\n"
                       "AMPL customer support.");
      static_cast<Converter&>(cvt).PropagateResult(
            GetConstraint(i), lb, ub, ctx);
    } catch (const std::exception& exc) {
      MP_RAISE(Converter::GetTypeName() +
               std::string(": propagating result for constraint ") +
               std::to_string(i) + "\nof type '" +
               Constraint::GetTypeName() +
               "':\n" + exc.what());
    }
  }

  /// Result variable of constraint \a i. Returns -1 if none
  int GetResultVar(int i) const override
  { assert(check_index(i)); return cons_[i].GetCon().GetResultVar(); }

  /// Conversion priority. Uses that from Converter
  double ConversionPriority() const
  { return Converter::ConstraintCvtPriority((Constraint*)nullptr); }

  /// Whether we consider conversion of this type.
  /// It can still not happen if not convertible, or asked to delay
  /// individual items.
  bool IfConsiderConversion() const {
    return
        (ConstraintAcceptanceLevel::NotAccepted == GetChosenAcceptanceLevel()
         || (ConstraintAcceptanceLevel::NotAccepted == GetModelAPIAcceptance()
             && (!GetConverter().IfWantNLOutput()
                 || ExpressionAcceptanceLevel::NotAccepted
                        == GetChosenAcceptanceLevelEXPR())
                && 2 != AccLevelCommon()));             // Not when acc:_all=2
  }

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

  /// Mark whether we could use result vars of functional constraints
  /// as expressions
  void MarkExprResultVars(BasicFlatConverter& cvt) override {
    assert(&cvt == &GetConverter());         // Using the same Converter
    DoMarkForResultVars();
  }

  /// Then, mark arguments of flat constraints as proper vars
  void MarkArguments(BasicFlatConverter& cvt) override {
    assert(&cvt == &GetConverter());         // Using the same Converter
    DoMarkForArguments();
  }

  /// Convert to use expressions
  void ConvertAllWithExpressions(BasicFlatConverter& cvt) override {
    assert(&cvt == &GetConverter());         // Using the same Converter
    DoCvtWithExprs();
  }

  /// acc:_expr==1 ?
  bool IfWantNLOutput() const override
  { return GetConverter().IfWantNLOutput(); }

  /// Converter's ability to convert the constraint type
  bool IfConverterConverts(
      BasicFlatConverter& cvt ) const override {
    return static_cast<Converter&>(cvt).
        IfHasConversion((const Constraint*)nullptr);
  }

  /// Acceptance level of this constraint type in the ModelAPI
  ConstraintAcceptanceLevel GetModelAPIAcceptance(
      const BasicFlatModelAPI&  = *(BasicFlatModelAPI*)nullptr)
      const override {
    return
        // static_cast<const Backend&>( ba ).
        Backend::
        AcceptanceLevel((Constraint*)nullptr);
  }

  /// Acceptance level of the corresponding expression type in the ModelAPI
  ExpressionAcceptanceLevel GetModelAPIAcceptanceEXPR(
      const BasicFlatModelAPI&  = *(BasicFlatModelAPI*)nullptr)
      const override {
    return
        // static_cast<const Backend&>( ba ).
        Backend::
        AcceptanceLevel((FlatExprType*)nullptr);
  }

  /// Acceptance level of the overall expression interface in the ModelAPI
  ExpressionAcceptanceLevel GetModelAPIAcceptance_EXPR_INTF(
      const BasicFlatModelAPI& ba) const override {
    return
        static_cast<const Backend&>( ba ).
        ExpressionInterfaceAcceptanceLevel();
  }

  /// @return acc:_all
  int AccLevelCommon() const override
  { return GetConverter().AcceptanceLevelCommon(); }

  /// Constraint type_info
  const std::type_info& GetTypeInfo() const override
  { return typeid(ConstraintType); }

  /// Report total number
  int Size() const {
    return (int)cons_.size();
  }

  /// Report how many will be added to Backend
  int GetNumberOfAddable() const override {
    // printf("   N ADDABLE '%s': %d - %d = %d\n", GetShortTypeName(),
    //       Size(), n_bridged_or_unused_, Size()-n_bridged_or_unused_);
    assert(Size() >= n_bridged_or_unused_);
    assert(0 <= n_bridged_or_unused_);
    return Size()-n_bridged_or_unused_;
  }

  /// Group number of this constraint type in the Backend.
  /// This is needed for pre- / postsolve to group solution values
  int GetConstraintGroup(const BasicFlatModelAPI& ba) const override {
    return static_cast<const Backend&>( ba ).
        GroupNumber((Constraint*)nullptr);
  }

  /// Add remaining constraints to Backend
  void AddUnbridgedToBackend(
      BasicFlatModelAPI& be,
      ItemNamer& vnam) override {
    if (ExpressionAcceptanceLevel::NotAccepted
        == GetChosenAcceptanceLevelEXPR()
        || !GetConverter().IfWantNLOutput()) {
      try {
        AddAllUnbridged(be, vnam);
      } catch (const std::exception& exc) {
        MP_RAISE(std::string("Adding constraint of type '") +
                 Constraint::GetTypeName() + "' to " +
                 Backend::GetTypeName() + std::string(": ") +
                 exc.what());
      }
    }
  }

  /// Store the solver's native expression for constraint \a i.
  /// Have to abandon type safety - an alternative would be to
  /// parameterize BasicConstraintKeeper by ConverterType
  /// and ModelAPI incl. ExprType.
  void StoreSolverExpression(
      BasicFlatModelAPI& be, int i, void* pexpr) override {
    if constexpr (ExpressionAcceptanceLevel::NotAccepted
        != Backend::ExpressionInterfaceAcceptanceLevel()) {
      assert(!cons_[i].IsRedundant());
      assert(!cons_[i].IsExprAdded());
      cons_[i].MarkExprAdded();
      *(typename Backend::Expr*)pexpr =
          static_cast<Backend&>(be).AddExpression(cons_[i].GetExpr());
    }
  }

  /// Log constraint group
  void LogConstraintGroup(
      BasicFlatModelAPI& be) override {
    auto cg = GetConstraintGroup(be);
    if (cg>=0 && GetLogger()) {
      fmt::MemoryWriter wrt;
      {
        MiniJSONWriter jw(wrt);
        jw["CON_TYPE"] = GetShortTypeName();
        jw["CON_GROUP"] = ConGroupName(cg);
        jw["CON_GROUP_index"] = cg;
      }
      wrt.write("\n");                     // EOL
      GetLogger()->Append(wrt);
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
  class Container {
  public:
    Container(int d, Constraint&& c) noexcept
      : con_(std::move(c)), depth_(d) { }

    /// Depth in redef tree
    int GetDepth() const { return depth_; }

    /// Is redundant: bridged or unused
    bool IsRedundant() const
    { return IsBridged() || IsUnused(); }

    /// Bridged (reformulated.)
    bool IsBridged() const { return is_bridged_; }
    /// Mark as bridged
    void MarkAsBridged() { is_bridged_=true; }

    /// Unused, e.g., inlined in another expression.
    bool IsUnused() const { return is_unused_; }
    /// Mark as unused
    void MarkAsUnused() { is_unused_=true; }
    /// Mark as used
    void MarkAsUsed() { is_unused_=false; }

    /// Has the expression been added to the backend?
    bool IsExprAdded() const { return is_expr_stored_; }
    /// Mark as added
    void MarkExprAdded() { is_expr_stored_=true; }

    /// Context used for redefinition.
    /// We redefine just once currently #248.
    Context GetRedefContext() const { return ctx_redef_; }
    /// Set redef context
    void SetRedefContext(Context ctx) {
      assert(!ctx.IsNone());
      ctx_redef_ = ctx;
    }

    /// Get the flat constraint, const &
    const Constraint& GetCon() const { return con_.GetFlatConstraint(); }
    /// Get the flat constraint &
    Constraint& GetCon() { return con_.GetFlatConstraint(); }

    /// Get the expression, const &
    const FlatExprType& GetExpr() const { return con_; }
    /// Get the expression &
    FlatExprType& GetExpr() { return con_; }

  private:
    // Storing in the ExprWrapper,
    // so we can send (wrapper &) to ModelAPI::AddExpression().
    FlatExprType con_;
    int depth_ = 0;
    Context ctx_redef_;    // Context used for redefinition, if any
    char is_bridged_ = false;
    char is_unused_ = true;
    char is_expr_stored_ = false;
  };

	/// Convert all new constraints of this type
  bool ConvertAllFrom(int& i_last) {
    int i=i_last;
    const auto acceptanceLevel =
        GetChosenAcceptanceLevel();
    if (IfConsiderConversion()) {
      if (!IfConverterConverts(GetConverter())) {
        i = (int)cons_.size();             // skip unconverted items
      } else {
        for ( ; ++i!=(int)cons_.size(); )
          if (!cons_[i].IsRedundant() &&
              !GetConverter().IfDelayConversion(cons_[i].GetCon(), i)) {
            ConvertConstraint(cons_[i], i);
            if (-2 == GetLowLevelAcc()) {
              GetConverter().GetEnv().warn_from_cb(GetShortTypeName());
            }
          }
      }
    }
    else if (ConstraintAcceptanceLevel::AcceptedButNotRecommended == acceptanceLevel) {
      if (!IfConverterConverts(GetConverter())) {
        i = (int)cons_.size();
      } else {
        for (; ++i != (int)cons_.size(); ) {
          if (!cons_[i].IsRedundant() &&
              !GetConverter().IfDelayConversion(cons_[i].GetCon(), i)) {
            try {       // Try to convert all but allow failure
              ConvertConstraint(cons_[i], i);
            } catch (const ConstraintConversionGracefulFailure& ) {
              /// nothing
            } catch (const ConstraintConversionFailure& ccf) {
              GetConverter().AddWarning( ccf.key(), ccf.message() );
            }
          }
        }
      }
    } else { // Recommended == acceptanceLevel &&
      for (; ++i != (int)cons_.size(); )
        if (!cons_[i].IsRedundant() &&
            GetConverter().IfNeedsConversion(cons_[i].GetCon(), i))
          ConvertConstraint(cons_[i], i);
    }
    bool any_converted = i_last!=i-1;
    i_last = i-1;
    return any_converted;
  }

  /// Mark final func cons as expressions if accepted so
  void DoMarkForResultVars() {
    const auto eal        // as expr only
        = GetChosenAcceptanceLevelEXPR();
    if (ExpressionAcceptanceLevel::NotAccepted!=eal) {    // accepted
      for (int i=0; i< (int)cons_.size(); ++i) {
        const auto& cnt = cons_[i];
        if (!cnt.IsRedundant()) {
          const auto& con = cnt.GetCon();
          // Delegate actual logic to Converter
          GetConverter().ConsiderMarkingResultVar(con, i, eal);
        }
      }
    }
  }

  void DoMarkForArguments() {
    const auto eal        // expr only
        = GetChosenAcceptanceLevelEXPR();
    for (int i=0; i< (int)cons_.size(); ++i) {
      const auto& cnt = cons_[i];
      if (!cnt.IsRedundant()) {      // Delegate actual logic to Converter
        const auto& con = cnt.GetCon();
        GetConverter().ConsiderMarkingArguments(con, i, eal);
      }
    }
  }

  void DoCvtWithExprs() {
    // For acc:_all=0 to work, we have to accept non-convertible con types
    auto cal = IfConverterConverts(GetConverter())
                   ? GetChosenAcceptanceLevel()
                   : ConstraintAcceptanceLevel::Recommended;
    auto eal = GetChosenAcceptanceLevelEXPR();
    ForEachActive(
        [this, cal, eal](const auto& con, int i) {
          return this->GetConverter().ConvertWithExpressions(con, i, cal, eal);
        });
  }

public:
  /// Call Converter's RunConversion() and mark as "bridged".
  /// @param i constraint index
  void ConvertConstraint(int i) {
    ConvertConstraint(cons_[i], i);
  }

protected:
  /// Call Converter's RunConversion() and mark as "bridged".
  ///
	/// @param cnt the constraint container -
  ///   actually redundant, as \a i is enough to find it. But for speed.
  /// @param i constraint index, needed for bridging
  void ConvertConstraint(Container& cnt, int i) {
    assert(!cnt.IsRedundant());
    assert(cnt.GetRedefContext().IsNone());
    auto ctx_redef
        = GetConverter().RunConversion(cnt.GetCon(), i, cnt.GetDepth());
    cnt.SetRedefContext(ctx_redef);
    assert(!cnt.GetCon().HasResultVar()
           || cnt.GetCon().GetContext().IsSubsetOf(ctx_redef));
    MarkAsBridged(cnt, i);
  }

  /// Mark item as reformulated
  void MarkAsBridged(Container& cnt, int i) {
    if (!cnt.IsBridged()) {  // can be called 2x,
      printf("MarkBRIDGED: %s [%d], resvar=%d\n",
             GetShortTypeName(), i, cnt.GetCon().GetResultVar());
      cnt.MarkAsBridged();   // e.g. IfThen: 1st by RedefineVariable(),
      GetConverter().UncountArgRefs(cnt.GetCon());
      ++n_bridged_or_unused_;  // then in ConvertConstraint()
    }
  }

  /// Mark item as unused
  void MarkAsUnused(Container& cnt, int i, bool recurs=true) {
    printf("MarkUNUSED: %s [%d], resvar=%d\n",
           GetShortTypeName(), i, cnt.GetCon().GetResultVar());
    assert(!cnt.IsUnused());
    if (!cnt.IsUnused()) {        // in Release
      cnt.MarkAsUnused();
      if (recurs)
        GetConverter().UncountArgRefs(cnt.GetCon());
      ++n_bridged_or_unused_;
    }
  }

  /// Mark item as used
  void MarkAsUsed(Container& cnt, int i) {
    printf("MarkUSED: %s [%d], resvar=%d\n",
           GetShortTypeName(), i, cnt.GetCon().GetResultVar());
    assert(cnt.IsUnused());
    if (cnt.IsUnused()) {
      cnt.MarkAsUsed();
      GetConverter().CountArgRefs(cnt.GetCon());
      assert(!cnt.IsBridged());
      --n_bridged_or_unused_;
    }
  }

protected:
  /// Export (last added) constraint
  void ExportConstraint(int i_con, const Container& cnt) {
    if (GetLogger()) {
      fmt::MemoryWriter wrt;
      WriteCon2JSON(wrt, i_con, cnt);
      wrt.write("\n");                     // EOL
      GetLogger()->Append(wrt);
    }
  }

  /// Write constraint in JSON format
  void WriteCon2JSON(
      fmt::MemoryWriter& wrt, int i_con, const Container& cnt) {
    MiniJSONWriter jw(wrt);
    jw["CON_TYPE"] = GetShortTypeName();
    jw["index"] = i_con;
    if (*cnt.GetCon().name())
      jw["name"] = cnt.GetCon().name();
    jw["depth"] = cnt.GetDepth();
    WriteJSON(jw["data"], cnt.GetCon());
  }

  /// Export constraint status.
  /// This is called in the end,
  /// so printing the readable form.
  void ExportConStatus(int i_con, const Container& cnt,
                       ItemNamer& vnam, bool add2final) {
    if (GetLogger()) {
      fmt::MemoryWriter wrt;
      {
        MiniJSONWriter jw(wrt);
        jw["CON_TYPE"] = GetShortTypeName();
        jw["index"] = i_con;
        if (*cnt.GetCon().name()) {
          jw["name"] = cnt.GetCon().name();
          {
            fmt::MemoryWriter pr;
            WriteFlatCon(pr, cnt.GetCon(), vnam);
            jw["printed"] = pr.c_str();
          }
        }
        jw["depth"] = cnt.GetDepth();
        jw["unused"] = (int)cnt.IsUnused();
        jw["bridged"] = (int)cnt.IsBridged();
        jw["final"] = (int)add2final;
      }
      wrt.write("\n");                     // EOL
      GetLogger()->Append(wrt);
    }
  }

public:
  /// Mark cons[\a i] as reformulated.
  /// Use index only.
  void MarkAsBridged(int i) override {
    MarkAsBridged(cons_.at(i), i);
	}

  /// Mark cons[\a i] as unused.
  /// Use index only.
  void MarkAsUnused(int i) override {
    MarkAsUnused(cons_.at(i), i);
  }

  /// Mark cons[\a i] as unused.
  /// Do not propagate to arguments.
  /// Use index only.
  void MarkAsUnused_ThisOnly(int i) override {
    MarkAsUnused(cons_.at(i), i, false);
  }

  /// Mark cons[\a i] as used.
  /// Use index only.
  void MarkAsUsed(int i) override {
    MarkAsUsed(cons_.at(i), i);
  }

  /// Is constraint \a i reformulated?
  bool IsBridged(int i) const override {
    return cons_.at(i).IsBridged();
  }

  /// Is bridging of constraint \a i
  /// to be considered yet?
  /// @todo repeated redefinition cycle
  bool IsBridgingToBeConsidered(int i) const override {
    return i > i_cvt_last_;
  }

  /// Is constraint \a i unused?
  bool IsUnused(int i) const override {
    return cons_.at(i).IsUnused();
  }

  /// Copy names from ValueNodes
  void CopyNamesFromValueNodes() override {
    const auto& vn = GetValueNode().GetStrVec();
    assert(vn.size()==cons_.size());
    for (auto i=vn.size(); i--; )
      cons_[i].GetCon().SetName(vn[i].MakeCurrentName());
  }

  /// Copy names to ValueNodes
  void CopyNames2ValueNodes() {
    auto& vn = GetValueNode().GetStrVec();
    assert(vn.size()==cons_.size());
    for (auto i=vn.size(); i--; )
      vn[i] = std::string(cons_[i].GetCon().name());
  }

  /// ForEachActive().
  /// Deletes every constraint where fn() returns true.
  /// @return number of deleted / decreased-usage items.
	template <class Fn>
  int ForEachActive(
      Fn fn, int i_start = 0, int i_end = INT_MAX) {
    int ndel = 0;
    for (int i=i_start;
         i<std::min((int)cons_.size(), i_end);
         ++i)
      if (!cons_[i].IsRedundant())
        if (fn(cons_[i].GetCon(), i)) {
          MarkAsBridged(cons_[i], i);
          ++ndel;
        }
    return ndel;
	}

  /// Compute result for constraint \a i
  /// (for functional constraints).
  double
  ComputeValue(int i, const VarInfoRecomp& vir) override {
    assert(cons_[i].GetCon().GetResultVar() >= 0);
    return mp::ComputeValue(cons_[i].GetCon(), vir);
  }

  /// Compute violations for this constraint type.
  /// We do it for redefined (intermediate) ones too.
  void ComputeViolations(SolCheck& chk) override {
    if (cons_.size()) {
      auto& conviolmap =
          cons_.front().GetCon().IsLogical() ?
            chk.ConViolLog() :
            chk.ConViolAlg();
      const auto& x = chk.x_ext();
      ViolSummArray<3>* conviolarray {nullptr};
      for (int i=(int)cons_.size(); i--; ) {
        if (!cons_[i].IsUnused()) {
          int c_class = 0;    // class of this constraint
          if (!cons_[i].IsBridged())
            c_class |= 8;     // solver-side constraint
          if (!cons_[i].GetDepth())
            c_class |= 2;     // top-level
          if (!c_class)
            c_class = 4;      // intermediate
          if (c_class & chk.check_mode()) {
            auto viol = cons_[i].GetCon().ComputeViolation(x);
            auto cr = viol.Check(
                  chk.GetFeasTol(), chk.GetFeasTolRel());
            if (cr.first) {
              if (!conviolarray)
                conviolarray =         // lazy map access
                    &conviolmap[GetShortTypeName()];
              /// index==0,1,2: original, interm, solver-side
              /// If both orig and solver, report as orig
              int index = (c_class & 2) ? 0
                                        : (c_class & 8)
                                          ? 2 : 1;
              assert(index < (int)conviolarray->size());
              (*conviolarray)[index].CountViol(
                    viol, cr.second, cons_[i].GetCon().name());
            }
          }
        }
      }
    }
  }


protected:
  /// Add all non-converted items to ModelAPI.
  /// Export all constraints if desired.
  void AddAllUnbridged(BasicFlatModelAPI& be,
                       ItemNamer& vnam) {
    auto con_group = GetConstraintGroup(be);
		for ( ; i_2add_next_ < (int)cons_.size(); ++i_2add_next_) {
      const auto& cont = cons_[i_2add_next_];
      bool adding = !cont.IsRedundant();            // includes 'unused'
      if (adding) {
        static_cast<Backend&>(be).AddConstraint(cont.GetCon());
        GetConverter().GetCopyLink().             // Linking to the "final" nodes
            AddEntry({
                       GetValueNode().Select(i_2add_next_),
                       GetConverter().GetValuePresolver().GetTargetNodes().
                         GetConValues()(con_group).Add()
                     });
      }
      ExportConStatus(i_2add_next_, cont, vnam, adding);
    }
  }


private:
  Converter& cvt_;
  std::deque<Container> cons_;
  int i_cvt_last_ = -1;               // Last converted constraint.
  int n_bridged_or_unused_ = 0;       // Number of converted items,
                                      // they won't go to Backend.
  int i_2add_next_ = 0;               // Next constraint to consider
                                      // for adding to Backend.
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
  GetConstraintKeeper(const Constraint* ) const { \
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
    Constraint, optionNames, prior) \
  STORE_CONSTRAINT_TYPE__INTERNAL( \
    Constraint, optionNames) \
  int MapFind__Impl(const Constraint& ) { return -1; } \
  bool MapInsert__Impl(const Constraint&, int ) \
    { return true; } \
  DEFINE_CONSTRAINT_PRIORITY(Constraint, prior)


/// Define a constraint keeper
/// with a subexpression map.
/// The Converter storing the Constraint
/// should define MapFind / MapInsert accessing
/// the GET_(CONST_)CONSTRAINT_MAP(Constraint)
#define STORE_CONSTRAINT_TYPE__WITH_MAP( \
    Constraint, optionNames, prior) \
  STORE_CONSTRAINT_TYPE__INTERNAL( \
    Constraint, optionNames) \
  STORE_CONSTRAINT_MAP(Constraint) \
  DEFINE_CONSTRAINT_PRIORITY(Constraint, prior)

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

/// Define constraint priority
#define DEFINE_CONSTRAINT_PRIORITY(Constraint, prior) \
  static constexpr double ConstraintCvtPriority( \
    const Constraint*) { return prior; }


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

} // namespace mp

#endif // CONSTRAINT_KEEPER_H
