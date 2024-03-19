#ifndef CONSTRAINT_KEEPER_H
#define CONSTRAINT_KEEPER_H

#include <deque>
#include <unordered_map>
#include <functional>
#include <cmath>

#include "mp/common.h"
#include "mp/format.h"
#include "mp/env.h"
#include "mp/utils-math.h"
#include "mp/utils-file.h"
#include "mp/util-json-write.hpp"

#include "mp/flat/model_api_base.h"
#include "mp/flat/constr_hash.h"
#include "mp/flat/redef/redef_base.h"
#include "mp/valcvt-node.h"
#include "mp/flat/constr_eval.h"

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
  /// Check if this violation should be counted.
  void CheckViol(
      Violation viol,
      double epsabs, double epsrel,
      const char* nm) {
    auto chk = viol.Check(epsabs, epsrel);
    if (chk.first)
      CountViol(viol, chk.second, nm);
  }
  /// Count violation
  void CountViol(
      Violation viol, double violRel, const char* nm) {
    ++N_;
    if (epsAbsMax_ < viol.viol_) {
      epsAbsMax_ = viol.viol_;
      nameAbs_ = nm;
    }
    if (epsRelMax_ < violRel) {
      epsRelMax_ = violRel;
      nameRel_ = nm;
    }
  }
  int N_ {0};
  double epsAbsMax_ {0.0};
  const char* nameAbs_ {nullptr};
  double epsRelMax_ {0.0};
  const char* nameRel_ {nullptr};
};


/// Array of violation summaries.
/// For different kinds, e.g., original / aux vars.
template <int Nkinds>
using ViolSummArray = std::array<ViolSummary, Nkinds>;

class VarInfoRecomp;

/// Function prototype to recompute
/// variable at index \a i.
using VarsRecomputeFn
    = std::function<
        double(int i, const VarInfoRecomp& x) >;

/// Var vector managing recomputation
class VarVecRecomp {
public:
  /// Construct
  VarVecRecomp(std::vector<double> x,
               VarsRecomputeFn rec_fn)
    : x_(std::move(x)), is_recomp_(x_.size()),
      recomp_fn_(rec_fn) { assert(recomp_fn_); }
  /// Set p_var_info_recomp_
  void set_p_var_info(const VarInfoRecomp* p) const
  { p_var_info_recomp_ = p; }
  /// Number of vars
  int size() const { return (int)x_.size(); }
  /// Access variable value.
  /// Recompute if not yet.
  double operator[]( int i ) const {
    assert(i>=0 && i<(int)x_.size());
    assert(p_var_info_recomp_ != nullptr);
    if (!is_recomp_[i]) {
      x_[i] = recomp_fn_(i, *p_var_info_recomp_);
      is_recomp_[i] = true;
    }
    return x_[i];
  }
  /// Expose begin()
  std::vector<double>::iterator begin() { return x_.begin(); }
  /// Expose end()
  std::vector<double>::iterator end() { return x_.end(); }
  /// Move out x
  std::vector<double>& get_x() const { return x_; }

private:
  mutable std::vector<double> x_;
  mutable std::vector<bool> is_recomp_;
  VarsRecomputeFn recomp_fn_;
  mutable const VarInfoRecomp* p_var_info_recomp_{nullptr};
};


/// Static var vector
using VarVecStatic = std::vector<double>;


/// Variable information used by solution check
template <class VarVec>
class VarInfoImpl {
public:
  /// Constructor
  VarInfoImpl(double ft, bool recomp_vals,
          VarVec x,
          ArrayRef<double> x_raw,
          ArrayRef<var::Type> type,
          ArrayRef<double> lb, ArrayRef<double> ub,
          int sol_rnd, int sol_prec)
    : feastol_(ft), recomp_vals_(recomp_vals),
      x_(std::move(x)), x_raw_(x_raw),
      type_(type), lb_(lb), ub_(ub) {
    assert((int)x_.size()>=(int)type_.size());  // feasrelax can add more
    assert(type_.size()==lb_.size());
    assert(type_.size()==ub_.size());
    apply_precision_options(sol_rnd, sol_prec); // after recomp?
  }
  /// Number of vars
  int size() const { return (int)x_.size(); }
  /// Access variable value
  double operator[]( int i ) const {
    assert(i>=0 && i<(int)x_.size());
    return x_[i];
  }
  /// Access VarVec
  const VarVec& get_x() const { return x_; }
  /// Access raw variables
  double raw(int i) const {
    assert(i < (int)x_raw_.size()
           && "Can only access raw solver values "
              "in idealistic mode and they should be available");
    return x_raw_[i];
  }

  /// Access integrality condition
  bool is_var_int(int i) const {
    assert(i>=0 && i<(int)type_.size());
    return var::INTEGER==type_[i];
  }
  /// Variable value nonzero?
  bool is_nonzero(int i) const {
    return
        std::fabs( (*this)[i] )
        >= (is_var_int(i) ? 0.5 : feastol());
  }
  /// Variable value positive?
  bool is_positive(int i) const {
    return
        (*this)[i]
        >= (is_var_int(i) ? 0.5 : feastol());
  }
  /// Is at lb?
  bool is_at_lb(int i) const
  { return (*this)[i] - lb_[i] <= feastol(); }
  /// Is at ub?
  bool is_at_ub(int i) const
  { return -((*this)[i] - ub_[i]) <= feastol(); }

  /// Bounds violation
  double bounds_viol(int i) const {
    assert(i>=0 && i<(int)type_.size());
    return std::max(lb_[i] - x_[i], x_[i] - ub_[i]);
  }

  /// Feasibility tolerance
  double feastol() const { return feastol_; }
  /// Using recomputed auxiliary vars?
  bool recomp_vals() const { return recomp_vals_; }
  /// Using idealistic checking of solution
  /// (without tolerances)?
  bool idealistic() const { return recomp_vals(); }
  /// sol_rnd as string
  std::string solution_round() const
  { return sol_rnd_ < 100 ? std::to_string(sol_rnd_) : ""; }
  /// sol_rnd as string
  std::string solution_precision() const
  { return sol_prec_ < 100 ? std::to_string(sol_prec_) : ""; }


protected:
  void apply_precision_options(
      int sol_rnd, int sol_prec) {
    try {                 // Apply sol_rnd
      if (sol_rnd<100) {
        sol_rnd_ = (sol_rnd);
        auto scale = std::pow(10, sol_rnd_);
        auto scale_rec = 1.0/scale;
        for (auto& x: x_)
          x = std::round(x * scale) * scale_rec;
      }
    } catch (...) { sol_rnd_=100; }     // Could add a warning
    try {                 // Apply sol_prec
      if (sol_prec<100) {
        sol_prec_ = (sol_prec);
        for (auto& x: x_)
          x = round_to_digits(x, sol_prec_);
      }
    } catch (...) { sol_prec_=100; }     // Could add a warning
  }

private:
  double feastol_;
  bool recomp_vals_;   // variables are recomputed

  VarVec x_;   // can be rounded, recomputed, etc.
  ArrayRef<double> x_raw_;     // solver values
  const ArrayRef<var::Type> type_;
  const ArrayRef<double> lb_;
  const ArrayRef<double> ub_;
  int sol_rnd_=100;    // AMPL option solution_round, if used
  int sol_prec_=100;   // AMPL option solution_precision, if used
};


/// VarInfoRecompTypedef
using VarInfoRecompTypedef = VarInfoImpl<VarVecRecomp>;

/// Define VarInfoRecomp
class VarInfoRecomp : public VarInfoRecompTypedef {
public:
  /// Inherit ctor's
  using VarInfoRecompTypedef::VarInfoRecompTypedef;
};

/// VarInfoStatic
using VarInfoStatic = VarInfoImpl<VarVecStatic>;


/// Solution check data
struct SolCheck {
  /// Construct.
  /// @param chk_mode: can be subset of 1+2+4+8+16
  SolCheck(ArrayRef<double> x,
           const pre::ValueMapDbl& , //duals,
           ArrayRef<double> obj,
           ArrayRef<double> x_raw,
           ArrayRef<var::Type> vtype,
           ArrayRef<double> lb,  ArrayRef<double> ub,
           double feastol, double feastolrel,
           int sol_rnd, int sol_prec,
           bool recomp_vals, int chk_mode)
    : x_(feastol, recomp_vals,
         x, x_raw, vtype, lb, ub, sol_rnd, sol_prec),
      obj_(obj),
      feastol_(feastol), feastolrel_(feastolrel),
      fRecomputedVals_(recomp_vals),
      check_mode_(chk_mode) { }
  /// Any violations?
  bool HasAnyViols() const
  { return HasAnyConViols() || HasAnyObjViols(); }
  /// Any constraint violations?
  bool HasAnyConViols() const {
    return viol_var_bnds_[0].N_ || viol_var_bnds_[1].N_
        || viol_var_int_[0].N_ || viol_var_int_[1].N_
        || viol_cons_alg_.size()
        || viol_cons_log_.size();
  }
  /// Any objective value violations?
  bool HasAnyObjViols() const
  { return viol_obj_.N_; }

  /// Summary
  const std::string& GetReport() const { return report_; }

  /// VarInfo, can be used like x() for templates
  const VarInfoStatic& x_ext() const { return x_; }
  /// x[i]
  double x(int i) const { return x_[i]; }
  /// objective values
  const ArrayRef<double>& obj_vals() const
  { return obj_; }

  /// Absolute feasibility tolerance
  double GetFeasTol() const { return feastol_; }
  /// Relative feasibility tolerance
  double GetFeasTolRel() const { return feastolrel_; }

  /// Using recomputed aux vars?
  bool if_recomputed() const { return fRecomputedVals_; }
  /// Check mode
  int check_mode() const { return check_mode_; }

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

  /// Obj viols
  ViolSummary& ObjViols() { return viol_obj_; }

  /// Set report
  void SetReport(std::string rep)
  { report_ = std::move(rep); }

private:
  VarInfoStatic x_;
  ArrayRef<double> obj_;
  double feastol_;
  double feastolrel_;
  bool fRecomputedVals_;
  int check_mode_;

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
  virtual void AddUnbridgedToBackend(
      BasicFlatModelAPI& be, const std::vector<std::string>* vnames) = 0;

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

  /// Constraint name
  const char* GetConstraintName() const { return constr_name_; }

  /// Acceptance option names
  virtual const char* GetAcceptanceOptionNames() const
  { return solver_opt_nm_; }

  /// Constraint type short name.
  /// Ideally should be in the constraint itself,
  /// but currently we derive it from acceptance options.
  virtual const char* GetShortTypeName() const;

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
    bool optadded = true;
    if (cancvt) {   // If can convert. But see if ModelAPI accepts too
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
        } else
          optadded = false;
    } else
      optadded = false;
    if (!optadded)         // Still add as hidden option
      env.AddStoredOption(GetAcceptanceOptionNames(),
                          "HIDDEN",
                          GetAccLevRef(), 0, 2);
    // Description table
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
  int& GetAccLevRef() { return acceptance_level_; }


private:
  pre::ValueNode value_node_;
  const char* const constr_name_;
  const char* const solver_opt_nm_;
  mutable std::string type_name_short_;
  int acceptance_level_ {-1};
  BasicLogger* exporter_{};
};

const char*
BasicConstraintKeeper::GetShortTypeName() const {
  if (type_name_short_.empty()) {
    std::string acc_opt = GetAcceptanceOptionNames();
    assert(acc_opt.size());
    auto word_end = std::min(acc_opt.find(' '),
                             acc_opt.size());
    auto colon_pos = acc_opt.find(':');
    if (colon_pos>word_end)
      colon_pos = 0;
    type_name_short_ = acc_opt.substr(
          colon_pos, word_end-colon_pos);
    for (auto& c: type_name_short_)
      if (':'==c)
        c = '_';                // Markdown
    assert(type_name_short_.size());
  }
  return type_name_short_.c_str();
}


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

  /// Constrint Keeper description
  const std::string& GetDescription() const override
  { return desc_; }

  /// Assume Converter has the Backend
  Backend& GetBackend(BasicFlatConverter& cvt)
  { return static_cast<Converter&>(cvt).GetModelAPI(); }

  /// Add a pre-constructed constraint (or just arguments)
  /// @return index of the new constraint
  template <class... Args>
  int AddConstraint(int d, Args&&... args)
  {
    cons_.emplace_back( d, std::move(args)... );
    ExportConstraint(cons_.size()-1, cons_.back());
    return cons_.size()-1;
  }

  /// Get const constraint \a i
  const Constraint& GetConstraint(int i) const
  { assert(check_index(i)); return cons_[i].con_; }

  /// Get constraint \a i
  Constraint& GetConstraint(int i)
  { assert(check_index(i)); return cons_[i].con_; }

  /// Get constraint depth in the reformulation tree
  int GetConstraintDepth(int i) const
  { assert(check_index(i)); return cons_[i].GetDepth(); }

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
    return (int)cons_.size()-n_bridged_or_unused_;
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
      const std::vector<std::string>* pvnam) override {
    try {
      AddAllUnbridged(be, pvnam);
    } catch (const std::exception& exc) {
      MP_RAISE(std::string("Adding constraint of type '") +
                             Constraint::GetTypeName() + "' to " +
                             Backend::GetTypeName() + std::string(": ") +
                             exc.what());
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
  struct Container {
    Container(int d, Constraint&& c) noexcept
      : con_(std::move(c)), depth_(d) { }

    /// Depth in redef tree
    int GetDepth() const { return depth_; }

    /// Bridged (reformulated or just unused.)
    /// If only reformulated, can still be checked
    /// for solution correctness.
    bool IsBridged() const { return is_bridged_; }
    /// Mark as bridged
    void MarkAsBridged() { is_bridged_=true; }

    /// Unused (should not be checked)
    bool IsUnused() const { return is_unused_; }
    /// Mark as unused
    void MarkAsUnused() {
      MarkAsBridged();
      is_unused_=true;
    }

    Constraint con_;
    int depth_ = 0;
    bool is_bridged_ = false;
    bool is_unused_ = false;
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
    GetConverter().RunConversion(cnt.con_, i, cnt.GetDepth());
    MarkAsBridged(cnt, i);
  }

  /// Mark item as reformulated
  void MarkAsBridged(Container& cnt, int ) {
		cnt.MarkAsBridged();
    ++n_bridged_or_unused_;
	}

  /// Mark item as unused
  void MarkAsUnused(Container& cnt, int ) {
    cnt.MarkAsUnused();
    ++n_bridged_or_unused_;
  }

protected:
  /// Export (last added) constraint
  void ExportConstraint(int i_con, const Container& cnt) {
    if (GetLogger()) {
      fmt::MemoryWriter wrt;
      {
        MiniJSONWriter jw(wrt);
        jw["CON_TYPE"] = GetShortTypeName();
        jw["index"] = i_con;
        if (*cnt.con_.name())
          jw["name"] = cnt.con_.name();
        jw["depth"] = cnt.GetDepth();
        WriteJSON(jw["data"], cnt.con_);
      }
      wrt.write("\n");                     // EOL
      GetLogger()->Append(wrt);
    }
  }
  /// Export constraint status.
  /// This is called in the end,
  /// so printing the readable form.
  void ExportConStatus(int i_con, const Container& cnt,
                       const std::vector<std::string>* pvnam) {
    if (GetLogger()) {
      fmt::MemoryWriter wrt;
      {
        MiniJSONWriter jw(wrt);
        jw["CON_TYPE"] = GetShortTypeName();
        jw["index"] = i_con;
        if (*cnt.con_.name()) {
          jw["name"] = cnt.con_.name();
          if (pvnam && pvnam->size()) {
            fmt::MemoryWriter pr;
            WriteFlatCon(pr, cnt.con_, *pvnam);
            jw["printed"] = pr.c_str();
          }
        }
        jw["depth"] = cnt.GetDepth();
        jw["unused"] = (int)cnt.IsUnused();
        jw["bridged"] = (int)cnt.IsBridged();
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

  /// Is constraint \a i unused?
  bool IsUnused(int i) const override {
    return cons_.at(i).IsUnused();
  }

  /// Copy names from ValueNodes
  void CopyNamesFromValueNodes() override {
    const auto& vn = GetValueNode().GetStrVec();
    assert(vn.size()==cons_.size());
    for (auto i=vn.size(); i--; )
      cons_[i].con_.SetName(vn[i].MakeCurrentName());
  }

  /// Copy names to ValueNodes
  void CopyNames2ValueNodes() {
    auto& vn = GetValueNode().GetStrVec();
    assert(vn.size()==cons_.size());
    for (auto i=vn.size(); i--; )
      vn[i] = std::string(cons_[i].con_.name());
  }

  /// ForEachActive().
  /// Deletes every constraint where fn() returns true.
	template <class Fn>
	void ForEachActive(Fn fn) {
		for (int i=0; i<(int)cons_.size(); ++i)
			if (!cons_[i].IsBridged())
				if (fn(cons_[i].con_, i))
          MarkAsBridged(cons_[i], i);
	}

  /// Compute result for constraint \a i
  /// (for functional constraints).
  double
  ComputeValue(int i, const VarInfoRecomp& vir) override {
    assert(cons_[i].con_.GetResultVar() >= 0);
    return mp::ComputeValue(cons_[i].con_, vir);
  }

  /// Compute violations for this constraint type.
  /// We do it for redefined (intermediate) ones too.
  void ComputeViolations(SolCheck& chk) override {
    if (cons_.size()) {
      auto& conviolmap =
          cons_.front().con_.IsLogical() ?
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
            auto viol = cons_[i].con_.ComputeViolation(x);
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
                    viol, cr.second, cons_[i].con_.name());
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
                       const std::vector<std::string>* pvnam) {
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
      ExportConStatus(con_index, cont, pvnam);
      ++con_index;                      // increment index
    }
  }


private:
  Converter& cvt_;
  std::deque<Container> cons_;
  int i_cvt_last_ = -1;               // Last converted constraint.
  int n_bridged_or_unused_ = 0;       // Number of converted items,
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

} // namespace mp

#endif // CONSTRAINT_KEEPER_H
