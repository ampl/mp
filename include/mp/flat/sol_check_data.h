#ifndef SOL_CHECK_DATA_H
#define SOL_CHECK_DATA_H

#include <map>
#include <cassert>

#include "mp/arrayref.h"
#include "mp/common.h"
#include "mp/valcvt-node.h"
#include "mp/utils-math.h"

namespace mp {

/// Constraint/obj violation
struct Violation {
  double viol_;    // abs violation: >0 if really violated
  double valX_;    // value compared to
  /// Compute whether violated + relative violation.
  /// Both absolute and relative should be violated
  /// (relative only if refvar!=0.)
  std::pair<bool, double> Check(
      double epsabs, double epsrel) const {
    double violRel {0.0};
    if (viol_ > epsabs
        && (!valX_
            || (violRel=std::fabs(viol_/valX_))>epsrel))
      return {true, violRel};
    return {false, 0.0};
  }
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
  /// @param is_final: if provided, which values are final
  VarVecRecomp(std::vector<double> x,
               VarsRecomputeFn rec_fn,
               std::vector<bool> is_final = {})
      : x_(std::move(x)), is_recomp_(std::move(is_final)),
      recomp_fn_(rec_fn) {
    assert(recomp_fn_);
    is_recomp_.resize(x_.size());   // x_ now
  }
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
  /// x can be moved out
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
/// @param VarVec should provide recomputation info
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
  /// sol_prec as string
  std::string solution_precision() const
  { return sol_prec_ < 100 ? std::to_string(sol_prec_) : ""; }


protected:
  /// AMPL options solution_round, solution_prec
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

}  // namespace mp

#endif // SOL_CHECK_DATA_H
