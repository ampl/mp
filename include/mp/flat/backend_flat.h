#ifndef BACKEND_FLAT_H
#define BACKEND_FLAT_H

#include "mp/common.h"
#include "mp/backend-with-valcvt.h"
#include "mp/backend-mip.h"
#include "mp/suffix.h"

namespace mp {

/// Store presolved sensitivities
struct SensRangesPresolved {
  pre::ModelValuesDbl
    varlblo, varlbhi, varublo, varubhi,
    varobjlo, varobj, varobjhi,
    conrhslo, conrhshi,                   // for rhs-constraints
    conlblo, conlbhi, conublo, conubhi;   // for range constraints
};


/// Solver Backends using ModelAPI+FlatConverter
/// should derive from this.
/// This is a wrapper providing convenient methods
/// to automatically pre- / popstsolve standard suffixes etc
template <class BaseBackend>
class FlatBackend :
    public BaseBackend,
    public BackendWithValuePresolver {
public:
  /// Default GetSolution() for flat backends.
  /// Invokes postsolver.
  /// @note requires solve code retrieved.
  Solution GetSolution() override {
    auto x = PrimalSolution();
    auto y = DualSolution();
    auto objs = GetObjectiveValues();
    auto fKnownInfeasOrUnb = BaseBackend::IsProblemInfeasible();
    auto mv = GetValuePresolver().PostsolveSolution(
          { x,
            y,
            objs,
         (void*)fKnownInfeasOrUnb } );
    auto x1 = std::move(mv.GetVarValues()());
    if (x.empty())
      x1.clear();      // don't send variable values
    auto y1 = std::move(mv.GetConValues()());
    if (y.Empty())
        y1.clear();      // don't send dual values
    auto objs1 = std::move(mv.GetObjValues()());
    if (objs.empty())    // e.g., MP2NL
      objs1.clear();
    return{ std::move(x1),
            std::move(y1),
            std::move(objs1) };
  }

  /// Redeclare GetObjectiveValues()
  ArrayRef<double> GetObjectiveValues() override { return {}; }
	/// Separate PrimalSolution() for flat backends.
	/// @return empty vector if no primal solution.
  virtual ArrayRef<double> PrimalSolution() { return {}; }
	/// Separate DualSolution() for flat backends.
	/// @return empty map if no dual solution.
  virtual pre::ValueMapDbl DualSolution() { return {}; }

  /// Obtain presolved sens info
  virtual SensRangesPresolved GetSensRangesPresolved()
  { MP_RAISE("SensRangesPresolved() not implemented"); }

  /// Postsolve sensitivity info for flat backends.
  /// Currenlty needs con... maps to be correctly populated
  SensRanges GetSensRanges() override {
    SensRangesPresolved senspre = GetSensRangesPresolved();

    /// Should postsolve to obtain CON suffixes .sens(lb/ub)(lo/hi)
    /// containing sens ranges for range constraints.
    auto mvlbhi = GetValuePresolver().PostsolveGenericDbl(senspre.varlbhi);
    auto mvlblo = GetValuePresolver().PostsolveGenericDbl(senspre.varlblo);
    auto mvubhi = GetValuePresolver().PostsolveGenericDbl(senspre.varubhi);
    auto mvublo = GetValuePresolver().PostsolveGenericDbl(senspre.varublo);
    auto mvobjhi = GetValuePresolver().PostsolveGenericDbl(senspre.varobjhi);
    auto mvobj = GetValuePresolver().PostsolveGenericDbl(senspre.varobj);
    auto mvobjlo = GetValuePresolver().PostsolveGenericDbl(senspre.varobjlo);
    auto mvconlbhi = GetValuePresolver().PostsolveGenericDbl(senspre.conlbhi);
    auto mvconlblo = GetValuePresolver().PostsolveGenericDbl(senspre.conlblo);
    auto mvconubhi = GetValuePresolver().PostsolveGenericDbl(senspre.conubhi);
    auto mvconublo = GetValuePresolver().PostsolveGenericDbl(senspre.conublo);
    auto mvrhshi = GetValuePresolver().PostsolveGenericDbl(senspre.conrhshi);
    auto mvrhslo = GetValuePresolver().PostsolveGenericDbl(senspre.conrhslo);

    SensRanges sensr;
    sensr.varlbhi = mvlbhi.GetVarValues()();
    sensr.varlblo = mvlblo.GetVarValues()();
    sensr.varubhi = mvubhi.GetVarValues()();
    sensr.varublo = mvublo.GetVarValues()();
    sensr.varobjhi = mvobjhi.GetVarValues()();
    sensr.varobj = mvobj.GetVarValues()();
    sensr.varobjlo = mvobjlo.GetVarValues()();
    sensr.conlbhi = mvconlbhi.GetConValues()();
    sensr.conlblo = mvconlblo.GetConValues()();
    sensr.conubhi = mvconubhi.GetConValues()();
    sensr.conublo = mvconublo.GetConValues()();
    sensr.conrhshi = mvrhshi.GetConValues()();
    sensr.conrhslo = mvrhslo.GetConValues()();

    return sensr;
  }


protected:
  /// Supplementary LP data in solver coordinates, published atomically only
  /// after an optimal fixed solve and successful extraction.
  struct FixedLPResult {
    std::vector<double> primal, dual;
    std::vector<int> varstt, constt;
    double objective = 0;
    int iterations = 0;
  };
  virtual FixedLPResult SolveFixedLP(ArrayRef<double>) {
    MP_RAISE("Fixed-model emulation not implemented");
  }
  bool FixedModelReady() const override {
    return !this->IsFixedModelEmulated() || fixed_ready_;
  }
  const FixedLPResult& FixedLP() const { return fixed_lp_; }
  bool HasFixedLP() const { return fixed_ready_; }

  /// Called after recording MIP status/messages/pool, before standard reports.
  void ConsiderEmulatedFixedModel() {
    fixed_ready_ = false;
    fixed_lp_ = {};
    if (!BaseBackend::need_fixed_MIP() || !this->IsMIP())
      return;
    try {
      if (!this->GetMM().IsLinearAlgebraic())
        MP_RAISE("unsupported model: requires linear algebraic constraints and objective, without SOS/PWL/logical constructs");
      auto x = PrimalSolution();
      if (x.empty())
        MP_RAISE("no integer-feasible incumbent");
      if (this->fixed_model_time_limit() == 0 || this->fixed_model_iter_limit() == 0)
        MP_RAISE("fixed-model time or iteration limit is zero");
      auto result = SolveFixedLP(x);
      fixed_lp_ = std::move(result);
      fixed_ready_ = true;
      this->AddToSolverMessage(fmt::format(
          "Fixed MIP for mip:basis (emulated): {} iterations\n", fixed_lp_.iterations));
    } catch (const std::exception& e) {
      this->AddToSolverMessage(fmt::format(
          "Fixed MIP for mip:basis (emulated) failed: {}.\n", e.what()));
    }
  }

  SolutionBasis FixedLPBasis(int group) {
    if (!fixed_ready_)
      return {};
    auto mv = GetValuePresolver().PostsolveBasis(
        {fixed_lp_.varstt, {{{group, fixed_lp_.constt}}}});
    return {mv.GetVarValues()(), mv.GetConValues()()};
  }

  static double FixedIntegerValue(double x, double lb, double ub, double tol) {
    double rounded = std::round(x);
    if (!std::isfinite(x) || std::fabs(x-rounded) > tol ||
        rounded < lb-tol || rounded > ub+tol)
      MP_RAISE("incumbent has an invalid integer value");
    return rounded;
  }

  /// Feasibility relaxation penalties in solver coordinates:
  /// bounds of the variables and right-hand sides of the constraints
  /// of group \a con_group (e.g., CG_Linear). An empty vector means
  /// that no element may be relaxed.
  struct FeasrelaxPenalties {
    std::vector<double> lb, ub, rhs;
  };

  /// The penalties from the suffixes and options (alg:lbpen etc.).
  /// Entries >= \a inf, and all variables and constraints that the
  /// model conversion added, get \a not_relaxable (the solver's
  /// marker): relaxing an auxiliary element at no cost would violate
  /// the original constraint it represents.
  FeasrelaxPenalties FeasrelaxPenaltiesInSolver(
      int n_vars, int con_group, double inf, double not_relaxable) {
    const auto& frx = this->feasrelax();
    return { PadPenalties(frx.lbpen(), n_vars, inf, not_relaxable),
             PadPenalties(frx.ubpen(), n_vars, inf, not_relaxable),
             PresolveRhsPenalties(frx.rhspen(), con_group,
                                  inf, not_relaxable) };
  }

private:
  /// Penalties of the original variables, padded for the
  /// auxiliary ones (which follow them)
  static std::vector<double> PadPenalties(ArrayRef<double> pen, int n,
      double inf, double not_relaxable) {
    std::vector<double> result(pen.begin(), pen.end());
    if (result.empty())
      return result;
    result.resize(std::max<size_t>(result.size(), n), not_relaxable);
    for (auto& p: result)
      if (p >= inf)
        p = not_relaxable;
    return result;
  }

  /// Constraint penalties mapped through the value presolver.
  /// A constraint without an original one (it receives 0
  /// from the presolver for a vector of ones) is not relaxable.
  std::vector<double> PresolveRhsPenalties(ArrayRef<double> pen,
      int con_group, double inf, double not_relaxable) {
    if (pen.empty())
      return {};
    auto mv = GetValuePresolver().PresolveSolution({ {}, pen });
    std::vector<double> ones(pen.size(), 1.0);
    auto mv1 = GetValuePresolver().PresolveSolution({ {}, ones });
    std::vector<double> result = mv.GetConValues()(con_group);
    const auto& orig = mv1.GetConValues()(con_group);
    for (size_t i = 0; i < result.size(); ++i)
      if (i >= orig.size() || 0.0 == orig[i] || result[i] >= inf)
        result[i] = not_relaxable;
    return result;
  }

  FixedLPResult fixed_lp_;
  bool fixed_ready_ = false;

protected:
  /// Convenience feature:
  /// Read int suffixes for several entities (var/con/obj).
  /// @return suffix value map
  /// (can be checked for non-empty with operator bool())
  pre::ModelValuesInt ReadModelSuffixInt(const ModelSuffixDef<int>& msd)
  { return ReadModelSuffix(msd); }

  /// Read double suffixes for several entities (var/con/obj).
  /// @return suffix value map
  /// (can be checked for non-empty with operator bool())
  pre::ModelValuesDbl ReadModelSuffixDbl(const ModelSuffixDef<double>& msd)
  { return ReadModelSuffix(msd); }

  /// Implementation of ReadModelSuffix<>
  template <class T>
  pre::MVOverEl<T> ReadModelSuffix(const ModelSuffixDef<T>& msd) {
    assert(msd.kind() &
           (suf::Kind::VAR_BIT | suf::Kind::CON_BIT | suf::Kind::OBJ_BIT));
    return {
      suf::Kind::VAR_BIT & msd.kind() ?
            BaseBackend::template ReadSuffix<T>( {msd.name(), suf::Kind::VAR} ) :
            ArrayRef<T>{},
          suf::Kind::CON_BIT & msd.kind() ?
                BaseBackend::template ReadSuffix<T>( {msd.name(), suf::Kind::CON} ) :
            ArrayRef<T>{},
          suf::Kind::OBJ_BIT & msd.kind() ?
                BaseBackend::template ReadSuffix<T>( {msd.name(), suf::Kind::OBJ} ) :
            ArrayRef<T>{},
    };
  }
};

}  // namespace mp

#endif // BACKEND_FLAT_H
