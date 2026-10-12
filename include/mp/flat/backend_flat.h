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

private:
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
