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
    varobjlo, varobjhi,
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
  Solution GetSolution() {
    auto mv = GetValuePresolver().PostsolveSolution(
          { PrimalSolution(),
            DualSolution(),
            GetObjectiveValues() } );
    return{ mv.GetVarValues()(),
          mv.GetConValues()(),
          mv.GetObjValues()() };
  }

  /// Redeclare GetObjectiveValues()
  ArrayRef<double> GetObjectiveValues() override { return {}; }
  /// Separate PrimalSolution() for flat backends
  virtual ArrayRef<double> PrimalSolution() { return {}; }
  /// Separate DualSolution() for flat backends
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
