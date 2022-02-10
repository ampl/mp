/*
 Flat model API.
 This is how Backend receives a flat model.
 TODO Separate NL entities from Flat objects.

 Copyright (C) 2021 AMPL Optimization Inc

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization Inc disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.
*/
#ifndef FLAT_MODEL_API_H_
#define FLAT_MODEL_API_H_

#include <type_traits>

#include "mp/flat/flat_model_api_basic.h"
#include "mp/flat/std_constr.h"
#include "mp/flat/std_obj.h"
#include "mp/presolve_base.h"

namespace mp {

/// FlatBackend: Backends receiving flat constraints
/// and handling Presolver (typycally, flat-API backends)
/// should derive from this
/// TODO rename into BackendWithPresolver
/// TODO move the model manip methods into FlatModelAPI
template <class Impl>
class FlatBackend : public BasicFlatBackend {
public:
  using BaseFlatBackend = BasicFlatBackend;

  ////////////////// Some standard items /////////////////
  void SetLinearObjective(int , const LinearObjective& ) {
    throw MakeUnsupportedError("FlatBackend::SetLinearObjective()");
  }

  void SetQuadraticObjective(int , const QuadraticObjective& ) {
    throw MakeUnsupportedError("FlatBackend::SetQuadraticObjective()");
  }

  USE_BASE_CONSTRAINT_HANDLERS(BaseFlatBackend)

  /// Optionally exclude LDCs from being posted,
  /// then all those are converted to LinearConstraint's first
  ACCEPT_CONSTRAINT(LinearDefiningConstraint, NotAccepted, CG_Linear)
  void AddConstraint(const LinearDefiningConstraint& ldc) {
    MP_DISPATCH( AddConstraint(ldc.to_linear_constraint()) );
  }

  /// FlatConverter should provide pPre before FlatBackend can run solving
  /// and request pre- / postsolving of suffix values etc
  void ProvidePresolver(pre::BasicPresolver* pPre) {
    pPresolver_ = pPre;
  }

protected:
  const pre::BasicPresolver& GetPresolver() const
  { assert(pPresolver_); return *pPresolver_; }
  pre::BasicPresolver& GetPresolver()
  { assert(pPresolver_); return *pPresolver_; }

private:
  pre::BasicPresolver* pPresolver_ = nullptr;
};

} // namespace mp

#endif // FLAT_MODEL_API_H_
