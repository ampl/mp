#ifndef BACKEND_WITH_PRE_H
#define BACKEND_WITH_PRE_H

#include "mp/presolve-base.h"

namespace mp {

/// Flat backends need presolver
/// for pre- / postsolving of suffix values etc
class BackendWithPresolver {
protected:
  void SetPresolver(pre::BasicPresolver* pPre) {
    pPresolver_ = pPre;
  }

  const pre::BasicPresolver& GetPresolver() const
  { assert(pPresolver_); return *pPresolver_; }
  pre::BasicPresolver& GetPresolver()
  { assert(pPresolver_); return *pPresolver_; }


private:
  pre::BasicPresolver* pPresolver_ = nullptr;
};

} // namespace mp

#endif // BACKEND_WITH_PRE_H
