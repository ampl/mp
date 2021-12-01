#ifndef BACKEND_BASE_H
#define BACKEND_BASE_H

#include "mp/flat/nlsolver_proxy_base.h"

namespace mp {

/// Abstract backend API
class BasicBackend {
public:
  virtual ~BasicBackend() = default;

  virtual void SolveAndReport() = 0;

  /// NLSolver (or whatever it is) should provide pCQ
  /// before Backend can run solving and query/provide values
  virtual void ProvideNLSolverProxyObject(NLSolverProxy* pCQ) = 0;
};

} // namespace mp

#endif // BACKEND_BASE_H
