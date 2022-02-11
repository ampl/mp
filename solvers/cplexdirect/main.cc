#include "mp/flat/backend_app.h"
#include "mp/flat/nlsolver.h"
#include "mp/flat/MIP/mp2mip.h"
#include "mp/flat/backend.h"

#include "cplexbackend.h"

extern "C" int main1(int, char **argv) {
  try {
    using CplexNLSolverWithFlatConversion =
        mp::NLSolverWithFlatBackend<mp::CplexBackend, mp::MIPFlatConverter>;
    using CplexBackendApp =
      mp::BackendApp<CplexNLSolverWithFlatConversion>;
    CplexBackendApp s;
    return s.Run(argv);
  } catch (const std::exception &e) {
    fmt::print(stderr, "Error: {}\n", e.what());
  }
  return 0;
}
