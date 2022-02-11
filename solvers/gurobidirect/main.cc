#include "mp/error.h"
#include "mp/flat/backend_app.h"
#include "mp/flat/nlsolver.h"
#include "mp/flat/MIP/mp2mip.h"
#include "mp/flat/backend.h"

#include "gurobibackend.h"

extern "C" int main1(int, char **argv) {
  try {
    using GurobiNLSolverWithFlatConversion =
        mp::NLSolverWithFlatBackend<mp::GurobiBackend, mp::MIPFlatConverter>;
    using GurobiBackendApp =
      mp::BackendApp<GurobiNLSolverWithFlatConversion>;
    GurobiBackendApp s;
    return s.Run(argv);
  } catch (const mp::Error &e) {
    fmt::print(stderr, "Error: {}\n", e.what());
    return e.exit_code();
  } catch (const std::exception &e) {
    fmt::print(stderr, "Error: {}\n", e.what());
  }
  return 0;
}
