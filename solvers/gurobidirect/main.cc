#include "mp/error.h"
#include "mp/flat/interface_app.h"
#include "mp/flat/expr_flattener.h"
#include "mp/flat/MIP/mp2mip.h"
#include "mp/flat/backend.h"

#include "gurobibackend.h"

extern "C" int main1(int, char **argv) {
  try {
    using GurobiNLSolverWithFlatConversion =
        mp::NLSolverWithFlatBackend<mp::GurobiBackend, mp::MPToMIPConverter>;
    using GurobiInterfaceApp =
      mp::InterfaceApp<GurobiNLSolverWithFlatConversion>;
    GurobiInterfaceApp s;
    return s.RunFromNLFile(argv);
  } catch (const mp::Error &e) {
    fmt::print(stderr, "Error: {}\n", e.what());
    return e.exit_code();
  } catch (const std::exception &e) {
    fmt::print(stderr, "Error: {}\n", e.what());
  }
  return 1;
}
