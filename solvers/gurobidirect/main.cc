#include "mp/error.h"
#include "mp/convert/interface_app.h"
#include "mp/convert/expr_flattener.h"
#include "mp/convert/MIP/mp2mip.h"
#include "mp/convert/backend.h"

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
