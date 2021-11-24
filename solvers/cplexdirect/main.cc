#include "mp/flat/interface_app.h"
#include "mp/flat/expr_flattener.h"
#include "mp/flat/MIP/mp2mip.h"
#include "mp/flat/backend.h"

#include "cplexbackend.h"

extern "C" int main1(int, char **argv) {
  try {
    using CplexNLSolverWithFlatConversion =
        mp::NLSolverWithFlatBackend<mp::CplexBackend, mp::MPToMIPConverter>;
    using CplexInterfaceApp =
      mp::InterfaceApp<CplexNLSolverWithFlatConversion>;
    CplexInterfaceApp s;
    return s.RunFromNLFile(argv);
  } catch (const std::exception &e) {
    fmt::print(stderr, "Error: {}\n", e.what());
  }
  return 1;
}
