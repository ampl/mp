#include "mp/convert/interface_app.h"
#include "mp/convert/expr_flattener.h"
#include "mp/convert/MIP/mp2mip.h"
#include "mp/convert/backend.h"

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
