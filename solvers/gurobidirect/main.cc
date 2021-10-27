#include "mp/error.h"
#include "mp/convert/interface_app.h"
#include "mp/convert/MIP/mp2mip.h"
#include "mp/convert/backend.h"

#include "gurobibackend.h"

extern "C" int main1(int, char **argv) {
  try {
    using GurobiInterface =
        mp::Interface<mp::MPToMIPConverter, mp::GurobiBackend>;
    using GurobiInterfaceApp = mp::InterfaceApp<GurobiInterface>;
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
