#include "mp/convert/interface_app.h"
#include "mp/convert/MIP/mp2mip.h"
#include "mp/convert/backend.h"

#include "cplexbackend.h"

extern "C" int main1(int, char **argv) {
  try {
    using CplexInterface =
        mp::Interface<mp::MPToMIPConverter, mp::CplexBackend>;
    using CplexInterfaceApp = mp::InterfaceApp<CplexInterface>;
    CplexInterfaceApp s;
    return s.RunFromNLFile(argv);
  } catch (const std::exception &e) {
    fmt::print(stderr, "Error: {}\n", e.what());
  }
  return 1;
}
