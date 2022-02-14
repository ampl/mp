#include "mp/error.h"
#include "mp/backend_app.h"

#include "gurobibackend.h"

extern "C" int main1(int, char **argv) {
  try {
    using GurobiBackendApp =
      mp::BackendApp<mp::GurobiBackend>;
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
