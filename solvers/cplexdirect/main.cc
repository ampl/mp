#include "mp/backend_app.h"

#include "cplexbackend.h"

extern "C" int main1(int, char **argv) {
  try {
    using CplexBackendApp =
      mp::BackendApp<mp::CplexBackend>;
    CplexBackendApp s;
    return s.Run(argv);
  } catch (const std::exception &e) {
    fmt::print(stderr, "Error: {}\n", e.what());
  }
  return 0;
}
