#include "mp/backend_app.h"

/// Declare a backend factory
std::unique_ptr<mp::BasicBackend> CreateCplexBackend();

extern "C" int main1(int, char **argv) {
  try {
    mp::BackendApp s(CreateCplexBackend());
    return s.Run(argv);
  } catch (const std::exception &e) {
    fmt::print(stderr, "Error: {}\n", e.what());
  }
  return 0;
}
