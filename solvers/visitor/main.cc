#include "mp/backend-app.h"
#include "visitorbackend.h"

/// Declare a backend factory
std::unique_ptr<mp::BasicBackend> CreateVisitorBackend();

extern "C" int main1(int, char **argv) {
  try {
    mp::BackendApp s(CreateVisitorBackend());
    return s.Run(argv);
  } catch (const std::exception &e) {
    fmt::print(stderr, "Error: {}\n", e.what());
  }
  return 0;
}
