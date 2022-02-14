#include "mp/error.h"
#include "mp/backend_app.h"

/// Declare a backend factory
std::unique_ptr<mp::BasicBackend> CreateGurobiBackend();

extern "C" int main1(int, char **argv) {
  try {
    mp::BackendApp s(CreateGurobiBackend());
    return s.Run(argv);
  } catch (const mp::Error &e) {
    fmt::print(stderr, "Error: {}\n", e.what());
    return e.exit_code();
  } catch (const std::exception &e) {
    fmt::print(stderr, "Error: {}\n", e.what());
  }
  return 0;
}
