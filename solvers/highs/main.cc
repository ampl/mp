#include "mp/backend-app.h"
#include "highsbackend.h"

/// Declare a backend factory
std::unique_ptr<mp::BasicBackend> CreateHighsBackend();

int main(int, char **argv) {
  try {
    mp::BackendApp s(CreateHighsBackend());
    return s.Run(argv);
  } catch (const std::exception &e) {
    fmt::print(stderr, "Error: {}\n", e.what());
  }
  return 0;
}
