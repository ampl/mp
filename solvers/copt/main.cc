#include "mp/backend-app.h"
#include "coptbackend.h"

/// Declare a backend factory
std::unique_ptr<mp::BasicBackend> CreateCoptBackend();

extern "C" int main1(int, char **argv) {
  try {
    mp::BackendApp s(CreateCoptBackend());
    return s.Run(argv);
  } catch (const std::exception &e) {
    fmt::print(stderr, "Error: {}\n", e.what());
  }
  return 0;
}


extern "C" int main2(int, char** argv,
  void* (*init)(), void (*check)(size_t, size_t, size_t)) {
  mp::BasicBackend::Callbacks callbacks = { init, check };
  return mp::RunBackendApp(argv, CreateCoptBackend, callbacks);
}
