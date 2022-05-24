#include "mp/backend-app.h"

/// Declare a backend factory
std::unique_ptr<mp::BasicBackend> CreateCplexBackend();

extern "C" int main1(int, char **argv) {
  return
      mp::RunBackendApp(argv, CreateCplexBackend);
}

extern "C" int main2(int, char** argv,
  void* (*init)(), void (*check)(size_t, size_t, size_t)) {
  mp::BasicBackend::Callbacks callbacks = { init, check };
  return mp::RunBackendApp(argv, CreateCplexBackend, callbacks);
}
