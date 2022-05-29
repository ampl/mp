#include "mp/backend-app.h"
#include "mp/ampls-c-api.h" // for CCallbacks

/// Declare a backend factory
std::unique_ptr<mp::BasicBackend> CreateGurobiBackend();

extern "C" int main1(int, char **argv) {
  return mp::RunBackendApp(argv, CreateGurobiBackend);
}

extern "C" int main2(int, char** argv, CCallbacks cb) {
  mp::BasicBackend::Callbacks callbacks = { cb.init,
    cb.check, cb.additionalText, cb.diagnostics };
  return mp::RunBackendApp(argv, CreateGurobiBackend, callbacks);
}
