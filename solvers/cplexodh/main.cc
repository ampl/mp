#include "mp/backend-app.h"

/// Declare a backend factory
std::unique_ptr<mp::BasicBackend> CreateCplexODHBackend();

extern "C" int main1(int, char **argv) {
  return mp::RunBackendApp(argv, CreateCplexODHBackend);
}

extern "C" int main2(int, char** argv, CCallbacks cb) {
  return mp::RunBackendApp(argv, CreateCplexODHBackend, cb);
}
