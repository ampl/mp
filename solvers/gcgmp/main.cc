#include "mp/backend-app.h"

std::unique_ptr<mp::BasicBackend> CreateGcgBackend();

extern "C" int main1(int, char** argv) {
  return mp::RunBackendApp(argv, CreateGcgBackend);
}

extern "C" int main2(int, char** argv, CCallbacks cb) {
  return mp::RunBackendApp(argv, CreateGcgBackend, cb);
}
