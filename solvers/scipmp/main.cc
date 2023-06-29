#include "mp/backend-app.h"

std::unique_ptr<mp::BasicBackend> CreateScipBackend();

#ifndef SOLVER_LICNAME
int main(int, char** argv) {
  return mp::RunBackendApp(argv, CreateScipBackend);
}
#endif

extern "C" int main2(int, char** argv, CCallbacks cb) {
  return mp::RunBackendApp(argv, CreateScipBackend, cb);
}
