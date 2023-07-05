#include "mp/backend-app.h"

std::unique_ptr<mp::BasicBackend> CreateGcgBackend();

#ifndef SOLVER_LICNAME
int main(int, char** argv) {
  return mp::RunBackendApp(argv, CreateGcgBackend);
}
#endif

extern "C" int main2(int, char** argv, CCallbacks cb) {
  return mp::RunBackendApp(argv, CreateGcgBackend, cb);
}
