#include "mp/backend-app.h"

/// Declare a backend factory
std::unique_ptr<mp::BasicBackend> CreateCuoptmpBackend();

extern "C" int main(int, char **argv) {
  return
      mp::RunBackendApp(argv, CreateCuoptmpBackend);
}
