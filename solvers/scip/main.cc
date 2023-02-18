#include "mp/backend-app.h"

/// Declare a backend factory
std::unique_ptr<mp::BasicBackend> CreateScipBackend();

extern "C" int main1(int, char **argv) {
  return
      mp::RunBackendApp(argv, CreateScipBackend);
}
