#include "mp/backend-app.h"

/// Declare a backend factory
std::unique_ptr<mp::BasicBackend> CreateCplexBackend();

extern "C" int main1(int, char **argv) {
  return
      mp::RunBackendApp(argv, CreateCplexBackend);
}
