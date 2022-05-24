#include "mp/backend-app.h"
#include "coptbackend.h"

/// Declare a backend factory
std::unique_ptr<mp::BasicBackend> CreateCoptBackend();

extern "C" int main1(int, char **argv) {
  return
      mp::RunBackendApp(argv, CreateCoptBackend);
}
