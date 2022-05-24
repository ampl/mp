#include "mp/backend-app.h"
#include "visitorbackend.h"

/// Declare a backend factory
std::unique_ptr<mp::BasicBackend> CreateVisitorBackend();

extern "C" int main1(int, char **argv) {
  return
      mp::RunBackendApp(argv, CreateVisitorBackend);
}
