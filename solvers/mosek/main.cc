#include "mp/backend-app.h"

/// Declare a backend factory
std::unique_ptr<mp::BasicBackend> CreateMosekBackend();

extern "C" int main1(int, char **argv) {
  return
      mp::RunBackendApp(argv, CreateMosekBackend);
}
