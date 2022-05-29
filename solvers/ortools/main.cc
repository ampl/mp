#include "mp/backend-app.h"
#include "ortoolsbackend.h"

/// Declare a backend factory
std::unique_ptr<mp::BasicBackend> CreateOrtoolsBackend();

int main(int, char** argv) {
  return
    mp::RunBackendApp(argv, CreateOrtoolsBackend);
}
