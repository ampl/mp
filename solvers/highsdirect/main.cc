#include "mp/backend-app.h"
#include "highsdirectbackend.h"

/// Declare a backend factory
std::unique_ptr<mp::BasicBackend> CreateHighsBackend();

int main(int, char **argv) {
  return
      mp::RunBackendApp(argv, CreateHighsBackend);
}
