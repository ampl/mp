#if defined(MP_LINK_WITH_SHARED_LIB) && !defined(MP_SOLVER_LIBRARY)

// If compiling the executable linked with the shared library,
// simply consume the entry point exported by the library.
// Do not define the backend here, so that main() cannot
// resolve back into itself.
extern "C" int mp2nl_main(int argc, char** argv);

int main(int argc, char** argv) {
  return mp2nl_main(argc, argv);
}

#else

#include "mp/backend-app.h"
#include "mp/ampls-c-api.h"

/// Declare a backend factory
std::unique_ptr<mp::BasicBackend> CreateMP2NLBackend();

#ifdef MP_SOLVER_LIBRARY

// If compiling the dynamic library, implement and export the entry point
// so that the executable can import it shallowly.
extern "C" AMPLS_C_EXPORT int mp2nl_main(int, char** argv) {
  return mp::RunBackendApp(argv, CreateMP2NLBackend);
}

#else

int main(int, char** argv) {
  return mp::RunBackendApp(argv, CreateMP2NLBackend);
}

#endif

#endif
