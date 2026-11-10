#if defined(MP_LINK_WITH_SHARED_LIB) && !defined(MP_SOLVER_LIBRARY)
extern "C" int mp2nl_main(int argc, char** argv);

int main(int argc, char** argv)
{
	return mp2nl_main(argc, argv);
}
#else

#include "mp/backend-app.h"

std::unique_ptr<mp::BasicBackend> CreateMP2NLBackend();

#ifndef SOLVER_LICNAME
int main(int, char** argv) {
  return mp::RunBackendApp(argv, CreateMP2NLBackend);
}
#endif

#ifdef MP_LINK_WITH_SHARED_LIB
// Entry point for the executable linked with the shared library.
// Call RunBackendApp() directly: forwarding to main() would, on Linux,
// bind to the executable's main(), which calls back into mp2nl_main().
extern "C" AMPLS_C_EXPORT int mp2nl_main(int, char** argv) {
  return mp::RunBackendApp(argv, CreateMP2NLBackend);
}
#endif

extern "C" int main2(int, char** argv, CCallbacks cb) {
  return mp::RunBackendApp(argv, CreateMP2NLBackend, cb);
}

#endif
