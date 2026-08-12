#ifdef MP_LINK_WITH_SHARED_LIB

#ifdef MP_SOLVER_LIBRARY

// If compiling the  dynamic library, implement the main function so that
// the executable can import it shallowly
#include "mp/backend-app.h"

/// Declare a backend factory
std::unique_ptr<mp::BasicBackend> CreateHighsBackend();

extern "C" AMPLS_C_EXPORT int highs_main(int, char** argv) {
	return mp::RunBackendApp(argv, CreateHighsBackend);
}

#else

// If compiling the executable, simply consume the function above
extern "C" int highs_main(int argc, char** argv);

int main(int argc, char** argv)
{
	return highs_main(argc, argv);
}
#endif

#else

int main(int, char **argv) {
  return
      mp::RunBackendApp(argv, CreateHighsBackend);
}

#endif