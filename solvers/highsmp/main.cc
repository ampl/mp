#if defined(MP_LINK_WITH_SHARED_LIB) && !defined(MP_SOLVER_LIBRARY)
extern "C" int highs_main(int argc, char** argv);

int main(int argc, char** argv)
{
	return highs_main(argc, argv);
}
#else

#include "mp/backend-app.h"

/// Declare a backend factory
std::unique_ptr<mp::BasicBackend> CreateHighsBackend();

int main(int, char **argv) {
  return
      mp::RunBackendApp(argv, CreateHighsBackend);
}

#endif