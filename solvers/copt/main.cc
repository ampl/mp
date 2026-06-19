#if defined(MP_LINK_WITH_SHARED_LIB) && !defined(MP_SOLVER_LIBRARY)
extern "C" int copt_main(int argc, char** argv);

int main(int argc, char** argv)
{
	return copt_main(argc, argv);
}

#else

#include "mp/backend-app.h"

/// Declare a backend factory
std::unique_ptr<mp::BasicBackend> CreateCoptBackend();

extern "C" int main1(int, char **argv) {
  return
      mp::RunBackendApp(argv, CreateCoptBackend);
}
extern "C" int main2(int, char** argv, CCallbacks cb) {
  return mp::RunBackendApp(argv, CreateCoptBackend, cb);
}

#endif