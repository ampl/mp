#if defined(MP_LINK_WITH_SHARED_LIB) && !defined(MP_SOLVER_LIBRARY)
extern "C" int grb_main(int argc, char** argv);

int main(int argc, char** argv)
{
	return grb_main(argc, argv);
}
#else

#include "mp/backend-app.h"

/// Declare a backend factory
std::unique_ptr<mp::BasicBackend> CreateGurobiBackend();

extern "C" int main1(int, char **argv) {
  return mp::RunBackendApp(argv, CreateGurobiBackend);
}

#ifndef SLV_MAIN_IN_MAIN_CC
extern "C" int grb_main(int, char **argv) {
  return mp::RunBackendApp(argv, CreateGurobiBackend);
}
#endif  // SLV_MAIN_IN_MAIN_CC

extern "C" int main2(int, char** argv, CCallbacks cb) {
  return mp::RunBackendApp(argv, CreateGurobiBackend, cb);
}
#endif
