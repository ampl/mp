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

extern "C" int main2(int, char** argv, CCallbacks cb) {
  return mp::RunBackendApp(argv, CreateMP2NLBackend, cb);
}

#endif