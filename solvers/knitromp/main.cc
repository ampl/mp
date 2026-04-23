#include "mp/backend-app.h"

/// Declare a backend factory
std::unique_ptr<mp::BasicBackend> CreateKnitrompBackend();
extern "C" int main1(int, char** argv) {
    return
        mp::RunBackendApp(argv, CreateKnitrompBackend);
}

extern "C" int main2(int, char** argv, CCallbacks cb) {
    return mp::RunBackendApp(argv, CreateKnitrompBackend, cb);
}
