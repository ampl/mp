#include "cplexodhbackend.h"
#include "cplexmp/cplexmpbackend.h"

std::unique_ptr<mp::BasicBackend> CreateCplexODHBackend() {
  return std::unique_ptr<mp::BasicBackend>{new mp::CplexODHBackend()};
}

namespace mp {
  void CplexODHBackend::InitOptionParsing() {
    OpenSolver();
    OpenODH(env());
  }

  void CplexODHBackend::Solve() {
    int status = HEURopt(HEURptr(), env());
    status = HEURgetstat(HEURptr());
    fmt::print("At end of optimization the return code is %{}: {}\n", 
      status, GetODHStatusMsg(status));
    WindupCPLEXSolve();
  }

  CplexODHBackend::~CplexODHBackend() {
    CloseODH();
  }
}