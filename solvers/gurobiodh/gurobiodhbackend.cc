#include "gurobiodhbackend.h"
#include "gurobi/gurobibackend.h"

std::unique_ptr<mp::BasicBackend> CreateGurobiODHBackend() {
  return std::unique_ptr<mp::BasicBackend>{new mp::GurobiODHBackend()};
}

namespace mp {
  void GurobiODHBackend::InitOptionParsing() {
    OpenGurobi();
    OpenODH(env());
  }

  void GurobiODHBackend::Solve() {
    PrepareGurobiSolve();
    int status = HEURopt(HEURptr(), model());
    status = HEURgetstat(HEURptr());
    fmt::format("At end of optimization the return code is {%d}:\n{}\n", status, GetODHStatusMsg(status));
    WindupGurobiSolve();
  }
  GurobiODHBackend::~GurobiODHBackend() {
    CloseODH();
  }
}