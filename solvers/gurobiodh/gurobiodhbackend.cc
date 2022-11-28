#include "gurobiodhbackend.h"
#include "gurobi/gurobibackend.h"

extern "C" {
  #include "gurobi/gurobi-ampls-c-api.h"    // Gurobi AMPLS C API
}
#include "mp/ampls-cpp-api.h"

std::unique_ptr<mp::BasicBackend> CreateGurobiODHBackend() {
  return std::unique_ptr<mp::BasicBackend>{new mp::GurobiODHBackend()};
}

namespace mp {
  void GurobiODHBackend::InitOptionParsing() {
    OpenGurobi();
    OpenODH();
  }
  void GurobiODHBackend::OpenODH() {
    openODH(env());
  }

  void GurobiODHBackend::Solve() {
    PrepareGurobiSolve();
    int status = HEURopt(HEURptr(), model());
    status = HEURgetstat(HEURptr());
    fmt::format("At end of optimization the return code is {%d}:\n{}\n", status, getStatusMsg(status));
    WindupGurobiSolve();
  }
}