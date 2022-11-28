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
    WindupGurobiSolve();
  }

  void  GurobiODHBackend::ReportResults()  {
    ReportGurobiResults();
    FlatBackend< MIPBackend<GurobiBackend> >::ReportResults();
  }

  void  GurobiODHBackend::ReportGurobiResults() {
    SetStatus(ConvertODHStatus());
    AddGurobiMessage();
    if (need_multiple_solutions())
      ReportGurobiPool();
    if (need_fixed_MIP())
      ConsiderGurobiFixedModel();
  }

  GurobiODHBackend::~GurobiODHBackend() {
    CloseODH();
  }

}