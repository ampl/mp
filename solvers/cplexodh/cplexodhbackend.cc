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
    WindupCPLEXSolve();
  }

  void CplexODHBackend::ReportResults() {
    ReportODHResults();
    FlatBackend< MIPBackend<CplexBackend> >::ReportResults();
  }

  void CplexODHBackend::ReportODHResults() {
    SetStatus(ConvertODHStatus());
    AddCPLEXMessages();
  }


  CplexODHBackend::~CplexODHBackend() {
    CloseODH();
  }
}