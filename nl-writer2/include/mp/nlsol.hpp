#ifndef NLSOL_HPP
#define NLSOL_HPP

#include "mp/nlsol.h"
#include "mp/nl-writer2.hpp"
#include "mp/sol-reader2.hpp"

namespace mp {

void NLSOL::SetFileStub(std::string stub) {
  if (stub.size()) {
    filestub_ = stub;
    filestubCustom_ = true;
  }
}

template <class NLFeeder2>
bool NLSOL::LoadModel(NLFeeder2& nlf) {
  if (GetFileStub().empty())
    return (err_msg_="WriteNL error: provide filestub.", false);
  auto result = mp::WriteNLFile(GetFileStub(), nlf, Utils());
  if (NLW2_WriteNL_OK != result.first)
    return (err_msg_ = "WriteNL error: " + result.second, false);
  return true;
}

bool NLSOL::Solve(const std::string& solver,
                  const std::string& solver_opts) {
  if (GetFileStub().empty())
    return (err_msg_="NLSOL: provide filestub.", false);
  if (solver.empty())
    return (err_msg_="NLSOL: provide solver.", false);
  auto call = solver
      + ' ' + GetFileStub()
      + " -AMPL "
      + solver_opts;
  if (auto status = std::system(call.c_str()))
    return (err_msg_="NLSOL: call \""
        + call + "\" failed (code "
        + std::to_string(status) + ").", false);
  return true;
}

template <class SOLHandler2>
bool NLSOL::ReadSolution(SOLHandler2& solh) {
  if (GetFileStub().empty())
    return (err_msg_="SOLReader: provide filename.", false);
  auto status = mp::ReadSOLFile(
        GetFileStub() + ".sol", solh, Utils());
  if (NLW2_SOLRead_OK != status.first)
    return (err_msg_="SOLReader error: "+status.second, false);
  return true;
}

}  // namespace mp

#endif // NLSOL_HPP
