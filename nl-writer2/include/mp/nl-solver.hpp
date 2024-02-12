#ifndef NLSOLVER_HPP
#define NLSOLVER_HPP

#include "mp/nl-solver.h"
#include "mp/nl-writer2.hpp"
#include "mp/sol-reader2.hpp"

namespace mp {

template <class NLFeeder2>
bool NLSolver::LoadModel(NLFeeder2& nlf) {
  if (GetFileStub().empty())
    return (err_msg_="WriteNL error: provide filestub.", false);
  auto result = mp::WriteNLFile(GetFileStub(), nlf, Utils());
  if (NLW2_WriteNL_OK != result.first)
    return (err_msg_ = "WriteNL error: " + result.second, false);
  return true;
}

template <class SOLHandler2>
bool NLSolver::ReadSolution(SOLHandler2& solh) {
  if (GetFileStub().empty())
    return (err_msg_="SOLReader: provide filename.", false);
  auto status = mp::ReadSOLFile(
        GetFileStub() + ".sol", solh, Utils());
  if (NLW2_SOLRead_OK != status.first)
    return (err_msg_="SOLReader error: "+status.second, false);
  return true;
}

}  // namespace mp

#endif // NLSOLVER_HPP
