/*
 Copyright (C) 2024 AMPL Optimization Inc.

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization Inc disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Gleb Belov
 */
#ifndef NLSOLVER_HPP
#define NLSOLVER_HPP

#include "mp/nl-solver.h"
#include "mp/nl-writer2.hpp"
#include "mp/sol-reader2.hpp"

namespace mp {

template <class NLFeeder>
bool NLSolver::LoadModel(NLFeeder& nlf) {
  if (GetFileStub().empty())
    InitAutoStub();
  auto result = mp::WriteNLFile(GetFileStub(), nlf, Utils());
  if (NLW2_WriteNL_OK != result.first)
    return (err_msg_ = "WriteNL error: " + result.second, false);
  return true;
}

template <class SOLHandler>
bool NLSolver::ReadSolution(SOLHandler& solh) {
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
