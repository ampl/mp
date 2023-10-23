/**
 mp::NLSOL.
 Manager for solving optimization models via NL files.
 It performs zero-overhead model/solution transmission.
 In particular, it does not store any intermediate
 model representation.

 NL is a format for representing optimization problems such as linear,
 quadratic, nonlinear, complementarity and constraint programming problems
 in discrete or continuous variables. It is described in the technical report
 "Writing .nl Files" (http://www.cs.sandia.gov/~dmgay/nlwrite.pdf).


 Copyright (C) 2023 AMPL Optimization, Inc.

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization Inc. disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Gleb Belov
 */
#ifndef NLSOL_H
#define NLSOL_H

#include <utility>
#include <cassert>
#include <cstdlib>

#include "mp/nl-writer2.h"
#include "mp/sol-reader2.h"

namespace mp {

/// Class NLSOL.
///
/// Manager for solving optimization models via NL files.
/// It performs zero-overhead model/solution transmission.
/// In particular, it does not store any intermediate
/// model representation.
///
/// Usage:
///
///   MyNLFeeder feeder;
///   MySOLHandler handler;
///   mp::NLUtils nlutils;
///   mp::NLSOL<MyNLFeeder, MySOLHandler>
///     nlsol(feeder, handler, nlutils);
///   nlsol.SetSolver("minos");
///   nlsol.SetSolverOptions("outlev=1");
///   if (!nlsol.Solve("filestub")) {
///     printf("%s\n", nlsol.GetErrorMessage());
///   }
///
/// If separate steps are needed, use a subset of
///
///   if (!nlsol.WriteNLFile(filestub)
///     || !nlsol.InvokeSolver(filestub)
///       || !nlsol.ReadSolution(filestub + ".sol")) ...
///
/// See mp::NLFeeder2 and mp::SOLHandler2 interfaces
/// for model/solution transmission.
///
/// @param MyNLFeeder: a class implementing
///   the mp::NLFeeder2 interface.
/// @param MySOLHandler: a class implementing
///   the mp::SOLHandler2 interface.
template <class MyNLFeeder, class MySOLHandler>
class NLSOL {
public:
  /// Construct
  NLSOL(MyNLFeeder& fd, MySOLHandler& sh, mp::NLUtils& ut)
    : feeder_(fd), handler_(sh), utils_(ut) { }

  /// NLFeederType
  using NLFeederType = MyNLFeeder;
  /// SOLHandlerType
  using SOLHandlerType = MySOLHandler;

  /// Set solver, such as "gurobi", "highs", "ipopt"
  void SetSolver(std::string solver);

  /// Set solver options, such as "outlev=1 lim:time=500"
  void SetSolverOptions(std::string sopts);

  /// Solve.
  /// @param filestub: filename stub to be used
  /// for input files (.nl, .col., .row, etc.),
  /// and output files (.sol).
  /// @return true if all ok.
  bool Solve(const std::string& filestub);

  /// Get error message.
  const char* GetErrorMessage() const { return err_msg_.c_str(); }

  /// Substep: write NL and any accompanying files.
  bool WriteNLFile(const std::string& filestub);

  /// Substep: invoke chosen solver for \a filestub.
  bool InvokeSolver(const std::string& filestub);

  /// Substep: read solution.
  /// @param filename: complete file name,
  /// normally (stub).sol.
  bool ReadSolution(const std::string& filename);


private:
  MyNLFeeder& feeder_;
  MySOLHandler& handler_;
  mp::NLUtils& utils_;

  std::string solver_;
  std::string solver_options_;

  std::string err_msg_;
};


template <class MyNLFeeder, class MySOLHandler>
void NLSOL<MyNLFeeder, MySOLHandler>::
SetSolver(std::string solver) { solver_ = std::move(solver); }

template <class MyNLFeeder, class MySOLHandler>
void NLSOL<MyNLFeeder, MySOLHandler>::
SetSolverOptions(std::string sopts)
{ solver_options_ = std::move(sopts); }

template <class MyNLFeeder, class MySOLHandler>
bool NLSOL<MyNLFeeder, MySOLHandler>::
Solve(const std::string& filestub) {
  return (WriteNLFile(filestub)
    && InvokeSolver(filestub)
      && ReadSolution(filestub + ".sol"));
}

template <class MyNLFeeder, class MySOLHandler>
bool NLSOL<MyNLFeeder, MySOLHandler>::
WriteNLFile(const std::string& filestub) {
  if (filestub.empty())
    return (err_msg_="WriteNL error: provide filestub.", false);
  auto result = mp::WriteNLFile(filestub, feeder_, utils_);
  if (mp::WriteNL_OK != result.first)
    return (err_msg_ = "WriteNL error: " + result.second, false);
  return true;
}

template <class MyNLFeeder, class MySOLHandler>
bool NLSOL<MyNLFeeder, MySOLHandler>::
InvokeSolver(const std::string& filestub) {
  if (filestub.empty())
    return (err_msg_="NLSOL: provide filestub.", false);
  if (solver_.empty())
    return (err_msg_="NLSOL: provide solver.", false);
  auto call = solver_
      + ' ' + filestub
      + " -AMPL "
      + solver_options_;
  if (auto status = std::system(call.c_str()))
    return (err_msg_="NLSOL: call \""
        + call + "\" failed (code "
        + std::to_string(status) + ").", false);
  return true;
}

template <class MyNLFeeder, class MySOLHandler>
bool NLSOL<MyNLFeeder, MySOLHandler>::
ReadSolution(const std::string& filename) {
  if (filename.empty())
    return (err_msg_="SOLReader: provide filename.", false);
  auto status = mp::ReadSOLFile(
        filename, handler_, utils_);
  if (mp::SOL_Read_OK != status.first)
    return (err_msg_="SOLReader error: "+status.second, false);
  return true;
}

}  // namespace mp

#endif // NLSOL_H
