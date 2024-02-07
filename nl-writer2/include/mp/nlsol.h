/**
 mp::NLSOL.
 Manager for solving optimization models via NL files.
 It performs zero-overhead model/solution transmission.
 In particular, it does not store any intermediate
 model/solution representation.

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
#include <cstdio>
#include <cstdlib>

#include "mp/nl-writer2.h"
#include "mp/sol-reader2.h"

namespace mp {

/// \rst
/// Class NLSOL.
///
/// Manager for solving optimization models via NL files.
/// It performs zero-overhead model/solution transmission.
/// In particular, it does not store any intermediate
/// model/solution representation.
///
/// This class offers full NL functionality.
/// For simplified interface for special model classes,
/// see `~NLSOL_Easy`.
///
/// Usage (see tests/examples):
///
///   ExampleModel emdl;
///   ExampleNLFeeder2 nlf(emdl, binary);
///   ExampleSOLHandler2 esolh(emdl);
///   mp::NLUtils utils;
///
///   mp::NLSOL nlsol(&utils);
///   nlsol.SetFileStub(stub);
///   if (!nlsol.LoadModel(nlf)
///       || !nlsol.Solve(solver, sopts)
///       || !nlsol.ReadSolution(esolh)) {
///     printf("%s\n", nlsol.GetErrorMessage());
///     return EXIT_FAILURE;
///   } else {
///     esolh.PrintSolution(stub);
///   }
///
/// See mp::NLFeeder2 and mp::SOLHandler2 interfaces
/// for model/solution transmission.
/// \endrst
class NLSOL {
public:
  /// Construct.
  ///
  /// @param put: pointer to NLUtils or a derived object
  ///   (optional).
  NLSOL(mp::NLUtils* put=nullptr)
    : p_ut_(put ? put : &utils_) { Init(); }

  /// Destruct.
  ~NLSOL() { Destroy(); }

  /// Set NLUtils [OPTIONAL].
  ///
  /// If not provided, default is used.
  void SetNLUtils(mp::NLUtils* p_ut) { p_ut_ = p_ut; }

  /// Retrieve NLUtils
  mp::NLUtils* GetNLUtils() const { return p_ut_; }

  /// Set file stub [OPTIONAL].
  ///
  /// Used for filename base of .nl, .col, row, etc. input files,
  /// as well as .sol output files.
  ///
  /// If not provided, a temporary filename is used;
  /// then, .nl is deleted upon object destruction.
  void SetFileStub(std::string stub);

  /// Retrieve file stub.
  const std::string& GetFileStub() const
  { return filestub_; }

  /// Get error message.
  /// Nonempty iff error occurred.
  const char* GetErrorMessage() const { return err_msg_.c_str(); }

  /// Get model load result code.
  WriteNLResultCode GetLoadModelResultCode() const
  { return nl_result_; }

  /// Get solution read result code.
  SOLReadResultCode GetSolReadResultCode() const
  { return sol_result_; }

  /// Write NL and any accompanying files.
  ///
  /// @param nlf: NL feeder.
  ///
  /// @return true if all ok, otherwise see
  ///   GetErrorMessage() and, possibly, GetLoadModelResultCode().
  template <class NLFeeder2>
  bool LoadModel(NLFeeder2& nlf);

  /// Solve.
  ///
  /// @param solver: solver executable, such as "gurobi".
  /// @param solver_opts: string of solver options,
  ///   such as "outlev=1 writeprob=model.lp".
  ///
  /// @return true if all ok.
  bool Solve(const std::string& solver,
             const std::string& solver_opts);

  /// Read solution.
  ///
  /// @param solh: solution handler.
  ///
  /// @return true if all ok, otherwise see
  ///   GetErrorMessage() and, possibly, GetSolReadResultCode().
  template <class SOLHandler2>
  bool ReadSolution(SOLHandler2& solh);

protected:
  void Init() {
    // init file stub
    char tmpn[L_tmpnam];
    tmpnam(tmpn);
    filestub_ = tmpn;
  }
  void Destroy() {
    // try & delete .nl
    if (!filestubCustom_)
      std::remove((filestub_ + ".nl").c_str());
  }
  mp::NLUtils& Utils() const { return *p_ut_; }

private:
  mp::NLUtils utils_;
  mp::NLUtils* p_ut_ = nullptr;

  std::string filestub_;
  bool filestubCustom_ = false;

  std::string err_msg_;
  mp::WriteNLResultCode nl_result_
  {WriteNL_Unset};
  mp::SOLReadResultCode sol_result_
  {SOLRead_Result_Not_Set};
};

}  // namespace mp

#endif // NLSOL_H
