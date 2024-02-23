/*
 mp::NLSolver.

 Manager for solving optimization models via NL files.
 It performs almost zero-overhead model/solution transmission.
 In particular, it does not store any intermediate
 model/solution representation.

 NL is a format for representing optimization problems such as linear,
 quadratic, nonlinear, complementarity and constraint programming problems
 in discrete or continuous variables. It is described in the technical report
 "Writing .nl Files" (http://www.cs.sandia.gov/~dmgay/nlwrite.pdf).


 Copyright (C) 2024 AMPL Optimization, Inc.

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
#ifndef NLSOLVER_H
#define NLSOLVER_H

#include <utility>
#include <cassert>
#include <cstdio>
#include <cstdlib>

#include "mp/nl-solver-basics-c.h"
#include "mp/nl-model.h"
#include "mp/nl-utils.h"

namespace mp {

struct NLHeader;

/**
  \rst
  Class `~mp::NLSolver`.

  Manager for solving optimization models via NL files.
  It performs zero-overhead model/solution transmission.
  In particular, it does not store any intermediate
  model/solution representation.

  This class offers both full NL functionality,
  as well as a simplified interface for special model classes,
  see `~NLModel`.

  **Usage with full API** (see
  `NLFeeder` and `~mp::SOLHandler` interfaces
  for model/solution transmission,
  as well as tests/examples):

  .. code-block:: c++

    ExampleModel emdl;
    ExampleNLFeeder nlf(emdl, binary);
    ExampleSOLHandler esolh(emdl);
    mp::NLUtils utils;

  .. code-block:: c++

    mp::NLSolver nlsol(&utils);
    nlsol.SetFileStub(stub);
    if (!nlsol.Solve(nlf, esolh, solver, sopts)) {
      printf("%s\n", nlsol.GetErrorMessage());
      return EXIT_FAILURE;
    } else {
      esolh.PrintSolution(stub);
    }

  \endrst
*/
class NLSolver {
public:
  /// Construct.
  NLSolver();

  /// Construct.
  ///
  /// @param put: pointer to NLUtils or a derived object
  ///   (optional).
  NLSolver(mp::NLUtils* put);

  /// Destruct.
  ~NLSolver();

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
  /// @note If not provided, a temporary filename is used;
  /// then, .nl is deleted upon object destruction.
  void SetFileStub(std::string stub);

  /// Retrieve file stub.
  const std::string& GetFileStub() const
  { return filestub_; }

  /// Set some NL options for NLModel output [OPTIONAL].
  ///
  /// If not provided, default is used.
  /// @note Not used for NLFeeder output.
  void SetNLOptions(NLW2_NLOptionsBasic_C nlo) { nl_opts_=nlo; }

  /// Get NLOptions for NLModel output.
  NLW2_NLOptionsBasic_C GetNLOptions() const { return nl_opts_; }


  /// Get error message.
  /// Nonempty iff error occurred.
  const char* GetErrorMessage() const { return err_msg_.c_str(); }

  /// Get model load result code.
  NLW2_WriteNLResultCode GetLoadModelResultCode() const
  { return nl_result_; }

  /// Get solution read result code.
  NLW2_SOLReadResultCode GetSolReadResultCode() const
  { return sol_result_; }

  /// Load and solve an NLModel instance.
  ///
  /// @return NLSolution object with computed obj_value_
  ///    (has operator bool() for checking
  ///     if any result was obtained.)
  ///
  /// See LoadModel(), Solve(), ReadSolution()
  /// for details.
  NLSolution Solve(const NLModel& mdl,
                   const std::string& solver,
                   const std::string& solver_opts) {
    NLSolution sol;
    if (LoadModel(mdl)
        && Solve(solver, solver_opts)) {
      sol = ReadSolution();
      if (sol.x_.size())
        sol.obj_val_ = mdl.ComputeObjValue(sol.x_.data());
    }
    return sol;
  }

  /// Load and solve model
  /// using NLFeeder and SOLHandler.
  ///
  /// @return true iff all ok.
  ///
  /// See LoadModel(), Solve(), ReadSolution()
  /// for details.
  template <class NLFeeder, class SOLHandler>
  bool Solve(NLFeeder& nlf, SOLHandler& solh,
             const std::string& solver,
             const std::string& solver_opts) {
    return LoadModel(nlf)
      && Solve(solver, solver_opts)
      && ReadSolution(solh);
  }

  /// Write NL and any accompanying files.
  /// NL file name base and some options
  /// can be provided, if non-defaults desired,
  /// via SetFileStub() and SetNLOptions().
  ///
  /// @return true if all ok, otherwise see
  ///   GetErrorMessage().
  bool LoadModel(const NLModel& mdl);

  /// Write NL and any accompanying files.
  ///
  /// @param nlf: NL feeder.
  ///
  /// @return true if all ok, otherwise see
  ///   GetErrorMessage() and GetLoadModelResultCode().
  template <class NLFeeder>
  bool LoadModel(NLFeeder& nlf);

  /// Solve after loading model.
  ///
  /// @param solver: solver executable, such as "gurobi".
  /// @param solver_opts: string of solver options,
  ///   such as "outlev=1 writeprob=model.lp".
  ///
  /// @return true if all ok.
  bool Solve(const std::string& solver,
             const std::string& solver_opts);

  /// Read solution after Solve()
  /// when the NL file was written from NLModel.
  ///
  /// @return NLSolution object
  ///    (has operator bool() for checking
  ///     if any result was obtained.)
  ///
  /// @note To compute objective value,
  ///   execute NLModel::ComputeObjValue()
  ///   if x_ available.
  NLSolution ReadSolution();

  /// Read solution after Solve() when NL file
  /// was written from NLFeeder.
  ///
  /// @param solh: solution handler.
  ///
  /// @return true if all ok, otherwise see
  ///   GetErrorMessage() and, possibly, GetSolReadResultCode().
  template <class SOLHandler>
  bool ReadSolution(SOLHandler& solh);

protected:
  void InitAutoStub();
  void DestroyAutoStub();
  mp::NLUtils& Utils() const { return *p_ut_; }

private:
  mp::NLUtils utils_;
  mp::NLUtils* p_ut_ = nullptr;

  std::string pathstr_;
  std::string filestub_;
  bool filestubCustom_ = false;
  NLW2_NLOptionsBasic_C nl_opts_;

  // NLModel stuff
  std::unique_ptr<NLHeader> p_nlheader_;
  NLModel::PreprocessData pd_;

  std::string err_msg_;
  NLW2_WriteNLResultCode nl_result_
  {NLW2_WriteNL_Unset};
  NLW2_SOLReadResultCode sol_result_
  {NLW2_SOLRead_Result_Not_Set};
};

}  // namespace mp

#endif // NLSOLVER_H
