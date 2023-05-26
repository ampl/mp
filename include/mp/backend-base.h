#ifndef BACKEND_BASE_H
#define BACKEND_BASE_H

#include <string>

extern "C" {
  #include "mp/ampls-c-api.h" // for CCallbacks
}

#include "mp/solver-base.h"

namespace mp {

/// Abstract backend API
class BasicBackend : public BasicSolver {
public:
  /// Virtual destructor
  virtual ~BasicBackend() = default;

  /// Initialize backend, incl. solver options.
  /// Not parse options yet, as we need to parse cmdline arguments first.
  /// @param argv: the command-line arguments, 0-terminated
  virtual void Init(char** argv) {
    set_exe_path(*argv);
  }

  /// Provide cmdline solver options list and flags
  virtual void SetOptionData(char** argv, int f) {
    argv_options_ = argv;
    flag_options_ = f;
  }

  /// Parse solver options such as "outlev=1" from env and argv.
  /// @param filename_no_ext: basname of the .nl file
  ///        (or whatever should be the basename of an output .sol file)
  /// @param argv: (remaining part of) vector of cmdline strings
  /// @param flags: 0 or \a Solver::NO_OPTION_ECHO
  virtual
  bool ParseSolverOptions(char **argv, unsigned flags = 0) {
    /// Chance e.g. for the Backend to init solver environment, etc
    InitOptionParsing();
    if (ParseOptions(argv, flags)) {
      /// Chance to consider options immediately (open cloud, etc)
      FinishOptionParsing();
      return true;
    }
    return false;
  }

  /// Runs Solver given the NL file name
  virtual
  void RunFromNLFile(const std::string& nl_filename,
                     const std::string& filename_no_ext) = 0;

  //////////////////////////// AMPLS //////////////////////////////
  /////////////// Detailed steps for AMPLS C API //////////////////

  /// Read NL.
  /// This is also used by the AMPLS C API.
  /// @param opts: 0-terminated list of extra options,
  /// to be read after the <solver>_options nev var.
  /// All model-related options should be here
  /// (obj:.../objno/multiobj, cvt:..., acc:...),
  /// they are not effective after NL input.
  virtual void ReadNL(const std::string& nl_filename,
                      const std::string& filename_no_ext,
                      char** opts = nullptr) = 0;

  /// Input warm start, suffixes, and all that can modify the model.
  /// This is used by the AMPLS C API
  virtual void InputExtras() = 0;

  /// Report results.
  /// This is used by the AMPLS C API
  virtual void ReportResults() = 0;

  /// Chance for the Backend to init solver environment, etc
  virtual void InitOptionParsing() { }
  /// Chance to consider options immediately (open cloud, etc)
  virtual void FinishOptionParsing() { }
  /// Having everything set up, solve the problem
  virtual void Solve() { }

  /// Callbacks typedef
  using Callbacks = CCallbacks;

  /// Obtain callbacks
  Callbacks& GetCallbacks() { return callbacks_; }
  void OverrideSolutionFile(const std::string& solFile) {
    solutionfileoverride_ = solFile;
  }
  std::string GetOverridenSolutionFile() { return solutionfileoverride_; }

  /// Can be overridden.
  /// For example, if we know the output name,
  /// this should write .sol file with the solve_result and msg.
  virtual void ReportError(int solve_result, fmt::CStringRef msg) = 0;


protected:
  /// 0-terminated list of custom options
  char** GetArgvOptions() const { return argv_options_; }

  /// option parser flags
  int GetOptionFlags() const { return flag_options_; }


private:
  Callbacks callbacks_;
  std::string solutionfileoverride_;

  char** argv_options_ {nullptr};   // 0-terminated list of options
  int flag_options_ {NO_OPTION_ECHO};
};

} // namespace mp

#endif // BACKEND_BASE_H
