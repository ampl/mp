#ifndef BACKEND_BASE_H
#define BACKEND_BASE_H

#include <string>

#include "mp/solver-base.h"

namespace mp {

/// Abstract backend API
class BasicBackend : public BasicSolver {
public:
  virtual ~BasicBackend() = default;

  /// Constructs a BasicBackend object.
  /// date:  The solver date in YYYYMMDD format.
  /// flags: Bitwise OR of zero or more of the following values
  ///          MULTIPLE_SOL
  ///          MULTIPLE_OBJ
  BasicBackend(fmt::CStringRef name, fmt::CStringRef long_name,
               long date, int flags) :
    BasicSolver(name, long_name, date , flags) { }

  /// Initialize backend, incl. solver options.
  /// Not parse options yet, as we need to parse cmdline arguments first.
  /// @param argv: the command-line arguments, 0-terminated
  virtual void Init(char** argv) {
    set_exe_path(*argv);
  }

  /// Parse solver options such as "outlev=1" from env and argv.
  /// @param filename_no_ext: basname of the .nl file
  ///        (or whatever should be the basename of an output .sol file)
  /// @param argv: (remaining part of) vector of cmdline strings
  /// @param flags: 0 or \a Solver::NO_OPTION_ECHO
  virtual
  bool ParseSolverOptions(const char* filename_no_ext,
                    char **argv, unsigned flags = 0) {
    /// Chance for the Backend to note base IO filename
    SetBasename(filename_no_ext);
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


protected:
  /// Chance for the Backend to note base IO filename
  virtual void SetBasename(const std::string& ) { }
  /// Chance for the Backend to init solver environment, etc
  virtual void InitOptionParsing() { }
  /// Chance to consider options immediately (open cloud, etc)
  virtual void FinishOptionParsing() { }
};

} // namespace mp

#endif // BACKEND_BASE_H
