/*
 A mathematical optimization app using abstract interfaces.

 Copyright (C) 2020 AMPL Optimization Inc

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

 Author: Gleb Belov <gleb.belov@monash.edu>
 */

#ifndef INTERFACE_APP_H_
#define INTERFACE_APP_H_

#include <string>
#include <memory>

#include "mp/solver.h"   // for namespace internal::

namespace mp {

/// NL solver app.
/// Reads command-line parameters like -AMPL and NL filename,
/// handles screen output and signals,
/// and calls solver class
template <typename NLSolver>
class NLSolverApp {
 protected:
  int GetResultCode() const { return result_code_; }
  const NLSolver& GetNLSolver() const { return nlsolver_; }
  NLSolver& GetNLSolver() { return nlsolver_; }

  /// TODO use Env
  using MPUtils = typename NLSolver::MPUtils;
  const MPUtils& GetMPUtils() const { return GetNLSolver().GetMPUtils(); }
  MPUtils& GetMPUtils() { return GetNLSolver().GetMPUtils(); }

private:
  int result_code_ = 0;
  unsigned banner_size_ = 0;

  NLSolver nlsolver_;

  std::string nl_filename_, filename_no_ext_;

  std::unique_ptr<internal::SignalHandler> p_sig_handler_;
  std::unique_ptr<internal::SolverAppOptionParser> p_option_parser_;

private:
  OutputHandler output_handler_;

 public:
  NLSolverApp() {
    p_sig_handler_ = std::unique_ptr<internal::SignalHandler>(
            new internal::SignalHandler(GetMPUtils()));
    p_option_parser_ = std::unique_ptr<internal::SolverAppOptionParser>(
        new internal::SolverAppOptionParser(GetMPUtils()));
    GetMPUtils().set_output_handler(&output_handler_);
  }

  /// Runs the application.
  /// It processes command-line arguments and, if the file name (stub) is
  /// specified, calls NLSolver class.
  /// @param argv: an array of command-line arguments terminated by a null pointer
  /// @param nl_reader_flags: currently can be 0 or \a READ_BOUNDS_FIRST, see nl.h
  /// @return app exit code
  int Run(char **argv, int nl_reader_flags = 0);

protected:
  bool Init(char** argv);
};

template <typename Interface>
int NLSolverApp<Interface>::Run(char **argv, int nl_reader_flags) {
  if (!Init(argv))
    return result_code_;
  GetNLSolver().RunFromNLFile(
        nl_filename_, filename_no_ext_, nl_reader_flags);
  return 0;
}

template <typename Interface>
bool NLSolverApp<Interface>::Init(char **argv) {
  /// Init solver/converter options
  GetNLSolver().InitOptions();

  // Parse command-line arguments.
  const char *filename = p_option_parser_->Parse(argv);
  if (!filename) return false;

  if (GetMPUtils().ampl_flag()) {
    fmt::MemoryWriter banner;
    banner.write("{}: ", GetMPUtils().long_name());
    std::fputs(banner.c_str(), stdout);
    std::fflush(stdout);
    output_handler_.banner_size = banner.size();
    output_handler_.has_output = false;
  }
  // TODO: test output

  // Add .nl extension if necessary.
  nl_filename_ = filename;
  filename_no_ext_ = nl_filename_;
  const char *ext = std::strrchr(filename, '.');
  if (!ext || std::strcmp(ext, ".nl") != 0)
    nl_filename_ += ".nl";
  else
    filename_no_ext_.resize(filename_no_ext_.size() - 3);
  internal::SetBasename(nlsolver_, &filename_no_ext_);

  // Parse solver options.
  unsigned flags =
      p_option_parser_->echo_solver_options() ? 0 : Solver::NO_OPTION_ECHO;
  if (!nlsolver_.ParseSolverOptions(filename_no_ext_.c_str(), argv, flags)) {
    result_code_ = 1;
    return false;
  }

  return true;
}

}  // namespace mp


#endif  // INTERFACE_APP_H_
