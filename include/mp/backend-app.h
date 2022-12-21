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
#include <functional>
#include <memory>

#include "mp/solver-app-base.h"
#include "mp/backend-base.h"
#include "mp/error.h"

namespace mp {

/// Backend app.
/// Reads command-line parameters like -AMPL and NL filename,
/// installs handlers for screen output and signals,
/// and calls backend class
class BackendApp {
public:
  BackendApp(std::unique_ptr<BasicBackend> pb)
    : pb_(std::move(pb)) { InitHandlers(); }

  /// Runs the application.
  /// It processes command-line arguments and, if the file name (stub) is
  /// specified, calls NLSolver class.
  /// @param argv: an array of command-line arguments terminated by a null pointer
  /// @return app exit code
  virtual int Run(char **argv);

  /// Get result code
  int GetResultCode() const { return result_code_; }

  /// Get backend, const
  const BasicBackend& GetBackend() const { return *pb_; }

  /// Get backend
  BasicBackend& GetBackend() { return *pb_; }


protected:
  virtual bool Init(char** argv);
  virtual void InitHandlers();


private:
  std::unique_ptr<BasicBackend> pb_;
  int result_code_ = 0;

  std::string nl_filename_, filename_no_ext_;

  std::unique_ptr<internal::SignalHandler> p_sig_handler_;
  std::unique_ptr<internal::SolverAppOptionParser> p_option_parser_;

  OutputHandler output_handler_;
};

int BackendApp::Run(char **argv) {
  if (!Init(argv))
    return result_code_;
  GetBackend().RunFromNLFile(
        nl_filename_, filename_no_ext_);
  return 0;
}

bool BackendApp::Init(char **argv) {
  /// Init solver/converter options
  GetBackend().Init(argv);

  // Parse command-line arguments.
  const char *filename = p_option_parser_->Parse(argv);
  if (!filename) return false;

  if (GetBackend().ampl_flag()) {
    fmt::MemoryWriter banner;
    banner.write("{}: ", GetBackend().long_name());
    std::fputs(banner.c_str(), stdout);
    std::fflush(stdout);
    output_handler_.banner_size = banner.size();
    output_handler_.has_output = false;
  }

  // Add .nl extension if necessary.
  nl_filename_ = filename;
  filename_no_ext_ = nl_filename_;
  const char *ext = std::strrchr(filename, '.');
  if (!ext || std::strcmp(ext, ".nl") != 0)
    nl_filename_ += ".nl";
  else
    filename_no_ext_.resize(filename_no_ext_.size() - 3);
  internal::SetBasename(GetBackend(), &filename_no_ext_);

  // Parse solver options.
  unsigned flags =
      p_option_parser_->echo_solver_options() ?
        0 : BasicSolver::NO_OPTION_ECHO;
  if (!GetBackend().ParseSolverOptions(
        filename_no_ext_.c_str(), argv, flags)) {
    result_code_ = 1;
    return false;
  }

  return true;
}

void BackendApp::InitHandlers() {
  p_sig_handler_ = std::unique_ptr<internal::SignalHandler>(
        new internal::SignalHandler(GetBackend()));
  p_option_parser_ = std::unique_ptr<internal::SolverAppOptionParser>(
        new internal::SolverAppOptionParser(GetBackend()));
  GetBackend().set_output_handler(&output_handler_);
}


/// Create and run a backend.
/// @param argv: command-line arguments, null-terminated vector
/// @param be_creator: backend factory
/// @param callbacks: special callbacks
///
/// @return application exit value
inline
int RunBackendApp(char** argv,
                        std::function<
                          std::unique_ptr<mp::BasicBackend>() > be_creator,
                        BasicBackend::Callbacks callbacks = {}) {
  try {
    mp::BackendApp s(be_creator());
    s.GetBackend().GetCallbacks() = callbacks;
    return s.Run(argv);
  } catch (const mp::Error &e) {   // Why does MSVC 2019 or 2022 go to std::exception and not here?
    fmt::print(stderr, "Error: {}\n", e.what());
    return e.exit_code();                // can deliver exit code
  } catch (const std::exception &e) {
    fmt::print(stderr, "Error: {}\n", e.what());
    return EXIT_FAILURE;
  }
  return 0;
}

}  // namespace mp

#endif  // INTERFACE_APP_H_
