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

template <typename Interface>
class InterfaceApp {
 private:
  int result_code_ = 0;
  unsigned banner_size = 0;

  Interface interface_;

  std::string nl_filename, filename_no_ext;
  typename Interface::Converter::NLReadResult nl_read_result;

  std::unique_ptr<internal::SignalHandler> p_sig_handler_;
  std::unique_ptr<internal::SolverAppOptionParser> p_option_parser_;

 protected:
  int GetResultCode() const { return result_code_; }
  const Interface& GetInterface() const { return interface_; }
  Interface& GetInterface() { return interface_; }

  using ProblemBuilder = typename Interface::ModelType;
  const ProblemBuilder& GetProblemBuilder() const
  { return GetInterface().GetInputModel(); }
  ProblemBuilder& GetProblemBuilder() { return GetInterface().GetInputModel(); }

  using OutputModel = typename Interface::OutputModelType;
  const OutputModel& GetOutputModel() const
  { return GetInterface().GetOutputModel(); }
  OutputModel& GetOutputModel() { return GetInterface().GetOutputModel(); }

  using Backend = typename Interface::BackendType;
  const Backend& GetBackend() const { return GetInterface().GetBackend(); }
  Backend& GetBackend() { return GetInterface().GetBackend(); }

  using MPUtils = typename Interface::MPUtils;
  const MPUtils& GetMPUtils() const { return GetInterface().GetMPUtils(); }
  MPUtils& GetMPUtils() { return GetInterface().GetMPUtils(); }

 private:
  struct AppOutputHandler : OutputHandler {
    bool has_output;
    AppOutputHandler() : has_output(false) {}
    void HandleOutput(fmt::CStringRef output) {
      has_output = true;
      std::fputs(output.c_str(), stdout);
    }
  };
  AppOutputHandler output_handler_;

 public:
  InterfaceApp() {
    p_sig_handler_ = std::unique_ptr<internal::SignalHandler>(
            new internal::SignalHandler(GetMPUtils()));
    p_option_parser_ = std::unique_ptr<internal::SolverAppOptionParser>(
        new internal::SolverAppOptionParser(GetMPUtils()));
    GetMPUtils().set_output_handler(&output_handler_);
  }


  // Runs the application.
  // It processes command-line arguments and, if the file name (stub) is
  // specified, reads a problem from an .nl file, parses solver options
  // from argv and environment variables, solves the problem and writes
  // solution(s).
  // argv: an array of command-line arguments terminated by a null pointer
  int RunFromNLFile(char **argv, int nl_reader_flags = 0);

protected:
  bool Init(char** argv, int nl_reader_flags);
  void ReadNL(int nl_reader_flags);
  void Solve();

  void Resolve();      // assuming the Solver has the corr. method
};

template <typename Interface>
int InterfaceApp<Interface>::RunFromNLFile(char **argv, int nl_reader_flags) {
  if (!Init(argv, nl_reader_flags))
    return result_code_;
  ReadNL(nl_reader_flags);
  Solve();
  return 0;
}

template <typename Interface>
bool InterfaceApp<Interface>::Init(char **argv, int nl_reader_flags) {
  internal::Unused(nl_reader_flags);

  // Parse command-line arguments.
  const char *filename = p_option_parser_->Parse(argv);
  if (!filename) return false;

  if (GetMPUtils().ampl_flag()) {
    fmt::MemoryWriter banner;
    banner.write("{}: ", GetMPUtils().long_name());
    std::fputs(banner.c_str(), stdout);
    std::fflush(stdout);
    banner_size = static_cast<unsigned>(banner.size());
    output_handler_.has_output = false;
  }
  // TODO: test output

  // Add .nl extension if necessary.
  nl_filename = filename;
  filename_no_ext = nl_filename;
  const char *ext = std::strrchr(filename, '.');
  if (!ext || std::strcmp(ext, ".nl") != 0)
    nl_filename += ".nl";
  else
    filename_no_ext.resize(filename_no_ext.size() - 3);
  internal::SetBasename(interface_, &filename_no_ext);

  // Parse solver options.
  unsigned flags =
      p_option_parser_->echo_solver_options() ? 0 : Solver::NO_OPTION_ECHO;
  if (!interface_.ParseOptions(argv, flags)) {
    result_code_ = 1;
    return false;
  }

  return true;
}

template <typename Interface>
void InterfaceApp<Interface>::ReadNL(int nl_reader_flags) {
  steady_clock::time_point start = steady_clock::now();

  nl_read_result = interface_.ReadNLFileAndUpdate(nl_filename, nl_reader_flags);

  double read_time = GetTimeAndReset(start);
  if (GetMPUtils().timing())      // TODO why print via backend? We are the app!
    GetMPUtils().Print("Input time = {:.6f}s\n", read_time);
}

template <typename Interface>
void InterfaceApp<Interface>::Solve() {
  ArrayRef<int> options(nl_read_result.handler_->options(),
                        nl_read_result.handler_->num_options());
  internal::AppSolutionHandler<MPUtils> sol_handler(
        filename_no_ext, GetMPUtils(), GetOutputModel(), options,
        output_handler_.has_output ? 0 : banner_size);
  GetInterface().Solve(sol_handler);
}

template <typename Interface>
void InterfaceApp<Interface>::Resolve() {
  ArrayRef<int> options(nl_read_result.handler_->options(),
                        nl_read_result.handler_->num_options());
  internal::AppSolutionHandler<Interface> sol_handler(
        filename_no_ext, GetInterface(), GetOutputModel(), options,
        output_handler_.has_output ? 0 : banner_size);
  throw 0;        // DON't USE THIS ONE
  GetInterface().Resolve(sol_handler);
}

}  // namespace mp


#endif  // INTERFACE_APP_H_
