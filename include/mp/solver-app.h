#ifndef SOLVERAPP_H
#define SOLVERAPP_H

#include "mp/solver-app-base.h"
#include "mp/solver-io.h"
#include "mp/clock.h"

namespace mp {

/// DEPRECATED.
///
/// A solver application with built-in NL file input.
/// Solver: optimization solver class; normally a subclass of SolverImpl
/// Reader: .nl reader
template <typename Solver, typename Reader = internal::NLFileReader<> >
class SolverApp : private Reader {
 private:
  int result_code_ = 0;
  unsigned banner_size = 0;

  Solver solver_;

  std::string nl_filename, filename_no_ext;
  typedef typename Solver::ProblemBuilder ProblemBuilder;
  std::unique_ptr<ProblemBuilder> builder_;
  std::unique_ptr< internal::SolverNLHandler<Solver> > handler_;

  internal::SignalHandler sig_handler;
  internal::SolverAppOptionParser option_parser_;

 protected:
  int GetResultCode() const { return result_code_; }
  Solver& GetSolver() { return solver_; }
  ProblemBuilder& GetProblemBuilder() { return *builder_; }

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
  SolverApp() :
    sig_handler(solver_),
    option_parser_(solver_) {
    solver_.set_output_handler(&output_handler_);
  }

  // Returns the list of command-line options.
  OptionList &options() { return option_parser_.options(); }

  // Returns the solver.
  Solver &solver() { return solver_; }

  // Returns the .nl reader.
  Reader &reader() { return *this; }

  // Runs the application.
  // It processes command-line arguments and, if the file name (stub) is
  // specified, reads a problem from an .nl file, parses solver options
  // from argv and environment variables, solves the problem and writes
  // solution(s).
  // argv: an array of command-line arguments terminated by a null pointer
  int Run(char **argv, int nl_reader_flags = 0);

protected:
  bool Init(char** argv, int nl_reader_flags);
  void ReadNL(int nl_reader_flags);
  void Solve();

  /// Methods for incremental interface
  void Resolve();      // assuming the Solver has the corr. method
};

template <typename Solver, typename Reader>
int SolverApp<Solver, Reader>::Run(char **argv, int nl_reader_flags) {
  if (!Init(argv, nl_reader_flags))
    return result_code_;
  ReadNL(nl_reader_flags);
  Solve();
  return 0;
}

template <typename Solver, typename Reader>
bool SolverApp<Solver, Reader>::Init(char **argv, int nl_reader_flags) {
  internal::Unused(nl_reader_flags);

  // Parse command-line arguments.
  const char *filename = option_parser_.Parse(argv);
  if (!filename) return false;

  if (solver_.ampl_flag()) {
    fmt::MemoryWriter banner;
    banner.write("{}: ", solver_.long_name());
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
  internal::SetBasename(solver_, &filename_no_ext);

  // Parse solver options.
  unsigned flags =
      option_parser_.echo_solver_options() ? 0 : Solver::NO_OPTION_ECHO;
  if (!solver_.ParseOptions(argv, flags)) {
    result_code_ = 1;
    return false;
  }

  return true;
}

template <typename Solver, typename Reader>
void SolverApp<Solver, Reader>::ReadNL(int nl_reader_flags) {
  steady_clock::time_point start = steady_clock::now();

  builder_.reset(new ProblemBuilder(solver_));
  handler_.reset(new internal::SolverNLHandler<Solver>(*builder_, solver_));
  this->Read(nl_filename, *handler_, nl_reader_flags);

  double read_time = GetTimeAndReset(start);
  if (solver_.timing())
    solver_.Print("Input time = {:.6f}s\n", read_time);
}

template <typename Solver, typename Reader>
void SolverApp<Solver, Reader>::Solve() {
  ArrayRef<int> options(handler_->options(), handler_->num_options());
  internal::AppSolutionHandler<Solver> sol_handler(
        filename_no_ext, solver_, *builder_, options,
        output_handler_.has_output ? 0 : banner_size);
  solver_.Solve(builder_->problem(), sol_handler);
}

template <typename Solver, typename Reader>
void SolverApp<Solver, Reader>::Resolve() {
  ArrayRef<int> options(handler_->options(), handler_->num_options());
  internal::AppSolutionHandler<Solver> sol_handler(
        filename_no_ext, solver_, *builder_, options,
        output_handler_.has_output ? 0 : banner_size);
  solver_.Resolve(builder_->problem(), sol_handler);
}

} // namespace mp

#endif // SOLVERAPP_H
