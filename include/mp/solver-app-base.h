#ifndef SOLVERAPPBASE_H
#define SOLVERAPPBASE_H

#include <string>
#include <csignal>

#if MP_USE_ATOMIC
# include <atomic>
#endif

#include "mp/option.h"
#include "mp/solver-base.h"


namespace mp {

namespace internal {

/// Command-line option parser for a solver application.
/// Not to be mistaken with solver option parser
/// built into the Solver class.
class SolverAppOptionParser {
 private:

  std::string currentOptionString;
  BasicSolver &solver_;

  /// Command-line options.
  OptionList options_;

  bool echo_solver_options_;

  /// Prints usage information and stops processing options.
  bool ShowUsage();

  /// Prints information about solver options.
  bool ShowSolverOptions();
  bool ShowSolverOptionsASL();


  bool WantSol() {
    solver_.set_wantsol(1);
    return true;
  }

  bool DontEchoSolverOptions() {
    echo_solver_options_ = false;
    return true;
  }

  /// Stops processing options.
  bool EndOptions() { return false; }

 public:
  explicit SolverAppOptionParser(BasicSolver &s);

  OptionList &options() { return options_; }

  /// Retruns true if assignments of solver options should be echoed.
  bool echo_solver_options() const { return echo_solver_options_; }

  /// Parses command-line options.
  const char *Parse(char **&argv);
};

#if MP_USE_ATOMIC
using std::atomic;
#else
// Dummy "atomic" for compatibility with legacy compilers.
template <typename T>
class atomic {
 private:
  T value_;

 public:
  atomic(T value = T()) : value_(value) {}
  operator T() const { return value_; }
};
#endif

#ifdef _WIN32
// Signal repeater used to pass signals across processes on Windows.
class SignalRepeater {
 private:
  fmt::ULongLong in_;
  fmt::ULongLong out_;

 public:
  // s: String in the "<int>,<int>" with integers representing the
  //    handles for the input and output ends of the pipe.
  explicit SignalRepeater(const char *s);

  fmt::ULongLong in() const { return in_; }
  fmt::ULongLong out() const { return out_; }
};
#else
struct SignalRepeater {
  explicit SignalRepeater(const char *) {}
};
#endif

// A SIGINT handler
class SignalHandler : public Interrupter {
 private:
  BasicSolver &solver_;
  std::string message_;
  static volatile std::sig_atomic_t stop_;

  static atomic<const char*> signal_message_ptr_;
  static atomic<unsigned> signal_message_size_;
  static atomic<InterruptHandler> handler_;
  static atomic<void*> data_;

  static void HandleSigInt(int sig);

  SignalRepeater repeater_;

 public:
  explicit SignalHandler(BasicSolver &s);

  ~SignalHandler();

  //static bool stop() { return stop_ != 0; }

  // Returns true if the execution should be stopped due to SIGINT.
  bool Stop() const { return stop_ != 0; }

  void SetHandler(InterruptHandler handler, void *data);
};

} // namespace internal


} // namespace mp

#endif // SOLVERAPPBASE_H
