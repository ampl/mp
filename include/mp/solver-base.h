#ifndef SOLVERBASE_H
#define SOLVERBASE_H

#include <unordered_map>

#include "solver-opt.h"

namespace mp {

/// An interface for receiving errors reported via Solver::ReportError.
class ErrorHandler {
 protected:
  virtual ~ErrorHandler() {}

 public:
  virtual void HandleError(fmt::CStringRef message) = 0;
};

/// A default interface for receiving solver output.
class OutputHandler {
 public:
  virtual ~OutputHandler() {}
  bool has_output = false;
  size_t banner_size = 0;
  virtual void HandleOutput(fmt::CStringRef output) {
    has_output = true;
    std::fputs(output.c_str(), stdout);
  }
};

/// An interface for receiving solutions.
class SolutionHandler {
 public:
  virtual ~SolutionHandler() {}

  // Receives a feasible solution.
  virtual void HandleFeasibleSolution(fmt::CStringRef message,
      const double *values, const double *dual_values, double obj_value) = 0;

  // Receives the final solution or a notification that the problem is
  // infeasible or unbounded.
  virtual void HandleSolution(int status, fmt::CStringRef message,
      const double *values, const double *dual_values, double obj_value) = 0;
};

/// "Silent" solution handler
class BasicSolutionHandler : public SolutionHandler {
 public:
  virtual void HandleFeasibleSolution(fmt::CStringRef,
      const double *, const double *, double) {}
  virtual void HandleSolution(int, fmt::CStringRef,
      const double *, const double *, double) {}
};


/// Interrupt handler.
/// Returns true if the solver was interrupted, false if it is not running.
typedef bool (*InterruptHandler)(void *);

/// An interface for interrupting solution process.
/// When a solver is run in a terminal it should respond to SIGINT (Ctrl-C)
/// by interrupting its execution and returning the best solution found.
/// This can be done in two ways. The first one is to check the return value
/// of Interrupter::Stop() periodically. The second is to register an object
/// that implements the Interruptible interface.
class Interrupter {
 public:
  virtual ~Interrupter() {}

  /// Returns true if the solution process should be stopped.
  virtual bool Stop() const = 0;

  /// Sets a handler function.
  /// The handler function must be safe to call from a signal handler.
  /// In particular, it must only call async-signal-safe library functions.
  virtual void SetHandler(InterruptHandler handler, void *data) = 0;
};


class ASLProblem;


/**
  A minimal set of standard features for AMPL solver logic.

  Can be used both for solver drivers, as well as for unit tests.
  Example:

  \rst
  .. code-block:: cpp

       class MySolver : public Solver {
        public:
         void GetTestOption(const char *name, int value, int info) {
           // Returns the option value; info is an arbitrary value passed as
           // the last argument to AddIntOption. It can be useful if the same
           // function handles multiple options.
           ...
         }
         void SetTestOption(const char *name, int value, int info) {
           // Set the option value; info is the same as in GetTestOption.
           ...
         }

         MySolver()  {
           AddIntOption("test", "This is a test option",
                        &MySolver::GetTestOption, &MySolver::SetTestOption, 42);
         }
       }
  \endrst
  */
class BasicSolver : private ErrorHandler,
    private OutputHandler, private Interrupter,
    public SolverOptionManager {
public:
  std::string GetOptionFile(const SolverOption &) const
  { return option_file_save_; }
  void UseOptionFile(const SolverOption &, fmt::StringRef value);

  int GetWantSol(const SolverOption &) const { return wantsol_; }
  void SetWantSol(const SolverOption &, int value) {
    if ((value & ~0xf) != 0)
      throw InvalidOptionValue("wantsol", value);
    wantsol_ = value;
  }

  std::string GetSolutionStub(const SolverOption &) const {
    return solution_stub_;
  }
  void SetSolutionStub(const SolverOption &, fmt::StringRef value) {
    solution_stub_ = value.to_string();
  }

  int GetObjNo(const SolverOption &) const { return std::abs(objno_); }
  void SetObjNo(const SolverOption &opt, int value) {
    if (value < 0)
      throw InvalidOptionValue(opt, value);
    objno_ = value;
  }

  /// Returns the solver name.
  /// This is used to extract solver options from
  /// the env variable <name>_options.
  const char *name() const { return name_.c_str(); }

  /// Returns the long solver name.
  /// This name is used in startup "banner".
  const char *long_name() const { return long_name_.c_str(); }

  /// Returns the executable path.
  const char *exe_path() const { return exe_path_.c_str(); }
  void set_exe_path(const char* p) { exe_path_ = p; }

  /// Returns the solver version.
  const char *version() const { return version_.c_str(); }

  /// Returns the solver date in YYYYMMDD format.
  long date() const { return date_; }

  /// Possible values for the wantsol option (can be combined with bitwise OR).
  enum {
    WRITE_SOL_FILE      = 1,  // write .sol file
    PRINT_SOLUTION      = 2,  // print primal variables to stdout
    PRINT_DUAL_SOLUTION = 4,  // print dual variables to stdout
    SUPPRESS_SOLVER_MSG = 8   // suppress solver message
  };

  /// Returns the value of the wantsol option which specifies what solution
  /// information to write in a stand-alone invocation (no -AMPL on the
  /// command line).
  int wantsol() const { return wantsol_; }
  void set_wantsol(int value) { wantsol_ = value; }

  // Returns true if -AMPL is specified.
  bool ampl_flag() { return (bool_options_ & AMPL_FLAG) != 0; }
  void set_ampl_flag(bool value = true) {
    if (value)
      bool_options_ |= AMPL_FLAG;
    else
      bool_options_ &= ~AMPL_FLAG;
  }

  /// True if verbose mode
  /// Should be set by the implementation,
  /// otherwise it's true
  bool verbose_mode() const { return verbose_; }
  /// Set verbosity, by the impl
  void set_verbose_mode(bool f) { verbose_=f; }

  /// True if need to debug
  /// Outputs test infos etc
  bool debug_mode() const { return debug_; }

  /// Returns the index of the objective to optimize starting from 1,
  /// 0 to not use objective.
  /// Both multiobj and objno are used in NLReader to select
  /// the objective(s), solvers should not use these options.
  int objno() const { return std::abs(objno_); }
  bool objno_specified() const { return objno_>=0; }

  /// Returns true if multiobjective optimization is enabled.
  /// Both multiobj and objno are used in NLReader to select
  /// the objective(s), solvers should not use these options.
  bool multiobj() const { return multiobj_ && objno_<0; }

  /// Returns true if the timing is enabled.
  bool timing() const { return timing_; }

  /// Returns the error handler.
  ErrorHandler *error_handler() { return error_handler_; }

  /// Sets the error handler.
  void set_error_handler(ErrorHandler *eh) { error_handler_ = eh; }

  /// Returns the output handler.
  OutputHandler *output_handler() { return output_handler_; }
  OutputHandler &get_output_handler() { return *output_handler_; }

  /// Sets the output handler.
  void set_output_handler(OutputHandler *oh) { output_handler_ = oh; }

  const Interrupter *interrupter() const { return interrupter_; }
  Interrupter *interrupter() { return interrupter_; }
  void set_interrupter(Interrupter *interrupter) {
    interrupter_ = interrupter ? interrupter : this;
  }

  const char *solution_stub() const { return solution_stub_.c_str(); }

  bool need_multiple_solutions() const {
    return count_solutions_ || !solution_stub_.empty();
  }


  /// Override methods from base service classes
  void HandleOutput(fmt::CStringRef output) override {
    std::fputs(output.c_str(), stdout);
  }

  void HandleError(fmt::CStringRef message) override {
    std::fputs(message.c_str(), stderr);
    std::fputc('\n', stderr);
  }

  /// The default implementation of Interrupter does nothing.
  bool Stop() const override { return false; }
  void SetHandler(InterruptHandler, void *) override {}

  virtual void HandleUnknownOption(const char *name) {
    ReportError("Unknown option or invalid key \"{}\"", name);
  }

  /// Reports an error printing the formatted error message to stderr.
  /// Usage: ReportError("File not found: {}") << filename;
  void ReportError(fmt::CStringRef format, const fmt::ArgList &args) {
    has_errors_ = true;
    fmt::MemoryWriter w;
    w.write(format, args);
    error_handler_->HandleError(w.c_str());
  }
  FMT_VARIADIC(void, ReportError, fmt::CStringRef)

  /// Formats a string and prints it to stdout or, if an output handler
  /// is registered, sends it to the output handler.
  void Print(fmt::CStringRef format, const fmt::ArgList &args) {
    fmt::MemoryWriter w;
    w.write(format, args);
    output_handler_->HandleOutput(w.c_str());
  }
  FMT_VARIADIC(void, Print, fmt::CStringRef)

  /// Add a warning.
  /// *key / *msg are not copied and should remain valid
  void AddWarning(const char* key, const char* msg);

  /// Print warnings
  void PrintWarnings();

  /// Prints version information.
  bool ShowVersion();

  /// Formats a double-prec value
  struct DoubleFormatter {
    double value;
    int precision;

    friend void format(fmt::BasicFormatter<char> &f,
                       const char *, DoubleFormatter df) {
      f.writer().write("{:.{}}", df.value, df.precision);
    }
  };

  /// The default precision used in FormatObjValue if the
  /// \a objective_precision env variable is not specified or 0.
  enum { DEFAULT_PRECISION = 15 };

  /// Returns a formatter that writes value using objective precision.
  /// Usage:
  ///   Print("objective {}", FormatObjValue(obj_value));
  DoubleFormatter FormatObjValue(double value);

  friend class ASLProblem;

  /// Parses solver options from the (solver)_options env variable
  /// and/or strings in argv. For example,
  /// solver_options='timelim=300' solver /tmp/diet.nl -AMPL "basis=0".
  /// @return true if there were no errors and
  /// false otherwise. It accepts a pointer to the problem because some
  /// options may depend on problem features.
  virtual bool ParseOptions(
      char **argv, unsigned flags = 0, const ASLProblem *p = 0);

  /// Solver flags
  enum {
    /// Multiple solutions support.
    /// Makes Solver register "countsolutions" and "solutionstub" options
    /// and write every solution passed to HandleFeastibleSolution to a file
    /// solutionstub & i & ".sol" where i is a solution number.
    MULTIPLE_SOL = 1,

    /// Multiple objectives support.
    /// Makes Solver register the "multiobj" option
    MULTIPLE_OBJ = 2
  };

  /// Testing constructor
  explicit BasicSolver();

protected:
  /// Constructs a BasicSolver object.
  /// date:  The solver date in YYYYMMDD format.
  /// flags: Bitwise OR of zero or more of the following values
  ///          MULTIPLE_SOL
  ///          MULTIPLE_OBJ
  BasicSolver(fmt::CStringRef name, fmt::CStringRef long_name,
              long date, int flags);

  void set_long_name(fmt::StringRef name)
  { long_name_ = name.to_string(); }
  void add_to_long_name(fmt::StringRef name)
  { long_name_ += name.to_string(); }
  void set_version(fmt::StringRef version)
  { version_ = version.to_string(); }
  void add_to_version(fmt::StringRef version)
  { version_ += version.to_string(); }

  /// Sets the flags for Problem::Read.
  void set_read_flags(unsigned flags) { read_flags_ = flags; }

  /// Parses a solver option string.
  void ParseOptionString(const char *s, unsigned flags);


  /// Map to count warnings by types.
  /// Stores char* to names / descriptions for speed
  /// Indexed by pointers to names, not values
  using WarningsMap =
    std::unordered_map< const char*,    // failure name
      std::pair<int, const char*> >;    // number, description
  /// Get warnings map
  WarningsMap& GetWarnings() { return warnings_; }


private:
  std::string name_;
  std::string long_name_;
  std::string exe_path_;
  std::string version_;
  std::string license_info_;
  long date_;
  int wantsol_;
  int obj_precision_;

  /// Index of the objective to optimize starting from 1, 0 to ignore
  /// objective, or -1 to use the first objective if there is one.
  int objno_;

  enum {SHOW_VERSION = 1, AMPL_FLAG = 2};
  int bool_options_;
  int option_flag_save_ = 0;
  std::string option_file_save_;

  // The filename stub for returning multiple solutions.
  std::string solution_stub_;

  // Specifies whether to return the number of solutions in the .nsol suffix.
  bool count_solutions_;

  unsigned read_flags_;  // flags passed to Problem::Read

  bool verbose_=true;
  bool debug_;
  bool timing_;
  bool multiobj_;

  bool has_errors_;
  OutputHandler *output_handler_;
  ErrorHandler *error_handler_;
  Interrupter *interrupter_;

  /// Warnings
  WarningsMap warnings_;
};


namespace internal {

template<typename T, T> struct Check {
  Check(int) {}
};

template <typename Solver>
inline void SetBasename(
    Solver &s, const std::string *basename,
    Check<void (Solver::*)(fmt::StringRef), &Solver::set_basename> = 0) {
  s.set_basename(*basename);
}

template <typename Solver>
inline void SetBasename(Solver &, ...) {}

} // namespace internal


} // namespace mp

#endif // SOLVERBASE_H
