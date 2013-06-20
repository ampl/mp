/*
 Utilities for writing AMPL solvers.

 Copyright (C) 2012 AMPL Optimization LLC

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization LLC disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Victor Zverovich
 */

#ifndef SOLVERS_UTIL_SOLVER_H_
#define SOLVERS_UTIL_SOLVER_H_

#include <csignal>
#include <cstring>

#include <limits>
#include <memory>

extern "C" {
#include "solvers/getstub.h"
#undef Char
}

#include "solvers/util/problem.h"

namespace ampl {

// Formats a double with objective precision.
// Usage: fmt::Format("{}") << ObjPrec(42.0);
class ObjPrec {
 private:
  double value_;

 public:
  explicit ObjPrec(double value) : value_(value) {}

  friend void Format(fmt::Writer &w, const fmt::FormatSpec &spec, ObjPrec op) {
    char buffer[32];
    g_fmtop(buffer, op.value_);
    w.Write(buffer, spec);
  }
};

namespace internal {
// Formats a string by indenting it and performing word wrap.
std::string Format(fmt::StringRef s, int indent = 0);
}

class BasicSolver;

// An interface for receiving errors reported via BasicSolver::ReportError.
class ErrorHandler {
 protected:
  ~ErrorHandler() {}

 public:
  virtual void HandleError(fmt::StringRef message) = 0;
};

// An interface for receiving solutions.
class SolutionHandler {
 protected:
  ~SolutionHandler() {}

 public:
  virtual void HandleSolution(BasicSolver &s, fmt::StringRef message,
      const double *values, const double *dual_values, double obj_value) = 0;
};

class Interruptible {
 protected:
  ~Interruptible() {}

 public:
  virtual void Interrupt() = 0;
};

// A signal handler.
// When a solver is run in a terminal it should respond to SIGINT (Ctrl-C)
// by interrupting its execution and returning the best solution found.
// This can be done in two ways. The first one is to check the return value
// of SignalHandler::stop() periodically. The second is to register an object
// that implements the Interruptible interface. In both cases you should
// create a SignalHandler object before solving.
class SignalHandler {
 private:
  static std::string signal_message_;
  static const char *signal_message_ptr_;
  static unsigned signal_message_size_;
  static volatile std::sig_atomic_t stop_;
  static Interruptible *interruptible_;

  static void HandleSigInt(int sig);

 public:
  explicit SignalHandler(const BasicSolver &s, Interruptible *i = 0);

  ~SignalHandler() {
    stop_ = 1;
    interruptible_ = 0;
  }

  // Returns true if the execution should be stopped due to SIGINT.
  static bool stop() { return stop_ != 0; }
};

// An option error.
class OptionError : public Error {
public:
  OptionError(fmt::StringRef message) : Error(message) {}
};

// An exception indicating an invalid option value.
class InvalidOptionValue : public OptionError {
 public:
  template <typename T>
  InvalidOptionValue(fmt::StringRef name, T value)
  : OptionError(fmt::Format("Invalid value {} for option {}")
      << value << name.c_str()) {}
};

// A solver option.
class SolverOption {
 private:
  std::string description_;
  bool is_keyword_;  // true if this is a keyword option not accepting values.

 public:
  SolverOption(const char *description, bool is_keyword = false)
  : description_(description), is_keyword_(is_keyword) {}
  virtual ~SolverOption() {}

  const char *description() const { return description_.c_str(); }
  bool is_keyword() const { return is_keyword_; }

  // Prints the current value of the option.
  // Throws OptionError in case of error.
  virtual void Print(fmt::StringRef name) = 0;

  // Parses and sets the value of the option. Throws InvalidOptionValue
  // if the value is invalid or OptionError in case of another error.
  virtual void Parse(fmt::StringRef name, const char *&str) = 0;
};

template <typename T>
class TypedSolverOption;

template <>
class TypedSolverOption<int> : public SolverOption {
 public:
  typedef int Arg;

  TypedSolverOption(const char *description) : SolverOption(description) {}

  void Print(fmt::StringRef name);
  void Parse(fmt::StringRef name, const char *&s);

  // Returns the option value.
  virtual int GetValue(fmt::StringRef name) = 0;

  // Sets the option value or throws InvalidOptionValue if the value is invalid.
  virtual void SetValue(fmt::StringRef name, int value) = 0;
};

template <>
class TypedSolverOption<double> : public SolverOption {
 public:
  typedef double Arg;

  TypedSolverOption(const char *description) : SolverOption(description) {}

  void Print(fmt::StringRef name);
  void Parse(fmt::StringRef name, const char *&s);

  // Returns the option value.
  virtual double GetValue(fmt::StringRef name) = 0;

  // Sets the option value or throws InvalidOptionValue if the value is invalid.
  virtual void SetValue(fmt::StringRef name, double value) = 0;
};

template <>
class TypedSolverOption<std::string> : public SolverOption {
 public:
  typedef const char *Arg;

  TypedSolverOption(const char *description) : SolverOption(description) {}

  void Print(fmt::StringRef name);
  void Parse(fmt::StringRef name, const char *&s);

  // Returns the option value.
  virtual std::string GetValue(fmt::StringRef name) = 0;

  // Sets the option value or throws InvalidOptionValue if the value is invalid.
  virtual void SetValue(fmt::StringRef name, const char *value) = 0;
};

// Base class for all solver classes.
class BasicSolver
  : private ErrorHandler, private SolutionHandler, private Option_Info {
 private:
  Problem problem_;

  std::string name_;
  std::string long_name_;
  std::string options_var_name_;
  std::string version_;

  bool has_errors_;
  ErrorHandler *error_handler_;
  SolutionHandler *sol_handler_;

  double read_time_;

  typedef std::map<std::string, SolverOption*> OptionMap;
  OptionMap options_;
  keyword cl_option_;  // command-line option '='

  static char *PrintOptionsAndExit(Option_Info *oi, keyword *kw, char *value);

  void HandleError(fmt::StringRef message) {
    std::fputs(message.c_str(), stderr);
    std::fputc('\n', stderr);
  }

  void HandleSolution(BasicSolver &, fmt::StringRef message,
        const double *values, const double *dual_values, double) {
    write_sol_ASL(reinterpret_cast<ASL*>(problem_.asl_),
        const_cast<char*>(message.c_str()), const_cast<double*>(values),
        const_cast<double*>(dual_values), this);
  }

  class ErrorReporter {
   private:
    ErrorHandler *handler_;

   public:
    ErrorReporter(ErrorHandler *h) : handler_(h) {}

    void operator()(const fmt::Formatter &f) const {
      handler_->HandleError(fmt::StringRef(f.c_str(), f.size()));
    }
  };

  // Parses an option string.
  void ParseOptionString(const char *s, unsigned flags);

 protected:
  // Constructs a BasicSolver object.
  // date: The solver date in YYYYMMDD format.
  BasicSolver(fmt::StringRef name, fmt::StringRef long_name, long date);

  void set_long_name(fmt::StringRef name) {
    long_name_ = name;
    bsname = const_cast<char*>(long_name_.c_str());
  }

  void set_version(fmt::StringRef version) {
    version_ = version;
    Option_Info::version = const_cast<char*>(version_.c_str());
  }

  double read_time() const { return read_time_; }

  virtual std::string GetOptionHeader() { return std::string(); }

  typedef std::auto_ptr<SolverOption> SolverOptionPtr;

  void AddOption(const char *name, SolverOptionPtr opt) {
    // First insert the option, then release a pointer to it. Doing the other
    // way around may lead to a memory leak if insertion throws.
    options_[name] = opt.get();
    opt.release();
  }

 public:
  virtual ~BasicSolver();

  // Flags for ParseOptions.
  enum {
    // Don't echo options during parsing.
    NO_OPTION_ECHO = 1
  };

  // Returns the current problem.
  Problem &problem() { return problem_; }

  // Returns the solver name.
  const char *name() const { return sname; }

  // Returns the long solver name.
  // This name is used in startup "banner".
  const char *long_name() const { return bsname; }

  // Returns the name of <solver>_options environment variable.
  const char *options_var_name() const { return opname; }

  // Returns the solver version.
  const char *version() const { return Option_Info::version; }

  // Returns the solver date in YYYYMMDD format.
  long date() const { return driver_date; }

  // Returns the solver flags.
  // Possible values that can be combined with bitwise OR:
  //   ASL_OI_want_funcadd
  //   ASL_OI_keep_underscores
  //   ASL_OI_show_version
  int flags() const { return Option_Info::flags; }

  // Returns the value of the wantsol option which specifies what solution
  // information to write in a stand-alone invocation (no -AMPL on the
  // command line).  Possible values that can be combined with bitwise OR:
  //   1 = write .sol file
  //   2 = primal variables to stdout
  //   4 = dual variables to stdout
  //   8 = suppress solution message
  int wantsol() const { return Option_Info::wantsol; }

  // Returns the error handler.
  ErrorHandler *error_handler() { return error_handler_; }

  // Sets the error handler.
  void set_error_handler(ErrorHandler *eh) { error_handler_ = eh; }

  // Returns the solution handler.
  SolutionHandler *solution_handler() { return sol_handler_; }

  // Sets the solution handler.
  void set_solution_handler(SolutionHandler *sh) { sol_handler_ = sh; }

  // Processes command-line arguments, reads a problem from an .nl file
  // if the file name (stub) is specified and parses solver options.
  // Returns true if the arguments contain the file name and options were
  // parsed successfully; false otherwise.
  // If there was an  error parsing arguments or reading the problem
  // ProcessArgs will print an error message and call std::exit (this is
  // likely to change in the future version).
  bool ProcessArgs(char **&argv, unsigned flags = 0);

  // Parses solver options and returns true if there were no errors and
  // false otherwise.
  virtual bool ParseOptions(char **argv, unsigned flags = 0);

  // Passes a solution to the solution handler.
  void HandleSolution(fmt::StringRef message,
      const double *values, const double *dual_values,
      double obj_value = std::numeric_limits<double>::quiet_NaN()) {
    sol_handler_->HandleSolution(*this, message, values, dual_values, obj_value);
  }

  // Reports an error printing the formatted error message to stderr.
  // Usage: ReportError("File not found: {}") << filename;
  fmt::TempFormatter<ErrorReporter> ReportError(fmt::StringRef format) {
    has_errors_ = true;
    return fmt::TempFormatter<ErrorReporter>(
        format.c_str(), ErrorReporter(error_handler_));
  }

  // Solves a problem.
  // The solutions are reported via the registered solution handler.
  virtual void Solve(Problem &p) = 0;

  // Runs the solver.
  int Run(char **argv);
};

// An AMPL solver.
// Impl should be a class derived from Solver that will receive notifications
// about parsed options.
//
// Example:
//
// class MySolver : public Solver<MySolver> {
//  public:
//   void SetTestOption(const char *name, int value, int info) {
//     // Set option - called when the test option is parsed;
//     // info is an arbitrary value passed as a third argument to
//     // AddIntOption. It can be useful if the same function handles
//     // multiple options.
//     ...
//   }
//
//   MySolver()  {
//     AddIntOption("test", "This is a test option",
//                  &MySolver::SetTestOption, 42);
//   }
// };
template <typename Impl>
class Solver : public BasicSolver {
 private:
  template <typename T>
  class ConcreteOption : public TypedSolverOption<T> {
   private:
    typedef typename TypedSolverOption<T>::Arg Arg;
    typedef T (Impl::*Getter)(const char *);
    typedef void (Impl::*Setter)(const char *, Arg);

    Impl &impl_;
    Getter getter_;
    Setter setter_;

   public:
    ConcreteOption(const char *description,
        Solver *s, Getter getter, Setter setter)
    : TypedSolverOption<T>(description), impl_(static_cast<Impl&>(*s)),
      getter_(getter), setter_(setter) {}

    T GetValue(fmt::StringRef name) {
      return (impl_.*getter_)(name.c_str());
    }

    void SetValue(fmt::StringRef name, Arg value) {
      (impl_.*setter_)(name.c_str(), value);
    }
  };

  template <typename T, typename Info, typename InfoArg = Info>
  class ConcreteOptionWithInfo : public TypedSolverOption<T> {
   private:
    typedef typename TypedSolverOption<T>::Arg Arg;
    typedef T (Impl::*Getter)(const char *, InfoArg);
    typedef void (Impl::*Setter)(const char *, Arg, InfoArg);

    Impl &impl_;
    Getter getter_;
    Setter setter_;
    Info info_;

   public:
    ConcreteOptionWithInfo(const char *description,
        Solver *s, Getter getter, Setter setter, InfoArg info)
    : TypedSolverOption<T>(description), impl_(static_cast<Impl&>(*s)),
      getter_(getter), setter_(setter), info_(info) {}

    T GetValue(fmt::StringRef name) {
      return (impl_.*getter_)(name.c_str(), info_);
    }

    void SetValue(fmt::StringRef name, Arg value) {
      (impl_.*setter_)(name.c_str(), value, info_);
    }
  };

 protected:
  // Adds an integer option. The arguments getter and setter should be
  // pointers to member functions in the solver class. They are used to
  // get and set an option value respectively.
  void AddIntOption(const char *name, const char *description,
      int (Impl::*getter)(const char *),
      void (Impl::*setter)(const char *, int)) {
    AddOption(name, SolverOptionPtr(
        new ConcreteOption<int>(description, this, getter, setter)));
  }

  // Adds an integer option with additional information. The arguments getter
  // and setter should be pointers to member functions in the solver class.
  // They are used to get and set an option value respectively.
  template <typename Info>
  void AddIntOption(const char *name, const char *description,
      int (Impl::*getter)(const char *, const Info &),
      void (Impl::*setter)(const char *, int, const Info &),
      const Info &info) {
    AddOption(name, SolverOptionPtr(
        new ConcreteOptionWithInfo<int, Info, const Info &>(
            description, this, getter, setter, info)));
  }

  // The same as above but with Info argument passed by value.
  template <typename Info>
  void AddIntOption(const char *name, const char *description,
      int (Impl::*getter)(const char *, Info),
      void (Impl::*setter)(const char *, int, Info), Info info) {
    AddOption(name, SolverOptionPtr(new ConcreteOptionWithInfo<int, Info>(
            description, this, getter, setter, info)));
  }

  // Adds a double option. The arguments getter and setter should be
  // pointers to member functions in the solver class. They are used to
  // get and set an option value respectively.
  void AddDblOption(const char *name, const char *description,
      double (Impl::*getter)(const char *),
      void (Impl::*setter)(const char *, double)) {
    AddOption(name, SolverOptionPtr(new ConcreteOption<double>(
        description, this, getter, setter)));
  }

  // Adds a double option with additional information. The arguments getter
  // and setter should be pointers to member functions in the solver class.
  // They are used to get and set an option value respectively.
  template <typename Info>
  void AddDblOption(const char *name, const char *description,
      double (Impl::*getter)(const char *, const Info &),
      void (Impl::*setter)(const char *, double, const Info &),
      const Info &info) {
    AddOption(name, SolverOptionPtr(
        new ConcreteOptionWithInfo<double, Info, const Info &>(
            description, this, getter, setter, info)));
  }

  // The same as above but with Info argument passed by value.
  template <typename Info>
  void AddDblOption(const char *name, const char *description,
      double (Impl::*getter)(const char *, Info),
      void (Impl::*setter)(const char *, double, Info), Info info) {
    AddOption(name, SolverOptionPtr(new ConcreteOptionWithInfo<double, Info>(
            description, this, getter, setter, info)));
  }

  // Adds a string option. The arguments getter and setter should be
  // pointers to member functions in the solver class. They are used to
  // get and set an option value respectively.
  void AddStrOption(const char *name, const char *description,
      std::string (Impl::*getter)(const char *),
      void (Impl::*setter)(const char *, const char *)) {
    AddOption(name, SolverOptionPtr(new ConcreteOption<std::string>(
        description, this, getter, setter)));
  }

  // Adds a string option with additional information. The arguments getter
  // and setter should be pointers to member functions in the solver class.
  // They are used to get and set an option value respectively.
  template <typename Info>
  void AddStrOption(const char *name, const char *description,
      std::string (Impl::*getter)(const char *, const Info &),
      void (Impl::*setter)(const char *, const char *, const Info &),
      const Info &info) {
    AddOption(name, SolverOptionPtr(
        new ConcreteOptionWithInfo<std::string, Info, const Info &>(
            description, this, getter, setter, info)));
  }

  // The same as above but with Info argument passed by value.
  template <typename Info>
  void AddStrOption(const char *name, const char *description,
      std::string (Impl::*getter)(const char *, Info),
      void (Impl::*setter)(const char *, const char *, Info),
      const Info &info) {
    typedef void (Impl::*Func)(const char *, const char *, Info);
    AddOption(name, SolverOptionPtr(
        new ConcreteOptionWithInfo<std::string, Info>(
            description, this, getter, setter, info)));
  }

 public:
  Solver(fmt::StringRef name, fmt::StringRef long_name = 0, long date = 0)
  : BasicSolver(name, long_name, date) {}
};
}

#endif  // SOLVERS_UTIL_SOLVER_H_
