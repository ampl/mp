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

namespace internal {

template <typename T>
struct OptionParser;

template <>
struct OptionParser<int> {
  int operator()(Option_Info *oi, keyword *kw, char *&s);
};

template <>
struct OptionParser<double> {
  double operator()(Option_Info *oi, keyword *kw, char *&s);
};

template <>
class OptionParser<const char*> {
 private:
  std::string value_;

 public:
  const char* operator()(Option_Info *, keyword *, char *&s);
};

// Formats a string by indenting it and performing word wrap.
std::string Format(fmt::StringRef s, int indent = 0);
}

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

// A solver option.
class SolverOption {
 private:
  std::string description_;
  bool is_keyword_;  // true if this is a keyword option not accepting values.

 public:
  explicit SolverOption(const char *description, bool is_keyword = false)
  : description_(description), is_keyword_(is_keyword) {}

  virtual ~SolverOption() {}

  const char *description() const { return description_.c_str(); }

  bool is_keyword() const { return is_keyword_; }

  // Handles the option and returns true if there were no errors;
  // false otherwise.
  virtual bool Handle(BasicSolver &s, keyword *kw, char *&value) = 0;
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

  std::vector<keyword> cl_options_;  // command-line options
  bool options_sorted_;

  bool has_errors_;
  ErrorHandler *error_handler_;
  SolutionHandler *sol_handler_;

  double read_time_;

  typedef std::map<std::string, SolverOption*> OptionMap;
  OptionMap options_;

  void SortOptions();

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

  typedef bool (BasicSolver::*HandleOption)(keyword *kw, char *&value);

  class BasicOption : public SolverOption {
   private:
    HandleOption handle_;

   public:
    BasicOption(const char *description, bool is_keyword, HandleOption h)
    : SolverOption(description, is_keyword), handle_(h) {}

    virtual bool Handle(BasicSolver &s, keyword *kw, char *&value) {
      return (s.*handle_)(kw, value);
    }
  };

  void AddBasicOption(const char *name,
      const char *description, bool is_keyword, HandleOption handle) {
    AddOption(name, std::auto_ptr<SolverOption>(
        new BasicOption(description, is_keyword, handle)));
  }

  bool PrintVersion(keyword *, char *&) {
    show_version_ASL(this);
    return true;
  }

  bool SetWantSol(keyword *kw, char *&value) {
    Option_Info oi = {};
    Option_Info::wantsol = internal::OptionParser<int>()(&oi, kw, value);
    return oi.n_badopts == 0;
  }

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

  void AddOption(const char *name, std::auto_ptr<SolverOption> opt) {
    // First insert the option, then release a pointer to it. Doing the other
    // way around may lead to a memory leak if insertion throws.
    options_[name] = opt.get();
    opt.release();
  }

  template <typename T>
  void ReportInvalidOptionValue(const char *name, T value) {
    ReportError("Invalid value {} for option {}") << value << name;
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
  // false otherwise. Note that handler functions can report errors with
  // BasicSolver::ReportError and ParseOptions will take them into account
  // as well, returning false if there was at least one such error.
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
  template <typename Value, typename Result = Value>
  class ConcreteOption : public SolverOption {
   private:
    typedef Result (Impl::*Getter)(const char *);
    typedef void (Impl::*Setter)(const char *, Value);

    Getter getter_;
    Setter setter_;

   public:
    ConcreteOption(const char *description, Getter getter, Setter setter)
    : SolverOption(description), getter_(getter), setter_(setter) {}

    bool Handle(BasicSolver &s, keyword *kw, char *&value) {
      Option_Info oi = {};
      oi.eqsign = "=";
      (static_cast<Impl&>(s).*setter_)(
          kw->name, internal::OptionParser<Value>()(&oi, kw, value));
      return oi.n_badopts == 0;
    }
  };

  template <typename Func, typename Info, typename Value>
  class ConcreteOptionWithInfo : public SolverOption {
   private:
    Func func_;
    Info info_;

   public:
    ConcreteOptionWithInfo(const char *description, Func func, const Info &info)
    : SolverOption(description), func_(func), info_(info) {}

    bool Handle(BasicSolver &s, keyword *kw, char *&value) {
      Option_Info oi = {};
      oi.eqsign = "=";
      (static_cast<Impl&>(s).*func_)(
          kw->name, internal::OptionParser<Value>()(&oi, kw, value), info_);
      return oi.n_badopts == 0;
    }
  };

 protected:
  // Adds an integer option. The arguments getter and setter should be
  // pointers to member functions in the solver class. They are used to
  // get and set an option value respectively.
  void AddIntOption(const char *name, const char *description,
      int (Impl::*getter)(const char *),
      void (Impl::*setter)(const char *, int)) {
    AddOption(name, std::auto_ptr<SolverOption>(
        new ConcreteOption<int>(description, getter, setter)));
  }

  // Adds an integer option with additional information. The argument f
  // should be a pointer to a member function in the solver class.
  // This function is called after the option is parsed:
  //   (solver.*f)(name, value, info);
  // where handler is a reference to a Handler object, name is an option name
  // of type "const char*", value is an option value of type "int" and info
  // is a reference to the Info object passed to AddIntOption.
  template <typename Info>
  void AddIntOption(const char *name, const char *description,
      void (Impl::*f)(const char *, int, const Info &), const Info &info) {
    typedef void (Impl::*Func)(const char *, int, const Info &);
    AddOption(name, std::auto_ptr<SolverOption>(
        new ConcreteOptionWithInfo<Func, Info, int>(description, f, info)));
  }

  // The same as above but with Info argument passed by value.
  template <typename Info>
  void AddIntOption(const char *name, const char *description,
      void (Impl::*f)(const char *, int, Info), Info info) {
    typedef void (Impl::*Func)(const char *, int, Info);
    AddOption(name, std::auto_ptr<SolverOption>(
        new ConcreteOptionWithInfo<Func, Info, int>(description, f, info)));
  }

  // Adds a double option. The arguments getter and setter should be
  // pointers to member functions in the solver class. They are used to
  // get and set an option value respectively.
  void AddDblOption(const char *name, const char *description,
      double (Impl::*getter)(const char *),
      void (Impl::*setter)(const char *, double)) {
    AddOption(name, std::auto_ptr<SolverOption>(
        new ConcreteOption<double>(description, getter, setter)));
  }

  // Adds a double option with additional information. The argument f
  // should be a pointer to a member function in the solver class.
  // This function is called after the option is parsed:
  //   (solver.*f)(name, value, info);
  // where handler is a reference to a Handler object, name is an option name
  // of type "const char*", value is an option value of type "double" and info
  // is a reference to the Info object passed to AddDblOption.
  template <typename Info>
  void AddDblOption(const char *name, const char *description,
      void (Impl::*f)(const char *, double, const Info &), const Info &info) {
    typedef void (Impl::*Func)(const char *, double, const Info &);
    AddOption(name, std::auto_ptr<SolverOption>(
        new ConcreteOptionWithInfo<Func, Info, double>(description, f, info)));
  }

  // The same as above but with Info argument passed by value.
  template <typename Info>
  void AddDblOption(const char *name, const char *description,
      void (Impl::*f)(const char *, double, Info), Info info) {
    typedef void (Impl::*Func)(const char *, double, Info);
    AddOption(name, std::auto_ptr<SolverOption>(
        new ConcreteOptionWithInfo<Func, Info, double>(description, f, info)));
  }

  // Adds a string option. The arguments getter and setter should be
  // pointers to member functions in the solver class. They are used to
  // get and set an option value respectively.
  void AddStrOption(const char *name, const char *description,
      std::string (Impl::*getter)(const char *),
      void (Impl::*setter)(const char *, const char *)) {
    AddOption(name, std::auto_ptr<SolverOption>(
        new ConcreteOption<const char*, std::string>(
            description, getter, setter)));
  }

  // Adds a string option with additional information. The argument f
  // should be a pointer to a member function in the Handler class.
  // This function is called after the option is parsed:
  //   (handler.*f)(name, value, info);
  // where handler is a reference to a Handler object, name is an option name
  // of type "const char*", value is an option value of type "const char*"
  // and info is a reference to the Info object passed to AddStrOption.
  template <typename Info>
  void AddStrOption(const char *name, const char *description,
      void (Impl::*f)(const char *, const char *, const Info &),
      const Info &info) {
    typedef void (Impl::*Func)(const char *, const char *, const Info &);
    AddOption(name, std::auto_ptr<SolverOption>(
        new ConcreteOptionWithInfo<Func, Info, const char*>(
            description, f, info)));
  }

  // The same as above but with Info argument passed by value.
  template <typename Info>
  void AddStrOption(const char *name, const char *description,
      void (Impl::*f)(const char *, const char *, Info),
      const Info &info) {
    typedef void (Impl::*Func)(const char *, const char *, Info);
    AddOption(name, std::auto_ptr<SolverOption>(
        new ConcreteOptionWithInfo<Func, Info, const char*>(
            description, f, info)));
  }

 public:
  Solver(fmt::StringRef name, fmt::StringRef long_name = 0, long date = 0)
  : BasicSolver(name, long_name, date) {}
};
}

#endif  // SOLVERS_UTIL_SOLVER_H_
