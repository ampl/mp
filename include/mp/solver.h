/*
 A mathematical optimization solver.

 Copyright (C) 2012 AMPL Optimization Inc

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

 Author: Victor Zverovich
 */

#ifndef MP_SOLVER_H_
#define MP_SOLVER_H_

#include <cctype>
#include <csignal>
#include <cstring>

#include <stdint.h>

#include <limits>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "mp/arrayref.h"
#include "mp/clock.h"
#include "mp/error.h"
#include "mp/format.h"
#include "mp/nl.h"
#include "mp/option.h"
#include "mp/sol.h"

namespace mp {

class Problem;
class Solver;

// Information about a possible option value.
struct OptionValueInfo {
  const char *value;
  const char *description;
  intptr_t data;  // Solver-specific data associated with this value.
};

// A reference to an array of OptionValueInfo objects.
class ValueArrayRef {
 private:
  const OptionValueInfo *values_;
  int size_;

 public:
  ValueArrayRef() : values_(), size_() {}

  template <int SIZE>
  ValueArrayRef(const OptionValueInfo (&values)[SIZE], int offset = 0)
  : values_(values + offset), size_(SIZE - offset) {
    assert(offset >= 0 && offset < SIZE);
  }

  int size() const { return size_; }

  typedef const OptionValueInfo *iterator;

  iterator begin() const { return values_; }
  iterator end() const { return values_ + size_; }
};

namespace internal {

// Formats the string s containing the reStructuredText (RST) markup
// and writes it to w.
// w:      stream writer used for output
// indent: indentation to use for the formatted text
// s:      string containing reStructuredText to format
// values: information about possible option values to be formatted by the
//         value-table directive
void FormatRST(fmt::Writer &w, fmt::StringRef s,
    int indent = 0, ValueArrayRef values = ValueArrayRef());

// A helper class for implementing an option of type T.
template <typename T>
struct OptionHelper;

template <>
struct OptionHelper<int> {
  typedef int Arg;
  static const char TYPE_NAME[];
  static void Write(fmt::Writer &w, Arg value) { w << value; }
  static int Parse(const char *&s);
  static int CastArg(int value) { return value; }
};

template <>
struct OptionHelper<double> {
  typedef double Arg;
  static const char TYPE_NAME[];
  static void Write(fmt::Writer &w, double value) { w << value; }
  static double Parse(const char *&s);
  static double CastArg(double value) { return value; }
};

template <>
struct OptionHelper<std::string> {
  typedef const char *Arg;
  static const char TYPE_NAME[];
  static void Write(fmt::Writer &w, const std::string &s) { w << s; }
  static std::string Parse(const char *&s);
  static const char *CastArg(const std::string &s) { return s.c_str(); }
};
}  // namespace internal

// An interface for receiving errors reported via Solver::ReportError.
class ErrorHandler {
 protected:
  ~ErrorHandler() {}

 public:
  virtual void HandleError(fmt::StringRef message) = 0;
};

// An interface for receiving solver output.
class OutputHandler {
 protected:
  ~OutputHandler() {}

 public:
  virtual void HandleOutput(fmt::StringRef output) = 0;
};

// An interface for receiving solutions.
class SolutionHandler {
 public:
  virtual ~SolutionHandler() {}

  // Receives a feasible solution.
  virtual void HandleFeasibleSolution(fmt::StringRef message,
      const double *values, const double *dual_values, double obj_value) = 0;

  // Receives the final solution or a notification that the problem is
  // infeasible or unbounded.
  virtual void HandleSolution(fmt::StringRef message,
      const double *values, const double *dual_values, double obj_value) = 0;
};

class BasicSolutionHandler : public SolutionHandler {
 public:
  virtual void HandleFeasibleSolution(fmt::StringRef,
      const double *, const double *, double) {}
  virtual void HandleSolution(fmt::StringRef,
      const double *, const double *, double) {}
};

// An interface for receiving interrupt notifications.
class Interruptible {
 protected:
  ~Interruptible() {}

 public:
  virtual void Interrupt() = 0;
};

class Solver;

// TODO:
// Don't use SignalHandler in the solvers' DoSolve methods, because DoSolve
// can be called from multiple threads & SignalHandler has static data.
// Instead, move SignalHandler to the implementation of Solve, use it
// from a separate method like RunFromMain & pass an interface that allows
// interrupting the Solver:
//   DoSolve(Interrupter *interrupter);  // better name?
// To make sure that SignalHandler is no used by mistake make it local to
// the RunFromMain method.

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
  explicit SignalHandler(const Solver &s, Interruptible *i = 0);

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
  const char *name_;
  const char *description_;
  ValueArrayRef values_;
  bool is_flag_;

 public:
  // Constructs a SolverOption object.
  //
  // The solver option stores pointers to the passed name and description and
  // doesn't copy the strings. Normally both the name and the description are
  // string literals and have static storage duration but if this is not the
  // case make sure that these strings' lifetimes are longer than that of the
  // option object.
  //
  // The description should be written in a subset of reStructuredText (RST).
  // Currently the following RST constructs are supported:
  //
  // * paragraphs
  // * bullet lists
  // * literal blocks
  // * line blocks
  // * the value-table directive (.. value-table::) which is replaced by a
  //   table of option values as given by the values array
  //
  // name:        option name
  // description: option description
  // values:      information about possible option values
  SolverOption(const char *name, const char *description,
      ValueArrayRef values = ValueArrayRef(), bool is_flag = false)
  : name_(name), description_(description),
    values_(values), is_flag_(is_flag) {}
  virtual ~SolverOption() {}

  // Returns the option name.
  const char *name() const { return name_; }

  // Returns the option description.
  const char *description() const { return description_; }

  // Returns the information about possible values.
  ValueArrayRef values() const { return values_; }

  // Returns true if this option is a flag, i.e. it doesn't take a value.
  bool is_flag() const { return is_flag_; }

  // Formats the option value. Throws OptionError in case of error.
  virtual void Write(fmt::Writer &w) = 0;

  // Parses a string and sets the option value. Throws InvalidOptionValue
  // if the value is invalid or OptionError in case of another error.
  virtual void Parse(const char *&s) = 0;
};

// An exception thrown when an invalid value is provided for an option.
class InvalidOptionValue : public OptionError {
 private:
  template <typename T>
  static std::string Format(fmt::StringRef name, T value) {
    return fmt::format("Invalid value \"{}\" for option \"{}\"", value, name);
  }

 public:
  template <typename T>
  InvalidOptionValue(fmt::StringRef name, T value)
  : OptionError(Format(name, value)) {}

  template <typename T>
  InvalidOptionValue(const SolverOption &opt, T value)
  : OptionError(Format(opt.name(), value)) {}
};

template <typename T>
class TypedSolverOption : public SolverOption {
 public:
  TypedSolverOption(const char *name, const char *description,
      ValueArrayRef values = ValueArrayRef())
  : SolverOption(name, description, values) {}

  void Write(fmt::Writer &w) {
    internal::OptionHelper<T>::Write(w, GetValue());
  }

  void Parse(const char *&s) {
    const char *start = s;
    T value = internal::OptionHelper<T>::Parse(s);
    if (*s && !std::isspace(*s)) {
      do ++s;
      while (*s && !std::isspace(*s));
      throw InvalidOptionValue(name(), std::string(start, s - start));
    }
    SetValue(internal::OptionHelper<T>::CastArg(value));
  }

  // Returns the option value.
  virtual T GetValue() const = 0;

  // Sets the option value or throws InvalidOptionValue if the value is invalid.
  virtual void SetValue(typename internal::OptionHelper<T>::Arg value) = 0;
};

// A mathematical optimization solver.
//
// Example:
//
// class MySolver : public Solver {
//  public:
//   void GetTestOption(const char *name, int value, int info) {
//     // Returns the option value; info is an arbitrary value passed as
//     // the last argument to AddIntOption. It can be useful if the same
//     // function handles multiple options.
//     ...
//   }
//   void SetTestOption(const char *name, int value, int info) {
//     // Set the option value; info is the same as in GetTestOption.
//     ...
//   }
//
//   MySolver()  {
//     AddIntOption("test", "This is a test option",
//                  &MySolver::GetTestOption, &MySolver::SetTestOption, 42);
//   }
// };
class Solver : private ErrorHandler, private OutputHandler {
 private:
  std::string name_;
  std::string long_name_;
  std::string version_;
  std::string license_info_;
  long date_;
  int wantsol_;
  int obj_precision_;

  enum {SHOW_VERSION = 1};
  int flags_;

  bool has_errors_;
  OutputHandler *output_handler_;
  ErrorHandler *error_handler_;

  // The filename stub for returning multiple solutions.
  std::string solution_stub_;

  // Specifies whether to return the number of solutions in the .nsol suffix.
  bool count_solutions_;

  unsigned read_flags_;  // flags passed to Problem::Read

  struct OptionNameLess {
    bool operator()(const SolverOption *lhs, const SolverOption *rhs) const;
  };

  std::string option_header_;
  typedef std::set<SolverOption*, OptionNameLess> OptionSet;
  OptionSet options_;

  bool timing_;

  class SuffixInfo {
   private:
    const char *name_;
    const char *table_;
    int kind_;
    int nextra_;

   public:
    SuffixInfo(const char *name, const char *table, int kind, int nextra)
      : name_(name), table_(table), kind_(kind), nextra_(nextra) {}

    const char *name() const { return name_; }
    const char *table() const { return table_; }
    int kind() const { return kind_; }
    int nextra() const { return nextra_; }
  };
  typedef std::vector<SuffixInfo> SuffixList;
  SuffixList suffixes_;

  friend class ASLSolver;

  void HandleOutput(fmt::StringRef output) {
    std::fputs(output.c_str(), stdout);
  }

  void HandleError(fmt::StringRef message) {
    std::fputs(message.c_str(), stderr);
    std::fputc('\n', stderr);
  }

  // Finds an option and returns a pointer to it if found or null otherwise.
  SolverOption *FindOption(const char *name) const;

  // Returns the option with specified name.
  SolverOption *GetOption(const char *name) const {
    SolverOption *opt = FindOption(name);
    if (!opt)
      throw OptionError(fmt::format("Unknown option \"{}\"", name));
    return opt;
  }

  template <typename T>
  T GetOptionValue(const char *name) const {
    const TypedSolverOption<T> *opt =
        dynamic_cast<TypedSolverOption<T> *>(GetOption(name));
    if (!opt)
      throw OptionTypeError(name, internal::OptionHelper<T>::TYPE_NAME);
    return opt->GetValue();
  }

  template <typename T>
  void SetOptionValue(const char *name,
      typename internal::OptionHelper<T>::Arg value) {
    TypedSolverOption<T> *opt =
        dynamic_cast<TypedSolverOption<T> *>(GetOption(name));
    if (!opt)
      throw OptionTypeError(name, internal::OptionHelper<T>::TYPE_NAME);
    opt->SetValue(value);
  }

  static OptionError OptionTypeError(fmt::StringRef name, fmt::StringRef type) {
    fmt::Writer msg;
    msg.write("Option \"{}\" is not of type \"{}\"", name, type);
    return OptionError(msg.c_str());
  }

  // Parses an option string.
  void ParseOptionString(const char *s, unsigned flags);

  // Handler should be a class derived from Solver that will receive
  // notifications about parsed options.
  template <typename Handler, typename T>
  class ConcreteOption : public TypedSolverOption<T> {
   private:
    typedef typename internal::OptionHelper<T>::Arg Arg;
    typedef T (Handler::*Get)(const SolverOption &) const;
    typedef void (Handler::*Set)(const SolverOption &, Arg);

    Handler &handler_;
    Get get_;
    Set set_;

   public:
    ConcreteOption(const char *name, const char *description,
        Solver *s, Get get, Set set, ValueArrayRef values = ValueArrayRef())
    : TypedSolverOption<T>(name, description, values),
      handler_(static_cast<Handler&>(*s)), get_(get), set_(set) {}

    T GetValue() const { return (handler_.*get_)(*this); }
    void SetValue(Arg value) { (handler_.*set_)(*this, value); }
  };

  template <typename Handler, typename T,
            typename Info, typename InfoArg = Info>
  class ConcreteOptionWithInfo : public TypedSolverOption<T> {
   private:
    typedef typename internal::OptionHelper<T>::Arg Arg;
    typedef T (Handler::*Get)(const SolverOption &, InfoArg) const;
    typedef void (Handler::*Set)(const SolverOption &, Arg, InfoArg);

    Handler &handler_;
    Get get_;
    Set set_;
    Info info_;

   public:
    ConcreteOptionWithInfo(const char *name, const char *description, Solver *s,
        Get get, Set set, InfoArg info, ValueArrayRef values = ValueArrayRef())
    : TypedSolverOption<T>(name, description, values),
      handler_(static_cast<Handler&>(*s)), get_(get), set_(set), info_(info) {}

    T GetValue() const { return (handler_.*get_)(*this, info_); }
    void SetValue(Arg value) { (handler_.*set_)(*this, value, info_); }
  };

 public:
#ifdef MP_USE_UNIQUE_PTR
  typedef std::unique_ptr<SolverOption> OptionPtr;
#else
  typedef std::auto_ptr<SolverOption> OptionPtr;
  static OptionPtr move(OptionPtr p) { return p; }
#endif

  struct DoubleFormatter {
    double value;
    int precision;

    friend void format(fmt::BasicFormatter<char> &f,
                       const char *, DoubleFormatter df) {
      f.writer().write("{:.{}}", df.value, df.precision);
    }
  };

 protected:
  // Flags for the Solver constructor.
  enum {
    // Multiple solution support.
    // Makes Solver register "countsolutions" and "solutionstub" options
    // and write every solution passed to HandleFeastibleSolution to a file
    // solutionstub & i & ".sol" where i is a solution number.
    MULTIPLE_SOL = 1
  };

  // Constructs a Solver object.
  // date: The solver date in YYYYMMDD format.
  Solver(fmt::StringRef name, fmt::StringRef long_name, long date, int flags);

  void set_long_name(fmt::StringRef name) { long_name_ = name; }
  void set_version(fmt::StringRef version) { version_ = version; }

  // Sets the flags for Problem::Read.
  void set_read_flags(unsigned flags) { read_flags_ = flags; }

  // Sets a text to be displayed before option descriptions.
  void set_option_header(const char *header) { option_header_ = header; }

  void AddOption(OptionPtr opt) {
    // First insert the option, then release a pointer to it. Doing the other
    // way around may lead to a memory leak if insertion throws an exception.
    options_.insert(opt.get());
    opt.release();
  }

  // Adds an integer option.
  // The option stores pointers to the name and the description so make
  // sure that these strings have sufficient lifetimes (normally these are
  // string literals).
  // The arguments get and set should be pointers to member functions in the
  // solver class. They are used to get and set an option value respectively.
  template <typename Handler>
  void AddIntOption(const char *name, const char *description,
      int (Handler::*get)(const SolverOption &) const,
      void (Handler::*set)(const SolverOption &, int)) {
    AddOption(OptionPtr(new ConcreteOption<Handler, int>(
        name, description, this, get, set)));
  }

  // Adds an integer option with additional information.
  // The option stores pointers to the name and the description so make
  // sure that these strings have sufficient lifetimes (normally these are
  // string literals).
  // The arguments get and set should be pointers to member functions in the
  // solver class. They are used to get and set an option value respectively.
  template <typename Handler, typename Info>
  void AddIntOption(const char *name, const char *description,
      int (Handler::*get)(const SolverOption &, const Info &) const,
      void (Handler::*set)(const SolverOption &, int, const Info &),
      const Info &info) {
    AddOption(OptionPtr(
        new ConcreteOptionWithInfo<Handler, int, Info, const Info &>(
            name, description, this, get, set, info)));
  }

  // The same as above but with Info argument passed by value.
  template <typename Handler, typename Info>
  void AddIntOption(const char *name, const char *description,
      int (Handler::*get)(const SolverOption &, Info) const,
      void (Handler::*set)(const SolverOption &, int, Info), Info info) {
    AddOption(OptionPtr(new ConcreteOptionWithInfo<Handler, int, Info>(
            name, description, this, get, set, info)));
  }

  // Adds a double option.
  // The option stores pointers to the name and the description so make
  // sure that these strings have sufficient lifetimes (normally these are
  // string literals).
  // The arguments get and set should be pointers to member functions in the
  // solver class. They are used to get and set an option value respectively.
  template <typename Handler>
  void AddDblOption(const char *name, const char *description,
      double (Handler::*get)(const SolverOption &) const,
      void (Handler::*set)(const SolverOption &, double)) {
    AddOption(OptionPtr(new ConcreteOption<Handler, double>(
        name, description, this, get, set)));
  }

  // Adds a double option with additional information.
  // The option stores pointers to the name and the description so make
  // sure that these strings have sufficient lifetimes (normally these are
  // string literals).
  // The arguments get and set should be pointers to member functions in the
  // solver class. They are used to get and set an option value respectively.
  template <typename Handler, typename Info>
  void AddDblOption(const char *name, const char *description,
      double (Handler::*get)(const SolverOption &, const Info &) const,
      void (Handler::*set)(const SolverOption &, double, const Info &),
      const Info &info) {
    AddOption(OptionPtr(
        new ConcreteOptionWithInfo<Handler, double, Info, const Info &>(
            name, description, this, get, set, info)));
  }

  // The same as above but with Info argument passed by value.
  template <typename Handler, typename Info>
  void AddDblOption(const char *name, const char *description,
      double (Handler::*get)(const SolverOption &, Info) const,
      void (Handler::*set)(const SolverOption &, double, Info), Info info) {
    AddOption(OptionPtr(new ConcreteOptionWithInfo<Handler, double, Info>(
            name, description, this, get, set, info)));
  }

  // Adds a string option.
  // The option stores pointers to the name and the description so make
  // sure that these strings have sufficient lifetimes (normally these are
  // string literals).
  // The arguments get and set should be pointers to member functions in the
  // solver class. They are used to get and set an option value respectively.
  template <typename Handler>
  void AddStrOption(const char *name, const char *description,
      std::string (Handler::*get)(const SolverOption &) const,
      void (Handler::*set)(const SolverOption &, const char *),
      ValueArrayRef values = ValueArrayRef()) {
    AddOption(OptionPtr(new ConcreteOption<Handler, std::string>(
        name, description, this, get, set, values)));
  }

  // Adds a string option with additional information.
  // The option stores pointers to the name and the description so make
  // sure that these strings have sufficient lifetimes (normally these are
  // string literals).
  // The arguments get and set should be pointers to member functions in the
  // solver class. They are used to get and set an option value respectively.
  template <typename Handler, typename Info>
  void AddStrOption(const char *name, const char *description,
      std::string (Handler::*get)(const SolverOption &, const Info &) const,
      void (Handler::*set)(const SolverOption &, const char *, const Info &),
      const Info &info, ValueArrayRef values = ValueArrayRef()) {
    AddOption(OptionPtr(
        new ConcreteOptionWithInfo<Handler, std::string, Info, const Info &>(
            name, description, this, get, set, info, values)));
  }

  // The same as above but with Info argument passed by value.
  template <typename Handler, typename Info>
  void AddStrOption(const char *name, const char *description,
      std::string (Handler::*get)(const SolverOption &, Info) const,
      void (Handler::*set)(const SolverOption &, const char *, Info),
      Info info, ValueArrayRef values = ValueArrayRef()) {
    AddOption(OptionPtr(new ConcreteOptionWithInfo<Handler, std::string, Info>(
            name, description, this, get, set, info, values)));
  }

  virtual void HandleUnknownOption(const char *name) {
    ReportError("Unknown option \"{}\"", name);
  }

  // Adds a suffix.
  void AddSuffix(const char *name, const char *table,
                 int kind, int nextra = 0) {
    suffixes_.push_back(SuffixInfo(name, table, kind, nextra));
  }

 public:
  virtual ~Solver();

  // Returns the solver name.
  const char *name() const { return name_.c_str(); }

  // Returns the long solver name.
  // This name is used in startup "banner".
  const char *long_name() const { return long_name_.c_str(); }

  // Returns the solver version.
  const char *version() const { return version_.c_str(); }

  // Returns the solver date in YYYYMMDD format.
  long date() const { return date_; }

  // Returns the solver flags.
  // Possible values that can be combined with bitwise OR:
  //   ASL_OI_want_funcadd
  //   ASL_OI_keep_underscores
  //   ASL_OI_show_version
  int flags() const { return flags_; }

  // Returns the value of the wantsol option which specifies what solution
  // information to write in a stand-alone invocation (no -AMPL on the
  // command line).  Possible values that can be combined with bitwise OR:
  //   1 = write .sol file
  //   2 = primal variables to stdout
  //   4 = dual variables to stdout
  //   8 = suppress solution message
  int wantsol() const { return wantsol_; }
  void set_wantsol(int value) { wantsol_ = value; }

  // Returns true if the timing is enabled.
  bool timing() const { return timing_; }

  // Returns the error handler.
  ErrorHandler *error_handler() { return error_handler_; }

  // Sets the error handler.
  void set_error_handler(ErrorHandler *eh) { error_handler_ = eh; }

  // Returns the output handler.
  OutputHandler *output_handler() { return output_handler_; }

  // Sets the output handler.
  void set_output_handler(OutputHandler *oh) { output_handler_ = oh; }

  const char *solution_stub() const { return solution_stub_.c_str(); }

  bool need_multiple_solutions() const {
    return count_solutions_ || !solution_stub_.empty();
  }

  // Returns the number of options.
  int num_options() const { return static_cast<int>(options_.size()); }

  // Option iterator.
  class option_iterator :
    public std::iterator<std::forward_iterator_tag, SolverOption> {
   private:
    OptionSet::const_iterator it_;

    friend class Solver;

    explicit option_iterator(OptionSet::const_iterator it) : it_(it) {}

   public:
    option_iterator() {}

    const SolverOption &operator*() const { return **it_; }
    const SolverOption *operator->() const { return *it_; }

    option_iterator &operator++() {
      ++it_;
      return *this;
    }

    option_iterator operator++(int ) {
      option_iterator it(*this);
      ++it_;
      return it;
    }

    bool operator==(option_iterator other) const { return it_ == other.it_; }
    bool operator!=(option_iterator other) const { return it_ != other.it_; }
  };

  option_iterator option_begin() const {
    return option_iterator(options_.begin());
  }
  option_iterator option_end() const {
    return option_iterator(options_.end());
  }

  // Returns the option header.
  const char *option_header() const { return option_header_.c_str(); }

  // Flags for ParseOptions.
  enum {
    // Don't echo options during parsing.
    NO_OPTION_ECHO = 1
  };

  // Parses solver options and returns true if there were no errors and
  // false otherwise. It accepts a pointer to the problem because some
  // options may depend on problem features.
  virtual bool ParseOptions(
      char **argv, unsigned flags = 0, const Problem *p = 0);

  // Returns the value of an integer option.
  // Throws OptionError if there is no such option or it has a different type.
  int GetIntOption(const char *name) const {
    return GetOptionValue<int>(name);
  }

  // Sets the value of an integer option.
  // Throws OptionError if there is no such option or it has a different type.
  void SetIntOption(const char *name, int value) {
    SetOptionValue<int>(name, value);
  }

  // Returns the value of a double option.
  // Throws OptionError if there is no such option or it has a different type.
  double GetDblOption(const char *name) const {
    return GetOptionValue<double>(name);
  }

  // Sets the value of a double option.
  // Throws OptionError if there is no such option or it has a different type.
  void SetDblOption(const char *name, double value) {
    SetOptionValue<double>(name, value);
  }

  // Returns the value of a string option.
  // Throws OptionError if there is no such option or it has a different type.
  std::string GetStrOption(const char *name) const {
    return GetOptionValue<std::string>(name);
  }

  // Sets the value of a string option.
  // Throws OptionError if there is no such option or it has a different type.
  void SetStrOption(const char *name, fmt::StringRef value) {
    SetOptionValue<std::string>(name, value.c_str());
  }

  const SuffixList &suffixes() const { return suffixes_; }

  // Reports an error printing the formatted error message to stderr.
  // Usage: ReportError("File not found: {}") << filename;
  void ReportError(fmt::StringRef format, const fmt::ArgList &args) {
    has_errors_ = true;
    fmt::Writer w;
    w.write(format, args);
    error_handler_->HandleError(fmt::StringRef(w.c_str(), w.size()));
  }
  FMT_VARIADIC(void, ReportError, fmt::StringRef)

  // Formats a string and prints it to stdout or, if an output handler
  // is registered, sends it to the output handler.
  void Print(fmt::StringRef format, const fmt::ArgList &args) {
    fmt::Writer w;
    w.write(format, args);
    output_handler_->HandleOutput(fmt::StringRef(w.c_str(), w.size()));
  }
  FMT_VARIADIC(void, Print, fmt::StringRef)

  // Prints version information.
  bool ShowVersion();

  // The default precision used in FormatObjValue if the objective_precision
  // option is not specified or 0.
  enum { DEFAULT_PRECISION = 15 };

  // Returns a formatter that writes value using objective precision.
  // Usage:
  //   Print("objective {}", FormatObjValue(obj_value));
  DoubleFormatter FormatObjValue(double value);
};

template <typename ProblemBuilderT>
class SolverImpl : public Solver {
 public:
  typedef ProblemBuilderT ProblemBuilder;

  SolverImpl(fmt::StringRef name, fmt::StringRef long_name = 0,
             long date = 0, int flags = 0)
    : Solver(name, long_name, date, flags) {}
};

// Adapts a solution for WriteSol.
template <typename ProblemBuilder>
class SolutionAdapter {
 private:
  ProblemBuilder *builder_;
  const char *message_;
  mp::ArrayRef<int> options_;
  mp::ArrayRef<double> values_;
  mp::ArrayRef<double> dual_values_;

 public:
  SolutionAdapter(ProblemBuilder *pb, const char *message,
                  mp::ArrayRef<int> options, mp::ArrayRef<double> values,
                  mp::ArrayRef<double> dual_values)
    : builder_(pb), message_(message), options_(options),
      values_(values), dual_values_(dual_values) {}

  const char *message() const { return message_; }

  int num_options() const { return options_.size(); }
  int option(int index) const { return options_[index]; }

  int num_values() const { return values_.size(); }
  int value(int index) const { return values_[index]; }

  int num_dual_values() const { return dual_values_.size(); }
  int dual_value(int index) const { return dual_values_[index]; }

  const typename ProblemBuilder::SuffixMap *suffixes(int kind) const {
    return builder_ ? &builder_->suffixes(kind) : 0;
  }
};

class NullSolutionHandler : public SolutionHandler {
 public:
  void HandleFeasibleSolution(
        fmt::StringRef, const double *, const double *, double) {}
  void HandleSolution(
        fmt::StringRef, const double *, const double *, double) {}
};

// The default .sol file writer.
class SolFileWriter {
 public:
  template <typename Solution>
  void Write(fmt::StringRef filename, const Solution &sol) {
    WriteSolFile(filename, sol);
  }
};

// A solution writer.
// Solver: optimization solver class
// Writer: .sol writer
template <typename Solver, typename Writer = SolFileWriter>
class SolutionWriter : private Writer, public SolutionHandler {
 private:
  std::string filename_;

  Solver &solver_;

  typedef typename Solver::ProblemBuilder ProblemBuilder;
  ProblemBuilder &builder_;

  // The number of feasible solutions found.
  int num_solutions_;

 public:
  SolutionWriter(fmt::StringRef stub, Solver &s, ProblemBuilder &b)
    : filename_(stub.c_str() + std::string(".sol")), solver_(s), builder_(b),
      num_solutions_(0) {}

  // Returns the .sol writer.
  Writer &sol_writer() { return *this; }

  void HandleFeasibleSolution(fmt::StringRef message,
        const double *values, const double *dual_values, double);

  // Writes the solution to a .sol file.
  void HandleSolution(fmt::StringRef message,
        const double *values, const double *dual_values, double);
};

template <typename Solver, typename Writer>
void SolutionWriter<Solver, Writer>::HandleFeasibleSolution(
    fmt::StringRef message, const double *values,
    const double *dual_values, double) {
  ++num_solutions_;
  const char *solution_stub = solver_.solution_stub();
  if (!*solution_stub)
    return;
  SolutionAdapter<ProblemBuilder> sol(
        0, message.c_str(), ArrayRef<int>(0, 0),
        MakeArrayRef(values, values ? builder_.num_vars() : 0),
        MakeArrayRef(dual_values, dual_values ? builder_.num_cons() : 0));
  fmt::Writer filename;
  filename << solution_stub << num_solutions_ << ".sol";
  this->Write(filename.c_str(), sol);
}

template <typename Solver, typename Writer>
void SolutionWriter<Solver, Writer>::HandleSolution(
    fmt::StringRef message, const double *values,
    const double *dual_values, double) {
  if (solver_.need_multiple_solutions()) {
    typedef typename ProblemBuilder::SuffixPtr SuffixPtr;
    SuffixPtr nsol_suffix = builder_.suffixes(suf::PROBLEM).Find("nsol");
    if (nsol_suffix)
      nsol_suffix->set_value(0, num_solutions_);
  }
  // TODO: pass to WriteSol
  //option_info.bsname = const_cast<char*>(solver_.long_name());
  //option_info.wantsol = solver_.wantsol();
  //const fint *options = problem_.asl_->i.ampl_options_;
  SolutionAdapter<ProblemBuilder> sol(
        &builder_, message.c_str(), ArrayRef<int>(0, 0),
        MakeArrayRef(values, values ? builder_.num_vars() : 0),
        MakeArrayRef(dual_values, dual_values ? builder_.num_cons() : 0));
  this->Write(filename_, sol);
}

// The default .nl file reader.
class NLFileReader {
 public:
  template <typename Handler>
  void Read(fmt::StringRef filename, Handler &handler) {
    ReadNLFile(filename, handler);
  }
};

namespace internal {

// Command-line option parser for a solver application.
// Not to be mistaken with solver option parser built into the Solver class.
class SolverAppOptionParser {
 private:
  Solver &solver_;

  // Command-line options.
  OptionList options_;

  bool echo_solver_options_;

  // Prints usage information and stops processing options.
  bool ShowUsage();

  // Prints information about solver options.
  bool ShowSolverOptions();

  bool WantSol() {
    solver_.set_wantsol(1);
    return true;
  }

  bool DontEchoSolverOptions() {
    echo_solver_options_ = false;
    return true;
  }

  // Stops processing options.
  bool EndOptions() { return false; }

 public:
  explicit SolverAppOptionParser(Solver &s);

  OptionList &options() { return options_; }

  // Retruns true if assignments of solver options should be echoed.
  bool echo_solver_options() const { return echo_solver_options_; }

  // Parses command-line options.
  const char *Parse(char **&argv);
};
}  // namespace internal

// A solver application.
// Solver: optimization solver class; normally a subclass of SolverImpl
// Reader: .nl reader
template <typename Solver, typename Reader = NLFileReader>
class SolverApp : private Reader {
 private:
  Solver solver_;
  internal::SolverAppOptionParser option_parser_;

 public:
  SolverApp() : option_parser_(solver_) {}

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
  int Run(char **argv);
};

template <typename Solver, typename Reader>
int SolverApp<Solver, Reader>::Run(char **argv) {
  // Parse command-line arguments.
  const char *filename = option_parser_.Parse(argv);
  if (!filename) return 0;

  // Add .nl extension if necessary.
  std::string nl_filename = filename, filename_no_ext = nl_filename;
  const char *ext = std::strrchr(filename, '.');
  if (!ext || std::strcmp(ext, ".nl") != 0)
    nl_filename += ".nl";
  else
    filename_no_ext.resize(filename_no_ext.size() - 3);

  // Read the problem.
  steady_clock::time_point start = steady_clock::now();
  typedef typename Solver::ProblemBuilder ProblemBuilder;
  // TODO: use name provider instead of passing filename to builder
  ProblemBuilder builder(solver_.GetProblemBuilder(filename_no_ext));
  ProblemBuilderToNLAdapter<ProblemBuilder> adapter(builder);
  this->Read(nl_filename, adapter);
  builder.EndBuild();
  double read_time = GetTimeAndReset(start);

  // Parse solver options.
  unsigned flags =
      option_parser_.echo_solver_options() ? 0 : Solver::NO_OPTION_ECHO;
  solver_.ParseOptions(argv, flags);
  if (solver_.timing())
    solver_.Print("Input time = {:.6f}s\n", read_time);

  // Solve the problem writing solution(s) if necessary.
  std::auto_ptr<SolutionHandler> sol_handler;
  if (solver_.wantsol() != 0) {
    sol_handler.reset(
          new SolutionWriter<Solver>(filename_no_ext, solver_, builder));
  } else {
    sol_handler.reset(new NullSolutionHandler());
  }
  solver_.Solve(builder, *sol_handler);
  return 0;
}

#ifdef MP_USE_UNIQUE_PTR
typedef std::unique_ptr<Solver> SolverPtr;
#else
typedef std::auto_ptr<Solver> SolverPtr;
inline SolverPtr move(SolverPtr p) { return p; }
#endif

// Implement this function in your code returning a new concrete solver object.
// options: Solver initialization options.
// Example:
//   SolverPtr CreateSolver(const char *) { return SolverPtr(new MySolver()); }
SolverPtr CreateSolver(const char *options);
}  // namespace mp

#endif  // MP_SOLVER_H_
