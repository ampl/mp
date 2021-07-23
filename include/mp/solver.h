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

#if MP_USE_ATOMIC
# include <atomic>
#endif

#include "mp/arrayref.h"
#include "mp/clock.h"
#include "mp/error.h"
#include "mp/format.h"
#include "mp/nl-reader.h"
#include "mp/option.h"
#include "mp/os.h"
#include "mp/problem-builder.h"
#include "mp/sol.h"
#include "mp/suffix.h"

namespace mp {

class ASLProblem;
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
void FormatRST(fmt::Writer &w, fmt::CStringRef s,
    int indent = 0, ValueArrayRef values = ValueArrayRef());

// A helper class for implementing an option of type T.
template <typename T>
struct OptionHelper;

template <>
struct OptionHelper<int> {
  typedef int Arg;
  static void Write(fmt::Writer &w, Arg value) { w << value; }
  static int Parse(const char *&s);
  static int CastArg(fmt::LongLong value) { return static_cast<int>(value); }
};

template <>
struct OptionHelper<fmt::LongLong> {
  typedef fmt::LongLong Arg;
  static void Write(fmt::Writer &w, Arg value) { w << value; }
  static fmt::LongLong Parse(const char *&s) {
    return OptionHelper<int>::Parse(s); // TODO
  }
  static fmt::LongLong CastArg(fmt::LongLong value) { return value; }
};

template <>
struct OptionHelper<double> {
  typedef double Arg;
  static void Write(fmt::Writer &w, double value) { w << value; }
  static double Parse(const char *&s);
  static double CastArg(double value) { return value; }
};

template <>
struct OptionHelper<std::string> {
  typedef fmt::StringRef Arg;
  static void Write(fmt::Writer &w, const std::string &s) { w << s; }
  static std::string Parse(const char *&s);
  static fmt::StringRef CastArg(fmt::StringRef s) { return s; }
};

inline OptionError OptionTypeError(fmt::StringRef name, fmt::StringRef type) {
  return OptionError(
        fmt::format("Option \"{}\" is not of type \"{}\"", name, type));
}

// Perform demo version checks if necessary.
void CheckDemoVersion(const NLHeader &h);
}  // namespace internal

// An interface for receiving errors reported via Solver::ReportError.
class ErrorHandler {
 protected:
  ~ErrorHandler() {}

 public:
  virtual void HandleError(fmt::CStringRef message) = 0;
};

// An interface for receiving solver output.
class OutputHandler {
 protected:
  ~OutputHandler() {}

 public:
  virtual void HandleOutput(fmt::CStringRef output) = 0;
};

// An interface for receiving solutions.
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

class BasicSolutionHandler : public SolutionHandler {
 public:
  virtual void HandleFeasibleSolution(fmt::CStringRef,
      const double *, const double *, double) {}
  virtual void HandleSolution(int, fmt::CStringRef,
      const double *, const double *, double) {}
};

// Interrupt handler.
// Returns true if the solver was interrupted, false if it is not running.
typedef bool (*InterruptHandler)(void *);

// An interface for interrupting solution process.
// When a solver is run in a terminal it should respond to SIGINT (Ctrl-C)
// by interrupting its execution and returning the best solution found.
// This can be done in two ways. The first one is to check the return value
// of Interrupter::Stop() periodically. The second is to register an object
// that implements the Interruptible interface.
class Interrupter {
 public:
  virtual ~Interrupter() {}

  // Returns true if the solution process should be stopped.
  virtual bool Stop() const = 0;

  // Sets a handler function.
  // The handler function must be safe to call from a signal handler.
  // In particular, it must only call async-signal-safe library functions.
  virtual void SetHandler(InterruptHandler handler, void *data) = 0;
};

// A solver option.
class SolverOption {
 private:
  const char *name_;
  const char* qualifiedname_;
  const char* description_;

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
  SolverOption(const char *name, const char* qualifiedName, const char *description,
      ValueArrayRef values = ValueArrayRef(), bool is_flag = false)
  : name_(name), qualifiedname_(qualifiedName), description_(description),
    values_(values), is_flag_(is_flag) {}

  virtual ~SolverOption() {}

  // Returns the option name.
  const char *name() const { return name_; }

  // Returns the option name.
  const char* qualifiedName() const { return qualifiedname_; }

  // Return/set the option description.
  const char *description() const { return description_; }
  void SetDesc(const char* d) { description_=d; }

  // Returns the information about possible values.
  ValueArrayRef values() const { return values_; }

  // Returns true if this option is a flag, i.e. it doesn't take a value.
  bool is_flag() const { return is_flag_; }

  // Returns the option value.
  virtual void GetValue(fmt::LongLong &) const {
    throw internal::OptionTypeError(name_, "int");
  }
  virtual void GetValue(double &) const {
    throw internal::OptionTypeError(name_, "double");
  }
  virtual void GetValue(std::string &) const {
    throw internal::OptionTypeError(name_, "string");
  }

  virtual void GetValue(int &int_value) const {
    fmt::LongLong value = 0;
    GetValue(value);
    if (value < std::numeric_limits<int>::min() ||
        value > std::numeric_limits<int>::max()) {
      throw Error("Value {} doesn't fit in int", value);
    }
    int_value = static_cast<int>(value);
  }

  template <typename T>
  T GetValue() const {
    T value = T();
    GetValue(value);
    return value;
  }

  // Sets the option value or throws InvalidOptionValue if the value is invalid.
  virtual void SetValue(fmt::LongLong) {
    throw internal::OptionTypeError(name_, "int");
  }
  virtual void SetValue(double) {
    throw internal::OptionTypeError(name_, "double");
  }
  virtual void SetValue(fmt::StringRef) {
    throw internal::OptionTypeError(name_, "string");
  }
  virtual void SetValue(int value) {
    fmt::LongLong long_value = value;
    SetValue(long_value);
  }

  // Formats the option value. Throws OptionError in case of error.
  virtual void Write(fmt::Writer &w) = 0;

  // Parses a string and sets the option value. Throws InvalidOptionValue
  // if the value is invalid or OptionError in case of another error.
  virtual void Parse(const char *&s) = 0;

  virtual std::string echo() {
    return name();
  }
};

/*
*/
class SolverOptionSynonym : public SolverOption
{
  SolverOption* real_;
  std::string desc_;
public:
  SolverOptionSynonym(const char* name, SolverOption& real) :
    SolverOption(name, name, NULL), real_(&real) {
    desc_ = fmt::sprintf("Synonym for %s.", real_->name());
    SetDesc( desc_.c_str() );
  }
  SolverOption* getRealOption() const { return real_; }

  virtual void Write(fmt::Writer& w) {
    real_->Write(w);
  }
  virtual void Parse(const char*& s) {
    real_->Parse(s);
  }

  virtual std::string echo() {
    return fmt::format("{} ({})", name(), real_->echo());
  }
};
// An exception thrown when an invalid value is provided for an option.
class InvalidOptionValue : public OptionError {
 private:
  template <typename T>
  static std::string Format(fmt::StringRef name, T value, fmt::StringRef msg) {
    if (0!=msg.size())
      return fmt::format("Invalid value \"{}\" for option \"{}\", {}",
                         value, name, msg);
    else
      return fmt::format("Invalid value \"{}\" for option \"{}\"", value, name);
  }

 public:
  template <typename T>
  InvalidOptionValue(fmt::StringRef name, T value, fmt::StringRef msg="")
  : OptionError(Format(name, value, msg)) {}

  template <typename T>
  InvalidOptionValue(const SolverOption &opt, T value, fmt::StringRef msg="")
  : OptionError(Format(opt.name(), value, msg)) {}
};

template <typename T>
class TypedSolverOption : public SolverOption {
 public:
  TypedSolverOption(const char *name, const char* qualifiedName, const char *description,
      ValueArrayRef values = ValueArrayRef())
  : SolverOption(name, qualifiedName, description, values) {}

  void Write(fmt::Writer &w) { w << GetValue<T>(); }

  void Parse(const char *&s) {
    const char *start = s;
    T value = internal::OptionHelper<T>::Parse(s);
    if (*s && !std::isspace(*s)) {
      do ++s;
      while (*s && !std::isspace(*s));
      throw InvalidOptionValue(name(), std::string(start, s - start));
    }
    SetValue(value);
  }
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
class Solver : private ErrorHandler,
    private OutputHandler, private Interrupter {
 private:
  std::string name_;
  std::string long_name_;
  std::string version_;
  std::string license_info_;
  long date_;
  int wantsol_;
  int obj_precision_;

  // Index of the objective to optimize starting from 1, 0 to ignore
  // objective, or -1 to use the first objective if there is one.
  int objno_;

  enum {SHOW_VERSION = 1, AMPL_FLAG = 2};
  int bool_options_;

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
  bool multiobj_;

  bool has_errors_;
  OutputHandler *output_handler_;
  ErrorHandler *error_handler_;
  Interrupter *interrupter_;

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

 public:
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

 private:
  SuffixList suffixes_;

  friend class ASLSolver;

  void HandleOutput(fmt::CStringRef output) {
    std::fputs(output.c_str(), stdout);
  }

  void HandleError(fmt::CStringRef message) {
    std::fputs(message.c_str(), stderr);
    std::fputc('\n', stderr);
  }

  // The default implementation of Interrupter does nothing.
  bool Stop() const { return false; }
  void SetHandler(InterruptHandler, void *) {}

  // Returns the option with specified name.
  SolverOption *GetOption(const char *name) const {
    SolverOption *opt = FindOption(name);
    if (!opt)
      throw OptionError(fmt::format("Unknown option \"{}\"", name));
    return opt;
  }

  // Parses an option string.
  void ParseOptionString(const char *s, unsigned flags);

  // Handler should be a class derived from Solver that will receive
  // notifications about parsed options.
  template <typename Handler, typename T, typename AccessorT = T>
  class ConcreteOption : public TypedSolverOption<T> {
   private:
    typedef AccessorT (Handler::*Get)(const SolverOption &) const;
    typedef void (Handler::*Set)(
        const SolverOption &, typename internal::OptionHelper<AccessorT>::Arg);

    Handler &handler_;
    Get get_;
    Set set_;

   public:
    ConcreteOption(const char *name, const char* qualifiedName, const char *description,
        Solver *s, Get get, Set set, ValueArrayRef values = ValueArrayRef())
    : TypedSolverOption<T>(name, qualifiedName, description, values),
      handler_(static_cast<Handler&>(*s)), get_(get), set_(set) {}

    void GetValue(T &value) const { value = (handler_.*get_)(*this); }
    void SetValue(typename internal::OptionHelper<T>::Arg value) {
      (handler_.*set_)(*this,
                       internal::OptionHelper<AccessorT>::CastArg(value));
    }
  };
  public:
  template <typename Handler, typename T,
            typename Info, typename InfoArg = Info, typename AccessorT = T>
  class ConcreteOptionWithInfo : public TypedSolverOption<T> {
   private:
    typedef AccessorT (Handler::*Get)(const SolverOption &, InfoArg) const;
    typedef void (Handler::*Set)(
        const SolverOption &,
        typename internal::OptionHelper<AccessorT>::Arg, InfoArg);

    Handler &handler_;
    Get get_;
    Set set_;
    Info info_;

   public:
    ConcreteOptionWithInfo(const char *name, const char* qualifiedName, 
        const char *description, Solver *s, Get get, Set set, InfoArg info, 
      ValueArrayRef values = ValueArrayRef())
    : TypedSolverOption<T>(name, qualifiedName, description, values),
      handler_(static_cast<Handler&>(*s)), get_(get), set_(set), info_(info) {}
    ConcreteOptionWithInfo(const char *name, const char* qualifiedName, const char *description, 
      Handler *s, Get get, Set set, InfoArg info, ValueArrayRef values = ValueArrayRef())
    : TypedSolverOption<T>(name, qualifiedName, description, values),
      handler_(*s), get_(get), set_(set), info_(info) {}

    void GetValue(T &value) const { value = (handler_.*get_)(*this, info_); }
    void SetValue(typename internal::OptionHelper<T>::Arg value) {
      (handler_.*set_)(*this, internal::OptionHelper<AccessorT>::CastArg(value),
                       info_);
    }
  };

  int GetObjNo(const SolverOption &) const { return std::abs(objno_); }
  void SetObjNo(const SolverOption &opt, int value) {
    if (value < 0)
      throw InvalidOptionValue(opt, value);
    objno_ = value;
  }


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

  // Solver flags
  enum {
    // Multiple solutions support.
    // Makes Solver register "countsolutions" and "solutionstub" options
    // and write every solution passed to HandleFeastibleSolution to a file
    // solutionstub & i & ".sol" where i is a solution number.
    MULTIPLE_SOL = 1,

    // Multiple objectives support.
    // Makes Solver register the "multiobj" option
    MULTIPLE_OBJ = 2
  };

 protected:
  // Constructs a Solver object.
  // date:  The solver date in YYYYMMDD format.
  // flags: Bitwise OR of zero or more of the following values
  //          MULTIPLE_SOL
  //          MULTIPLE_OBJ
  Solver(fmt::CStringRef name, fmt::CStringRef long_name, long date, int flags);

  void set_long_name(fmt::StringRef name) { long_name_ = name.to_string(); }
  void add_to_long_name(fmt::StringRef name) { long_name_ += name.to_string(); }
  void set_version(fmt::StringRef version) { version_ = version.to_string(); }
  void add_to_version(fmt::StringRef version) { version_ += version.to_string(); }

  // Sets the flags for Problem::Read.
  void set_read_flags(unsigned flags) { read_flags_ = flags; }

  // Sets a text to be displayed before option descriptions.
  void set_option_header(const char *header) { option_header_ = header; }

  // Add more text to be displayed before option descriptions.
  void add_to_option_header(const char *header_more) { option_header_ += header_more; }

  void AddOption(OptionPtr opt) {
    // First insert the option, then release a pointer to it. Doing the other
    // way around may lead to a memory leak if insertion throws an exception.
    if (!options_.insert(opt.get()).second)
      throw std::logic_error(
          fmt::format("Option {} already defined", opt.get()->name()));
    opt.release();
  }

  void AddOptionSynonym(const char* name, const char* realName)
  {
    SolverOption* real = FindOption(realName);
    if (!real)
      throw std::logic_error(
          fmt::format("Option {} referred to by {} is unknown", realName, name));
    OptionPtr opt = OptionPtr(new SolverOptionSynonym(name, *real));
    options_.insert(opt.get());
    opt.release();
  }

  // Adds an integer option.
  // The option stores pointers to the name and the description so make
  // sure that these strings have sufficient lifetimes (normally these are
  // string literals).
  // The arguments get and set should be pointers to member functions in the
  // solver class. They are used to get and set an option value respectively.
  template <typename Handler, typename Int>
  void AddIntOption(const char *name, const char* qualifiedName, 
    const char *description, Int (Handler::*get)(const SolverOption &) const,
      void (Handler::*set)(const SolverOption &, Int)) {
    AddOption(OptionPtr(new ConcreteOption<Handler, fmt::LongLong, Int>(
        name, qualifiedName, description, this, get, set)));
  }

  // Adds an integer option with additional information.
  // The option stores pointers to the name and the description so make
  // sure that these strings have sufficient lifetimes (normally these are
  // string literals).
  // The arguments get and set should be pointers to member functions in the
  // solver class. They are used to get and set an option value respectively.
  template <typename Handler, typename Info>
  void AddIntOption(const char *name, const char* qualifiedName, const char *description,
      int (Handler::*get)(const SolverOption &, const Info &) const,
      void (Handler::*set)(const SolverOption &, int, const Info &),
      const Info &info) {
    AddOption(OptionPtr(new ConcreteOptionWithInfo<
                        Handler, fmt::LongLong, Info, const Info &, int>(
                          name, qualifiedName, description, this, get, set, info)));
  }

  // The same as above but with Info argument passed by value.
  template <typename Handler, typename Info>
  void AddIntOption(const char *name, const char* qualifiedName, const char *description,
      int (Handler::*get)(const SolverOption &, Info) const,
      void (Handler::*set)(const SolverOption &, int, Info), Info info) {
    AddOption(OptionPtr(new ConcreteOptionWithInfo<
                        Handler, fmt::LongLong, Info, Info, int>(
                          name, qualifiedName, description, this, get, set, info)));
  }

  // Adds a double option.
  // The option stores pointers to the name and the description so make
  // sure that these strings have sufficient lifetimes (normally these are
  // string literals).
  // The arguments get and set should be pointers to member functions in the
  // solver class. They are used to get and set an option value respectively.
  template <typename Handler>
  void AddDblOption(const char *name, const char* qualifiedName, const char *description,
      double (Handler::*get)(const SolverOption &) const,
      void (Handler::*set)(const SolverOption &, double)) {
    AddOption(OptionPtr(new ConcreteOption<Handler, double>(
        name, qualifiedName, description, this, get, set)));
  }

  // Adds a double option with additional information.
  // The option stores pointers to the name and the description so make
  // sure that these strings have sufficient lifetimes (normally these are
  // string literals).
  // The arguments get and set should be pointers to member functions in the
  // solver class. They are used to get and set an option value respectively.
  template <typename Handler, typename Info>
  void AddDblOption(const char *name, const char* qualifiedName, const char *description,
      double (Handler::*get)(const SolverOption &, const Info &) const,
      void (Handler::*set)(const SolverOption &, double, const Info &),
      const Info &info) {
    AddOption(OptionPtr(
        new ConcreteOptionWithInfo<Handler, double, Info, const Info &>(
            name, qualifiedName, description, this, get, set, info)));
  }

  // The same as above but with Info argument passed by value.
  template <typename Handler, typename Info>
  void AddDblOption(const char *name, const char* qualifiedName, const char *description,
      double (Handler::*get)(const SolverOption &, Info) const,
      void (Handler::*set)(const SolverOption &, double, Info), Info info) {
    AddOption(OptionPtr(new ConcreteOptionWithInfo<Handler, double, Info>(
            name, qualifiedName, description, this, get, set, info)));
  }

  // Adds a string option.
  // The option stores pointers to the name and the description so make
  // sure that these strings have sufficient lifetimes (normally these are
  // string literals).
  // The arguments get and set should be pointers to member functions in the
  // solver class. They are used to get and set an option value respectively.
  template <typename Handler>
  void AddStrOption(const char *name, const char* qualifiedName, const char *description,
      std::string (Handler::*get)(const SolverOption &) const,
      void (Handler::*set)(const SolverOption &, fmt::StringRef),
      ValueArrayRef values = ValueArrayRef()) {
    AddOption(OptionPtr(new ConcreteOption<Handler, std::string>(
        name, qualifiedName, description, this, get, set, values)));
  }

  // Adds a string option with additional information.
  // The option stores pointers to the name and the description so make
  // sure that these strings have sufficient lifetimes (normally these are
  // string literals).
  // The arguments get and set should be pointers to member functions in the
  // solver class. They are used to get and set an option value respectively.
  template <typename Handler, typename Info>
  void AddStrOption(const char *name, const char* qualifiedName, const char *description,
      std::string (Handler::*get)(const SolverOption &, const Info &) const,
      void (Handler::*set)(const SolverOption &, fmt::StringRef, const Info &),
      const Info &info, ValueArrayRef values = ValueArrayRef()) {
    AddOption(OptionPtr(
        new ConcreteOptionWithInfo<Handler, std::string, Info, const Info &>(
            name, qualifiedName, description, this, get, set, info, values)));
  }

  // The same as above but with Info argument passed by value.
  template <typename Handler, typename Info>
  void AddStrOption(const char *name, const char* qualifiedName, const char *description,
      std::string (Handler::*get)(const SolverOption &, Info) const,
      void (Handler::*set)(const SolverOption &, fmt::StringRef, Info),
      Info info, ValueArrayRef values = ValueArrayRef()) {
    AddOption(OptionPtr(new ConcreteOptionWithInfo<Handler, std::string, Info>(
            name, qualifiedName, description, this, get, set, info, values)));
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

  // Possible values for the wantsol option (can be combined with bitwise OR).
  enum {
    WRITE_SOL_FILE      = 1,  // write .sol file
    PRINT_SOLUTION      = 2,  // print primal variables to stdout
    PRINT_DUAL_SOLUTION = 4,  // print dual variables to stdout
    SUPPRESS_SOLVER_MSG = 8   // suppress solver message
  };

  // Returns the value of the wantsol option which specifies what solution
  // information to write in a stand-alone invocation (no -AMPL on the
  // command line).
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

  /// Returns the index of the objective to optimize starting from 1,
  /// 0 to not use objective.
  /// Both multiobj and objno are used in NLReader to select
  /// the objective(s), solvers should not use these options.
  int objno() const { return std::abs(objno_); }

  /// Returns true if multiobjective optimization is enabled.
  /// Both multiobj and objno are used in NLReader to select
  /// the objective(s), solvers should not use these options.
  bool multiobj() const { return multiobj_ && objno_<0; }

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

  Interrupter *interrupter() { return interrupter_; }
  void set_interrupter(Interrupter *interrupter) {
    interrupter_ = interrupter ? interrupter : this;
  }

  const char *solution_stub() const { return solution_stub_.c_str(); }

  bool need_multiple_solutions() const {
    return count_solutions_ || !solution_stub_.empty();
  }

  // Returns the number of options.
  int num_options() const { return static_cast<int>(options_.size()); }

  // Finds an option and returns a pointer to it if found or null otherwise.
  SolverOption *FindOption(const char *name) const;

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
      char **argv, unsigned flags = 0, const ASLProblem *p = 0);

  // Returns the value of an integer option.
  // Throws OptionError if there is no such option or it has a different type.
  fmt::LongLong GetIntOption(const char *name) const {
    return GetOption(name)->GetValue<fmt::LongLong>();
  }

  // Sets the value of an integer option.
  // Throws OptionError if there is no such option or it has a different type.
  void SetIntOption(const char *name, fmt::LongLong value) {
    GetOption(name)->SetValue(value);
  }

  // Returns the value of a double option.
  // Throws OptionError if there is no such option or it has a different type.
  double GetDblOption(const char *name) const {
    return GetOption(name)->GetValue<double>();
  }

  // Sets the value of a double option.
  // Throws OptionError if there is no such option or it has a different type.
  void SetDblOption(const char *name, double value) {
    GetOption(name)->SetValue(value);
  }

  // Returns the value of a string option.
  // Throws OptionError if there is no such option or it has a different type.
  std::string GetStrOption(const char *name) const {
    return GetOption(name)->GetValue<std::string>();
  }

  // Sets the value of a string option.
  // Throws OptionError if there is no such option or it has a different type.
  void SetStrOption(const char *name, fmt::StringRef value) {
    GetOption(name)->SetValue(value);
  }

  const SuffixList &suffixes() const { return suffixes_; }

  // Reports an error printing the formatted error message to stderr.
  // Usage: ReportError("File not found: {}") << filename;
  void ReportError(fmt::CStringRef format, const fmt::ArgList &args) {
    has_errors_ = true;
    fmt::MemoryWriter w;
    w.write(format, args);
    error_handler_->HandleError(w.c_str());
  }
  FMT_VARIADIC(void, ReportError, fmt::CStringRef)

  // Formats a string and prints it to stdout or, if an output handler
  // is registered, sends it to the output handler.
  void Print(fmt::CStringRef format, const fmt::ArgList &args) {
    fmt::MemoryWriter w;
    w.write(format, args);
    output_handler_->HandleOutput(w.c_str());
  }
  FMT_VARIADIC(void, Print, fmt::CStringRef)

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
  typedef internal::NLProblemBuilder<ProblemBuilder> NLProblemBuilder;

  SolverImpl(fmt::CStringRef name, fmt::CStringRef long_name = 0,
             long date = 0, int flags = 0)
    : Solver(name, long_name, date, flags) {}
};

// Adapts a solution for WriteSol.
template <typename ProblemBuilder>
class SolutionAdapter {
 private:
  int status_;
  ProblemBuilder *builder_;
  const char *message_;
  mp::ArrayRef<int> options_;
  mp::ArrayRef<double> values_;
  mp::ArrayRef<double> dual_values_;
  int objno_;

 public:
  SolutionAdapter(int status, ProblemBuilder *pb, const char *message,
                  mp::ArrayRef<int> options, mp::ArrayRef<double> values,
                  mp::ArrayRef<double> dual_values,
                  int on)
    : status_(status), builder_(pb), message_(message), options_(options),
      values_(values), dual_values_(dual_values), objno_(on) {}

  int status() const { return status_; }

  const char *message() const { return message_; }

  int num_options() const { return static_cast<int>(options_.size()); }
  int option(int index) const { return options_[index]; }

  int num_values() const { return static_cast<int>(values_.size()); }
  double value(int index) const { return values_[index]; }

  int num_dual_values() const { return static_cast<int>(dual_values_.size()); }
  double dual_value(int index) const { return dual_values_[index]; }

  int objno() const { return objno_; }

  const typename ProblemBuilder::SuffixSet *suffixes(suf::Kind kind) const {
    return builder_ ? &builder_->suffixes(kind) : 0;
  }
};

class NullSolutionHandler : public SolutionHandler {
 public:
  void HandleFeasibleSolution(
        fmt::CStringRef, const double *, const double *, double) {}
  void HandleSolution(
        int, fmt::CStringRef, const double *, const double *, double) {}
};

// The default .sol file writer.
class SolFileWriter {
 public:
  template <typename Solution>
  void Write(fmt::CStringRef filename, const Solution &sol) {
    WriteSolFile(filename, sol);
  }
};

// A solution writer.
// Solver: optimization solver class
// Writer: .sol writer
template <typename Solver, typename Writer = SolFileWriter>
class SolutionWriter : private Writer, public SolutionHandler {
 private:
  std::string stub_;
  Solver &solver_;

  typedef typename Solver::ProblemBuilder ProblemBuilder;
  ProblemBuilder &builder_;

  ArrayRef<int> options_;

  // The number of feasible solutions found.
  int num_solutions_;

 protected:
  Solver &solver() { return solver_; }
  ProblemBuilder &builder() { return builder_; }
  const std::string &stub() const { return stub_; }

 public:
  SolutionWriter(fmt::StringRef stub, Solver &s, ProblemBuilder &b,
                 ArrayRef<int> options = mp::ArrayRef<int>(0, 0))
    : stub_(stub.to_string()), solver_(s), builder_(b),
      options_(options), num_solutions_(0) {}

  // Returns the .sol writer.
  Writer &sol_writer() { return *this; }

  void HandleFeasibleSolution(fmt::CStringRef message,
        const double *values, const double *dual_values, double);

  // Writes the solution to a .sol file.
  void HandleSolution(int status, fmt::CStringRef message,
        const double *values, const double *dual_values, double);
};

template <typename Solver, typename Writer>
void SolutionWriter<Solver, Writer>::HandleFeasibleSolution(
    fmt::CStringRef message, const double *values,
    const double *dual_values, double) {
  ++num_solutions_;
  const char *solution_stub = solver_.solution_stub();
  if (!*solution_stub)
    return;
  SolutionAdapter<ProblemBuilder> sol(
        sol::UNCERTAIN, 0, message.c_str(), ArrayRef<int>(0, 0),
        MakeArrayRef(values, values ? builder_.num_vars() : 0),
        MakeArrayRef(dual_values,
                     dual_values ? builder_.num_algebraic_cons() : 0),
        solver_.objno());
  fmt::MemoryWriter filename;
  filename << solution_stub << num_solutions_ << ".sol";
  this->Write(filename.c_str(), sol);
}

template <typename Solver, typename Writer>
void SolutionWriter<Solver, Writer>::HandleSolution(
    int status, fmt::CStringRef message, const double *values,
    const double *dual_values, double) {
  typedef typename ProblemBuilder::IntSuffix IntSuffix;
  if (solver_.need_multiple_solutions()) {
    builder_.AddIntSuffix("nsol", suf::PROBLEM | suf::OUTPUT | suf::OUTONLY, 0);

    IntSuffix nsol_suffix = Cast<IntSuffix>(
          builder_.suffixes(suf::PROBLEM).Find("nsol"));
    nsol_suffix.set_value(0, num_solutions_);
  }
  SolutionAdapter<ProblemBuilder> sol(
        status, &builder_, message.c_str(), options_,
        MakeArrayRef(values, values ? builder_.num_vars() : 0),
        MakeArrayRef(dual_values,
                     dual_values ? builder_.num_algebraic_cons() : 0),
        solver_.objno());
  this->Write(stub_ + ".sol", sol);
}

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
  Solver &solver_;
  std::string message_;
  static volatile std::sig_atomic_t stop_;

  static atomic<const char*> signal_message_ptr_;
  static atomic<unsigned> signal_message_size_;
  static atomic<InterruptHandler> handler_;
  static atomic<void*> data_;

  static void HandleSigInt(int sig);

  SignalRepeater repeater_;

 public:
  explicit SignalHandler(Solver &s);

  ~SignalHandler();

  //static bool stop() { return stop_ != 0; }

  // Returns true if the execution should be stopped due to SIGINT.
  bool Stop() const { return stop_ != 0; }

  void SetHandler(InterruptHandler handler, void *data);
};

// An .nl handler for SolverApp.
template <typename Solver>
class SolverNLHandler : public Solver::NLProblemBuilder {
 private:
  Solver &solver_;
  int num_options_;
  int options_[MAX_AMPL_OPTIONS];

  typedef typename Solver::NLProblemBuilder Base;

 public:
  SolverNLHandler(typename Solver::ProblemBuilder &pb, Solver &s)
    : Base(pb), solver_(s), num_options_(0) {}

  int num_options() const { return num_options_; }
  const int *options() const { return options_; }

  int objno() const override { return solver_.objno(); }
  bool multiobj() const override { return solver_.multiobj(); }

  void OnHeader(const NLHeader &h);
};

template <typename Solver>
void SolverNLHandler<Solver>::OnHeader(const NLHeader &h) {
  int objno = solver_.objno();
  if (objno > h.num_objs)
    throw InvalidOptionValue("objno", objno,
                             fmt::format("expected value between 0 and {}", h.num_objs));
  num_options_ = h.num_ampl_options;
  std::copy(h.ampl_options, h.ampl_options + num_options_, options_);
  Base::OnHeader(h);
#ifndef MP_DATE
  internal::CheckDemoVersion(h);
#endif
}

// A variable or constraint name provider.
class NameProvider {
 private:
  std::vector<const char *> names_;
  std::string gen_name_;
  NameReader reader_;
  fmt::MemoryWriter writer_;

 public:
  NameProvider(fmt::CStringRef filename, fmt::CStringRef gen_name,
               std::size_t num_items);

  // Returns the name of the item at specified index.
  fmt::StringRef name(std::size_t index);
};

// Prints a solution to stdout.
void PrintSolution(const double *values, int num_values, const char *name_col,
                   const char *value_col, NameProvider &np);

// Solution handler for a solver application.
template <typename Solver, typename Writer = SolFileWriter>
class AppSolutionHandler : public SolutionWriter<Solver, Writer> {
 private:
  unsigned banner_size_;

 public:
  AppSolutionHandler(fmt::StringRef stub, Solver &s,
                    typename Solver::ProblemBuilder &b,
                    ArrayRef<int> options, unsigned banner_size)
  : SolutionWriter<Solver, Writer>(stub, s, b, options),
    banner_size_(banner_size) {}

  void HandleSolution(int status, fmt::CStringRef message, const double *values,
                      const double *dual_values, double obj_value);
};

template <typename Solver, typename Writer>
void AppSolutionHandler<Solver, Writer>::HandleSolution(
    int status, fmt::CStringRef message, const double *values,
    const double *dual_values, double obj_value) {
  Solver &solver = this->solver();
  int wantsol = solver.wantsol();
  if (solver.ampl_flag() || (wantsol & Solver::WRITE_SOL_FILE) != 0) {
    // "Erase" the banner so that it is not duplicated when printing
    // the solver message.
    if (solver.ampl_flag() && banner_size_ != 0) {
      fmt::MemoryWriter w;
      w << fmt::pad("", banner_size_, '\b');
      solver.Print("{}", w.c_str());
    }
    SolutionWriter<Solver, Writer>::HandleSolution(
          status, message, values, dual_values, obj_value);
  }
  if (solver.ampl_flag())
    return;
  if ((wantsol & Solver::SUPPRESS_SOLVER_MSG) == 0)
    solver.Print("{}\n", message.c_str() + banner_size_);
  using internal::PrintSolution;
  if ((wantsol & Solver::PRINT_SOLUTION) != 0) {
    int num_vars = this->builder().num_vars();
    internal::NameProvider np(this->stub() + ".col", "_svar", num_vars);
    PrintSolution(values, num_vars, "variable", "value", np);
  }
  if ((wantsol & Solver::PRINT_DUAL_SOLUTION) != 0) {
    int num_cons = this->builder().num_algebraic_cons();
    internal::NameProvider np(this->stub() + ".row", "_scon", num_cons);
    PrintSolution(dual_values, num_cons, "constraint", "dual value", np);
  }
}

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
}  // namespace internal

// A solver application.
// Solver: optimization solver class; normally a subclass of SolverImpl
// Reader: .nl reader
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

#ifdef MP_USE_UNIQUE_PTR
typedef std::unique_ptr<Solver> SolverPtr;
#else
typedef std::auto_ptr<Solver> SolverPtr;
inline SolverPtr move(SolverPtr p) { return p; }
#endif
}  // namespace mp

#endif  // MP_SOLVER_H_
