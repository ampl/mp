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

#include <algorithm>
#include <deque>
#include <limits>
#include <memory>

extern "C" {
#include "solvers/getstub.h"
#undef Char
}

#include "solvers/util/expr.h"

namespace ampl {

template <typename Grad>
class LinearExpr;

// A linear term.
template <typename GradT>
class LinearTerm {
 private:
  typedef GradT Grad;
  Grad *grad_;

  friend class LinearExpr< LinearTerm<Grad> >;

  explicit LinearTerm(Grad *g) : grad_(g) {}

 public:
  // Returns the coefficient.
  double coef() const { return grad_->coef; }

  // Returns the variable index.
  int var_index() const { return grad_->varno; }
};

typedef LinearTerm<ograd> LinearObjTerm;
typedef LinearTerm<cgrad> LinearConTerm;

class Problem;

// A linear expression.
template <typename Term>
class LinearExpr : public std::iterator<std::forward_iterator_tag, Term> {
 private:
  Term first_term_;

  friend class Problem;

  explicit LinearExpr(typename Term::Grad *first_term)
  : first_term_(Term(first_term)) {}

 public:
  class iterator {
   private:
    Term term_;

   public:
    explicit iterator(Term t) : term_(t) {}

    Term operator*() const { return term_; }
    const Term *operator->() const { return &term_; }

    iterator &operator++() {
      term_ = Term(term_.grad_->next);
      return *this;
    }

    iterator operator++(int ) {
      iterator it(*this);
      term_ = Term(term_.grad_->next);
      return it;
    }

    bool operator==(iterator other) const {
      return term_.grad_ == other.term_.grad_;
    }
    bool operator!=(iterator other) const {
      return term_.grad_ != other.term_.grad_;
    }
  };

  iterator begin() { return iterator(first_term_); }
  iterator end() { return iterator(Term(0)); }
};

typedef LinearExpr<LinearObjTerm> LinearObjExpr;
typedef LinearExpr<LinearConTerm> LinearConExpr;

// An objective type.
enum ObjType { MIN = 0, MAX = 1 };

class Noncopyable {
 private:
  // Do not implement!
  Noncopyable(const Noncopyable &);
  void operator=(const Noncopyable &);

 public:
  Noncopyable() {}
};

// A solution of an optimization problem.
class Solution : Noncopyable {
 private:
  int solve_code_;
  int num_vars_;
  int num_cons_;
  double *values_;
  double *dual_values_;

 public:
  // Solution status.
  enum Status {
    UNKNOWN,
    SOLVED,
    SOLVED_MAYBE,
    INFEASIBLE,
    UNBOUNDED,
    LIMIT,
    FAILURE
  };

  // Constructs a solution with zero variables and constraints and the
  // solve code -1.
  Solution();

  ~Solution();

  // Swaps this solution with other.
  void Swap(Solution &other);

  // Returns the solution status.
  Status status() const {
    return solve_code_ < 0 || solve_code_ >= 600 ?
        UNKNOWN : static_cast<Status>(solve_code_ / 100 + 1);
  }

  // Returns the solve code.
  int solve_code() const { return solve_code_; }

  // Returns the number of variables.
  int num_vars() const { return num_vars_; }

  // Returns the number of constraints.
  int num_cons() const { return num_cons_; }

  // Returns the values of all variables.
  const double *values() const { return values_; }

  // Returns the values of all dual variables.
  const double *dual_values() const { return dual_values_; }

  // Returns the value of a variable.
  double value(int var) const {
    assert(var >= 0 && var < num_vars_);
    return values_[var];
  }

  // Returns the value of a dual variable corresponding to constraint con.
  double dual_value(int con) const {
    assert(con >= 0 && con < num_cons_);
    return dual_values_[con];
  }

  // Reads a solution from the file <stub>.sol.
  void Read(const char *stub, int num_vars, int num_cons);
};

class ProblemChanges;

// An optimization problem.
class Problem : Noncopyable {
 private:
  ASL_fg *asl_;

  friend class BasicSolver;

 public:
  Problem();
  ~Problem();

  // Returns the number of variables.
  int num_vars() const { return asl_->i.n_var_; }

  // Returns the number of objectives.
  int num_objs() const { return asl_->i.n_obj_; }

  // Returns the number of constraints excluding logical constraints.
  int num_cons() const { return asl_->i.n_con_; }

  // Returns the number of integer variables including binary.
  int num_integer_vars() const {
    return asl_->i.nbv_ + asl_->i.niv_ + asl_->i.nlvbi_ +
        asl_->i.nlvci_ + asl_->i.nlvoi_;
  }

  // Returns the number of continuous variables.
  int num_continuous_vars() const {
    return num_vars() - num_integer_vars();
  }

  // Returns the number of nonlinear objectives.
  int num_nonlinear_objs() const { return asl_->i.nlo_; }

  // Returns the number of nonlinear constraints.
  int num_nonlinear_cons() const { return asl_->i.nlc_; }

  // Returns the number of logical constraints.
  int num_logical_cons() const { return asl_->i.n_lcon_; }

  // Returns the lower bounds for the variables.
  const double *var_lb() const { return asl_->i.LUv_; }

  // Returns the lower bound for the variable.
  double var_lb(int var_index) const {
    assert(var_index >= 0 && var_index < num_vars());
    return asl_->i.LUv_[var_index];
  }

  // Returns the upper bounds for the variables.
  const double *var_ub() const { return asl_->i.Uvx_; }

  // Returns the upper bound for the variable.
  double var_ub(int var_index) const {
    assert(var_index >= 0 && var_index < num_vars());
    return asl_->i.Uvx_[var_index];
  }

  // Returns the lower bounds for the constraints.
  const double *con_lb() const { return asl_->i.LUrhs_; }

  // Returns the lower bound for the constraint.
  double con_lb(int con_index) const {
    assert(con_index >= 0 && con_index < num_cons());
    return asl_->i.LUrhs_[con_index];
  }

  // Returns the upper bounds for the constraints.
  const double *con_ub() const { return asl_->i.Urhsx_; }

  // Returns the upper bound for the constraint.
  double con_ub(int con_index) const {
    assert(con_index >= 0 && con_index < num_cons());
    return asl_->i.Urhsx_[con_index];
  }

  // Returns the objective type.
  ObjType obj_type(int obj_index) const {
    assert(obj_index >= 0 && obj_index < num_objs());
    return static_cast<ObjType>(asl_->i.objtype_[obj_index]);
  }

  // Returns the linear part of an objective expression.
  LinearObjExpr linear_obj_expr(int obj_index) const {
    assert(obj_index >= 0 && obj_index < num_objs());
    return LinearObjExpr(asl_->i.Ograd_[obj_index]);
  }

  // Returns the linear part of a constraint expression.
  LinearConExpr linear_con_expr(int con_index) const {
    assert(con_index >= 0 && con_index < num_cons());
    return LinearConExpr(asl_->i.Cgrad_[con_index]);
  }

  // Returns the nonlinear part of an objective expression.
  NumericExpr nonlinear_obj_expr(int obj_index) const {
    assert(obj_index >= 0 && obj_index < num_objs());
    return Expr::Create<NumericExpr>(asl_->I.obj_de_[obj_index].e);
  }

  // Returns the nonlinear part of a constraint expression.
  NumericExpr nonlinear_con_expr(int con_index) const {
    assert(con_index >= 0 && con_index < num_cons());
    return Expr::Create<NumericExpr>(asl_->I.con_de_[con_index].e);
  }

  // Returns a logical constraint expression.
  LogicalExpr logical_con_expr(int lcon_index) const {
    assert(lcon_index >= 0 && lcon_index < num_logical_cons());
    return Expr::Create<LogicalExpr>(asl_->I.lcon_de_[lcon_index].e);
  }

  // Returns the solve code.
  int solve_code() const { return asl_->p.solve_code_; }

  // Sets the solve code.
  void set_solve_code(int value) {
    asl_->p.solve_code_ = value;
  }

  // Flags for Solve.
  enum { IGNORE_FUNCTIONS = 1 };

  // Solves the problem.
  void Solve(const char *solver_name, Solution &sol,
      ProblemChanges *pc = 0, unsigned flags = 0);
};

// Writes a problem in the AMPL format.
fmt::Writer &operator<<(fmt::Writer &w, const Problem &p);

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

// Changes (additions) to an optimization problem.
class ProblemChanges {
 private:
  const Problem *problem_;
  std::vector<double> var_lb_;
  std::vector<double> var_ub_;
  std::vector<double> con_lb_;
  std::vector<double> con_ub_;
  std::deque<ograd> con_terms_;
  std::deque<ograd> obj_terms_;
  std::vector<ograd*> cons_;
  std::vector<ograd*> objs_;
  std::vector<char> obj_types_;
  NewVCO vco_;

  friend class Problem;

  NewVCO *vco();

 public:
  explicit ProblemChanges(const Problem &p) : problem_(&p), vco_() {}

  // Returns the number of additional variables.
  int num_vars() const { return static_cast<int>(var_lb_.size()); }

  // Returns the number of additional constraints.
  int num_cons() const { return static_cast<int>(cons_.size()); }

  // Adds a variable.
  int AddVar(double lb, double ub) {
    var_lb_.push_back(lb);
    var_ub_.push_back(ub);
    return static_cast<int>(problem_->num_vars() + var_lb_.size() - 1);
  }

  // Adds an objective.
  void AddObj(ObjType type,
      unsigned size, const double *coefs, const int *vars);

  // Adds a constraint.
  void AddCon(const double *coefs, double lb, double ub);
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

// Base class for all solver classes.
class BasicSolver
  : private ErrorHandler, private SolutionHandler, private Option_Info {
 private:
  Problem problem_;

  std::string name_;
  std::string long_name_;
  std::string options_var_name_;
  std::string version_;
  std::string version_desc_;
  std::string wantsol_desc_;

  std::vector<keyword> keywords_;
  bool options_sorted_;

  bool has_errors_;
  ErrorHandler *error_handler_;
  SolutionHandler *sol_handler_;

  void SortOptions();

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

 protected:
  // Make Option_Info accessible in subclasses despite private inheritance.
  typedef Option_Info Option_Info;

  // Constructs a BasicSolver object.
  // date: The solver date in YYYYMMDD format.
  BasicSolver(fmt::StringRef name, fmt::StringRef long_name, long date);
  virtual ~BasicSolver();

  void set_long_name(fmt::StringRef name) {
    long_name_ = name;
    bsname = const_cast<char*>(long_name_.c_str());
  }

  void set_version(fmt::StringRef version) {
    version_ = version;
    Option_Info::version = const_cast<char*>(version_.c_str());
  }

  template <typename SolverT>
  static SolverT *GetSolver(Option_Info *oi) {  // throw()
    return static_cast<SolverT*>(oi);
  }

  // Parses solver options.
  // Returns true if there were no errors and false otherwise.
  bool DoParseOptions(char **argv, unsigned flags) {
    has_errors_ = false;
    SortOptions();
    ASL *asl = reinterpret_cast<ASL*>(problem_.asl_);
    if ((flags & NO_OPTION_ECHO) != 0) {
      option_echo |= ASL_OI_echothis;
      option_echo &= ~ASL_OI_echo;
    } else {
      option_echo |= ASL_OI_echo;
    }
    return getopts_ASL(asl, argv, this) == 0 && !has_errors_;
  }

  void AddKeyword(const char *name,
      const char *description, Kwfunc func, const void *info);

  // Formats an option description by indenting it and performing word wrap.
  static std::string FormatDescription(const char *description);

 public:
  // Flags for DoParseOptions.
  enum {
    // Don't echo options during parsing.
    NO_OPTION_ECHO = 1
  };

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

  // Parses command-line options starting with '-' and reads a problem from
  // an .nl file if the file name (stub) is specified. Returns true if the
  // arguments contain the file name and false otherwise. If there was an
  // error parsing arguments or reading the problem ReadProblem will print an
  // error message and call std::exit (this is likely to change in the future
  // version).
  bool ReadProblem(char **&argv);

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
};

// An AMPL solver.
// OptionHandler is a class that will receive notification about
// parsed options. Often OptionHandler is a class derived from Solver,
// but this is not required.
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
template <typename OptionHandler>
class Solver : public BasicSolver {
 private:
  class Option {
   private:
    std::string description_;

   public:
    Option(const char *description)
    : description_(FormatDescription(description)) {}

    virtual ~Option() {}

    const char *description() const { return description_.c_str(); }

    virtual void Handle(
        OptionHandler &h, Option_Info *oi, keyword *kw, char *&value) = 0;
  };

  template <typename Func, typename Value>
  class ConcreteOption : public Option {
   private:
    Func func_;

   public:
    ConcreteOption(const char *description, Func func)
    : Option(description), func_(func) {}

    void Handle(OptionHandler &h, Option_Info *oi, keyword *kw, char *&s) {
      (h.*func_)(kw->name, internal::OptionParser<Value>()(oi, kw, s));
    }
  };

  template <typename Func, typename Info, typename Value>
  class ConcreteOptionWithInfo : public Option {
   private:
    Func func_;
    Info info_;

   public:
    ConcreteOptionWithInfo(const char *description, Func func, const Info &info)
    : Option(description), func_(func), info_(info) {}

    void Handle(OptionHandler &h, Option_Info *oi, keyword *kw, char *&s) {
      (h.*func_)(kw->name, internal::OptionParser<Value>()(oi, kw, s), info_);
    }
  };

  OptionHandler *handler_;
  std::vector<Option*> options_;

  static char *HandleOption(Option_Info *oi, keyword *kw, char *value) {
    Solver *self = GetSolver<Solver>(oi);
    try {
      Option *opt = self->options_[reinterpret_cast<size_t>(kw->info)];
      opt->Handle(*self->handler_, oi, kw, value);
    } catch (const std::exception &e) {
      self->ReportError(e.what());
    } catch (...) {
      self->ReportError("Unknown exception in option handler");
    }
    return value;
  }

  void AddOption(const char *name, std::auto_ptr<Option> opt) {
    AddKeyword(name, opt->description(), HandleOption,
        reinterpret_cast<void*>(options_.size()));
    options_.push_back(0);
    options_.back() = opt.release();
  }

  struct Deleter {
    void operator()(Option *opt) { delete opt; }
  };

 protected:
  // Adds an integer option. The argument f should be a pointer to a member
  // function in the Handler class. This function is called after the option
  // is parsed:
  //   (handler.*f)(name, value);
  // where handler is a reference to a Handler object, name is an option name
  // of type "const char*" and value is an option value of type "int".
  template <typename Func>
  void AddIntOption(const char *name, const char *description, Func f) {
    AddOption(name, std::auto_ptr<Option>(
        new ConcreteOption<Func, int>(description, f)));
  }

  // Adds an integer option with additional information. The argument f
  // should be a pointer to a member function in the Handler class.
  // This function is called after the option is parsed:
  //   (handler.*f)(name, value, info);
  // where handler is a reference to a Handler object, name is an option name
  // of type "const char*", value is an option value of type "int" and info
  // is a reference to the Info object passed to AddIntOption.
  template <typename Info, typename Func>
  void AddIntOption(const char *name,
      const char *description, Func f, const Info &info) {
    AddOption(name, std::auto_ptr<Option>(
        new ConcreteOptionWithInfo<Func, Info, int>(description, f, info)));
  }

  // Adds a double option. The argument f should be a pointer to a member
  // function in the Handler class. This function is called after the option
  // is parsed:
  //   (handler.*f)(name, value);
  // where handler is a reference to a Handler object, name is an option name
  // of type "const char*" and value is an option value of type "double".
  template <typename Func>
  void AddDblOption(const char *name, const char *description, Func f) {
    AddOption(name, std::auto_ptr<Option>(
        new ConcreteOption<Func, double>(description, f)));
  }

  // Adds a double option with additional information. The argument f
  // should be a pointer to a member function in the Handler class.
  // This function is called after the option is parsed:
  //   (handler.*f)(name, value, info);
  // where handler is a reference to a Handler object, name is an option name
  // of type "const char*", value is an option value of type "double" and info
  // is a reference to the Info object passed to AddDblOption.
  template <typename Info, typename Func>
  void AddDblOption(const char *name,
      const char *description, Func f, const Info &info) {
    AddOption(name, std::auto_ptr<Option>(
        new ConcreteOptionWithInfo<Func, Info, double>(description, f, info)));
  }

  // Adds a string option. The argument f should be a pointer to a member
  // function in the Handler class. This function is called after the option
  // is parsed:
  //   (handler.*f)(name, value);
  // where handler is a reference to a Handler object, name is an option name
  // of type "const char*" and value is an option value of type "const char*".
  template <typename Func>
  void AddStrOption(const char *name, const char *description, Func f) {
    AddOption(name, std::auto_ptr<Option>(
        new ConcreteOption<Func, const char*>(description, f)));
  }

  // Adds a string option with additional information. The argument f
  // should be a pointer to a member function in the Handler class.
  // This function is called after the option is parsed:
  //   (handler.*f)(name, value, info);
  // where handler is a reference to a Handler object, name is an option name
  // of type "const char*", value is an option value of type "const char*"
  // and info is a reference to the Info object passed to AddStrOption.
  template <typename Info, typename Func>
  void AddStrOption(const char *name,
      const char *description, Func f, const Info &info) {
    AddOption(name, std::auto_ptr<Option>(
        new ConcreteOptionWithInfo<Func, Info, const char*>(
            description, f, info)));
  }

 public:
  Solver(fmt::StringRef name, fmt::StringRef long_name = 0, long date = 0)
  : BasicSolver(name, long_name, date), handler_(0) {}

  ~Solver() {
    std::for_each(options_.begin(), options_.end(), Deleter());
  }

  // Parses solver options and returns true if there were no errors and
  // false otherwise. Note that handler functions can report errors with
  // BasicSolver::ReportError and ParseOptions will take them into account
  // as well, returning false if there was at least one such error.
  bool ParseOptions(char **argv, OptionHandler &h, unsigned flags = 0) {
    handler_ = &h;
    return DoParseOptions(argv, flags);
  }

  // Processes command-line arguments, reads a problem from an .nl file
  // if the file name (stub) is specified and parses solver options.
  // Returns true if the arguments contain the file name and options were
  // parsed successfully; false otherwise.
  // If there was an  error parsing arguments or reading the problem
  // ProcessArgs will print an error message and call std::exit (this is
  // likely to change in the future version).
  bool ProcessArgs(char **&argv, OptionHandler &h, unsigned flags = 0) {
    return ReadProblem(argv) && ParseOptions(argv, h, flags);
  }
};
}

#endif  // SOLVERS_UTIL_SOLVER_H_
