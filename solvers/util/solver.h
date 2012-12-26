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

#include <cstring>
#include <memory>
#include <sstream>

extern "C" {
#include "solvers/getstub.h"
}

#include "solvers/util/expr.h"
#include "solvers/util/format.h"

namespace ampl {

// An AMPL problem.
class Problem {
 private:
  ASL_fg *asl_;

  // Do not implement.
  Problem(const Problem&);
  Problem& operator=(const Problem&);

  friend class SolverBase;

 public:
  Problem();
  virtual ~Problem();

  // Returns the number of variables.
  int num_vars() const { return asl_->i.n_var_; }

  // Returns the number of objectives.
  int num_objs() const { return asl_->i.n_obj_; }

  // Returns the number of constraints.
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

  // Returns the variable lower bound.
  double GetVarLB(int var_index) const {
    assert(var_index >= 0 && var_index < num_vars());
    return asl_->i.LUv_[var_index];
  }

  // Returns the variable lower bound.
  double GetVarUB(int var_index) const {
    assert(var_index >= 0 && var_index < num_vars());
    return asl_->i.Uvx_[var_index];
  }

  // Returns the constraint lower bound.
  double GetConLB(int con_index) const {
    assert(con_index >= 0 && con_index < num_cons());
    return asl_->i.LUrhs_[con_index];
  }

  // Returns the constraint lower bound.
  double GetConUB(int con_index) const {
    assert(con_index >= 0 && con_index < num_cons());
    return asl_->i.Urhsx_[con_index];
  }

  enum ObjType { MIN = 0, MAX = 1 };

  // Returns the objective type.
  ObjType GetObjType(int obj_index) const {
    assert(obj_index >= 0 && obj_index < num_objs());
    return static_cast<ObjType>(asl_->i.objtype_[obj_index]);
  }

  // Returns the linear part of an objective expression.
  ograd *GetLinearObjExpr(int obj_index) const {
    assert(obj_index >= 0 && obj_index < num_objs());
    return asl_->i.Ograd_[obj_index];
  }

  // Returns the linear part of a constraint expression.
  cgrad *GetLinearConExpr(int con_index) const {
    assert(con_index >= 0 && con_index < num_cons());
    return asl_->i.Cgrad_[con_index];
  }

  // Returns the nonlinear part of an objective expression.
  NumericExpr GetNonlinearObjExpr(int obj_index) const {
    assert(obj_index >= 0 && obj_index < num_objs());
    return Expr::Create<NumericExpr>(asl_->I.obj_de_[obj_index].e);
  }

  // Returns the nonlinear part of a constraint expression.
  NumericExpr GetNonlinearConExpr(int con_index) const {
    assert(con_index >= 0 && con_index < num_cons());
    return Expr::Create<NumericExpr>(asl_->I.con_de_[con_index].e);
  }

  // Returns a logical constraint expression.
  LogicalExpr GetLogicalConExpr(int lcon_index) const {
    assert(lcon_index >= 0 && lcon_index < num_logical_cons());
    return Expr::Create<LogicalExpr>(asl_->I.lcon_de_[lcon_index].e);
  }

  // Returns the solve code.
  int solve_code() const { return asl_->p.solve_code_; }

  // Sets the solve code.
  void SetSolveCode(int value) {
    asl_->p.solve_code_ = value;
  }
};

// Formats a double with objective precision.
// Usage: fmt::Format("{0}") << ObjPrec(42.0);
class ObjPrec {
 private:
  double value_;

 public:
  explicit ObjPrec(double value) : value_(value) {}

  friend void Format(
      fmt::ArgFormatter &af, const fmt::FormatSpec &spec, ObjPrec op) {
    char buffer[32];
    g_fmtop(buffer, op.value_);
    af.Write(buffer, spec);
  }
};

class SolverBase;

class SolutionHandler {
 protected:
  ~SolutionHandler() {}

 public:
  virtual void HandleSolution(SolverBase &d, fmt::StringRef message,
      const double *primal, const double *dual, double obj_value) = 0;
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

class DummyOptionHandler {};

// Base class for all solver classes.
class SolverBase : private SolutionHandler, private Option_Info {
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

  SolutionHandler *sol_handler_;
  bool has_errors_;

  void SortOptions();

  void HandleSolution(SolverBase &, fmt::StringRef message,
        const double *primal, const double *dual, double) {
    write_sol_ASL(reinterpret_cast<ASL*>(problem_.asl_),
        const_cast<char*>(message.c_str()), const_cast<double*>(primal),
        const_cast<double*>(dual), this);
  }

  struct PrintError {
    void operator()(const fmt::Formatter &f) const {
      std::fputs(f.c_str(), stderr);
      std::fputc('\n', stderr);
    }
  };

 protected:
  static DummyOptionHandler dummy_option_handler;

  // Constructs a SolverBase object.
  // date: The solver date in YYYYMMDD format.
  explicit SolverBase(fmt::StringRef name,
      fmt::StringRef long_name = 0, long date = 0);
  virtual ~SolverBase();

  void set_long_name(fmt::StringRef name) {
    long_name_ = name;
    bsname = const_cast<char*>(long_name_.c_str());
  }

  void set_version(fmt::StringRef version) {
    version_ = version;
    Option_Info::version = const_cast<char*>(version_.c_str());
  }

  void EnableOptionEcho(int mask) { option_echo |= mask; }
  void DisableOptionEcho(int mask) { option_echo &= ~mask; }

  template <typename SolverT>
  static SolverT *GetSolver(Option_Info *oi) {
    return static_cast<SolverT*>(oi);
  }

  bool ParseOptions(char **argv) {
    has_errors_ = false;
    SortOptions();
    ASL *asl = reinterpret_cast<ASL*>(problem_.asl_);
    return getopts_ASL(asl, argv, this) == 0 && !has_errors_;
  }

  void AddKeyword(const char *name,
      const char *description, Kwfunc func, const void *info);

  // Formats an option description by indenting it and performing word wrap.
  static std::string FormatDescription(const char *description);

 public:
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

  // Returns the solution handler.
  SolutionHandler *solution_handler() { return sol_handler_; }

  // Sets the solution handler.
  void set_solution_handler(SolutionHandler *sh) { sol_handler_ = sh; }

  // Reads the problem from an .nl file.
  bool ReadProblem(char **&argv);

  // Passes a solution to the solution handler.
  void HandleSolution(fmt::StringRef message,
      const double *primal, const double *dual,
      double obj_value = std::numeric_limits<double>::quiet_NaN()) {
    sol_handler_->HandleSolution(*this, message, primal, dual, obj_value);
  }

  // Reports an error.
  fmt::TempFormatter<PrintError> ReportError(const char *format) {
    has_errors_ = true;
    return fmt::TempFormatter<PrintError>(format);
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
template <typename OptionHandler = DummyOptionHandler>
class Solver : public SolverBase {
 private:
  class Option {
   private:
    std::string description_;

   public:
    Option(const char *description)
    : description_(FormatDescription(description)) {}

    virtual ~Option() {}

    const char *description() const { return description_.c_str(); }

    virtual char *Handle(
        OptionHandler &h, Option_Info *oi, keyword *kw, char *value) = 0;
  };

  template <typename Func, typename Value>
  class ConcreteOption : public Option {
   private:
    Func func_;

   public:
    ConcreteOption(const char *description, Func func)
    : Option(description), func_(func) {}

    char *Handle(OptionHandler &h, Option_Info *oi, keyword *kw, char *s) {
      (h.*func_)(kw->name, internal::OptionParser<Value>()(oi, kw, s));
      return s;
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

    char *Handle(OptionHandler &h, Option_Info *oi, keyword *kw, char *s) {
      (h.*func_)(kw->name, internal::OptionParser<Value>()(oi, kw, s), info_);
      return s;
    }
  };

  OptionHandler *handler_;
  std::vector<Option*> options_;

  static char *HandleOption(Option_Info *oi, keyword *kw, char *value) {
    Solver *self = GetSolver<Solver>(oi);
    Option *opt = self->options_[reinterpret_cast<size_t>(kw->info)];
    return opt->Handle(*self->handler_, oi, kw, value);
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
  : SolverBase(name, long_name, date), handler_(0) {}
  ~Solver() {
    std::for_each(options_.begin(), options_.end(), Deleter());
  }

  // Parses options and returns true if there were no errors, false otherwise.
  // Note that handler functions can report errors with SolverBase::ReportError
  // and ParseOptions will take them into account as well, returning false if
  // there was at least one such error.
  bool ParseOptions(char **argv, OptionHandler &h = dummy_option_handler) {
    handler_ = &h;
    return SolverBase::ParseOptions(argv);
  }
};
}

#endif  // SOLVERS_UTIL_SOLVER_H_

