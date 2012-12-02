/*
 Utilities for writing AMPL solver drivers.

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

#ifndef SOLVERS_UTIL_DRIVER_H_
#define SOLVERS_UTIL_DRIVER_H_

#include <memory>

#include "solvers/getstub.h"
#include "solvers/util/expr.h"

namespace ampl {

// An AMPL problem.
class Problem {
 private:
  ASL_fg *asl_;

  // Do not implement.
  Problem(const Problem&);
  Problem& operator=(const Problem&);

  friend class Driver;

 public:
  Problem();
  virtual ~Problem();

  // Reads the problem form a .nl file.
  bool Read(char **&argv, Option_Info *oi);

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

class BaseOptionInfo : public Option_Info {
 private:
  std::vector<keyword> keywords_;
  bool sorted_;

  void Sort();

  friend class Driver;

 protected:
  BaseOptionInfo();

  void AddKeyword(const char *name,
      const char *description, Kwfunc func, const void *info);
};

template <typename Handler>
class OptionInfo : public BaseOptionInfo {
 private:
  class Option {
   public:
    virtual ~Option() {}

    virtual char *Handle(
        Handler &h, Option_Info *oi, keyword *kw, char *value) = 0;
  };

  template <typename Func, typename Value>
  class ConcreteOption : public Option {
   private:
    Func func_;

   public:
    ConcreteOption(Func func) : func_(func) {}

    char *Handle(Handler &h, Option_Info *oi, keyword *kw, char *s) {
      (h.*func_)(kw->name, OptionParser<Value>()(oi, kw, s));
      return s;
    }
  };

  template <typename Func, typename Info, typename Value>
  class ConcreteOptionWithInfo : public Option {
   private:
    Func func_;
    Info info_;

   public:
    ConcreteOptionWithInfo(Func func, const Info &info)
    : func_(func), info_(info) {}

    char *Handle(Handler &h, Option_Info *oi, keyword *kw, char *s) {
      (h.*func_)(kw->name, OptionParser<Value>()(oi, kw, s), info_);
      return s;
    }
  };

  Handler &handler_;
  std::vector<Option*> options_;

  static char *HandleOption(Option_Info *oi, keyword *kw, char *value) {
    OptionInfo *self = static_cast<OptionInfo*>(oi);
    Option *opt = self->options_[reinterpret_cast<size_t>(kw->info)];
    return opt->Handle(self->handler_, oi, kw, value);
  }

  void AddOption(const char *name,
      const char *description, std::auto_ptr<Option> opt) {
    AddKeyword(name, description, HandleOption,
        reinterpret_cast<void*>(options_.size()));
    options_.push_back(0);
    options_.back() = opt.release();
  }

  struct Deleter {
    void operator()(Option *opt) { delete opt; }
  };

 public:
  OptionInfo(Handler &h): handler_(h) {}

  ~OptionInfo() {
    std::for_each(options_.begin(), options_.end(), Deleter());
  }

  template <typename Func>
  void AddIntOption(const char *name, const char *description, Func f) {
    AddOption(name, description, std::auto_ptr<Option>(
        new ConcreteOption<Func, int>(f)));
  }

  template <typename Info, typename Func>
  void AddIntOption(const char *name,
      const char *description, Func f, const Info &info) {
    AddOption(name, description, std::auto_ptr<Option>(
        new ConcreteOptionWithInfo<Func, Info, int>(f, info)));
  }

  template <typename Func>
  void AddDblOption(const char *name, const char *description, Func f) {
    AddOption(name, description, std::auto_ptr<Option>(
        new ConcreteOption<Func, double>(f)));
  }

  template <typename Info, typename Func>
  void AddDblOption(const char *name,
      const char *description, Func f, const Info &info) {
    AddOption(name, description, std::auto_ptr<Option>(
        new ConcreteOptionWithInfo<Func, Info, double>(f, info)));
  }

  template <typename Func>
  void AddStrOption(const char *name, const char *description, Func f) {
    AddOption(name, description, std::auto_ptr<Option>(
        new ConcreteOption<Func, const char*>(f)));
  }

  template <typename Info, typename Func>
  void AddStrOption(const char *name,
      const char *description, Func f, const Info &info) {
    AddOption(name, description, std::auto_ptr<Option>(
        new ConcreteOptionWithInfo<Func, Info, const char*>(f, info)));
  }
};

#undef printf

// An AMPL solver driver.
class Driver {
 private:
  Problem problem_;
  bool has_errors_;

 public:
  Driver() : has_errors_(false) {}

  Problem &problem() { return problem_; }

  bool has_errors() const { return has_errors_; }

  // Reports an error.
  void ReportError(const char *format, ...)
    __attribute__((format(printf, 2, 3)));

  // Gets the options.
  bool GetOptions(char **argv, BaseOptionInfo &oi);

  // Writes the solution.
  void WriteSolution(char *msg, double *x, double *y, Option_Info* oi) {
    write_sol_ASL(reinterpret_cast<ASL*>(problem_.asl_), msg, x, y, oi);
  }
};
}

#endif  // SOLVERS_UTIL_DRIVER_H_

