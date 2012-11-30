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

#include <cstring>

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

// An AMPL solver driver.
class Driver {
 private:
  Problem problem_;

 public:
  Problem &problem() { return problem_; }

  // Gets the options.
  int GetOptions(char **argv, Option_Info *oi);

  // Writes the solution.
  void WriteSolution(char *msg, double *x, double *y, Option_Info* oi) {
    write_sol_ASL(reinterpret_cast<ASL*>(problem_.asl_), msg, x, y, oi);
  }
};

template <typename Handler>
class OptionInfo : public Option_Info {
 public:
  struct Option {
    char *(Handler::*handler)(Option_Info *oi, keyword *kw, char *value);
    const void *info;
  };

 private:
  Handler &handler_;
  std::vector<Option> options_;
  std::vector<keyword> keywords_;

  struct KeywordNameLess {
    bool operator()(const keyword &lhs, const keyword &rhs) const {
      return std::strcmp(lhs.name, rhs.name) < 0;
    }
  };

  static char *HandleOption(Option_Info *oi, keyword *kw, char *value) {
    OptionInfo *self = static_cast<OptionInfo*>(oi);
    Option &opt = self->options_[reinterpret_cast<size_t>(kw->info)];
    keyword thiskw(*kw);
    thiskw.info = const_cast<void*>(opt.info);
    return (self->handler_.*opt.handler)(oi, &thiskw, value);
  }

  void AddKeyword(const char *name,
      const char *description, Kwfunc func, const void *info) {
    keywords_.push_back(keyword());
    keyword &kw = keywords_.back();
    kw.name = const_cast<char*>(name);
    kw.desc = const_cast<char*>(description);
    kw.kf = func;
    kw.info = const_cast<void*>(info);

    // TODO: set once
    std::sort(keywords_.begin(), keywords_.end(), KeywordNameLess());
    keywds = &keywords_[0];
    n_keywds = keywords_.size();
  }

 public:
  OptionInfo(Handler &h): Option_Info(), handler_(h) {
    // TODO: align text
    AddKeyword("version",
        "Single-word phrase:  report version details\n"
        "before solving the problem.\n", Ver_val, 0);
    AddKeyword("wantsol",
        "In a stand-alone invocation (no -AMPL on the\n"
        "command line), what solution information to\n"
        "write.  Sum of\n"
        "      1 = write .sol file\n"
        "      2 = primal variables to stdout\n"
        "      4 = dual variables to stdout\n"
        "      8 = suppress solution message\n", WS_val, 0);
  }

  void AddOption(const char *name, const char *description,
      char *(Handler::*handler)(Option_Info *oi, keyword *kw, char *value),
      const void* info) {
    AddKeyword(name, description, HandleOption,
        reinterpret_cast<void*>(options_.size()));
    Option opt = {handler, info};
    options_.push_back(opt);
  }
};
}

#endif  // SOLVERS_UTIL_DRIVER_H_

