/*
 IBM/ILOG CP solver for AMPL.

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

#ifndef SOLVERS_ILOGCP_ILOGCP_H_
#define SOLVERS_ILOGCP_ILOGCP_H_

#if __clang__
# pragma clang diagnostic push
# pragma clang diagnostic ignored "-Wunused-parameter"
# pragma clang diagnostic ignored "-Wunused-private-field"
#elif _MSC_VER
# pragma warning(push)
# pragma warning(disable: 4244)
#endif

#include <ilcplex/ilocplex.h>
#include <ilcp/cp.h>

#if __clang__
# pragma clang diagnostic pop
#elif _MSC_VER
# pragma warning(pop)
#endif

#include <string.h> /* This and -fpermissive seem to be needed for MacOSX, */
                    /* at least with g++ 4.6.  Otherwise there are errors */
                    /* with iloconcert/iloenv.h . */
#include <limits.h> /* Needed for g++ -m32 on MacOSX. */
#include <memory>
#include <string>
#include <vector>

#include "solvers/util/solver.h"
#include "solvers/util/format.h"

namespace ampl {

class Optimizer : public Interruptible {
 private:
  IloRangeArray cons_;

 public:
  Optimizer(IloEnv env) : cons_(env) {}
  virtual ~Optimizer();

  void AllocateCons(int num_cons) { cons_.setSize(num_cons); }
  IloRangeArray cons() const { return cons_; }

  virtual IloAlgorithm algorithm() const = 0;

  virtual void StartSearch() = 0;
  virtual void EndSearch() = 0;
  virtual bool FindNextSolution() = 0;

  virtual void GetSolutionInfo(
      fmt::Writer &w, std::vector<double> &dual_values) const = 0;
};

class CPLEXOptimizer : public Optimizer {
 private:
  IloCplex cplex_;
  IloCplex::Aborter aborter_;
  bool started_;

 public:
  CPLEXOptimizer(IloEnv env);

  IloCplex cplex() const { return cplex_; }
  IloAlgorithm algorithm() const { return cplex_; }

  void StartSearch() { started_ = true; }
  void EndSearch() { started_ = false; }

  bool FindNextSolution();

  void GetSolutionInfo(fmt::Writer &w, std::vector<double> &dual_values) const;

  void Interrupt() { aborter_.abort(); }
};

class CPOptimizer : public Optimizer {
 private:
  IloCP cp_;

 public:
  CPOptimizer(IloEnv env, const Problem *p);

  IloCP solver() const { return cp_; }
  IloAlgorithm algorithm() const { return cp_; }

  void StartSearch() { cp_.startNewSearch(); }
  void EndSearch() { cp_.endSearch(); }
  bool FindNextSolution() { return cp_.next() != IloFalse; }

  void GetSolutionInfo(fmt::Writer &w, std::vector<double> &dual_values) const;

  void Interrupt() { cp_.abortSearch(); }
};

// IlogCP solver.
class IlogCPSolver : public Solver<IlogCPSolver> {
 private:
  IloEnv env_;
  std::auto_ptr<Optimizer> optimizer_;
  bool gotopttype_;

  // Do not implement.
  IlogCPSolver(const IlogCPSolver&);
  IlogCPSolver &operator=(const IlogCPSolver&);

 public:
  // Options accessible from AMPL.
  enum Option {
    DEBUGEXPR,
    OPTIMIZER,
    USENUMBEROF,
    NUM_OPTIONS
  };

  // Values for the OPTIMIZER option.
  enum {
    AUTO  = -1,
    CP    =  0,
    CPLEX =  1
  };

 private:
  int options_[NUM_OPTIONS];

  CPOptimizer *GetCPForOption(fmt::StringRef option_name) const;
  CPLEXOptimizer *GetCPLEXForOption(fmt::StringRef option_name) const;

  std::string GetOptimizer(const char *name) const;
  void SetOptimizer(const char *name, const char *value);

  int GetBoolOption(const char *, Option opt) const { return options_[opt]; }
  void SetBoolOption(const char *name, int value, Option opt);

  // An enumerated option of the constraint programming optimizer.
  class EnumCPOption : public TypedSolverOption<std::string> {
   private:
    IlogCPSolver &solver_;
    IloCP::IntParam param_;
    int start_;           // start value for the enumerated options
    const char **values_; // string values for enum options
    bool accepts_auto_;   // true if the option accepts IloCP::Auto value

   public:
    EnumCPOption(const char *name, const char *description,
        IlogCPSolver *s, IloCP::IntParam p, int start,
        const char **values, bool accepts_auto = false)
    : TypedSolverOption<std::string>(name, description),
      solver_(*s), param_(p), start_(start), values_(values),
      accepts_auto_(accepts_auto) {
    }

    std::string GetValue() const;
    void SetValue(const char *value);
  };

  // An integer option of the constraint programming optimizer.
  class IntCPOption : public TypedSolverOption<int> {
   private:
    IlogCPSolver &solver_;
    IloCP::IntParam param_;

   public:
    IntCPOption(const char *name, const char *description,
        IlogCPSolver *s, IloCP::IntParam p)
    : TypedSolverOption<int>(name, description), solver_(*s), param_(p) {}

    int GetValue() const;
    void SetValue(int value);
  };

  // Returns a double option of the constraint programming optimizer.
  double GetCPDblOption(const char *name, IloCP::NumParam param) const;

  // Sets a double option of the constraint programming optimizer.
  void SetCPDblOption(const char *name, double value, IloCP::NumParam param);

  // Returns an integer option of the CPLEX optimizer.
  int GetCPLEXIntOption(const char *name, int param) const;

  // Sets an integer option of the CPLEX optimizer.
  void SetCPLEXIntOption(const char *name, int value, int param);

  void CreateOptimizer(const Problem *p);

 protected:

  std::string GetOptionHeader();

 public:
  IlogCPSolver();
  virtual ~IlogCPSolver();

  IloEnv env() const { return env_; }

  IloAlgorithm alg() const {
    return optimizer_.get() ? optimizer_->algorithm() : IloAlgorithm();
  }

  Optimizer *optimizer() const { return optimizer_.get(); }

  bool ParseOptions(char **argv, unsigned flags = 0, const Problem *p = 0);

  int GetOption(Option opt) const {
    assert(opt >= 0 && opt < NUM_OPTIONS);
    return options_[opt];
  }

  void use_numberof(bool use = true) { options_[USENUMBEROF] = use; }

  void Solve(Problem &p);
};
}

#endif  // SOLVERS_ILOGCP_ILOGCP_H_
