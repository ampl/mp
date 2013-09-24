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
#include <string>
#include <vector>

#include "solvers/util/solver.h"
#include "solvers/util/format.h"

namespace ampl {

// IlogCP solver.
class IlogCPSolver : private Interruptible, public Solver<IlogCPSolver> {
 private:
  IloEnv env_;
  IloCP cp_;

  // Do not implement.
  IlogCPSolver(const IlogCPSolver&);
  IlogCPSolver &operator=(const IlogCPSolver&);

  void Interrupt() { cp_.abortSearch(); }

 public:
  // Options accessible from AMPL.
  enum Option {
    DEBUGEXPR,
    USENUMBEROF,
    NUM_OPTIONS
  };

 private:
  int options_[NUM_OPTIONS];

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

  void GetSolution(IloNumVarArray vars,
      std::vector<double> &solution, std::vector<double> &dual_solution);

 protected:

  std::string GetOptionHeader();

 public:
  IlogCPSolver();
  virtual ~IlogCPSolver();

  IloEnv env() const { return env_; }
  IloCP optimizer() const { return cp_; }

  int GetOption(Option opt) const {
    assert(opt >= 0 && opt < NUM_OPTIONS);
    return options_[opt];
  }

  void use_numberof(bool use = true) { options_[USENUMBEROF] = use; }

  void Solve(Problem &p);
};
}

#endif  // SOLVERS_ILOGCP_ILOGCP_H_
