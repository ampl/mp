/*
 SMPS writer implemented as an AMPL solver.

 Copyright (C) 2013 AMPL Optimization Inc

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

#ifndef MP_SOLVERS_SMPSWRITER_H_
#define MP_SOLVERS_SMPSWRITER_H_

#include <vector>

#include "asl/aslsolver.h"

namespace mp {

class FileWriter;

class SMPSWriter : public ASLSolver {
 private:
  // Information about a variable or constraint.
  struct VarConInfo {
    int core_index;  // index of this variable in the core problem
    int scenario_index;
    VarConInfo() : core_index(), scenario_index() {}
  };

  struct CoreConInfo {
    char type;
    double rhs;
    CoreConInfo() : type(), rhs() {}
  };

  struct CoreVarInfo {
    double lb;
    double ub;
    CoreVarInfo() : lb(), ub() {}
  };

  class Scenario {
   public:
    // A constraint expression term.
    struct ConTerm {
      int con_index;
      int var_index;
      double coef;
      ConTerm(int con_index, int var_index, double coef)
      : con_index(con_index), var_index(var_index), coef(coef) {}
    };

    struct RHS {
      int con_index;
      double rhs;
      RHS(int con_index, double rhs) : con_index(con_index), rhs(rhs) {}
    };

    struct Bound {
      int var_index;
      double bound;
      Bound(int var_index, double bound) : var_index(var_index), bound(bound) {}
    };

   private:
    std::vector<ConTerm> con_terms_;
    std::vector<RHS> rhs_;
    std::vector<Bound> lb_;
    std::vector<Bound> ub_;

   public:
    Scenario() {}

    void AddConTerm(int con_index, int var_index, double coef) {
      con_terms_.push_back(ConTerm(con_index, var_index, coef));
    }

    void AddRHS(int con_index, double rhs) {
      rhs_.push_back(RHS(con_index, rhs));
    }

    void AddLB(int var_index, double lb) {
      lb_.push_back(Bound(var_index, lb));
    }

    void AddUB(int var_index, double ub) {
      ub_.push_back(Bound(var_index, ub));
    }

    typedef std::vector<ConTerm>::const_iterator ConTermIterator;

    ConTermIterator con_term_begin() const { return con_terms_.begin(); }
    ConTermIterator con_term_end() const { return con_terms_.end(); }

    typedef std::vector<RHS>::const_iterator RHSIterator;

    RHSIterator rhs_begin() const { return rhs_.begin(); }
    RHSIterator rhs_end() const { return rhs_.end(); }

    typedef std::vector<Bound>::const_iterator BoundIterator;

    BoundIterator lb_begin() const { return lb_.begin(); }
    BoundIterator lb_end() const { return lb_.end(); }

    BoundIterator ub_begin() const { return ub_.begin(); }
    BoundIterator ub_end() const { return ub_.end(); }
  };

  std::vector<VarConInfo> var_info;
  std::vector<VarConInfo> con_info;
  std::vector<Scenario> scenarios;

  void SplitConRHSIntoScenarios(
      const Problem &p, std::vector<CoreConInfo> &core_cons);

  void SplitVarBoundsIntoScenarios(
      const Problem &p, std::vector<CoreVarInfo> &core_vars);

  void WriteColumns(FileWriter &writer, const Problem &p, int num_stages,
      int num_core_cons, const std::vector<double> &core_obj_coefs);

 protected:
  int DoSolve(Problem &p, SolutionHandler &sh);

 public:
  SMPSWriter();
};
}

#endif  // MP_SOLVERS_SMPSWRITER_H_
