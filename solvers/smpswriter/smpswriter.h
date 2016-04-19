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

#include "mp/problem.h"
#include "mp/solver.h"

namespace mp {

class FileWriter;
class SMPSNameReader;

class SMPSWriter : public SolverImpl<ColProblem> {
 private:
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

   private:
    std::vector<ConTerm> con_terms_;
    std::vector<RHS> rhs_;

   public:
    Scenario() {}

    void AddConTerm(int con_index, int var_index, double coef) {
      con_terms_.push_back(ConTerm(con_index, var_index, coef));
    }

    void AddRHS(int con_index, double rhs) {
      rhs_.push_back(RHS(con_index, rhs));
    }

    typedef std::vector<ConTerm>::const_iterator ConTermIterator;

    ConTermIterator con_term_begin() const { return con_terms_.begin(); }
    ConTermIterator con_term_end() const { return con_terms_.end(); }

    typedef std::vector<RHS>::const_iterator RHSIterator;

    RHSIterator rhs_begin() const { return rhs_.begin(); }
    RHSIterator rhs_end() const { return rhs_.end(); }
  };

  // core_var_indices_[i] is the index of variable i in the core problem.
  std::vector<int> core_var_indices_;

  int num_stage1_cons_;

  // con_indices_[i] is the index of core constraint i in the original problem.
  std::vector<int> con_indices_;

  // core_con_indices_[i] is the index of constraint i in the core problem.
  std::vector<int> core_con_indices_;

  std::vector<Scenario> scenarios_;

  void GetScenario(ColProblem &p, int scenario, std::vector<double> &coefs);

  void WriteColumns(FileWriter &writer, const ColProblem &p,
                    int num_stages, int num_core_cons,
                    const std::vector<double> &core_obj_coefs,
                    const std::vector<double> &coefs);

 public:
  SMPSWriter();

  typedef ColProblemBuilder NLProblemBuilder;

  void Solve(ColProblem &p, SolutionHandler &sh);
};
}

#endif  // MP_SOLVERS_SMPSWRITER_H_
