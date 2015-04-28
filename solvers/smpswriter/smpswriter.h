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

// An optimization problem with a column-wise constraint matrix.
class ColProblem : public Problem {
 private:
  // Column-wise constraint matrix.
  std::vector<int> col_starts_;
  std::vector<int> row_indices_;
  std::vector<double> coefs_;

  std::string name_;

  friend class ColProblemBuilder;

 public:
  ColProblem(const Solver &, fmt::StringRef name) { name_ = name; }

  int col_start(int col_index) const { return col_starts_[col_index]; }
  int row_index(int elt_index) const { return row_indices_[elt_index]; }
  double value(int elt_index) const { return coefs_[elt_index]; }

  const std::string &name() const { return name_; }
};

// An NL handler that builds a problem with a column-wise constraint matrix.
class ColProblemBuilder : public ProblemBuilderToNLAdapter<ColProblem> {
 private:
  typedef ProblemBuilderToNLAdapter<ColProblem> Base;

 public:
  explicit ColProblemBuilder(ColProblem &p) : ProblemBuilderToNLAdapter(p) {}

  void OnHeader(const NLHeader &h) {
    Base::OnHeader(h);
    builder_.col_starts_.reserve(h.num_vars + 1);
    builder_.col_starts_.resize(2);
    builder_.row_indices_.resize(h.num_con_nonzeros);
    builder_.coefs_.resize(h.num_con_nonzeros);
  }

  class ColumnSizeHandler {
   private:
    ColProblem *problem_;

    explicit ColumnSizeHandler(ColProblem &p) : problem_(&p) {}

    friend class ColProblemBuilder;

   public:
    void Add(int size) {
      // Convert column size to column offset.
      problem_->col_starts_.push_back(problem_->col_starts_.back() + size);
    }
  };

  ColumnSizeHandler OnColumnSizes() {
    return ColumnSizeHandler(builder_);
  }

  class LinearConHandler {
   private:
    ColProblem *problem_;
    int con_index_;

    friend class ColProblemBuilder;

    LinearConHandler(ColProblem &p, int con_index)
      : problem_(&p), con_index_(con_index) {}

   public:
    void AddTerm(int var_index, double coef) {
      int index = problem_->col_starts_[var_index + 1]++;
      problem_->row_indices_[index] = con_index_;
      problem_->coefs_[index] = coef;
    }
  };

  LinearConHandler OnLinearConExpr(int con_index, int) {
    // Pass zero as the number of linear terms as we store constraints
    // column-wise rather than row-wise.
    Base::OnLinearConExpr(con_index, 0);
    return LinearConHandler(builder_, con_index);
  }
};

class SMPSWriter : public SolverImpl<ColProblem> {
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
      const Problem &p, std::vector<CoreConInfo> &core_cons,
      SMPSNameReader &con_names);

  void SplitVarBoundsIntoScenarios(
      const Problem &p, std::vector<CoreVarInfo> &core_vars);

  void WriteColumns(FileWriter &writer, const ColProblem &p, int num_stages,
      int num_core_cons, const std::vector<double> &core_obj_coefs);

 public:
  SMPSWriter();

  typedef ColProblemBuilder NLProblemBuilder;

  void Solve(ColProblem &p, SolutionHandler &sh);
};
}

#endif  // MP_SOLVERS_SMPSWRITER_H_
