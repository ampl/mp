/*
 SMPS writer implemented as an AMPL solver.

 Copyright (C) 2013 - 2016 AMPL Optimization Inc

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

class SMPSWriter : public SolverImpl<ColProblem> {
 private:
  class RandomVector {
   private:
    std::vector<double> probabilities_;
    // A matrix with rows representing elements and columns representing
    // realizatons.
    std::vector<double> data_;

   public:
    void AddProbability(double prob) {
      probabilities_.push_back(prob);
    }

    void Add(double value) {
      data_.push_back(value);
    }

    int num_realizations() const {
      return static_cast<int>(probabilities_.size());
    }

    int num_elements() const {
      return static_cast<int>(data_.size() / probabilities_.size());
    }

    void set_num_realizations(int num_realizations) {
      probabilities_.resize(num_realizations, 1.0 / num_realizations);
    }

    double probability(int realization) const {
      return probabilities_[realization];
    }

    double value(int element, int realization) const {
      return data_[element * num_realizations() + realization];
    }
  };

  std::vector<RandomVector> rvs_;

  struct RVInfo {
    int var_index;      // Index of the variable in the original problem.
    int rv_index;       // Index of a random vector in rvs_.
    int element_index;  // Index of an element in the random vector.

    RVInfo(int var_index, int rv_index, int element_index)
      : var_index(var_index), rv_index(rv_index),
        element_index(element_index) {}
  };
  std::vector<RVInfo> rv_info_;

  // var_core2orig_[i] is the index of core variable i in the original problem.
  std::vector<int> var_core2orig_;

  // If var_orig2core_[i] >= 0 then it gives the index of variable i in the
  // core problem. Otherwise, variable i represents a random variable/parameter.
  std::vector<int> var_orig2core_;

  int num_stage1_cons_;

  // con_core2orig_[i] is the index of core constraint i in the original
  // problem.
  std::vector<int> con_core2orig_;

  // con_orig2core_[i] is the index of constraint i in the core problem.
  std::vector<int> con_orig2core_;

  mp::Function random_;

  // Add an element of a random vector.
  void AddRVElement(Expr arg, int rv_index, int element_index);

  // Extract random vectors from logical constraints with expressions of the
  // form random(x, a_1, ..., a_n).
  void GetRandomVectors(const Problem &p);

  void GetScenario(ColProblem &p, int scenario, std::vector<double> &coefs,
                   std::vector<double> &rhs);

  void WriteColumns(FileWriter &writer, const ColProblem &p, int num_core_cons,
                    const std::vector<double> &core_obj_coefs,
                    const std::vector<double> &coefs);

 public:
  SMPSWriter();

  typedef ColProblemBuilder NLProblemBuilder;

  void Solve(ColProblem &p, SolutionHandler &sh);

  int GetRVIndex(int var_index) const {
    int core_var_index = var_orig2core_[var_index];
    return core_var_index < 0 ? -(core_var_index + 1) : -1;
  }

  double GetRealization(int rv_index, int scenario) const {
    const auto &info = rv_info_[rv_index];
    return rvs_[info.rv_index].value(info.element_index, scenario);
  }
};
}

#endif  // MP_SOLVERS_SMPSWRITER_H_
