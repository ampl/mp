/*
 Stochastic programming support

 Copyright (C) 2016 AMPL Optimization Inc

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

#ifndef MP_SP_H_
#define MP_SP_H_

#include "mp/problem-builder.h"

#include <iterator>
#include <vector>

namespace mp {

class SparseMatrix {
 private:
  std::vector<int> starts_;
  std::vector<int> indices_;
  std::vector<double> coefs_;

 public:
  SparseMatrix() {}

  void resize_major(int major_size) { starts_.resize(major_size + 1); }

  void resize_elements(int num_elements) {
    indices_.resize(num_elements);
    coefs_.resize(num_elements);
  }

  int start(int major_index) const { return starts_[major_index]; }
  int &start(int major_index) { return starts_[major_index]; }

  int index(int element_index) const { return indices_[element_index]; }
  int &index(int element_index) { return indices_[element_index]; }

  double coef(int element_index) const { return coefs_[element_index]; }
  double &coef(int element_index) { return coefs_[element_index]; }
};

// Adapts ColProblem to stochastic programming problem API.
class SPAdapter {
 private:
  const ColProblem &problem_;
  Function random_;
  SparseMatrix linear_random_;

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
      return probabilities_.empty() ?
            0 : static_cast<int>(data_.size() / probabilities_.size());
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

  struct RandomVar {
    int var_index;      // Index of the variable in the original problem.
    int rv_index;       // Index of a random vector in rvs_.
    int element_index;  // Index of an element in the random vector.

    RandomVar(int var_index, int rv_index, int element_index)
      : var_index(var_index), rv_index(rv_index),
        element_index(element_index) {}
  };
  std::vector<RandomVar> random_vars_;

  // var_core2orig_[i] is the index of core variable i in the original problem.
  std::vector<int> var_core2orig_;

  // If var_orig2core_[i] >= 0 then it gives the index of variable i in the
  // core problem. Otherwise, variable i represents a random variable/parameter.
  std::vector<int> var_orig2core_;

  int num_stages_;
  std::vector<int> num_stage_vars_;
  std::vector<int> num_stage_cons_;

  // con_core2orig_[i] is the index of core constraint i in the original
  // problem.
  std::vector<int> con_core2orig_;

  // con_orig2core_[i] is the index of constraint i in the core problem.
  std::vector<int> con_orig2core_;

  LinearExpr linear_obj_;

  // Nonlinear parts of objective expressions.
  // The array can be empty if the problem is linear.
  std::vector<NumericExpr> nonlinear_objs_;

  // Extract realizations of a random variable from a call to random(...).
  void GetRealizations(int con_index, CallExpr random, int &arg_index);

  // Extract random vectors from logical constraints with expressions of the
  // form random(x, a_1, ..., a_n).
  void GetRandomVectors(const Problem &p);

  // Get the information about variable stages.
  template <typename Suffix>
  void GetVarStages(Suffix stage);

  void ProcessObjs();
  int ProcessCons();

  void ExtractRandomTerms();

 public:
  SPAdapter(const ColProblem &p);

  // Returns the number of variables in the core problem.
  int num_vars() const { return static_cast<int>(var_core2orig_.size()); }

  // Returns the number of constraints in the core problem.
  int num_cons() const { return problem_.num_algebraic_cons(); }

  // Returns the core variable with the specified index.
  Problem::Variable var(int index) const {
    return problem_.var(var_core2orig_[index]);
  }

  // Returns the core constraint with the specified index.
  Problem::AlgebraicCon con(int index) const {
    return problem_.algebraic_con(con_core2orig_[index]);
  }

  class Objective {
   private:
    const SPAdapter *adapter_;
    int index_;

    friend class SPAdapter;

    Objective(const SPAdapter *adapter, int index)
      : adapter_(adapter), index_(index) {}

   public:
    const LinearExpr &linear_expr() const {
      assert(index_ == 0);
      return adapter_->linear_obj_;
    }

    NumericExpr nonlinear_expr() const {
      return adapter_->nonlinear_objs_[index_];
    }
  };

  // Returns the core objective with the specified index.
  Objective obj(int index) const { return Objective(this, index); }

  // Returns the number of stages.
  int num_stages() const { return num_stages_; }

  class Stage {
   private:
    const SPAdapter *adapter_;
    int index_;

    friend class SPAdapter;

    Stage(const SPAdapter *adapter, int index)
      : adapter_(adapter), index_(index) {}

   public:
    // Returns the number of variables in the stage.
    int num_vars() const { return adapter_->num_stage_vars_[index_]; }

    // Returns the number of constraints in the stage.
    int num_cons() const { return adapter_->num_stage_cons_[index_]; }
  };

  // Returns the stage with the specified index.
  Stage stage(int index) const { return Stage(this, index); }

  class Column {
   private:
    const SPAdapter *adapter_;
    int var_index_;

    friend class SPAdapter;

    Column(const SPAdapter *adapter, int var_index)
      : adapter_(adapter), var_index_(var_index) {}

   public:
    class Term {
     private:
      int con_index_;
      double coef_;

     public:
      Term(int con_index, double coef) : con_index_(con_index), coef_(coef) {}

      int con_index() const { return con_index_; }
      double coef() const { return coef_; }
    };

    class Iterator : public std::iterator<std::forward_iterator_tag, Term> {
     private:
      const SPAdapter *adapter_;
      int coef_index_;

      friend class Column;

      Iterator(const SPAdapter *adapter, int coef_index)
        : adapter_(adapter), coef_index_(coef_index) {}

      class Proxy {
       private:
        Term term_;

        friend class Iterator;

        Proxy(const Term &term) : term_(term) {}

       public:
        const Term *operator->() const { return &term_; }
      };

     public:
      Term operator*() const {
        int orig_con_index = adapter_->problem_.row_index(coef_index_);
        return Term(adapter_->con_orig2core_[orig_con_index],
                    adapter_->problem_.value(coef_index_));
      }

      Proxy operator->() const { return Proxy(**this); }

      Iterator &operator++() {
        ++coef_index_;
        return *this;
      }

      Iterator operator++(int ) {
        Iterator it(*this);
        ++coef_index_;
        return it;
      }

      bool operator==(Iterator other) const {
        assert(adapter_ == other.adapter_);
        return coef_index_ == other.coef_index_;
      }

      bool operator!=(Iterator other) const {
        assert(adapter_ == other.adapter_);
        return coef_index_ != other.coef_index_;
      }
    };

    Iterator begin() const {
      return Iterator(adapter_, adapter_->problem_.col_start(var_index_));
    }
    Iterator end() const {
      return Iterator(adapter_, adapter_->problem_.col_start(var_index_ + 1));
    }
  };

  // Returns a constraint matrix column.
  Column column(int var_index) const {
    return Column(this, var_core2orig_[var_index]);
  }

  // Returns the number of random vectors.
  int num_rvs() const { return rvs_.size(); }

  // Returns the random vector with the specified index.
  const RandomVector &rv(int index) const { return rvs_[index]; }

  class Scenario {
   private:
    std::vector<double> rhs_offsets_;

    friend class SPAdapter;

   public:
    const std::vector<double> &rhs_offsets() const { return rhs_offsets_; }
  };

  void GetScenario(Scenario &s, int scenario_index) const;

  // If var_index refers to a random variable, returns its index in
  // random_vars_. Otherwise returns -1.
  int GetRandVarIndex(int var_index) const {
    int core_var_index = var_orig2core_[var_index];
    return core_var_index < 0 ? -(core_var_index + 1) : -1;
  }

  double GetRealization(int rv_index, int scenario) const {
    const auto &info = random_vars_[rv_index];
    return rvs_[info.rv_index].value(info.element_index, scenario);
  }
};
}  // namespace mp

#endif  // MP_SP_H_
