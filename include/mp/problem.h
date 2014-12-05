/*
 Optimization problem

 Copyright (C) 2014 AMPL Optimization Inc

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

#ifndef MP_PROBLEM_H_
#define MP_PROBLEM_H_

#include <cstddef>  // for std::size_t
#include <limits>
#include <vector>

#include "mp/expr.h"
#include "mp/suffix.h"

// Maximum index of a variable, objective or constraint.
#ifndef MP_MAX_PROBLEM_ITEMS
# define MP_MAX_PROBLEM_ITEMS \
  static_cast<std::size_t>(std::numeric_limits<int>::max())
#endif

namespace mp {

class LinearExpr {
 private:
  class Term {
   private:
    int var_index_;
    double coef_;

    friend class LinearExpr;

    Term(int var_index, double coef) : var_index_(var_index), coef_(coef) {}

   public:
    int var_index() const { return var_index_; }
    double coef() const { return coef_; }
  };
  std::vector<Term> terms_;

 public:
  int num_terms() const { return static_cast<int>(terms_.size()); }

  typedef std::vector<Term>::const_iterator iterator;

  iterator begin() const { return terms_.begin(); }
  iterator end() const { return terms_.end(); }

  void AddTerm(int var_index, double coef) {
    terms_.push_back(Term(var_index, coef));
  }

  void Reserve(int num_terms) {
    terms_.reserve(num_terms);
  }
};

// An optimization problem.
template <typename Alloc>
class BasicProblem : public ExprFactory, public SuffixManager {
 private:
  // A variable.
  struct Var {
    double lb;
    double ub;
    Var(double lb, double ub) : lb(lb), ub(ub) {}
  };
  std::vector<Var> vars_;

  // Packed variable type information.
  // is_var_int_[i] specifies whether variable i is integer.
  std::vector<bool> is_var_int_;

  class LinearExprBuilder {
   private:
    LinearExpr *expr_;

   public:
    explicit LinearExprBuilder(LinearExpr *expr) : expr_(expr) {}

    void AddTerm(int var_index, double coef) {
      expr_->AddTerm(var_index, coef);
    }
  };

  // Packed objective type information.
  // is_obj_max_[i] specifies whether objective i is maximization.
  std::vector<bool> is_obj_max_;

  // Linear parts of objective expessions.
  std::vector<LinearExpr> linear_objs_;

  // Nonlinear parts of objective expressions.
  // The array can be empty if the problem is linear.
  std::vector<NumericExpr> nonlinear_objs_;

  // Algebraic constraint information.
  struct AlgebraicConInfo {
    // Linear part of an algebraic constraint expression.
    // Nonlinear parts are stored in nonlinear_cons_ to avoid overhead
    // for linear problems.
    LinearExpr linear_expr;
    double lb;
    double ub;
    AlgebraicConInfo(double lb, double ub) : lb(lb), ub(ub) {}
  };
  std::vector<AlgebraicConInfo> algebraic_cons_;

  // Information about complementarity conditions.
  // compl_vars_[i] > 0 means constraint i complements variable
  // compl_vars_[i] - 1. The array can be empty if there are no
  // complementarity conditions.
  std::vector<unsigned> compl_vars_;

  // Nonlinear parts of algebraic constraint expressions.
  // The array can be empty if the problem is linear.
  std::vector<NumericExpr> nonlinear_cons_;

  // Logical constraint expressions.
  std::vector<LogicalExpr> logical_cons_;

  // Linear parts of common expressions.
  std::vector<LinearExpr> linear_exprs_;

  // Nonlinear parts of common expressions.
  std::vector<NumericExpr> nonlinear_exprs_;

  // Initial values for variables.
  std::vector<double> initial_values_;

  // Initial values for dual variables.
  std::vector<double> initial_dual_values_;

  std::vector<Function> funcs_;

  // Checks if index is in the range [0, size).
  void CheckIndex(int index, int size) const {
    MP_ASSERT(0 <= index && index < size, "invalid index");
  }

  // A list of problem elements.
  template <typename T>
  class List {
   private:
    const BasicProblem *problem_;

    friend class BasicProblem;

    explicit List(const BasicProblem *p) : problem_(p) {}

   public:
    class iterator : std::iterator<std::forward_iterator_tag, T> {
     private:
      T item_;

      friend class List<T>;

      iterator(const BasicProblem *p, int index) : item_(p, index) {}

     public:
      const T *operator->() const {
        MP_ASSERT(0 <= item_.index_ &&
                  item_.index_ < T::num_items(*item_.problem_),
                  "invalid access");
        return &item_;
      }

      T operator*() const { return *this->operator->(); }

      iterator &operator++() {
        ++item_.index_;
        return *this;
      }

      iterator operator++(int ) {
        iterator it(*this);
        ++item_.index_;
        return it;
      }

      bool operator==(iterator other) const {
        return item_ == other.item_;
      }
      bool operator!=(iterator other) const {
        return item_ != other.item_;
      }
    };

    // Returns an iterator to the first element in the list.
    iterator begin() const {
      return iterator(problem_, 0);
    }

    // Returns an iterator to the element following the last element
    // in the list. An attempt to access this element will result in
    // assertion failure if assertions are enabled and undefined behavoir
    // otherwise.
    iterator end() const {
      return iterator(problem_, T::num_items(*problem_));
    }
  };

  template <typename T>
  class SuffixHandler {
   private:
    Suffix *suffix_;

   public:
    explicit SuffixHandler(Suffix *s) : suffix_(s) {}

    // Sets the suffix value.
    void SetValue(int index, T value) {
      suffix_->set_value(index, value);
    }
  };

  int GetSuffixSize(int suffix_type);

  template <typename T>
  SuffixHandler<T> AddSuffix(fmt::StringRef name, int kind) {
    int type = kind & suf::MASK;
    SuffixSet::Set &set = suffixes(type).set_;
    Suffix &suffix = const_cast<Suffix&>(*set.insert(Suffix(name, kind)).first);
    suffix.InitValues(GetSuffixSize(type));
    return SuffixHandler<T>(&suffix);
  }

  struct ProblemItem {
    const BasicProblem *problem_;
    int index_;

    ProblemItem(const BasicProblem *p, int index)
      : problem_(p), index_(index) {}
  };

 public:
  // Returns the number of variables.
  int num_vars() const { return static_cast<int>(vars_.size()); }

  // Returns the number of objectives.
  int num_objs() const { return static_cast<int>(linear_objs_.size()); }

  // Returns the number of algebraic constraints.
  int num_algebraic_cons() const {
    return static_cast<int>(algebraic_cons_.size());
  }

  // Returns the number of logical constraints.
  int num_logical_cons() const {
    return static_cast<int>(logical_cons_.size());
  }

  // An optimization variable.
  class Variable : private ProblemItem {
   private:
    friend class BasicProblem;

    Variable(const BasicProblem *p, int index) : ProblemItem(p, index) {}

    static int num_items(const BasicProblem &p) {
      return p.num_vars();
    }

   public:
    // Returns the lower bound on the variable.
    double lb() const {
      return this->problem_->vars_[this->index_].lb;
    }

    // Returns the upper bound on the variable.
    double ub() const {
      return this->problem_->vars_[this->index_].ub;
    }

    // Returns the type of the variable.
    var::Type type() const {
      return this->problem_->is_var_int_[this->index_] ?
          var::INTEGER : var::CONTINUOUS;
    }

    bool operator==(Variable other) const {
      MP_ASSERT(this->problem_ == other.problem_,
                "comparing variables from different problems");
      return this->index_ == other.index_;
    }
    bool operator!=(Variable other) const {
      return !(*this == other);
    }
  };

  // A list of variables.
  typedef List<Variable> VarList;

  // Returns the list of problem variables.
  // It can be used to iterate over all variables in a problem:
  //   for (auto var: problem.vars()) {
  //     ...
  //   }
  VarList vars() const { return VarList(this); }

  // Returns the variable at the specified index.
  Variable var(int index) const {
    CheckIndex(index, num_vars());
    return Variable(this, index);
  }

  // Adds a variable.
  void AddVar(double lb, double ub, var::Type type = var::CONTINUOUS) {
    MP_ASSERT(vars_.size() < MP_MAX_PROBLEM_ITEMS, "too many variables");
    vars_.push_back(Var(lb, ub));
    is_var_int_.push_back(type != var::CONTINUOUS);
  }

  // An objective.
  class Objective : private ProblemItem {
   private:
    friend class BasicProblem;

    Objective(const BasicProblem *p, int index) : ProblemItem(p, index) {}

    static int num_items(const BasicProblem &p) {
      return p.num_objs();
    }

   public:
    // Returns the type of the objective.
    obj::Type type() const {
      return this->problem_->is_obj_max_[this->index_] ? obj::MAX : obj::MIN;
    }

    // Returns the linear part of an objective expression.
    const LinearExpr &linear_expr() const {
      return this->problem_->linear_objs_[this->index_];
    }

    // Returns the nonlinear part of an objective expression.
    NumericExpr nonlinear_expr() const {
      return this->problem_->nonlinear_objs_.empty() ?
            NumericExpr() : this->problem_->nonlinear_objs_[this->index_];
    }

    bool operator==(Objective other) const {
      MP_ASSERT(this->problem_ == other.problem_,
                "comparing objectives from different problems");
      return this->index_ == other.index_;
    }
    bool operator!=(Objective other) const {
      return !(*this == other);
    }
  };

  // A list of objectives.
  typedef List<Objective> ObjList;

  // Returns the list of problem objectives.
  // It can be used to iterate over all objectives in a problem:
  //   for (auto obj: problem.objs()) {
  //     ...
  //   }
  ObjList objs() const { return ObjList(this); }

  // Returns the objective at the specified index.
  Objective obj(int index) const {
    CheckIndex(index, num_objs());
    return Objective(this, index);
  }

  typedef LinearExprBuilder LinearObjBuilder;

  // Adds an objective.
  // Returns a handler for receiving linear terms in the objective.
  LinearObjBuilder AddObj(obj::Type type, NumericExpr expr,
                          int num_linear_terms = 0);

  LinearObjBuilder AddObj(obj::Type type, int num_linear_terms = 0) {
    return AddObj(type, NumericExpr(), num_linear_terms);
  }

  // An algebraic constraint.
  class AlgebraicCon : private ProblemItem {
   private:
    friend class BasicProblem;

    AlgebraicCon(const BasicProblem *p, int index) : ProblemItem(p, index) {}

    static int num_items(const BasicProblem &p) {
      return p.num_algebraic_cons();
    }

   public:
    // Returns the lower bound on the constraint.
    double lb() const {
      return this->problem_->algebraic_cons_[this->index_].lb;
    }

    // Returns the upper bound on the constraint.
    double ub() const {
      return this->problem_->algebraic_cons_[this->index_].ub;
    }

    // Returns the linear part of a constraint expression.
    const LinearExpr &linear_expr() const {
      return this->problem_->algebraic_cons_[this->index_].linear_expr;
    }

    // Returns the nonlinear part of a constraint expression.
    NumericExpr nonlinear_expr() const {
      return this->problem_->nonlinear_cons_.empty() ?
            NumericExpr() : this->problem_->nonlinear_cons_[this->index_];
    }

    bool operator==(AlgebraicCon other) const {
      MP_ASSERT(this->problem_ == other.problem_,
                "comparing constraints from different problems");
      return this->index_ == other.index_;
    }
    bool operator!=(AlgebraicCon other) const {
      return !(*this == other);
    }
  };

  // A list of algebraic constraints.
  typedef List<AlgebraicCon> AlgebraicConList;

  // Returns the list of algebraic constraints.
  // It can be used to iterate over all algebraic constraints in a problem:
  //   for (auto con: problem.algebraic_cons()) {
  //     ...
  //   }
  AlgebraicConList algebraic_cons() const { return AlgebraicConList(this); }

  // Returns the algebraic constraint at the specified index.
  AlgebraicCon algebraic_con(int index) const {
    CheckIndex(index, num_algebraic_cons());
    return AlgebraicCon(this, index);
  }

  typedef LinearExprBuilder LinearConBuilder;

  // Adds an algebraic constraint.
  // Returns a handler for receiving linear terms in the constraint.
  LinearConBuilder AddCon(double lb, double ub, NumericExpr expr,
                          int num_linear_terms = 0);

  LinearConBuilder AddCon(double lb, double ub, int num_linear_terms = 0) {
    return AddCon(lb, ub, NumericExpr(), num_linear_terms);
  }

  // Adds a logical constraint.
  void AddCon(LogicalExpr expr) {
    MP_ASSERT(logical_cons_.size() < MP_MAX_PROBLEM_ITEMS,
              "too many logical constraints");
    logical_cons_.push_back(expr);
  }

  // Begins building a common expression (defined variable).
  // Returns a handler for receiving linear terms in the common expression.
  LinearExprBuilder BeginCommonExpr(int num_linear_terms) {
    linear_exprs_.push_back(LinearExpr());
    LinearExpr &linear = linear_exprs_.back();
    linear.Reserve(num_linear_terms);
    return LinearExprBuilder(&linear);
  }

  // Ends building a common expression.
  void EndCommonExpr(LinearExprBuilder, NumericExpr expr, int) {
    nonlinear_exprs_.push_back(expr);
  }

  // Sets a complementarity relation.
  void SetComplement(int con_index, int var_index, int flags);

  // Sets the initial value for a variable.
  void SetInitialValue(int var_index, double value) {
    CheckIndex(var_index, num_vars());
    if (initial_values_.size() <= static_cast<unsigned>(var_index)) {
      initial_values_.reserve(vars_.capacity());
      initial_values_.resize(num_vars());
    }
    initial_values_[var_index] = value;
  }

  // Sets the initial value for a dual variable.
  void SetInitialDualValue(int con_index, double value) {
    MP_ASSERT(0 <= con_index && con_index <= num_algebraic_cons(),
              "invalid index");
    if (initial_dual_values_.size() <= static_cast<unsigned>(con_index)) {
      initial_dual_values_.reserve(algebraic_cons_.capacity());
      initial_dual_values_.resize(num_algebraic_cons());
    }
    initial_dual_values_[con_index] = value;
  }

  struct ColumnSizeHandler {
    void Add(int) {
      // Ignore column sizes as the constraints are stored row-wise.
    }
  };

  // Returns a handler that receives column sizes in Jacobian.
  ColumnSizeHandler GetColumnSizeHandler() {
    return ColumnSizeHandler();
  }

  typedef Suffix *SuffixPtr;
  typedef mp::SuffixSet SuffixSet;

  typedef SuffixHandler<int> IntSuffixHandler;

  // Adds an integer suffix.
  // name: Suffix name that may not be null-terminated.
  IntSuffixHandler AddIntSuffix(fmt::StringRef name, int kind, int) {
    return AddSuffix<int>(name, kind);
  }

  typedef SuffixHandler<double> DblSuffixHandler;

  // Adds an double suffix.
  // name: Suffix name that may not be null-terminated.
  DblSuffixHandler AddDblSuffix(fmt::StringRef name, int kind, int) {
    return AddSuffix<double>(name, kind);
  }

  // Sets problem information and reserves memory for problem elements.
  void SetInfo(const ProblemInfo &info);
};

template <typename Alloc>
int BasicProblem<Alloc>::GetSuffixSize(int suffix_type) {
  switch (suffix_type) {
  default:
    MP_ASSERT(false, "invalid suffix type");
    // Fall through.
  case suf::VAR:
    return vars_.capacity();
  case suf::CON:
    return algebraic_cons_.capacity();
  case suf::OBJ:
    return linear_objs_.capacity();
  case suf::PROBLEM:
    return 1;
  }
}

template <typename Alloc>
typename BasicProblem<Alloc>::LinearObjBuilder BasicProblem<Alloc>::AddObj(
    obj::Type type, NumericExpr expr, int num_linear_terms) {
  MP_ASSERT(linear_objs_.size() < MP_MAX_PROBLEM_ITEMS, "too many objectives");
  is_obj_max_.push_back(type != obj::MIN);
  linear_objs_.push_back(LinearExpr());
  LinearExpr &linear_expr = linear_objs_.back();
  linear_expr.Reserve(num_linear_terms);
  if (expr) {
    if (nonlinear_objs_.empty()) {
      nonlinear_objs_.reserve(linear_objs_.capacity());
      nonlinear_objs_.resize(linear_objs_.size() - 1);
    }
    nonlinear_objs_.push_back(expr);
  }
  return LinearObjBuilder(&linear_expr);
}

template <typename Alloc>
typename BasicProblem<Alloc>::LinearConBuilder BasicProblem<Alloc>::AddCon(
    double lb, double ub, NumericExpr expr, int num_linear_terms) {
  MP_ASSERT(algebraic_cons_.size() < MP_MAX_PROBLEM_ITEMS,
            "too many algebraic constraints");
  algebraic_cons_.push_back(AlgebraicConInfo(lb, ub));
  AlgebraicConInfo &con = algebraic_cons_.back();
  con.linear_expr.Reserve(num_linear_terms);
  if (expr) {
    if (nonlinear_cons_.empty()) {
      nonlinear_cons_.reserve(algebraic_cons_.capacity());
      nonlinear_cons_.resize(algebraic_cons_.size() - 1);
    }
    nonlinear_cons_.push_back(expr);
  }
  return LinearConBuilder(&con.linear_expr);
}

template <typename Alloc>
void BasicProblem<Alloc>::SetComplement(
    int con_index, int var_index, int flags) {
  MP_ASSERT(0 <= con_index && con_index <= num_algebraic_cons(),
            "invalid index");
  if (compl_vars_.size() <= con_index) {
    compl_vars_.reserve(algebraic_cons_.capacity());
    compl_vars_.resize(algebraic_cons_.size());
  }
  compl_vars_[con_index] = var_index + 1u;
  double inf = std::numeric_limits<double>::infinity();
  AlgebraicConInfo &con = algebraic_cons_[con_index];
  con.lb = (flags & comp::INF_LB) != 0 ? -inf : 0;
  con.ub = (flags & comp::INF_UB) != 0 ?  inf : 0;
}

template <typename Alloc>
void BasicProblem<Alloc>::SetInfo(const ProblemInfo &info) {
  vars_.reserve(info.num_vars);
  is_var_int_.reserve(info.num_vars);
  is_obj_max_.reserve(info.num_objs);
  linear_objs_.reserve(info.num_objs);
  if (info.num_nl_objs != 0)
    nonlinear_objs_.reserve(info.num_objs);
  algebraic_cons_.reserve(info.num_algebraic_cons);
  if (info.num_compl_conds != 0)
    compl_vars_.reserve(info.num_algebraic_cons);
  if (info.num_nl_cons != 0)
    nonlinear_cons_.reserve(info.num_algebraic_cons);
  logical_cons_.reserve(info.num_logical_cons);
  int num_common_exprs = info.num_common_exprs();
  linear_exprs_.reserve(num_common_exprs);
  nonlinear_exprs_.reserve(num_common_exprs);
  funcs_.reserve(info.num_funcs);
}

typedef BasicProblem< std::allocator<char> > Problem;
}  // namespace mp

#endif  // MP_PROBLEM_H_
