/*
 Optimization problem

 Copyright (C) 2014 - 2016 AMPL Optimization Inc

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
#include <cmath>
#include <vector>

#include "mp/expr.h"
#include "mp/suffix.h"

/// Maximum index of a variable, objective or constraint.
#ifndef MP_MAX_PROBLEM_ITEMS
# define MP_MAX_PROBLEM_ITEMS \
  static_cast<std::size_t>(std::numeric_limits<int>::max())
#endif

namespace mp {

/// Linear expression (not affine: no constant term)
/// used in `mp::BasicProblem<>`
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
    void set_coef(double c) { coef_=c; }
    void operator*=(double n) { coef_*=n; }
  };
  std::vector<Term> terms_;

 public:
  LinearExpr() { }

  template <class CoefVec=std::vector<double>, class VarVec=std::vector<int> >
  LinearExpr(CoefVec&& c, VarVec&& v) {
    ConstructFrom(std::forward<CoefVec>(c), std::forward<VarVec>(v));
  }
  template <int N>
  LinearExpr(const std::array<double, N>& c, const std::array<int, N>& v) {
    ConstructFrom(c, v);
  }

  int num_terms() const { return static_cast<int>(terms_.size()); }
  int capacity() const { return static_cast<int>(terms_.capacity()); }

  int var_index(int i) const { return terms_[i].var_index(); }
  double coef(int i) const { return terms_[i].coef(); }
  void set_coef(int i, double c) { terms_[i].set_coef(c); }

  typedef std::vector<Term>::const_iterator const_iterator;

  const_iterator begin() const { return terms_.begin(); }
  const_iterator end() const { return terms_.end(); }

  typedef std::vector<Term>::iterator iterator;

  iterator begin() { return terms_.begin(); }
  iterator end() { return terms_.end(); }

  void AddTerm(int var_index, double coef) {
    terms_.push_back(Term(var_index, coef));
  }

  void AddTerms(const LinearExpr& li) {
    terms_.insert(end(), li.begin(), li.end());
  }

  template <class CoefVec=std::vector<double>, class VarVec=std::vector<int> >
  void ConstructFrom(CoefVec&& c, VarVec&& v) {
    assert(c.size()==v.size());
    assert(0==num_terms());
    Reserve(c.size());
    for (size_t i=0; i<c.size(); ++i)
      AddTerm(v[i], c[i]);
  }

  void Reserve(std::size_t num_terms) {
    terms_.reserve(num_terms);
  }

  void SortTerms();
};

/// Placeholder for whatever params of BasicProblem<>
template < class A = int >
struct BasicProblemParams {
  using Alloc = A;
};

/// An optimization problem
template <typename ProblemParams = BasicProblemParams<> >
class BasicProblem : public ExprFactory, public SuffixManager {
 public:
  typedef mp::Function Function;
  typedef mp::Expr Expr;
  typedef mp::NumericExpr NumericExpr;
  typedef mp::LogicalExpr LogicalExpr;
  typedef mp::CountExpr CountExpr;
  typedef mp::Reference Reference;
  typedef internal::ExprTypes ExprTypes;

 private:
  /// Names
  std::vector<std::string> var_names_;
  std::vector<std::string> con_names_;
  std::vector<std::string> obj_names_;

  /// A variable.
  struct Var {
    double lb;
    double ub;
    Var(double lb, double ub) : lb(lb), ub(ub) { assert(lb<=ub); }
  };
  std::vector<Var> vars_;

  /// Packed variable type information.
  /// is_var_int_[i] specifies whether variable i is integer.
  std::vector<bool> is_var_int_;

  /// Variable deleted flags.
  /// Only markers
  std::vector<bool> is_var_deleted_;

  /// Packed objective type information.
  /// is_obj_max_[i] specifies whether objective i is maximization.
  std::vector<bool> is_obj_max_;

  /// Linear parts of objective expessions.
  std::vector<LinearExpr> linear_objs_;

  /// Nonlinear parts of objective expressions.
  /// The array can be empty if the problem is linear.
  std::vector<NumericExpr> nonlinear_objs_;

  /// Algebraic constraint information.
  struct AlgebraicConInfo {
    /// Linear part of an algebraic constraint expression.
    /// Nonlinear parts are stored in nonlinear_cons_ to avoid overhead
    /// for linear problems.
    LinearExpr linear_expr;
    double lb;
    double ub;
    AlgebraicConInfo() : lb(0), ub(0) {}
    AlgebraicConInfo(double lb, double ub) : lb(lb), ub(ub) {}
  };
  std::vector<AlgebraicConInfo> algebraic_cons_;

  /// Nonlinear parts of algebraic constraint expressions.
  /// The array can be empty if the problem is linear.
  std::vector<NumericExpr> nonlinear_cons_;

  /// Algebraic constraint deletion marker.
  std::vector<bool> is_alg_con_deleted_;

  /// Information about complementarity conditions.
  /// compl_vars_[i] > 0 means constraint i complements variable
  /// compl_vars_[i] - 1. The array can be empty if there are no
  /// complementarity conditions.
  std::vector<unsigned> compl_vars_;

  /// Logical constraint expressions.
  std::vector<LogicalExpr> logical_cons_;

  /// Linear parts of common expressions.
  std::vector<LinearExpr> linear_exprs_;

  /// Nonlinear parts of common expressions.
  std::vector<NumericExpr> nonlinear_exprs_;

  /// Initial values for variables.
  std::vector<double> initial_values_;
  /// Initial values for variables: sparsity.
  std::vector<int> iv_set_;

  /// Initial values for dual variables.
  std::vector<double> initial_dual_values_;
  /// Initial values for dual variables: sparsity.
  std::vector<int> idv_set_;

  /// Set nonlinear objective expression.
  void SetNonlinearObjExpr(int obj_index, NumericExpr expr) {
    internal::CheckIndex(obj_index, linear_objs_.size());
    if (nonlinear_objs_.size() <= static_cast<std::size_t>(obj_index))
      nonlinear_objs_.resize(obj_index + 1);
    nonlinear_objs_[obj_index] = expr;
  }

  /// Set nonlinear algebraic constraint expression.
  void SetNonlinearConExpr(int con_index, NumericExpr expr) {
    internal::CheckIndex(con_index, algebraic_cons_.size());
    if (nonlinear_cons_.size() <= static_cast<std::size_t>(con_index))
      nonlinear_cons_.resize(con_index + 1);
    nonlinear_cons_[con_index] = expr;
  }

  /// Mark algebraic constraint as deleted.
  void MarkAlgConDeleted(int con_index) {
    internal::CheckIndex(con_index, algebraic_cons_.size());
    if (is_alg_con_deleted_.size() <= static_cast<std::size_t>(con_index))
      is_alg_con_deleted_.resize(num_algebraic_cons());
    is_alg_con_deleted_[con_index] = true;
  }

  /// Mark variable as deleted.
  void MarkVarDeleted(int var_index) {
    if (is_var_deleted_.size() <= static_cast<unsigned>(var_index)) {
      is_var_deleted_.reserve(vars_.capacity());
      is_var_deleted_.resize(num_vars());
    }
    is_var_deleted_[var_index] = true;
  }

  /// Sets the initial value for a variable.
  void SetInitialValue(int var_index, double value) {
    if (initial_values_.size() <= static_cast<unsigned>(var_index)) {
      initial_values_.reserve(vars_.capacity());
      initial_values_.resize(num_vars());
      iv_set_.reserve(vars_.capacity());
      iv_set_.resize(num_vars());
    }
    initial_values_[var_index] = value;
    iv_set_[var_index] = 1;
  }

  /// Sets the initial value for a dual variable.
  void SetInitialDualValue(int con_index, double value) {
    MP_ASSERT(0 <= con_index && con_index < num_algebraic_cons(),
              "invalid index");
    if (initial_dual_values_.size() <= static_cast<unsigned>(con_index)) {
      initial_dual_values_.reserve(algebraic_cons_.capacity());
      initial_dual_values_.resize(num_algebraic_cons());
      idv_set_.reserve(algebraic_cons_.capacity());
      idv_set_.resize(num_algebraic_cons());
    }
    initial_dual_values_[con_index] = value;
    idv_set_[con_index] = 1;
  }

public:

  ////////////////////////////////////////////////////////////////////
  /// BasicProblemItem
  ////////////////////////////////////////////////////////////////////
  template <typename ProblemType>
  struct BasicProblemItem {
    ProblemType *problem_;
    ProblemType *problem() const { return problem_; }
    int index_;
    int index() const { return index_; }

    typedef ProblemType Problem;

    BasicProblemItem(ProblemType *p, int index)
      : problem_(p), index_(index) {}
  };

  typedef BasicProblemItem<const BasicProblem> ProblemItem;
  typedef BasicProblemItem<BasicProblem> MutProblemItem;

  /// Deprecated
  const std::vector<bool>& IsVarInt() const {
    return is_var_int_;
  }

  /// An optimization variable.
  template <typename Item>
  class BasicVariable : private Item {
   private:
    friend class BasicProblem;

    BasicVariable(typename Item::Problem *p, int index) : Item(p, index) {}

    static int num_items(const BasicProblem &p) {
      return p.num_vars();
    }

   public:
    using Item::index;

    /// Whether the variable is marked as deleted.
    bool is_marked_deleted() const {
        return index() < this->problem_->is_var_deleted_.size() ?
              this->problem_->is_var_deleted_[index()] : false;
    }

    /// Returns the lower bound on the variable.
    double lb() const {
      return this->problem_->vars_[this->index_].lb;
    }

    /// Returns the upper bound on the variable.
    double ub() const {
      return this->problem_->vars_[this->index_].ub;
    }

    /// Returns the type of the variable.
    var::Type type() const {
      return this->problem_->is_var_int_[this->index_] ?
          var::INTEGER : var::CONTINUOUS;
    }

    /// Returns the value of the variable.
    double value() const {
      std::size_t index = this->index_;
      return index < this->problem_->initial_values_.size() ?
            this->problem_->initial_values_[index] : 0;
    }

    template <typename OtherItem>
    bool operator==(BasicVariable<OtherItem> other) const {
      MP_ASSERT(this->problem_ == other.problem_,
                "comparing variables from different problems");
      return this->index_ == other.index_;
    }

    template <typename OtherItem>
    bool operator!=(BasicVariable<OtherItem> other) const {
      return !(*this == other);
    }
  };

  /// An objective.
  template <typename Item>
  class BasicObjective : private Item {
   private:
    friend class BasicProblem;
    friend class MutObjective;

    BasicObjective(typename Item::Problem *p, int index) : Item(p, index) {}

    static int num_items(const BasicProblem &p) {
      return p.num_objs();
    }

   public:
    /// Returns the objective index
    using Item::index;

    /// Returns the type of the objective.
    obj::Type type() const {
      return this->problem_->is_obj_max_[this->index_] ? obj::MAX : obj::MIN;
    }

    /// Returns the linear part of the objective expression.
    const LinearExpr &linear_expr() const {
      return this->problem_->linear_objs_[this->index_];
    }

    /// Returns the nonlinear part of the objective expression.
    NumericExpr nonlinear_expr() const {
      std::size_t index = this->index_;
      return index < this->problem_->nonlinear_objs_.size() ?
            this->problem_->nonlinear_objs_[index] : NumericExpr();
    }

    template <typename OtherItem>
    bool operator==(BasicObjective<OtherItem> other) const {
      MP_ASSERT(this->problem_ == other.problem_,
                "comparing objectives from different problems");
      return this->index_ == other.index_;
    }

    template <typename OtherItem>
    bool operator!=(BasicObjective<OtherItem> other) const {
      return !(*this == other);
    }
  };

  /// An algebraic constraint.
  /// This is a constraint of the form lb <= expr <= ub.
  template <typename Item>
  class BasicAlgebraicCon : private Item {
   private:
    friend class BasicProblem;

    BasicAlgebraicCon(typename Item::Problem *p, int index) : Item(p, index) {}

    static int num_items(const BasicProblem &p) {
      return p.num_algebraic_cons();
    }

   public:
    /// index(): constraint index
    using Item::index;

    /// Whether the alg con is marked as deleted.
    bool is_marked_deleted() const {
        return index() < this->problem_->is_alg_con_deleted_.size() ?
              this->problem_->is_alg_con_deleted_[index()] : false;
    }

    /// Returns the lower bound on the constraint.
    double lb() const {
      return this->problem_->algebraic_cons_[this->index_].lb;
    }

    /// Returns the upper bound on the constraint.
    double ub() const {
      return this->problem_->algebraic_cons_[this->index_].ub;
    }

    /// Returns the dual value.
    double dual() const {
      std::size_t index = this->index_;
      return index < this->problem_->initial_dual_values_.size() ?
            this->problem_->initial_dual_values_[index] : 0;
    }

    /// Returns the linear part of a constraint expression.
    const LinearExpr &linear_expr() const {
      return this->problem_->algebraic_cons_[this->index_].linear_expr;
    }

    /// Returns the nonlinear part of a constraint expression.
    NumericExpr nonlinear_expr() const {
      std::size_t index = this->index_;
      return index < this->problem_->nonlinear_cons_.size() ?
            this->problem_->nonlinear_cons_[index] : NumericExpr();
    }

    template <typename OtherItem>
    bool operator==(BasicAlgebraicCon<OtherItem> other) const {
      MP_ASSERT(this->problem_ == other.problem_,
                "comparing constraints from different problems");
      return this->index_ == other.index_;
    }

    template <typename OtherItem>
    bool operator!=(BasicAlgebraicCon<OtherItem> other) const {
      return !(*this == other);
    }
  };

 public:
  /** Constructs an empty optimization problem. */
  BasicProblem() {}

  /** Placeholder, some APIs assume
   *  a ProblemBuilder to use a Solver */
  template <class Solver>
  explicit BasicProblem(const Solver &) {}

  /** Returns the number of variables. */
  int num_vars() const { return static_cast<int>(vars_.size()); }

  /// Normal variable name
  const std::string& var_name(int i) {
    assert(0<=i && i<num_vars());
    return item_name(i, var_names_, num_vars(), "_x[");
  }
  /// Defined variable name
  const std::string& dvar_name(int i) {
    assert(0<=i && i<num_common_exprs());
    i += num_vars();
    return item_name(i, var_names_,
                     num_vars()+num_common_exprs(), "_x[",
                     num_vars(), "_sdvar[");
  }
  /// Constraint name
  const std::string& con_name(int i) {
    assert(0<=i && i<num_cons());
    return item_name(i, con_names_,
                     num_cons(),  "_CON",
                     num_algebraic_cons(), "_LCON");
  }
  /// Objective name
  const std::string& obj_name(int i) {
    assert(0<=i && i<num_objs());
    return item_name(i, obj_names_, num_objs(), "_OBJ");
  }

protected:
  /// Return names[i].
  /// Generate names[old_size...n-1] if needed,
  /// so that the whole vector is valid,
  /// even if only some variables/cons/objs were queried.
  /// From [n2], name 'stub2'[i-n2] is given.
  /// If stub/stub2 ends with a '[', ']' is added after the index.
  static const std::string& item_name(
      int i, std::vector<std::string>& names,
      int n, const char* stub,
      int n2=std::numeric_limits<int>::max(),
      const char* stub2=nullptr);

public:
  /** Returns the variable names (if present).
   *  After normal variables follow defined variables.
   */
  const std::vector<std::string>& var_names() { return var_names_; }
  /** Returns the constraint names (if present). */
  const std::vector<std::string>& con_names() { return con_names_; }
  /** Returns the objective names (if present). */
  const std::vector<std::string>& obj_names() { return obj_names_; }

  /// Variable namer
  class VarNamer {
  public:
    /// Construct
    VarNamer(BasicProblem& p) : p_(p) { }
    /// Normal var name
    const std::string& vname(int i) const
    { return p_.var_name(i); }
    /// Defined var name
    const std::string& dvname(int i) const
    { return p_.dvar_name(i); }
  private:
    BasicProblem& p_;
  };

  /// Obtain variable  namer
  VarNamer GetVarNamer() { return VarNamer(*this); }

  /** Returns the number of objectives. */
  int num_objs() const { return static_cast<int>(linear_objs_.size()); }

  /** Returns total number of constraints from the NL file. */
  int num_cons() const {
    return num_algebraic_cons() + num_logical_cons();
  }

  /** Returns the number of algebraic constraints. */
  int num_algebraic_cons() const {
    return static_cast<int>(algebraic_cons_.size());
  }

  /** Returns the number of logical constraints. */
  int num_logical_cons() const {
    return static_cast<int>(logical_cons_.size());
  }

  /// Return true if the problem has nonlinear constraints.
  bool has_nonlinear_cons() const {
    return !nonlinear_cons_.empty();
  }

  /** Returns the number of common expressions. */
  int num_common_exprs() const {
    return static_cast<int>(linear_exprs_.size());
  }

  /// An optimization variable.
  typedef BasicVariable<ProblemItem> Variable;

  /// A mutable variable.
  class MutVariable : public BasicVariable<MutProblemItem> {
   private:
    friend class BasicProblem;

    MutVariable(BasicProblem *p, int index)
      : BasicVariable<MutProblemItem>(p, index) {}

   public:
    operator Variable() const {
      return Variable(this->problem_, this->index_);
    }

    void set_lb(double lb) const {
      this->problem_->vars_[this->index_].lb = lb;
    }
    void set_ub(double ub) const {
      this->problem_->vars_[this->index_].ub = ub;
    }

    /// Sets the initial value.
    void set_value(double value) const {
      this->problem_->SetInitialValue(this->index_, value);
    }
  };

  /** A pair of iterators to problem elements. */
  template <typename T>
  class Range {
   private:
    const BasicProblem *problem_;

    friend class BasicProblem;

    explicit Range(const BasicProblem *p) : problem_(p) {}

   public:
    class iterator : public std::iterator<std::forward_iterator_tag, T> {
     private:
      T item_;

      friend class Range<T>;

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

    /** Returns an iterator to the first element in the range. */
    iterator begin() const {
      return iterator(problem_, 0);
    }

    /**
      Returns an iterator to the element following the last element
      in the range. An attempt to access this element will result in
      assertion failure if assertions are enabled and undefined behavior
      otherwise.
     */
    iterator end() const {
      return iterator(problem_, T::num_items(*problem_));
    }
  };

  /** A range of variables. */
  typedef Range<Variable> VarRange;

  /**
    \rst
    Returns a range representing all variables in this optimization problem.
    It can be used for iterating over variables::

      for (auto var: problem.vars()) {
        ...
      }
    \endrst
   */
  VarRange vars() const { return VarRange(this); }

  /// Returns the variable at the specified index.
  Variable var(int index) const {
    internal::CheckIndex(index, num_vars());
    return Variable(this, index);
  }
  MutVariable var(int index) {
    internal::CheckIndex(index, num_vars());
    return MutVariable(this, index);
  }
  
  /// Adds a variable.
  Variable AddVar(double lb, double ub, var::Type type = var::CONTINUOUS) {
    std::size_t index = vars_.size();
    MP_ASSERT(index < MP_MAX_PROBLEM_ITEMS, "too many variables");
    vars_.push_back(Var(lb, ub));
    is_var_int_.push_back(type != var::CONTINUOUS);
    return Variable(this, static_cast<int>(index));
  }

  void AddVars(int num_vars, var::Type type) {
    MP_ASSERT(num_vars >= 0, "invalid size");
    std::size_t new_size = val(SafeInt<int>(vars_.size()) + num_vars);
    vars_.resize(new_size, Var(0, 0));
    is_var_int_.resize(new_size, type != var::CONTINUOUS);
  }

  void AddVars(int num_vars,
               const double* lb, const double* ub, const var::Type* type) {
    MP_ASSERT(num_vars >= 0, "invalid size");
    std::size_t new_size = val(SafeInt<int>(vars_.size()) + num_vars);
    vars_.reserve(new_size);
    is_var_int_.reserve(new_size);
    for (int i=0; i!=num_vars; ++i) {
      vars_.push_back(Var(lb[i], ub[i]));
      is_var_int_.push_back(type[i] != var::CONTINUOUS);
    }
  }

  /// Add vector of variables. Type: var::CONTINUOUS by default
  std::vector<int> AddVars(std::size_t nvars,
                           double lb=-INFINITY, double ub=INFINITY,
                           var::Type type = var::CONTINUOUS) {
    std::vector<int> newVars(nvars);
    for (std::size_t  i=0; i<nvars; ++i)
      newVars[i] = AddVar(lb, ub, type).index();
    return newVars;
  }

  /// Set name vectors
  void SetVarNames(std::vector<std::string> names) {
    assert((size_t)(num_vars() + num_common_exprs()) == names.size());
    var_names_ = std::move( names );
  }
  void SetConNames(std::vector<std::string> names) {
    assert((size_t)num_cons() == names.size());
    con_names_ = std::move( names );
  }
  void SetObjNames(std::vector<std::string> names) {
    assert((size_t)num_objs() == names.size());
    obj_names_ = std::move( names );
  }

  class LinearExprBuilder {
   private:
    LinearExpr *expr_;

   public:
    explicit LinearExprBuilder(LinearExpr *expr) : expr_(expr) {}

    void AddTerm(int var_index, double coef) {
      expr_->AddTerm(var_index, coef);
    }
  };

  typedef LinearExprBuilder LinearObjBuilder;

  /// An objective.
  typedef BasicObjective<ProblemItem> Objective;

  /// A mutable objective.
  class MutObjective : public BasicObjective<MutProblemItem> {
   private:
    friend class BasicProblem;

    MutObjective(BasicProblem *p, int index)
      : BasicObjective<MutProblemItem>(p, index) {}

   public:
    operator Objective() const {
      return Objective(this->problem_, this->index_);
    }

    void set_type(obj::Type type) const {
      this->problem_->is_obj_max_[this->index_] = (type == obj::MAX);
    }

    /// Returns the linear part of the objective expression.
    LinearExpr &linear_expr() const {
      return this->problem_->linear_objs_[this->index_];
    }

    /// Sets the linear part of the objective expression.
    LinearObjBuilder set_linear_expr(int num_linear_terms) const {
      LinearExpr &expr = linear_expr();
      expr.Reserve(num_linear_terms);
      return LinearObjBuilder(&expr);
    }

    /// Sets the extra info of the objective expression.
    template <class EI>
    void set_extra_info(EI ei) const {
      this->problem_->SetObjExtraInfo(this->index_, std::forward<EI>(ei));
    }

    /// Sets the nonlinear part of the objective expression.
    void set_nonlinear_expr(NumericExpr expr) const {
      this->problem_->SetNonlinearObjExpr(this->index_, expr);
    }

    /// Unsets the nonlinear part of the objective.
    void unset_nonlinear_expr() {
      this->problem_->SetNonlinearObjExpr(this->index_, NumericExpr());
      assert(!this->nonlinear_expr());
    }

  };

  /** A range of objectives. */
  typedef Range<Objective> ObjRange;

  /**
    \rst
    Returns a range representing all objectives in this optimization problem.
    It can be used for iterating over objectives::

      for (auto obj: problem.objs()) {
        ...
      }
    \endrst
   */
  ObjRange objs() const { return ObjRange(this); }

  /// Returns the objective at the specified index.
  Objective obj(int index) const {
    internal::CheckIndex(index, num_objs());
    return Objective(this, index);
  }

  /// Returns the mutable objective at the specified index.
  MutObjective obj(int index) {
    internal::CheckIndex(index, num_objs());
    return MutObjective(this, index);
  }

  /// Adds an objective.
  /// Returns a builder for the linear part of an objective expression.
  LinearObjBuilder AddObj(obj::Type type, NumericExpr expr,
                          int num_linear_terms = 0);

  LinearObjBuilder AddObj(obj::Type type, int num_linear_terms = 0) {
    return AddObj(type, NumericExpr(), num_linear_terms);
  }

  void AddObjs(int num_objs) {
    linear_objs_.resize(num_objs);
    is_obj_max_.resize(num_objs);
  }

  typedef LinearExprBuilder LinearConBuilder;

  /// An algebraic constraint.
  typedef BasicAlgebraicCon<ProblemItem> AlgebraicCon;

  /// A mutable algebraic constraint.
  class MutAlgebraicCon : public BasicAlgebraicCon<MutProblemItem> {
   private:
    friend class BasicProblem;

    MutAlgebraicCon(BasicProblem *p, int index)
      : BasicAlgebraicCon<MutProblemItem>(p, index) {}

   public:
    operator AlgebraicCon() const {
      return AlgebraicCon(this->problem_, this->index_);
    }

    /// Sets the lower bound on the constraint.
    void set_lb(double lb) const {
      this->problem_->algebraic_cons_[this->index_].lb = lb;
    }

    /// Sets the upper bound on the constraint.
    void set_ub(double ub) const {
      this->problem_->algebraic_cons_[this->index_].ub = ub;
    }

    /// Sets the initial dual value.
    void set_dual(double value) const {
      this->problem_->SetInitialDualValue(this->index_, value);
    }

    /// Returns the linear part of the constraint expression.
    LinearExpr &linear_expr() const {
      return this->problem_->algebraic_cons_[this->index_].linear_expr;
    }

    /// Sets the linear part of the objective expression.
    LinearConBuilder set_linear_expr(int num_linear_terms) const {
      LinearExpr &expr = linear_expr();
      expr.Reserve(num_linear_terms);
      return LinearConBuilder(&expr);
    }

    /// Sets the nonlinear part of the constraint expression.
    void set_nonlinear_expr(NumericExpr expr) const {
      if (expr)
        this->problem_->SetNonlinearConExpr(this->index_, expr);
    }

    /// Unsets the nonlinear part of the constraint expression.
    void unset_nonlinear_expr() {
      this->problem_->SetNonlinearConExpr(this->index_, NumericExpr());
      assert(!this->nonlinear_expr());
    }

  };

  /** A range of algebraic constraints. */
  typedef Range<AlgebraicCon> AlgebraicConRange;

  /**
    \rst
    Returns a range representing all algebraic constraints in this
    optimization problem. It can be used for iterating over algebraic
    constraints::

      for (auto con: problem.algebraic_cons()) {
        ...
      }
    \endrst
   */
  AlgebraicConRange algebraic_cons() const { return AlgebraicConRange(this); }

  /// Returns the algebraic constraint at the specified index.
  AlgebraicCon algebraic_con(int index) const {
    internal::CheckIndex(index, num_algebraic_cons());
    return AlgebraicCon(this, index);
  }

  /// Returns the mutable algebraic constraint at the specified index.
  MutAlgebraicCon algebraic_con(int index) {
    internal::CheckIndex(index, num_algebraic_cons());
    return MutAlgebraicCon(this, index);
  }

  /// Adds an algebraic constraint.
  /// Returns a builder for the linear part of a constraint expression.
  MutAlgebraicCon AddCon(double lb, double ub) {
    std::size_t num_cons = algebraic_cons_.size();
    MP_ASSERT(num_cons < MP_MAX_PROBLEM_ITEMS,
              "too many algebraic constraints");
    algebraic_cons_.push_back(AlgebraicConInfo(lb, ub));
    return MutAlgebraicCon(this, static_cast<int>(num_cons));
  }

  void AddAlgebraicCons(int num_cons) {
    algebraic_cons_.resize(num_cons);
  }

  /// A logical constraint.
  template <typename Item>
  class BasicLogicalCon : private Item {
   private:
    friend class BasicProblem;

    BasicLogicalCon(typename Item::Problem *p, int index) : Item(p, index) {}

    static int num_items(const BasicProblem &p) {
      return p.num_logical_cons();
    }

   public:
    /// Returns the constraint expression.
    LogicalExpr expr() const {
      return this->problem_->logical_cons_[this->index_];
    }

    template <typename OtherItem>
    bool operator==(BasicLogicalCon<OtherItem> rhs) const {
      MP_ASSERT(this->problem_ == rhs.problem_,
                "comparing constraints from different problems");
      return this->index_ == rhs.index_;
    }

    template <typename OtherItem>
    bool operator!=(BasicLogicalCon<OtherItem> rhs) const {
      return !(*this == rhs);
    }
  };

  typedef BasicLogicalCon<ProblemItem> LogicalCon;

  class MutLogicalCon : public BasicLogicalCon<MutProblemItem> {
   private:
    friend class BasicProblem;

    MutLogicalCon(BasicProblem *p, int index)
      : BasicLogicalCon<MutProblemItem>(p, index) {}

   public:
    operator LogicalCon() const {
      return LogicalCon(this->problem_, this->index_);
    }

    void set_expr(LogicalExpr expr) {
      this->problem_->logical_cons_[this->index_] = expr;
    }
  };

  /** A range of logical constraints. */
  typedef Range<LogicalCon> LogicalConRange;

  /**
    \rst
    Returns a range representing all logical constraints in this
    optimization problem. It can be used for iterating over logical
    constraints::

      for (auto con: problem.logical_cons()) {
        ...
      }
    \endrst
   */
  LogicalConRange logical_cons() const { return LogicalConRange(this); }

  /// Returns the logical constraint at the specified index.
  LogicalCon logical_con(int index) const {
    internal::CheckIndex(index, num_logical_cons());
    return LogicalCon(this, index);
  }

  /// Returns the mutable logical constraint at the specified index.
  MutLogicalCon logical_con(int index) {
    internal::CheckIndex(index, num_logical_cons());
    return MutLogicalCon(this, index);
  }

  /// Adds a logical constraint.
  void AddCon(LogicalExpr expr) {
    MP_ASSERT(logical_cons_.size() < MP_MAX_PROBLEM_ITEMS,
              "too many logical constraints");
    logical_cons_.push_back(expr);
  }

  void AddLogicalCons(int num_cons) {
    logical_cons_.resize(num_cons);
  }

  /// A common expression.
  template <typename Item>
  class BasicCommonExpr : private Item {
   private:
    friend class BasicProblem;

    BasicCommonExpr(typename Item::Problem *p, int index) : Item(p, index) {}

   public:
    /// Returns the linear part of the common expression.
    const LinearExpr &linear_expr() const {
      return this->problem_->linear_exprs_[this->index_];
    }

    /// Returns the nonlinear part of the common expression.
    NumericExpr nonlinear_expr() const {
      std::size_t index = this->index_;
      return index < this->problem_->nonlinear_exprs_.size() ?
            this->problem_->nonlinear_exprs_[index] : NumericExpr();
    }

    template <typename OtherItem>
    bool operator==(BasicCommonExpr<OtherItem> other) const {
      MP_ASSERT(this->problem_ == other.problem_,
                "comparing expressions from different problems");
      return this->index_ == other.index_;
    }

    template <typename OtherItem>
    bool operator!=(BasicCommonExpr<OtherItem> other) const {
      return !(*this == other);
    }
  };

  typedef BasicCommonExpr<ProblemItem> CommonExpr;

  class MutCommonExpr : public BasicCommonExpr<MutProblemItem> {
   private:
    friend class BasicProblem;

    MutCommonExpr(BasicProblem *p, int index)
      : BasicCommonExpr<MutProblemItem>(p, index) {}

   public:
    LinearExprBuilder set_linear_expr(int num_linear_terms) const {
      LinearExpr &linear = this->problem_->linear_exprs_[this->index_];
      linear.Reserve(num_linear_terms);
      return LinearExprBuilder(&linear);
    }

    void set_nonlinear_expr(NumericExpr expr) const {
      this->problem_->nonlinear_exprs_[this->index_] = expr;
    }

    void set_position(int) const {}
  };

  /// Returns the common expression at the specified index.
  CommonExpr common_expr(int index) const {
    internal::CheckIndex(index, num_common_exprs());
    return CommonExpr(this, index);
  }
  MutCommonExpr common_expr(int index) {
    internal::CheckIndex(index, num_common_exprs());
    return MutCommonExpr(this, index);
  }

  /// Adds a common expression (defined variable).
  MutCommonExpr AddCommonExpr(NumericExpr expr) {
    std::size_t num_exprs = linear_exprs_.size();
    MP_ASSERT(num_exprs < MP_MAX_PROBLEM_ITEMS, "too many expressions");
    linear_exprs_.push_back(LinearExpr());
    nonlinear_exprs_.push_back(expr);
    return MutCommonExpr(this, static_cast<int>(num_exprs));
  }

  void AddCommonExprs(int num_exprs) {
    MP_ASSERT(num_exprs >= 0, "invalid size");
    std::size_t new_size = val(SafeInt<int>(linear_exprs_.size()) + num_exprs);
    linear_exprs_.resize(new_size, LinearExpr());
    nonlinear_exprs_.resize(new_size, NumericExpr());
  }

  /// Sets a complementarity condition.
  void SetComplementarity(int con_index, int var_index, ComplInfo info);

  /// Returns true if the problem has complementarity conditions.
  bool HasComplementarity() const { return !compl_vars_.empty(); }

  /// Returns complementarity variable for algebraic constraint \a i.
  /// Result > 0 means constraint \a i complements variable
  /// Result - 1.
  int GetComplementarityVariable(int i) const {
    return compl_vars_.size()>(size_t)i ? compl_vars_.at(i) : -1;
  }

  /// Variables' initial values
  ArrayRef<double> InitialValues() const { return initial_values_; }
  /// Variables' initial values: sparsity pattern
  ArrayRef<int> InitialValuesSparsity() const { return iv_set_; }

  /// Initial dual values
  ArrayRef<double> InitialDualValues() const { return initial_dual_values_; }
  /// Dual initial values: sparsity pattern
  ArrayRef<int> InitialDualValuesSparsity() const { return idv_set_; }

  /////////////////////////////////////////////////////////////////////////
  /// Suffixes
  /////////////////////////////////////////////////////////////////////////
  template <typename T>
  class SuffixHandler {
   private:
    BasicMutSuffix<T> suffix_;

   public:
    explicit SuffixHandler(BasicMutSuffix<T> s) : suffix_(s) {}

    /// Safe bool type.
    typedef void (internal::SuffixBase::*SafeBool)() const;
    operator SafeBool() const { return suffix_; }

    /// Sets the suffix value.
    void SetValue(int index, T value) {
      suffix_.set_value(index, value);
    }
  };

  int GetSuffixSize(suf::Kind kind);

  template <typename T>
  SuffixHandler<T> AddSuffix(fmt::StringRef name, int kind) {
    auto main_kind = (suf::Kind)(kind & suf::KIND_MASK);
    return SuffixHandler<T>(
          suffixes(main_kind).template Add<T>(name, kind, GetSuffixSize(main_kind)));
  }

  typedef SuffixHandler<int> IntSuffixHandler;

  /// Adds an integer suffix.
  /// name: Suffix name that may not be null-terminated.
  IntSuffixHandler AddIntSuffix(fmt::StringRef name, int kind, int=0) {
    return AddSuffix<int>(name, kind);
  }

  typedef SuffixHandler<double> DblSuffixHandler;

  /// Adds a double suffix.
  /// name: Suffix name that may not be null-terminated.
  DblSuffixHandler AddDblSuffix(fmt::StringRef name, suf::Kind kind, int) {
    return AddSuffix<double>(name, kind);
  }

  ////////////////////////// HIGH-LEVEL SUFFIX I/O //////////////////////////////
  template <class T>
  void ReportSuffix(const SuffixDef<T>& sufdef,
                    ArrayRef<T> values) {
    if (values.empty())
      return;
    auto suf = FindOrCreateSuffix(sufdef);
    auto suf_size = suf.num_values();
    /// Check this because Converter or solver can add more variables
    assert(suf_size <= (int)values.size());
    for (auto i=suf_size; i--; ) {
      suf.set_value(i, values[i]);
    }
  }

  template <class T>
  ArrayRef<T> ReadSuffix_OneTypeOnly(const SuffixDef<T>& sufdef) {
    auto suf = FindSuffix(sufdef);
    if (!suf)
      return {};
    return suf.get_values();
  }

  /// Read integer suffix
  ArrayRef<int> ReadIntSuffix(const SuffixDef<int>& sufdef)
  { return ReadSuffix_OneTypeOnly(sufdef); }

  /// Read double suffix.
  /// If absent but an integer suffix with the same name exists,
  /// take that
  ArrayRef<double> ReadDblSuffix(const SuffixDef<double>& sufdef) {
    auto suf_dbl = ReadSuffix_OneTypeOnly(sufdef);
    if (!suf_dbl) {
      auto suf_int = ReadSuffix_OneTypeOnly(sufdef.to_type<int>());
      if (suf_int)
        return std::vector<double>(suf_int.begin(), suf_int.end());
    }
    return suf_dbl;
  }

  template <class T>
  BasicMutSuffix<T> FindSuffix(const SuffixDef<T>& sufdef) {
    auto main_kind = (suf::Kind)(sufdef.kind() & suf::KIND_MASK);
    auto suf_raw = suffixes(main_kind).Find(sufdef.name());
    if (suf_raw)
      return Cast< BasicMutSuffix<T> >( suf_raw );
    return BasicMutSuffix<T>();
  }

  template <class T>
  BasicMutSuffix<T> FindOrCreateSuffix(const SuffixDef<T>& sufdef) {
    auto main_kind = (suf::Kind)(sufdef.kind() & suf::KIND_MASK);
    auto suf_raw = FindSuffix(sufdef);
    auto suf_size = GetSuffixSize(main_kind);    // can be < values.size()
    if (suf_raw) {
      suf_raw.or_kind(suf::OUTPUT);
      return suf_raw;
    }
    return suffixes(main_kind).template
            Add<T>(sufdef.name(), sufdef.kind() | suf::OUTPUT,
                   suf_size, sufdef.table());
  }

  ///////////////////////////////////////////////////////////////////////
  /// Sets problem information and reserves memory for problem elements.
  void SetInfo(const NLProblemInfo &info);

  /// Pushing the whole instance to a backend or converter.
  /// A responsible backend should handle all essential items
  template <class Backend>
  void PushModelTo(Backend& backend) const {
    InitProblemModificationPhase(backend);
    PushStandardMPItemsTo(backend);
    FinishProblemModificationPhase(backend);
  }

protected:
  template <class Backend>
  void PushStandardMPItemsTo(Backend& backend) const {
    PushVariablesTo(backend);
    PushCommonExprTo(backend);
    PushObjectivesTo(backend);
    PushAlgebraicConstraintsTo(backend);
    PushLogicalConstraintsTo(backend);
    PushComplementarityConstraintsTo(backend);
  }

  template <class Backend>
  void InitProblemModificationPhase(Backend& backend) const {
    backend.InitProblemModificationPhase();
  }

  template <class Backend>
  void PushVariablesTo(Backend& backend) const {
    const int nv = num_vars();
    for (int j = 0; j < nv; ++j) {
      backend.AddVariable(var(j));
    }
  }

  template <class Backend>
  void PushCommonExprTo(Backend& backend) const {
    const int nce = num_common_exprs();
    for (int i = 0; i < nce; ++i) {
      backend.AddCommonExpression(common_expr(i));
    }
  }

  template <class Backend>
  void PushObjectivesTo(Backend& backend) const {
    if (const int no = num_objs()) {
      for (int i = 0; i < no; ++i) {
        backend.AddObjective(obj(i));
      }
    }
  }

  template <class Backend>
  void PushAlgebraicConstraintsTo(Backend& backend) const {
    if (const int n_cons = num_algebraic_cons()) {
      for (int i = 0; i < n_cons; ++i) {
        backend.AddAlgebraicConstraint(algebraic_con(i));
      }
    }
  }

  template <class Backend>
  void PushLogicalConstraintsTo(Backend& backend) const {
    if (const int n_lcons = num_logical_cons()) {
      for (int i = 0; i < n_lcons; ++i) {
        backend.AddLogicalConstraint(logical_con(i));
      }
    }
  }

  template <class Backend>
  void PushComplementarityConstraintsTo(Backend& ) const {
    if (HasComplementarity()) {
      throw std::logic_error("mp::Problem cannot push complementarity to a backend yet.");
    }
  }

  template <class Backend>
  void FinishProblemModificationPhase(Backend& backend) const {
    backend.FinishProblemModificationPhase();
  }

public:
  typedef BasicProblem Builder;

  /// Returns the built problem. This is used for compatibility with the problem
  /// builder API.
  BasicProblem &problem() { return *this; }
};

/// A BasicProblem<> with default parameters
typedef BasicProblem< > Problem;
}  // namespace mp

#endif  // MP_PROBLEM_H_
