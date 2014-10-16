/*
 An ASL problem builder.

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

#ifndef MP_ASLBUILDER_H_
#define MP_ASLBUILDER_H_

#include "mp/arrayref.h"
#include "mp/format.h"
#include "mp/safeint.h"
#include "expr.h"
#include "problem.h"

struct Static;

namespace mp {

struct NLHeader;

namespace internal {

// An exception representing an ASL error.
class ASLError: public mp::Error {
private:
  int error_code_;

public:
  ASLError(int error_code, fmt::StringRef message) :
      Error(message), error_code_(error_code) {}

  int error_code() const { return error_code_; }
};

// Use standard opcodes when building expressions.
enum { ASL_STANDARD_OPCODES = 0x1000000 };

// Provides methods for building an ASL problem object.
class ASLBuilder {
 private:
  ASL *asl_;
  bool own_asl_;
  efunc **r_ops_;
  efunc *standard_opcodes_[expr::MAX_OPCODE + 1];
  int flags_;
  int nz_;
  int nderp_;
  static const double DVALUE[];
  SuffixView suffixes_[suf::NUM_KINDS];
  int var_index_;
  int obj_index_;
  int con_index_;
  int lcon_index_;

  // "Static" data for the functions in fg_read.
  Static *static_;

  FMT_DISALLOW_COPY_AND_ASSIGN(ASLBuilder);

  void SetBounds(double *lbs, double *&ubs, int index, double lb, double ub) {
    if (ubs) {
      lbs[index] = lb;
      ubs[index] = ub;
    } else {
      lbs[2 * index] = lb;
      lbs[2 * index + 1] = ub;
    }
  }

  template <typename ExprT>
  static void CheckKind(expr::Kind kind, const char *expr_name) {
    if (!internal::Is<ExprT>(kind))
      throw Error("invalid {} expression kind {}", expr_name, kind);
  }

  template <typename T>
  T *Allocate(SafeInt<int> size = SafeInt<int>(sizeof(T))) {
    if (size.value() > std::numeric_limits<int>::max())
      throw std::bad_alloc();
    return reinterpret_cast<T*>(mem_ASL(asl_, size.value()));
  }

  template <typename T>
  T *ZapAllocate(std::size_t size);

  template <typename T>
  T *AllocateSuffixValues(T *&values, int num_values, int nx, int nx1);

  // Sets objective or constraint expression; adapted from co_read.
  void SetObjOrCon(int index, cde *d, int *cexp1_end, ::expr *e, int **z);

  ::expr *MakeConstant(double value);

  ::expr *DoMakeUnary(expr::Kind kind, Expr arg);

  template <typename ExprT>
  ExprT MakeUnary(expr::Kind kind, Expr arg) {
    return Expr::Create<ExprT>(DoMakeUnary(kind, arg));
  }

  ::expr *DoMakeBinary(expr::Kind kind, Expr lhs, Expr rhs);

  template <typename ExprT>
  ExprT MakeBinary(expr::Kind kind, Expr lhs, Expr rhs, const char *name) {
    CheckKind<ExprT>(kind, name);
    return Expr::Create<ExprT>(DoMakeBinary(kind, lhs, rhs));
  }

  ::expr *MakeIf(expr::Kind kind,
      LogicalExpr condition, Expr true_expr, Expr false_expr);
  ::expr *MakeIterated(expr::Kind kind, ArrayRef<Expr> args);

  template <typename IteratedExpr>
  IteratedExpr MakeIterated(expr::Kind kind, ArrayRef<Expr> args) {
    return Expr::Create<IteratedExpr>(MakeIterated(kind, args));
  }

 public:
  int num_vars() const { return asl_->i.n_var_; }
  int num_cons() const { return asl_->i.n_con_; }

  class CallArgHandler {
   private:
    expr_f *expr_;
    int arg_index_, num_args_, num_constants_, num_symbolic_args_, num_ifsyms_;

    friend class ASLBuilder;

    explicit CallArgHandler(expr_f *expr, int num_args)
      : expr_(expr), arg_index_(0), num_args_(num_args),
        num_constants_(0), num_symbolic_args_(0), num_ifsyms_(0) {}

    void DoAddArg(Expr arg) {
      expr_->args[arg_index_] = arg.expr_;
      ++arg_index_;
    }

   public:
    void AddArg(NumericExpr arg) {
      DoAddArg(arg);
      if (Is<NumericConstant>(arg))
        ++num_constants_;
      // TODO
      //if (args[i].kind() == expr::IFSYM)
      //  ++num_ifsyms;
    }
    void AddArg(Expr arg) {
      DoAddArg(arg);
      ++num_symbolic_args_;
    }
  };

 private:
  CallArgHandler DoBeginCall(Function f, int num_args);

  void Init(ASL *asl);

 public:
  typedef mp::Expr Expr;
  typedef mp::NumericExpr NumericExpr;
  typedef mp::LogicalExpr LogicalExpr;
  typedef mp::Variable Variable;
  typedef mp::CountExpr CountExpr;

  explicit ASLBuilder(ASL *asl = 0) { Init(asl); }
  explicit ASLBuilder(Problem::Proxy proxy) {
    Init(proxy.asl_);
    proxy.asl_ = 0;
    own_asl_ = true;
    flags_ = proxy.flags_ | ASL_STANDARD_OPCODES | ASL_allow_CLP;
  }
  ~ASLBuilder();

  // Returns a built problem via proxy. No builder methods other than
  // the destructor should be called after calling this method.
  Problem::Proxy GetProblem() {
    if (!own_asl_)
      throw Error("ASL problem is not transferable");
    ASL *asl = asl_;
    own_asl_ = false;
    return Problem::Proxy(asl);
  }

  typedef ASLSuffixPtr SuffixPtr;
  typedef SuffixView SuffixSet;

  SuffixView &suffixes(int kind) {
    assert(kind < suf::NUM_KINDS);
    return suffixes_[kind];
  }

  void set_flags(int flags) { flags_ = flags; }
  void set_stub(const char *stub);

  // Initializes the ASL object in a similar way to jac0dim, but
  // doesn't read the .nl file as it is the responsibility of NLReader.
  // Instead it uses the information provided in NLHeader.
  void InitASL(const NLHeader &h);

  // Sets problem information.
  // Throws ASLError on error.
  void SetInfo(const ProblemInfo &pi);

  // Ends building the ASL object.
  void EndBuild();

  void AddVar(double lb, double ub, var::Type) {
    // TODO: check type
    SetBounds(asl_->i.LUv_, asl_->i.Uvx_, var_index_++, lb, ub);
  }

  template <typename Grad>
  class LinearExprBuilder {
   protected:
    ASLBuilder *builder_;
    Grad **term_;

   public:
    LinearExprBuilder(ASLBuilder *b, Grad **term) : builder_(b), term_(term) {}

    void AddTerm(int var_index, double coef) {
      Grad *og = builder_->Allocate<Grad>();
      *term_ = og;
      term_ = &og->next;
      og->next = 0;
      og->varno = var_index;
      og->coef = coef;
    }
  };

  typedef LinearExprBuilder<ograd> LinearObjBuilder;

  class LinearConBuilder : private LinearExprBuilder<cgrad> {
   private:
    int con_index_;
    double *a_vals_;
    int *a_rownos_;
    int *a_colstarts_;

   public:
    LinearConBuilder(ASLBuilder *b, int con_index);

    void AddTerm(int var_index, double coef) {
      if (!a_vals_)
        return LinearExprBuilder<cgrad>::AddTerm(var_index, coef);
      std::size_t elt_index = a_colstarts_[var_index]++;
      a_vals_[elt_index] = coef;
      a_rownos_[elt_index] = con_index_;
    }
  };

  class LinearVarBuilder {
   private:
    linpart *linpart_;

   public:
    explicit LinearVarBuilder(linpart *lp) : linpart_(lp) {}

    void AddTerm(int var_index, double coef) {
      linpart_->v.i = var_index;
      linpart_->fac = coef;
      ++linpart_;
    }
  };

  // Adds and objective.
  LinearObjBuilder AddObj(obj::Type type, NumericExpr expr, int);

  // Adds an algebraic constraint.
  LinearConBuilder AddCon(NumericExpr expr, double lb, double ub, int);

  // Adds a logical constraint.
  void AddCon(LogicalExpr expr);

  LinearVarBuilder GetLinearVarBuilder(int index, int num_terms);

  void SetComplement(int, int, int) {
    // TODO
  }

  void SetCommonExpr(int, NumericExpr, int) {
    // TODO
  }

  class ColumnSizeHandler {
   private:
    int *colstarts_;
    std::size_t index_;

   public:
    explicit ColumnSizeHandler(int *colstarts)
      : colstarts_(colstarts), index_(0) {}

    void Add(int size) {
      colstarts_[index_ + 1] = colstarts_[index_] + size;
      ++index_;
    }
  };

  ColumnSizeHandler GetColumnSizeHandler();

  void SetInitialValue(int, double) {
    // TODO
  }
  void SetInitialDualValue(int, double) {
    // TODO
  }

  Function AddFunction(const char *name, ufunc f, int num_args,
                       func::Type type = func::NUMERIC, void *info = 0);

  // Sets a function at the given index.
  // If the function with the specified name doesn't exist and the flag
  // ASL_allow_missing_funcs is not set, SetFunction throws ASLError.
  Function SetFunction(int index, fmt::StringRef name, int num_args,
                       func::Type type = func::NUMERIC);

  class SuffixHandler {
   private:
    int *int_values_;
    double *dbl_values_;
    int num_items_;

    void SetIntValue(int index, int value) {
      assert(0 <= index && index < num_items_);
      int_values_[index] = value;
    }
    void SetDblValue(int index, double value) {
      assert(0 <= index && index < num_items_);
      dbl_values_[index] = value;
    }

   public:
    explicit SuffixHandler(int *values = 0, int num_items = 0)
      : int_values_(values), dbl_values_(0), num_items_(num_items) {}
    explicit SuffixHandler(double *values, int num_items)
      : int_values_(0), dbl_values_(values), num_items_(num_items) {}

    void SetValue(int index, int value) {
      if (int_values_)
        SetIntValue(index, value);
      else if (dbl_values_)
        SetDblValue(index, value);
    }

    void SetValue(int index, double value) {
      if (int_values_)
        SetIntValue(index, static_cast<int>(value + 0.5));
      else if (dbl_values_)
        SetDblValue(index, value);
    }
  };

  SuffixHandler AddSuffix(int kind, int num_values, fmt::StringRef name);

  // The Make* methods construct expression objects. These objects are
  // local to the currently built ASL problem and shouldn't be used with
  // other problems. The expression objects are not accessible via the
  // problem API until they are added as a part of objective or constraint
  // expression. For this reason the methods are called Make* rather than Add*.

  NumericConstant MakeNumericConstant(double value) {
    return Expr::Create<NumericConstant>(MakeConstant(value));
  }

  Variable MakeVariable(int var_index);

  Variable MakeCommonExprRef(int index);

  UnaryExpr MakeUnary(expr::Kind kind, NumericExpr arg);

  BinaryExpr MakeBinary(expr::Kind kind, NumericExpr lhs, NumericExpr rhs) {
    return MakeBinary<BinaryExpr>(kind, lhs, rhs, "binary");
  }

  IfExpr MakeIf(LogicalExpr condition,
      NumericExpr true_expr, NumericExpr false_expr) {
    return Expr::Create<IfExpr>(
        MakeIf(expr::IF, condition, true_expr, false_expr));
  }

  class PLTermHandler {
   private:
    ::expr *expr_;
    double *data_;

    friend class ASLBuilder;

    explicit PLTermHandler(::expr *e) : expr_(e), data_(e->L.p->bs) {}

   public:

    void AddSlope(double slope) { *data_++ = slope; }
    void AddBreakpoint(double breakpoint) { *data_++ = breakpoint; }
  };

  PLTermHandler BeginPLTerm(int num_breakpoints);
  PiecewiseLinearExpr EndPLTerm(PLTermHandler h, NumericExpr var) {
    h.expr_->R.e = var.expr_;
    return Expr::Create<PiecewiseLinearExpr>(h.expr_);
  }

  PiecewiseLinearExpr MakePiecewiseLinear(int num_breakpoints,
      const double *breakpoints, const double *slopes, Variable var);

  CallArgHandler BeginCall(int func_index, int num_args) {
    Function f(asl_->i.funcs_[func_index]);
    if (!f)
      throw Error("function {} undefined", func_index);
    return DoBeginCall(f, num_args);
  }
  CallExpr EndCall(CallArgHandler h);

  CallExpr MakeCall(Function f, ArrayRef<Expr> args) {
    int num_args = SafeInt<int>(args.size()).value();
    CallArgHandler handler = DoBeginCall(f, num_args);
    for (int i = 0; i < num_args; ++i) {
      if (NumericExpr num = Cast<NumericExpr>(args[i]))
        handler.AddArg(num);
      else
        handler.AddArg(Cast<StringLiteral>(args[i]));
    }
    return EndCall(handler);
  }

  template <typename Arg>
  class ArgHandler {
   private:
    ::expr *expr_;
    int arg_index_;

    friend class ASLBuilder;

    ArgHandler(::expr *e) : expr_(e), arg_index_(0) {}

   public:
    void AddArg(Arg arg) { expr_->L.ep[arg_index_++] = arg.expr_; }
  };

  typedef ArgHandler<LogicalExpr> LogicalArgHandler;
  typedef ArgHandler<NumericExpr> NumericArgHandler;

  NumericArgHandler BeginVarArg(expr::Kind kind, int num_args) {
    return NumericArgHandler(
          MakeIterated(kind, ArrayRef<NumericExpr>(0, num_args)));
  }

  VarArgExpr EndVarArg(NumericArgHandler handler) {
    return Expr::Create<VarArgExpr>(handler.expr_);
  }

  VarArgExpr MakeVarArg(expr::Kind kind, ArrayRef<NumericExpr> args);

  NumericArgHandler BeginSum(int num_args) {
    return NumericArgHandler(
          MakeIterated(expr::SUM, ArrayRef<NumericExpr>(0, num_args)));
  }
  SumExpr EndSum(NumericArgHandler handler) {
    return Expr::Create<SumExpr>(handler.expr_);
  }

  SumExpr MakeSum(ArrayRef<NumericExpr> args) {
    return MakeIterated<SumExpr>(expr::SUM, args);
  }

  LogicalArgHandler BeginCount(int num_args) {
    return LogicalArgHandler(
          MakeIterated(expr::COUNT, ArrayRef<LogicalExpr>(0, num_args)));
  }
  CountExpr EndCount(LogicalArgHandler handler) {
    return Expr::Create<CountExpr>(handler.expr_);
  }

  CountExpr MakeCount(ArrayRef<LogicalExpr> args) {
    return MakeIterated<CountExpr>(expr::COUNT, args);
  }

  typedef NumericArgHandler NumberOfArgHandler;

  NumberOfArgHandler BeginNumberOf(NumericExpr value, int num_args) {
    NumericArgHandler handler(
          MakeIterated(expr::NUMBEROF, ArrayRef<NumericExpr>(0, num_args)));
    handler.AddArg(value);
    return handler;
  }
  NumberOfExpr EndNumberOf(NumberOfArgHandler handler) {
    return Expr::Create<NumberOfExpr>(handler.expr_);
  }

  NumberOfExpr MakeNumberOf(ArrayRef<NumericExpr> args) {
    assert(args.size() >= 1);
    return MakeIterated<NumberOfExpr>(expr::NUMBEROF, args);
  }

  LogicalConstant MakeLogicalConstant(bool value) {
    return Expr::Create<LogicalConstant>(MakeConstant(value));
  }

  NotExpr MakeNot(LogicalExpr arg) {
    return MakeUnary<NotExpr>(expr::NOT, arg);
  }

  BinaryLogicalExpr MakeBinaryLogical(
      expr::Kind kind, LogicalExpr lhs, LogicalExpr rhs) {
    return MakeBinary<BinaryLogicalExpr>(kind, lhs, rhs, "binary logical");
  }

  RelationalExpr MakeRelational(
      expr::Kind kind, NumericExpr lhs, NumericExpr rhs) {
    return MakeBinary<RelationalExpr>(kind, lhs, rhs, "relational");
  }

  LogicalCountExpr MakeLogicalCount(
      expr::Kind kind, NumericExpr lhs, CountExpr rhs) {
    return MakeBinary<LogicalCountExpr>(kind, lhs, rhs, "logical count");
  }

  ImplicationExpr MakeImplication(
      LogicalExpr condition, LogicalExpr true_expr, LogicalExpr false_expr) {
    return Expr::Create<ImplicationExpr>(
        MakeIf(expr::IMPLICATION, condition, true_expr, false_expr));
  }

  LogicalArgHandler BeginIteratedLogical(expr::Kind kind, int num_args) {
    return LogicalArgHandler(
          MakeIterated(kind, ArrayRef<LogicalExpr>(0, num_args)));
  }
  IteratedLogicalExpr EndIteratedLogical(LogicalArgHandler handler) {
    return Expr::Create<IteratedLogicalExpr>(handler.expr_);
  }

  IteratedLogicalExpr MakeIteratedLogical(
      expr::Kind kind, ArrayRef<LogicalExpr> args);

  typedef NumericArgHandler AllDiffArgHandler;

  AllDiffArgHandler BeginAllDiff(int num_args) {
    return NumericArgHandler(
          MakeIterated(expr::ALLDIFF, ArrayRef<NumericExpr>(0, num_args)));
  }
  AllDiffExpr EndAllDiff(AllDiffArgHandler handler) {
    return Expr::Create<AllDiffExpr>(handler.expr_);
  }

  AllDiffExpr MakeAllDiff(ArrayRef<NumericExpr> args) {
    return MakeIterated<AllDiffExpr>(expr::ALLDIFF, args);
  }

  // Constructs a StringLiteral object.
  // value: string value which may not be null-terminated.
  StringLiteral MakeStringLiteral(fmt::StringRef value);
};
}
}  // namespace mp

#endif  // MP_ASLBUILDER_H_
