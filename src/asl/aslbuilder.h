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

#include <cstring>

#include "mp/arrayref.h"
#include "mp/format.h"
#include "mp/safeint.h"
#include "aslexpr.h"
#include "aslproblem.h"

struct Static;

namespace mp {
struct NLHeader;

namespace asl {

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
  int func_index_;
  int expr_index_;

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

  template <typename ExprType>
  static void CheckKind(expr::Kind kind, const char *expr_name) {
    if (!internal::Is<ExprType>(kind))
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
  class SuffixHandler {
   private:
    T *values_;

#ifndef NDEBUG
    int num_items_;
    void set_num_items(int num_items) { num_items_ = num_items; }
#else
    void set_num_items(int num_items) { MP_UNUSED(num_items); }
#endif

   public:
    explicit SuffixHandler(T *values = 0, int num_items = 0) : values_(values) {
      set_num_items(num_items);
    }

    void SetValue(int index, T value) {
      assert(0 <= index && index < num_items_);
      if (values_)
        values_[index] = value;
    }
  };

  struct SuffixInfo {
    SufDesc *desc;
    int nx;
    int nx1;
    SuffixInfo() : desc(0), nx(0), nx1(0) {}
    SuffixInfo(SufDesc *d, int nx, int nx1) : desc(d), nx(nx), nx1(nx1) {}
  };

  SuffixInfo AddSuffix(fmt::StringRef name, int kind, int num_values);

  template <typename T>
  void AllocateSuffixValues(T *&values, int num_values, int nx, int nx1) {
    if (!values)
      values = Allocate<T>(nx1 * sizeof(T));
    if (num_values < nx)
      std::memset(values, 0, nx * sizeof(T));
    if (nx < nx1)
      std::memset(values + nx, 0, (nx1 - nx) * sizeof(T));
  }

  // Sets objective or constraint expression; adapted from co_read.
  void SetObjOrCon(int index, cde *d, int *cexp1_end, ::expr *e, int **z);

  ::expr *MakeConstant(double value);

  ::expr *DoMakeUnary(expr::Kind kind, Expr arg);

  template <typename ExprType>
  ExprType MakeUnary(expr::Kind kind, Expr arg) {
    return Expr::Create<ExprType>(DoMakeUnary(kind, arg));
  }

  ::expr *DoMakeBinary(expr::Kind kind, Expr lhs, Expr rhs);

  template <typename ExprType>
  ExprType MakeBinary(expr::Kind kind, Expr lhs, Expr rhs, const char *name) {
    CheckKind<ExprType>(kind, name);
    return Expr::Create<ExprType>(DoMakeBinary(kind, lhs, rhs));
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
  int num_algebraic_cons() const { return asl_->i.n_con_; }

  class CallExprBuilder {
   private:
    expr_f *expr_;
    int arg_index_, num_args_, num_constants_, num_symbolic_args_, num_ifsyms_;

    friend class ASLBuilder;

    explicit CallExprBuilder(expr_f *expr, int num_args)
      : expr_(expr), arg_index_(0), num_args_(num_args),
        num_constants_(0), num_symbolic_args_(0), num_ifsyms_(0) {}

    void DoAddArg(Expr arg) {
      expr_->args[arg_index_] = arg.impl_;
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
  CallExprBuilder DoBeginCall(Function f, int num_args);

  void Init(ASL *asl);

 public:
  typedef asl::Function Function;
  typedef asl::Expr Expr;
  typedef asl::NumericExpr NumericExpr;
  typedef asl::LogicalExpr LogicalExpr;
  typedef asl::Variable Variable;
  typedef asl::CountExpr CountExpr;

  explicit ASLBuilder(ASL *asl = 0) { Init(asl); }
  explicit ASLBuilder(ASLProblem::Proxy proxy) {
    Init(proxy.asl_);
    proxy.asl_ = 0;
    own_asl_ = true;
    flags_ = proxy.flags_ | ASL_STANDARD_OPCODES | ASL_allow_CLP;
  }
  ~ASLBuilder();

  // Returns a built problem via proxy. No builder methods other than
  // the destructor should be called after calling this method.
  ASLProblem::Proxy GetProblem() {
    if (!own_asl_)
      throw Error("ASL problem is not transferable");
    ASL *asl = asl_;
    own_asl_ = false;
    return ASLProblem::Proxy(asl);
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
  class BasicLinearExprBuilder {
   protected:
    ASLBuilder *builder_;
    Grad **term_;

   public:
    BasicLinearExprBuilder(ASLBuilder *b, Grad **term)
      : builder_(b), term_(term) {}

    void AddTerm(int var_index, double coef) {
      Grad *og = builder_->Allocate<Grad>();
      *term_ = og;
      term_ = &og->next;
      og->next = 0;
      og->varno = var_index;
      og->coef = coef;
    }
  };

  typedef BasicLinearExprBuilder<ograd> LinearObjBuilder;

  class LinearConBuilder : private BasicLinearExprBuilder<cgrad> {
   private:
    int con_index_;
    double *a_vals_;
    int *a_rownos_;
    int *a_colstarts_;

   public:
    LinearConBuilder(ASLBuilder *b, int con_index);

    void AddTerm(int var_index, double coef) {
      if (!a_vals_)
        return BasicLinearExprBuilder<cgrad>::AddTerm(var_index, coef);
      std::size_t elt_index = a_colstarts_[var_index]++;
      a_vals_[elt_index] = coef;
      a_rownos_[elt_index] = con_index_;
    }
  };

  class LinearExprBuilder {
   private:
    linpart *linpart_;
    int num_terms_;
    int index_;

   public:
    LinearExprBuilder(linpart *lp, int num_terms)
      : linpart_(lp), num_terms_(num_terms), index_(0) {}

    linpart *get() const { return linpart_; }
    int num_terms() const { return num_terms_; }

    void AddTerm(int var_index, double coef) {
      linpart_[index_].v.i = var_index;
      linpart_[index_].fac = coef;
      ++index_;
    }
  };

  // Adds and objective.
  LinearObjBuilder AddObj(obj::Type type, NumericExpr expr, int);

  // Adds an algebraic constraint.
  LinearConBuilder AddCon(double lb, double ub, NumericExpr expr, int);

  // Adds a logical constraint.
  void AddCon(LogicalExpr expr);

  LinearExprBuilder BeginCommonExpr(int num_terms);

  NumericExpr EndCommonExpr(LinearExprBuilder builder,
                            NumericExpr expr, int position);

  void SetComplement(int, int, int) {
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

  Function RegisterFunction(const char *name, ufunc f, int num_args,
                            func::Type type = func::NUMERIC, void *info = 0);

  // Adds a function.
  // If the function with the specified name doesn't exist and the flag
  // ASL_allow_missing_funcs is not set, SetFunction throws ASLError.
  Function AddFunction(fmt::StringRef name, int num_args,
                       func::Type type = func::NUMERIC);

  Function FindFunction(fmt::StringRef name) const {
    return Function(func_lookup_ASL(asl_, name.c_str(), 0));
  }

  typedef SuffixHandler<int> IntSuffixHandler;

  IntSuffixHandler AddIntSuffix(
      fmt::StringRef name, int kind, int num_values) {
    SuffixInfo info = AddSuffix(name, kind, num_values);
    AllocateSuffixValues(info.desc->u.i, num_values, info.nx, info.nx1);
    return IntSuffixHandler(info.desc->u.i, info.nx1);
  }

  typedef SuffixHandler<double> DblSuffixHandler;

  DblSuffixHandler AddDblSuffix(
      fmt::StringRef name, int kind, int num_values) {
    SuffixInfo info = AddSuffix(name, kind, num_values);
    AllocateSuffixValues(info.desc->u.r, num_values, info.nx, info.nx1);
    return DblSuffixHandler(info.desc->u.r, info.nx1);
  }

  // The Make* methods construct expression objects. These objects are
  // local to the currently built ASL problem and shouldn't be used with
  // other problems. The expression objects are not accessible via the
  // problem API until they are added as a part of objective or constraint
  // expression. For this reason the methods are called Make* rather than Add*.

  NumericConstant MakeNumericConstant(double value) {
    return Expr::Create<NumericConstant>(MakeConstant(value));
  }

  Variable MakeVariable(int var_index);
  Variable MakeCommonExpr(int index) {
    return MakeVariable(index);
  }

  UnaryExpr MakeUnary(expr::Kind kind, NumericExpr arg);

  BinaryExpr MakeBinary(expr::Kind kind, NumericExpr lhs, NumericExpr rhs) {
    return MakeBinary<BinaryExpr>(kind, lhs, rhs, "binary");
  }

  IfExpr MakeIf(LogicalExpr condition,
      NumericExpr true_expr, NumericExpr false_expr) {
    return Expr::Create<IfExpr>(
        MakeIf(expr::IF, condition, true_expr, false_expr));
  }

  class PLTermBuilder {
   private:
    ::expr *expr_;
    double *data_;

    friend class ASLBuilder;

    explicit PLTermBuilder(::expr *e) : expr_(e), data_(e->L.p->bs) {}

   public:

    void AddSlope(double slope) { *data_++ = slope; }
    void AddBreakpoint(double breakpoint) { *data_++ = breakpoint; }
  };

  PLTermBuilder BeginPLTerm(int num_breakpoints);
  PiecewiseLinearExpr EndPLTerm(PLTermBuilder b, NumericExpr var) {
    b.expr_->R.e = var.impl_;
    return Expr::Create<PiecewiseLinearExpr>(b.expr_);
  }

  PiecewiseLinearExpr MakePiecewiseLinear(int num_breakpoints,
      const double *breakpoints, const double *slopes, Variable var);

  CallExprBuilder BeginCall(Function f, int num_args) {
    if (!f)
      throw Error("undefined function");
    return DoBeginCall(f, num_args);
  }
  CallExpr EndCall(CallExprBuilder b);

  CallExpr MakeCall(Function f, ArrayRef<Expr> args) {
    int num_args = SafeInt<int>(args.size()).value();
    CallExprBuilder builder = DoBeginCall(f, num_args);
    for (int i = 0; i < num_args; ++i) {
      if (NumericExpr num = asl::Cast<NumericExpr>(args[i]))
        builder.AddArg(num);
      else
        builder.AddArg(asl::Cast<StringLiteral>(args[i]));
    }
    return EndCall(builder);
  }

  template <typename Arg>
  class ExprBuilder {
   private:
    ::expr *expr_;
    int arg_index_;

    friend class ASLBuilder;

    ExprBuilder(::expr *e) : expr_(e), arg_index_(0) {}

   public:
    void AddArg(Arg arg) { expr_->L.ep[arg_index_++] = arg.impl_; }
  };

  typedef ExprBuilder<NumericExpr> NumericExprBuilder;
  typedef ExprBuilder<LogicalExpr> LogicalExprBuilder;

  class VarArgExprBuilder {
   private:
    expr_va *expr_;
    int arg_index_;

    friend class ASLBuilder;

    VarArgExprBuilder(expr_va *e) : expr_(e), arg_index_(0) {}

   public:
    void AddArg(NumericExpr arg) { expr_->L.d[arg_index_++].e = arg.impl_; }
  };

  VarArgExprBuilder BeginVarArg(expr::Kind kind, int num_args);

  VarArgExpr EndVarArg(VarArgExprBuilder builder) {
    return Expr::Create<VarArgExpr>(reinterpret_cast< ::expr*>(builder.expr_));
  }

  VarArgExpr MakeVarArg(expr::Kind kind, ArrayRef<NumericExpr> args);

  NumericExprBuilder BeginSum(int num_args) {
    return NumericExprBuilder(
          MakeIterated(expr::SUM, ArrayRef<NumericExpr>(0, num_args)));
  }
  SumExpr EndSum(NumericExprBuilder builder) {
    return Expr::Create<SumExpr>(builder.expr_);
  }

  SumExpr MakeSum(ArrayRef<NumericExpr> args) {
    return MakeIterated<SumExpr>(expr::SUM, args);
  }

  LogicalExprBuilder BeginCount(int num_args) {
    return LogicalExprBuilder(
          MakeIterated(expr::COUNT, ArrayRef<LogicalExpr>(0, num_args)));
  }
  CountExpr EndCount(LogicalExprBuilder builder) {
    return Expr::Create<CountExpr>(builder.expr_);
  }

  CountExpr MakeCount(ArrayRef<LogicalExpr> args) {
    return MakeIterated<CountExpr>(expr::COUNT, args);
  }

  typedef NumericExprBuilder NumberOfExprBuilder;

  NumericExprBuilder BeginNumberOf(NumericExpr value, int num_args) {
    NumericExprBuilder builder(
          MakeIterated(expr::NUMBEROF, ArrayRef<NumericExpr>(0, num_args + 1)));
    builder.AddArg(value);
    return builder;
  }
  NumberOfExpr EndNumberOf(NumberOfExprBuilder builder) {
    return Expr::Create<NumberOfExpr>(builder.expr_);
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

  LogicalExprBuilder BeginIteratedLogical(expr::Kind kind, int num_args) {
    return LogicalExprBuilder(
          MakeIterated(kind, ArrayRef<LogicalExpr>(0, num_args)));
  }
  IteratedLogicalExpr EndIteratedLogical(LogicalExprBuilder builder) {
    return Expr::Create<IteratedLogicalExpr>(builder.expr_);
  }

  IteratedLogicalExpr MakeIteratedLogical(
      expr::Kind kind, ArrayRef<LogicalExpr> args);

  typedef NumericExprBuilder PairwiseExprBuilder;

  PairwiseExprBuilder BeginPairwise(expr::Kind kind, int num_args) {
    return PairwiseExprBuilder(
          MakeIterated(kind, ArrayRef<NumericExpr>(0, num_args)));
  }
  PairwiseExpr EndPairwise(PairwiseExprBuilder builder) {
    return Expr::Create<PairwiseExpr>(builder.expr_);
  }

  PairwiseExpr MakeAllDiff(ArrayRef<NumericExpr> args) {
    return MakeIterated<PairwiseExpr>(expr::ALLDIFF, args);
  }

  // Constructs a StringLiteral object.
  // value: string value which may not be null-terminated.
  StringLiteral MakeStringLiteral(fmt::StringRef value);
};
}
}  // namespace asl
}  // namespace mp

#endif  // MP_ASLBUILDER_H_
