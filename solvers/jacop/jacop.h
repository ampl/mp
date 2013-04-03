/*
 AMPL solver interface to JaCoP.

 Copyright (C) 2012 AMPL Optimization LLC

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization LLC disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Victor Zverovich
 */

#ifndef AMPL_SOLVERS_JACOP_H
#define AMPL_SOLVERS_JACOP_H

#include "solvers/util/solver.h"

#include <jni.h>
#include <cstdarg>

namespace ampl {

// Java Native Interface environment.
class Env {
 private:
  JNIEnv *env_;

  class String : Noncopyable {
   private:
    JNIEnv *env_;
    jstring str_;
    const char *utf_chars_;

   public:
    String(JNIEnv *env, jstring s) :
        env_(env), str_(s), utf_chars_(env->GetStringUTFChars(s, 0)) {
    }

    ~String() {
      env_->ReleaseStringUTFChars(str_, utf_chars_);
    }

    const char *c_str() const {
      return utf_chars_;
    }
  };

  void Throw(jthrowable exception, const char *method_name);

  // Checks the result of a method call returning pointer and throws an
  // exception in the case of a failure.
  template<typename T>
  T Check(T result, const char *method_name);

  void Check(const char *method_name) {
    if (jthrowable exception = env_->ExceptionOccurred())
      Throw(exception, method_name);
  }

 public:
  explicit Env(JNIEnv *env = 0) : env_(env) {}

  jclass FindClass(const char *name) {
    return Check(env_->FindClass(name), "FindClass");
  }

  jmethodID GetMethod(jclass cls, const char *name, const char *sig) {
    return Check(env_->GetMethodID(cls, name, sig), "GetMethodID");
  }

  jfieldID GetStaticFieldID(jclass cls, const char *name, const char *sig) {
    return Check(env_->GetStaticFieldID(cls, name, sig), "GetStaticFieldID");
  }

  jint GetStaticIntField(jclass cls, jfieldID field) {
    jint value = env_->GetStaticIntField(cls, field);
    Check("GetStaticIntField");
    return value;
  }

  jobject NewObject(jclass cls, jmethodID ctor, ...);

  jobject NewObjectV(jclass cls, jmethodID ctor, std::va_list args) {
    return Check(env_->NewObjectV(cls, ctor, args), "NewObjectV");
  }

  jobject NewObject(const char *class_name, const char *ctor_sig, ...);

  void CallVoidMethod(jobject obj, jmethodID method, ...);
  jboolean CallBooleanMethod(jobject obj, jmethodID method, ...);
  jint CallIntMethod(jobject obj, jmethodID method, ...);

  jobjectArray NewObjectArray(jsize length,
      jclass elementClass, jobject initialElement) {
    return Check(env_->NewObjectArray(length, elementClass, initialElement),
        "NewObjectArray");
  }

  void SetObjectArrayElement(jobjectArray array, jsize index, jobject value) {
    env_->SetObjectArrayElement(array, index, value);
    Check("SetObjectArrayElement");
  }

  jintArray NewIntArray(jsize length) {
    return Check(env_->NewIntArray(length), "NewIntArray");
  }

  void SetIntArrayRegion(jintArray array,
      jsize start, jsize length, const int *values) {
    env_->SetIntArrayRegion(array, start, length, values);
    Check("SetIntArrayRegion");
  }
};

template <typename T>
T Env::Check(T result, const char *method_name) {
  if (!result) {
    jthrowable exception = env_->ExceptionOccurred();
    if (!exception)
      throw Error(std::string(method_name) + " failed");
    Throw(exception, method_name);
  }
  return result;
}

// Java Virtual Machine.
// A process can have at most one JVM, therefore there is a single static JVM
// object which is initialized when JVM::env() is called for the first time.
class JVM {
 private:
  JavaVM *jvm_;
  Env env_;
  static JVM instance_;

  JVM() : jvm_() {}
  ~JVM();

 public:
  static Env env();
};

class ClassBase {
 protected:
  jclass class_;
  jmethodID ctor_;

  void True() const {}
  typedef void (ClassBase::*SafeBool)() const;

  virtual void DoInit(Env env) = 0;

 public:
  ClassBase() : class_(), ctor_() {}
  virtual ~ClassBase();

  void Init(Env env) {
    if (!class_)
      DoInit(env);
  }

  operator SafeBool() const { return class_ ? &ClassBase::True : 0; }

  jclass get() const { return class_; }

  jobject NewObject(Env env, ...);
};

// A reference to a Java class and one of its constructor.
template <typename Info>
class Class : public ClassBase {
 protected:
  void DoInit(Env env) {
    class_ = env.FindClass(Info::name());
    ctor_ = env.GetMethod(class_, "<init>", Info::ctor_sig());
  }
};

#define CLASS_INFO(class_name, jvm_class_name, jvm_ctor_sig) \
struct class_name { \
  static const char *name() { return jvm_class_name; } \
  static const char *ctor_sig() { return jvm_ctor_sig; } \
};

CLASS_INFO(IntVar, "JaCoP/core/IntVar", "(LJaCoP/core/Store;II)V")
CLASS_INFO(Sum, "JaCoP/constraints/Sum",
    "([LJaCoP/core/IntVar;LJaCoP/core/IntVar;)V")
CLASS_INFO(SumWeight, "JaCoP/constraints/SumWeight",
    "([LJaCoP/core/IntVar;[ILJaCoP/core/IntVar;)V")
CLASS_INFO(XplusYeqZ, "JaCoP/constraints/XplusYeqZ",
    "(LJaCoP/core/IntVar;LJaCoP/core/IntVar;LJaCoP/core/IntVar;)V")
CLASS_INFO(XplusCeqZ, "JaCoP/constraints/XplusCeqZ",
    "(LJaCoP/core/IntVar;ILJaCoP/core/IntVar;)V")
CLASS_INFO(XmulYeqZ, "JaCoP/constraints/XmulYeqZ",
    "(LJaCoP/core/IntVar;LJaCoP/core/IntVar;LJaCoP/core/IntVar;)V")
CLASS_INFO(XmulCeqZ, "JaCoP/constraints/XmulCeqZ",
    "(LJaCoP/core/IntVar;ILJaCoP/core/IntVar;)V")
CLASS_INFO(XdivYeqZ, "JaCoP/constraints/XdivYeqZ",
    "(LJaCoP/core/IntVar;LJaCoP/core/IntVar;LJaCoP/core/IntVar;)V")
CLASS_INFO(XmodYeqZ, "JaCoP/constraints/XmodYeqZ",
    "(LJaCoP/core/IntVar;LJaCoP/core/IntVar;LJaCoP/core/IntVar;)V")
CLASS_INFO(XexpYeqZ, "JaCoP/constraints/XexpYeqZ",
    "(LJaCoP/core/IntVar;LJaCoP/core/IntVar;LJaCoP/core/IntVar;)V")
CLASS_INFO(XeqY, "JaCoP/constraints/XeqY",
    "(LJaCoP/core/IntVar;LJaCoP/core/IntVar;)V")
CLASS_INFO(XeqC, "JaCoP/constraints/XeqC", "(LJaCoP/core/IntVar;I)V")
CLASS_INFO(XneqY, "JaCoP/constraints/XneqY",
    "(LJaCoP/core/IntVar;LJaCoP/core/IntVar;)V")
CLASS_INFO(XltY, "JaCoP/constraints/XltY",
    "(LJaCoP/core/IntVar;LJaCoP/core/IntVar;)V")
CLASS_INFO(XlteqY, "JaCoP/constraints/XlteqY",
    "(LJaCoP/core/IntVar;LJaCoP/core/IntVar;)V")
CLASS_INFO(XgtY, "JaCoP/constraints/XgtY",
    "(LJaCoP/core/IntVar;LJaCoP/core/IntVar;)V")
CLASS_INFO(XgteqY, "JaCoP/constraints/XgteqY",
    "(LJaCoP/core/IntVar;LJaCoP/core/IntVar;)V")
CLASS_INFO(XgteqC, "JaCoP/constraints/XgteqC",
    "(LJaCoP/core/IntVar;LJaCoP/core/IntVar;)V")
CLASS_INFO(AbsXeqY, "JaCoP/constraints/AbsXeqY",
    "(LJaCoP/core/IntVar;LJaCoP/core/IntVar;)V")
CLASS_INFO(Min, "JaCoP/constraints/Min",
    "([LJaCoP/core/IntVar;LJaCoP/core/IntVar;)V")
CLASS_INFO(Max, "JaCoP/constraints/Max",
    "([LJaCoP/core/IntVar;LJaCoP/core/IntVar;)V")
CLASS_INFO(Count, "JaCoP/constraints/Count",
    "([LJaCoP/core/IntVar;LJaCoP/core/IntVar;I)V")
CLASS_INFO(IfThen, "JaCoP/constraints/IfThen",
    "(LJaCoP/constraints/PrimitiveConstraint;"
    "LJaCoP/constraints/PrimitiveConstraint;)V")
CLASS_INFO(IfThenElse, "JaCoP/constraints/IfThenElse",
    "(LJaCoP/constraints/PrimitiveConstraint;"
    "LJaCoP/constraints/PrimitiveConstraint;"
    "LJaCoP/constraints/PrimitiveConstraint;)V")
CLASS_INFO(Or, "JaCoP/constraints/Or",
    "(LJaCoP/constraints/PrimitiveConstraint;"
    "LJaCoP/constraints/PrimitiveConstraint;)V")
CLASS_INFO(And, "JaCoP/constraints/And",
    "(LJaCoP/constraints/PrimitiveConstraint;"
    "LJaCoP/constraints/PrimitiveConstraint;)V")
CLASS_INFO(Not, "JaCoP/constraints/Not",
    "(LJaCoP/constraints/PrimitiveConstraint;)V")
CLASS_INFO(Eq, "JaCoP/constraints/Eq",
    "(LJaCoP/constraints/PrimitiveConstraint;"
    "LJaCoP/constraints/PrimitiveConstraint;)V")
CLASS_INFO(Alldiff, "JaCoP/constraints/Alldiff", "([LJaCoP/core/IntVar;)V")
CLASS_INFO(DepthFirstSearch, "JaCoP/search/DepthFirstSearch", "()V")

// Converter of constraint programming problems from NL to JaCoP format.
class NLToJaCoPConverter :
   public ExprVisitor<NLToJaCoPConverter, jobject, jobject> {
 private:
  Env env_;
  jobject store_;
  jmethodID impose_;
  jobjectArray var_array_;
  std::vector<jobject> vars_;
  Class<IntVar> var_class_;
  Class<Sum> sum_class_;
  Class<SumWeight> sum_weight_class_;
  Class<XplusYeqZ> plus_class_;
  Class<XplusCeqZ> plus_const_class_;
  Class<XmulYeqZ> mul_class_;
  Class<XmulCeqZ> mul_const_class_;
  Class<XdivYeqZ> div_class_;
  Class<XmodYeqZ> mod_class_;
  Class<XexpYeqZ> exp_class_;
  Class<XeqY> eq_class_;
  Class<XeqC> eq_const_class_;
  Class<XltY> lt_class_;
  Class<XlteqY> le_class_;
  Class<XgtY> gt_class_;
  Class<XgteqY> ge_class_;
  Class<XneqY> ne_class_;
  Class<AbsXeqY> abs_class_;
  Class<Min> min_class_;
  Class<Max> max_class_;
  Class<Count> count_class_;
  Class<IfThen> if_class_;
  Class<IfThenElse> if_else_class_;
  Class<Or> or_class_;
  Class<And> and_class_;
  Class<Not> not_class_;
  Class<Eq> eq_con_class_;
  Class<Alldiff> alldiff_class_;
  jclass constraint_class_;
  jmethodID or_array_ctor_;
  jmethodID and_array_ctor_;
  jobject one_var_;
  jint min_int_;
  jint max_int_;

  jint CastToInt(double value) const {
    jint int_value = static_cast<int>(value);
    if (int_value != value) {
      throw Error(str(
          fmt::Format("value {} can't be represented as int") << value));
    }
    if (int_value < min_int_ || int_value > max_int_)
      throw Error(str(fmt::Format("value {} is out of bounds") << value));
    return int_value;
  }

  jobject CreateVar() {
    return var_class_.NewObject(env_, store_, min_int_, max_int_);
  }

  jobjectArray CreateVarArray(jsize length) {
    return env_.NewObjectArray(length, var_class_.get(), 0);
  }

  // Creates a variable equal to a constant value.
  jobject CreateConst(int value) {
    jobject result_var = CreateVar();
    Impose(eq_const_class_.NewObject(env_, result_var, value));
    return result_var;
  }

  // Creates a constraint with an argument and a variable to hold the result.
  jobject CreateCon(ClassBase &cls, jobject arg) {
    jobject result_var = CreateVar();
    Impose(cls.NewObject(env_, arg, result_var));
    return result_var;
  }

  // Creates a constraint with two arguments and a variable to hold the result.
  template <typename Arg1, typename Arg2>
  jobject CreateCon(ClassBase &cls, Arg1 arg1, Arg2 arg2) {
    jobject result_var = CreateVar();
    Impose(cls.NewObject(env_, arg1, arg2, result_var));
    return result_var;
  }

  jobject CreateMinus(jobject lhs, jobject rhs) {
    return CreateCon(plus_class_, lhs, CreateCon(mul_const_class_, rhs, -1));
  }

  jobject Convert(VarArgExpr e, ClassBase &cls) {
    jobjectArray args = CreateVarArray(std::distance(e.begin(), e.end()));
    int index = 0;
    for (VarArgExpr::iterator i = e.begin(), end = e.end(); i != end; ++i)
      env_.SetObjectArrayElement(args, index++, Visit(*i));
    return CreateCon(cls, args);
  }

  // Converts a binary logical expression to a JaCoP constraint of class cls.
  template <typename ExprT>
  jobject Convert(ExprT e, ClassBase &cls) {
    return cls.NewObject(env_, Visit(e.lhs()), Visit(e.rhs()));
  }

  // Converts a logical count expression.
  jobject Convert(LogicalCountExpr e, ClassBase &cls) {
    return cls.NewObject(env_, Visit(e.value()), VisitCount(e.count()));
  }

  // Converts an iterated logical expression.
  jobject Convert(IteratedLogicalExpr e, ClassBase &cls, jmethodID &ctor);

  void Impose(jobject constraint) {
    env_.CallVoidMethod(store_, impose_, constraint);
  }

  static void RequireNonzeroConstRHS(
      BinaryExpr e, const std::string &func_name);

  template<typename Term>
  void ConvertExpr(LinearExpr<Term> linear,
      NumericExpr nonlinear, jobject result_var);

 public:
  explicit NLToJaCoPConverter();

  // Converts a logical constraint.
  void ConvertLogicalCon(LogicalExpr e, bool post = true);

  void Convert(const Problem &p);

  jobject store() const { return store_; }
  jobjectArray var_array() const { return var_array_; }
  const std::vector<jobject> &vars() const { return vars_; }
  Class<IntVar> &var_class() { return var_class_; }

  // The methods below perform conversion of AMPL NL expressions into
  // equivalent JaCoP expressions. JaCoP doesn't support the following
  // expressions/functions:
  // * division other than integer one
  // * trigonometric functions
  // * log, log10, exp, sqrt

  jobject VisitPlus(BinaryExpr e);

  jobject VisitMinus(BinaryExpr e) {
    return CreateMinus(Visit(e.lhs()), Visit(e.rhs()));
  }

  jobject VisitMult(BinaryExpr e) {
    return CreateCon(mul_class_, Visit(e.lhs()), Visit(e.rhs()));
  }

  jobject VisitRem(BinaryExpr e) {
    return CreateCon(mod_class_, Visit(e.lhs()), Visit(e.rhs()));
  }

  jobject VisitPow(BinaryExpr e) {
    return CreateCon(exp_class_, Visit(e.lhs()), Visit(e.rhs()));
  }

  jobject VisitNumericLess(BinaryExpr e);

  jobject VisitMin(VarArgExpr e) {
    return Convert(e, min_class_);
  }

  jobject VisitMax(VarArgExpr e) {
    return Convert(e, max_class_);
  }

  jobject VisitFloor(UnaryExpr e) {
    // floor does nothing because JaCoP supports only integer expressions.
    return Visit(e.arg());
  }

  jobject VisitCeil(UnaryExpr e) {
    // ceil does nothing because JaCoP supports only integer expressions.
    return Visit(e.arg());
  }

  jobject VisitAbs(UnaryExpr e) {
    return CreateCon(abs_class_, Visit(e.arg()));
  }

  jobject VisitUnaryMinus(UnaryExpr e) {
    return CreateCon(mul_const_class_, Visit(e.arg()), -1);
  }

  jobject VisitIf(IfExpr e);

  jobject VisitSum(SumExpr e);

  jobject VisitIntDiv(BinaryExpr e) {
    return CreateCon(div_class_, Visit(e.lhs()), Visit(e.rhs()));
  }

  jobject VisitRound(BinaryExpr e) {
    // round does nothing because JaCoP supports only integer expressions.
    RequireNonzeroConstRHS(e, "round");
    return Visit(e.lhs());
  }

  jobject VisitTrunc(BinaryExpr e) {
    // trunc does nothing because JaCoP supports only integer expressions.
    RequireNonzeroConstRHS(e, "trunc");
    return Visit(e.lhs());
  }

  jobject VisitCount(CountExpr e);

  jobject VisitNumberOf(NumberOfExpr e);

  jobject VisitPow2(UnaryExpr e) {
    jobject arg = Visit(e.arg());
    return CreateCon(mul_class_, arg, arg);
  }

  jobject VisitNumericConstant(NumericConstant c) {
    return CreateConst(CastToInt(c.value()));
  }

  jobject VisitVariable(Variable v) {
    return vars_[v.index()];
  }

  jobject VisitOr(BinaryLogicalExpr e) {
    return Convert(e, or_class_);
  }

  jobject VisitAnd(BinaryLogicalExpr e) {
    return Convert(e, and_class_);
  }

  jobject VisitLess(RelationalExpr e) {
    return Convert(e, lt_class_);
  }

  jobject VisitLessEqual(RelationalExpr e) {
    return Convert(e, le_class_);
  }

  jobject VisitEqual(RelationalExpr e) {
    return Convert(e, eq_class_);
  }

  jobject VisitGreaterEqual(RelationalExpr e) {
    return Convert(e, ge_class_);
  }

  jobject VisitGreater(RelationalExpr e) {
    return Convert(e, gt_class_);
  }

  jobject VisitNotEqual(RelationalExpr e) {
    return Convert(e, ne_class_);
  }

  jobject VisitNot(NotExpr e) {
    return not_class_.NewObject(env_, Visit(e.arg()));
  }

  jobject VisitAtLeast(LogicalCountExpr e) {
    return Convert(e, le_class_);
  }

  jobject VisitAtMost(LogicalCountExpr e) {
    return Convert(e, ge_class_);
  }

  jobject VisitExactly(LogicalCountExpr e) {
    return Convert(e, eq_class_);
  }

  jobject VisitNotAtLeast(LogicalCountExpr e) {
    return Convert(e, gt_class_);
  }

  jobject VisitNotAtMost(LogicalCountExpr e) {
    return Convert(e, lt_class_);
  }

  jobject VisitNotExactly(LogicalCountExpr e) {
    return Convert(e, ne_class_);
  }

  jobject VisitForAll(IteratedLogicalExpr e) {
    return Convert(e, and_class_, and_array_ctor_);
  }

  jobject VisitExists(IteratedLogicalExpr e) {
    return Convert(e, or_class_, or_array_ctor_);
  }

  jobject VisitImplication(ImplicationExpr e);

  jobject VisitIff(BinaryLogicalExpr e) {
    return eq_con_class_.NewObject(env_, Visit(e.lhs()), Visit(e.rhs()));
  }

  jobject VisitAllDiff(AllDiffExpr) {
    throw UnsupportedExprError("nested 'alldiff'");
    return jobject();
  }

  jobject VisitLogicalConstant(LogicalConstant c) {
    if (!one_var_)
      one_var_ = var_class_.NewObject(env_, store_, 1, 1);
    return eq_const_class_.NewObject(env_, one_var_, c.value());
  }
};

// TODO
/*template <typename Value>
struct OptionValue {
  const char *name;
  Value value;
};

template <typename T>
struct OptionInfo {
  const OptionValue<T> *values;
  T &value;

  OptionInfo(const OptionValue<T> *values, T &value)
  : values(values), value(value) {}
};*/

// JaCoP solver.
class JaCoPSolver : public Solver<JaCoPSolver> {
 private:
  bool debug_;
  /*bool output_;
  double output_frequency_;
  unsigned output_count_;
  std::string header_;

  Gecode::IntVarBranch var_branching_;
  Gecode::IntValBranch val_branching_;
  Gecode::Search::Options options_;
  double time_limit_; // Time limit in seconds.
  unsigned long node_limit_;
  unsigned long fail_limit_;
  std::size_t memory_limit_;*/

  void SetBoolOption(const char *name, int value, bool *option);
  /*void SetOutputFrequency(const char *name, int value);

  template <typename T>
  void SetStrOption(const char *name, const char *value,
      const OptionInfo<T> &info);

  template <typename T, typename OptionT>
  void SetOption(const char *name, T value, OptionT *option);

  void SetDblOption(const char *, double value, double *option) {
    *option = value;
  }

  fmt::TempFormatter<fmt::Write> Output(fmt::StringRef format);
*/
 public:
  JaCoPSolver();

  // Run the solver.
  int Run(char **argv);

  double var_min() const {
    Env env = JVM::env();
    jclass domain_class = env.FindClass("JaCoP/core/IntDomain");
    return env.GetStaticIntField(
        domain_class, env.GetStaticFieldID(domain_class, "MinInt", "I"));
  }

  double var_max() const {
    Env env = JVM::env();
    jclass domain_class = env.FindClass("JaCoP/core/IntDomain");
    return env.GetStaticIntField(
        domain_class, env.GetStaticFieldID(domain_class, "MaxInt", "I"));
  }

  void Solve(Problem &p);

  /*Gecode::IntVarBranch var_branching() const { return var_branching_; }
  Gecode::IntValBranch val_branching() const { return val_branching_; }
  const Gecode::Search::Options &options() const { return options_; }*/
};
}

#endif // AMPL_SOLVERS_JACOP_H
