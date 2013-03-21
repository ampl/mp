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
class JVM : private Noncopyable, public Env {
 private:
  JavaVM *jvm_;

 public:
  JVM();
  ~JVM() {
    jvm_->DestroyJavaVM();
  }
};

class ClassBase {
 protected:
  jclass class_;
  jmethodID ctor_;

  void True() const {}
  typedef void (ClassBase::*SafeBool)() const;

 public:
  ClassBase() : class_(), ctor_() {}

  operator SafeBool() const { return class_ ? &ClassBase::True : 0; }

  jclass get() const { return class_; }

  jobject NewObject(Env env, ...);
};

// A reference to a Java class and one of its constructor.
template <typename Info>
class Class : public ClassBase {
 public:
  void Init(Env env) {
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
CLASS_INFO(Sum, "JaCoP/constraints/Sum", "([LJaCoP/core/IntVar;)V")
CLASS_INFO(SumWeight, "JaCoP/constraints/SumWeight",
    "([LJaCoP/core/IntVar;[ILJaCoP/core/IntVar;)V")
CLASS_INFO(XplusYeqZ, "JaCoP/constraints/XplusYeqZ",
    "(LJaCoP/core/IntVar;LJaCoP/core/IntVar;LJaCoP/core/IntVar;)V")
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
CLASS_INFO(AbsXeqY, "JaCoP/constraints/AbsXeqY",
    "(LJaCoP/core/IntVar;LJaCoP/core/IntVar;)V")
CLASS_INFO(Min, "JaCoP/constraints/Min",
    "([LJaCoP/core/IntVar;LJaCoP/core/IntVar;)V")
CLASS_INFO(Max, "JaCoP/constraints/Max",
    "([LJaCoP/core/IntVar;LJaCoP/core/IntVar;)V")
CLASS_INFO(IfThenElse, "JaCoP/constraints/IfThenElse",
    "(LJaCoP/constraints/PrimitiveConstraint;"
    "LJaCoP/constraints/PrimitiveConstraint;"
    "LJaCoP/constraints/PrimitiveConstraint;)V")

class NLToJaCoPConverter :
   private ExprVisitor<NLToJaCoPConverter, jobject, jobject> {
 private:
  JVM &jvm_;
  jobject store_;
  jmethodID impose_;
  jobjectArray var_array_;
  std::vector<jobject> vars_;
  Class<IntVar> var_class_;
  Class<SumWeight> sum_class_;
  Class<SumWeight> sum_weight_class_;
  Class<XplusYeqZ> plus_class_;
  Class<XmulYeqZ> mul_class_;
  Class<XmulCeqZ> mul_const_class_;
  Class<XdivYeqZ> div_class_;
  Class<XmodYeqZ> mod_class_;
  Class<XexpYeqZ> exp_class_;
  Class<XeqY> eq_class_;
  Class<XeqC> eq_const_class_;
  Class<AbsXeqY> abs_class_;
  Class<Min> min_class_;
  Class<Max> max_class_;
  Class<IfThenElse> if_class_;
  jint min_int_;
  jint max_int_;

  friend class ExprVisitor<NLToJaCoPConverter, jobject, jobject>;

  static jint CastToInt(double value) {
    jint int_value = static_cast<int>(value);
    if (int_value != value) {
      throw Error(str(
          fmt::Format("value {} can't be represented as int") << value));
    }
    return int_value;
  }

  jobject CreateVar() {
    return var_class_.NewObject(jvm_, store_, min_int_, max_int_);
  }

  // Creates a variable equal to a constant value.
  jobject CreateConst(int value) {
    if (!eq_const_class_)
      eq_const_class_.Init(jvm_);
    jobject result_var = CreateVar();
    Impose(eq_const_class_.NewObject(jvm_, result_var, value));
    return result_var;
  }

  // Creates a constraint with a binary expression.
  template <typename ClassT, typename Arg>
  jobject CreateCon(ClassT &cls, Arg arg) {
    if (!cls)
      cls.Init(jvm_);
    jobject result_var = CreateVar();
    Impose(cls.NewObject(jvm_, arg, result_var));
    return result_var;
  }

  // Creates a constraint with a binary expression.
  template <typename ClassT, typename LHS, typename RHS>
  jobject CreateCon(ClassT &cls, LHS lhs, RHS rhs) {
    if (!cls)
      cls.Init(jvm_);
    jobject result_var = CreateVar();
    Impose(cls.NewObject(jvm_, lhs, rhs, result_var));
    return result_var;
  }

  jobject CreateMinus(jobject lhs, jobject rhs) {
    return CreateCon(plus_class_, lhs, CreateCon(mul_const_class_, rhs, -1));
  }

  jobjectArray ConvertVarArgs(VarArgExpr e) {
    jobjectArray args = jvm_.NewObjectArray(
        std::distance(e.begin(), e.end()), var_class_.get(), 0);
    int index = 0;
    for (VarArgExpr::iterator i = e.begin(), end = e.end(); i != end; ++i)
      jvm_.SetObjectArrayElement(args, index++, Visit(*i));
    return args;
  }

  void Impose(jobject constraint) {
    jvm_.CallVoidMethod(store_, impose_, constraint);
  }

  // TODO
  //jobject Convert(Gecode::BoolOpType op, IteratedLogicalExpr e);

  static void RequireNonzeroConstRHS(
      BinaryExpr e, const std::string &func_name);

  template<typename Term>
  void ConvertExpr(LinearExpr<Term> linear,
      NumericExpr nonlinear, jobject result_var);

  // The methods below perform conversion of AMPL NL expressions into
  // equivalent JaCoP expressions. JaCoP doesn't support the following
  // expressions/functions:
  // * division other than integer one
  // * trigonometric functions
  // * log, log10, exp, sqrt

  jobject VisitPlus(BinaryExpr e) {
    return CreateCon(plus_class_, Visit(e.lhs()), Visit(e.rhs()));
  }

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
    return CreateCon(min_class_, ConvertVarArgs(e));
  }

  jobject VisitMax(VarArgExpr e) {
    return CreateCon(max_class_, ConvertVarArgs(e));
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

  // TODO
  //jobject VisitNumberOf(NumberOfExpr e);

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

  // TODO
  /*jobject VisitOr(BinaryLogicalExpr e) {
    return Visit(e.lhs()) || Visit(e.rhs());
  }

  jobject VisitAnd(BinaryLogicalExpr e) {
    return Visit(e.lhs()) && Visit(e.rhs());
  }

  jobject VisitLess(RelationalExpr e) {
    return Visit(e.lhs()) < Visit(e.rhs());
  }

  jobject VisitLessEqual(RelationalExpr e) {
    return Visit(e.lhs()) <= Visit(e.rhs());
  }

  jobject VisitEqual(RelationalExpr e) {
    return Visit(e.lhs()) == Visit(e.rhs());
  }

  jobject VisitGreaterEqual(RelationalExpr e) {
    return Visit(e.lhs()) >= Visit(e.rhs());
  }

  jobject VisitGreater(RelationalExpr e) {
    return Visit(e.lhs()) > Visit(e.rhs());
  }

  jobject VisitNotEqual(RelationalExpr e) {
    return Visit(e.lhs()) != Visit(e.rhs());
  }

  jobject VisitNot(NotExpr e) {
    return !Visit(e.arg());
  }

  jobject VisitAtLeast(LogicalCountExpr e) {
    return Visit(e.value()) <= Visit(e.count());
  }

  jobject VisitAtMost(LogicalCountExpr e) {
    return Visit(e.value()) >= Visit(e.count());
  }

  jobject VisitExactly(LogicalCountExpr e) {
    return Visit(e.value()) == Visit(e.count());
  }

  jobject VisitNotAtLeast(LogicalCountExpr e) {
    return Visit(e.value()) > Visit(e.count());
  }

  jobject VisitNotAtMost(LogicalCountExpr e) {
    return Visit(e.value()) < Visit(e.count());
  }

  jobject VisitNotExactly(LogicalCountExpr e) {
    return Visit(e.value()) != Visit(e.count());
  }

  jobject VisitForAll(IteratedLogicalExpr e) {
    return Convert(Gecode::BOT_AND, e);
  }

  jobject VisitExists(IteratedLogicalExpr e) {
    return Convert(Gecode::BOT_OR, e);
  }

  jobject VisitImplication(ImplicationExpr e);

  jobject VisitIff(BinaryLogicalExpr e) {
    return Visit(e.lhs()) == Visit(e.rhs());
  }

  jobject VisitAllDiff(AllDiffExpr e);*/

  jobject VisitLogicalConstant(LogicalConstant c) {
    return CreateConst(c.value());
  }

 public:
  explicit NLToJaCoPConverter(JVM &jvm);

  // TODO
  //jobject ConvertFullExpr(NumericExpr e) { return Visit(e); }
  //jobject ConvertFullExpr(LogicalExpr e, bool post = true);
  void Convert(const Problem &p);

  jobject store() const { return store_; }
  jobjectArray var_array() const { return var_array_; }
  const std::vector<jobject> &vars() const { return vars_; }
  jclass var_class() const { return var_class_.get(); }
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
  JVM jvm_;
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
  std::size_t memory_limit_;

  void EnableOutput(const char *name, int value);
  void SetOutputFrequency(const char *name, int value);

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

  /*Gecode::IntVarBranch var_branching() const { return var_branching_; }
  Gecode::IntValBranch val_branching() const { return val_branching_; }
  const Gecode::Search::Options &options() const { return options_; }*/
};
}

#endif // AMPL_SOLVERS_JACOP_H
