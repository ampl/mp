/*
 AMPL solver interface to JaCoP.

 Copyright (C) 2012 AMPL Optimization Inc

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

#ifndef MP_SOLVERS_JACOP_H_
#define MP_SOLVERS_JACOP_H_

#include <vector>

#include "mp/clock.h"
#include "asl/aslexpr-visitor.h"
#include "asl/aslproblem.h"
#include "asl/aslsolver.h"
#include "jacop/java.h"

namespace mp {

CLASS_INFO(IntVar, "org/jacop/core/IntVar", "(Lorg/jacop/core/Store;II)V")
CLASS_INFO(Sum, "org/jacop/constraints/Sum",
    "([Lorg/jacop/core/IntVar;Lorg/jacop/core/IntVar;)V")
CLASS_INFO(SumWeight, "org/jacop/constraints/SumWeight",
    "([Lorg/jacop/core/IntVar;[ILorg/jacop/core/IntVar;)V")
CLASS_INFO(XplusYeqZ, "org/jacop/constraints/XplusYeqZ",
    "(Lorg/jacop/core/IntVar;Lorg/jacop/core/IntVar;Lorg/jacop/core/IntVar;)V")
CLASS_INFO(XplusCeqZ, "org/jacop/constraints/XplusCeqZ",
    "(Lorg/jacop/core/IntVar;ILorg/jacop/core/IntVar;)V")
CLASS_INFO(XmulYeqZ, "org/jacop/constraints/XmulYeqZ",
    "(Lorg/jacop/core/IntVar;Lorg/jacop/core/IntVar;Lorg/jacop/core/IntVar;)V")
CLASS_INFO(XmulCeqZ, "org/jacop/constraints/XmulCeqZ",
    "(Lorg/jacop/core/IntVar;ILorg/jacop/core/IntVar;)V")
CLASS_INFO(XdivYeqZ, "org/jacop/constraints/XdivYeqZ",
    "(Lorg/jacop/core/IntVar;Lorg/jacop/core/IntVar;Lorg/jacop/core/IntVar;)V")
CLASS_INFO(XmodYeqZ, "org/jacop/constraints/XmodYeqZ",
    "(Lorg/jacop/core/IntVar;Lorg/jacop/core/IntVar;Lorg/jacop/core/IntVar;)V")
CLASS_INFO(XexpYeqZ, "org/jacop/constraints/XexpYeqZ",
    "(Lorg/jacop/core/IntVar;Lorg/jacop/core/IntVar;Lorg/jacop/core/IntVar;)V")
CLASS_INFO(XeqY, "org/jacop/constraints/XeqY",
    "(Lorg/jacop/core/IntVar;Lorg/jacop/core/IntVar;)V")
CLASS_INFO(XeqC, "org/jacop/constraints/XeqC", "(Lorg/jacop/core/IntVar;I)V")
CLASS_INFO(XneqY, "org/jacop/constraints/XneqY",
    "(Lorg/jacop/core/IntVar;Lorg/jacop/core/IntVar;)V")
CLASS_INFO(XltY, "org/jacop/constraints/XltY",
    "(Lorg/jacop/core/IntVar;Lorg/jacop/core/IntVar;)V")
CLASS_INFO(XlteqY, "org/jacop/constraints/XlteqY",
    "(Lorg/jacop/core/IntVar;Lorg/jacop/core/IntVar;)V")
CLASS_INFO(XgtY, "org/jacop/constraints/XgtY",
    "(Lorg/jacop/core/IntVar;Lorg/jacop/core/IntVar;)V")
CLASS_INFO(XgteqY, "org/jacop/constraints/XgteqY",
    "(Lorg/jacop/core/IntVar;Lorg/jacop/core/IntVar;)V")
CLASS_INFO(XgteqC, "org/jacop/constraints/XgteqC",
    "(Lorg/jacop/core/IntVar;Lorg/jacop/core/IntVar;)V")
CLASS_INFO(AbsXeqY, "org/jacop/constraints/AbsXeqY",
    "(Lorg/jacop/core/IntVar;Lorg/jacop/core/IntVar;)V")
CLASS_INFO(Min, "org/jacop/constraints/Min",
    "([Lorg/jacop/core/IntVar;Lorg/jacop/core/IntVar;)V")
CLASS_INFO(Max, "org/jacop/constraints/Max",
    "([Lorg/jacop/core/IntVar;Lorg/jacop/core/IntVar;)V")
CLASS_INFO(Count, "org/jacop/constraints/Count",
    "([Lorg/jacop/core/IntVar;Lorg/jacop/core/IntVar;I)V")
CLASS_INFO(IfThen, "org/jacop/constraints/IfThen",
    "(Lorg/jacop/constraints/PrimitiveConstraint;"
    "Lorg/jacop/constraints/PrimitiveConstraint;)V")
CLASS_INFO(IfThenElse, "org/jacop/constraints/IfThenElse",
    "(Lorg/jacop/constraints/PrimitiveConstraint;"
    "Lorg/jacop/constraints/PrimitiveConstraint;"
    "Lorg/jacop/constraints/PrimitiveConstraint;)V")
CLASS_INFO(Or, "org/jacop/constraints/Or",
    "(Lorg/jacop/constraints/PrimitiveConstraint;"
    "Lorg/jacop/constraints/PrimitiveConstraint;)V")
CLASS_INFO(And, "org/jacop/constraints/And",
    "(Lorg/jacop/constraints/PrimitiveConstraint;"
    "Lorg/jacop/constraints/PrimitiveConstraint;)V")
CLASS_INFO(Not, "org/jacop/constraints/Not",
    "(Lorg/jacop/constraints/PrimitiveConstraint;)V")
CLASS_INFO(Eq, "org/jacop/constraints/Eq",
    "(Lorg/jacop/constraints/PrimitiveConstraint;"
    "Lorg/jacop/constraints/PrimitiveConstraint;)V")
CLASS_INFO(Alldiff, "org/jacop/constraints/Alldiff", "([Lorg/jacop/core/IntVar;)V")
CLASS_INFO(DepthFirstSearch, "org/jacop/search/DepthFirstSearch", "()V")
CLASS_INFO(SimpleTimeOut, "org/jacop/search/SimpleTimeOut", "()V")
CLASS_INFO(InterruptSearch, "InterruptSearch", "()V")
CLASS_INFO(InterruptingListener, "InterruptingListener", "(J)V")
CLASS_INFO(SolutionListener, "SolutionListener", "(J)V")

// Converter of constraint programming problems from NL to JaCoP format.
class NLToJaCoPConverter :
  public asl::ExprConverter<NLToJaCoPConverter, jobject> {
 private:
  Env env_;
  jobject store_;
  jmethodID impose_;
  jobjectArray var_array_;
  std::vector<jobject> vars_;
  jobject obj_;
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

  jint CastToInt(double value) const;

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

  jobject Convert(asl::VarArgExpr e, ClassBase &cls) {
    jobjectArray args = CreateVarArray(
        static_cast<jsize>(std::distance(e.begin(), e.end())));
    int index = 0;
    for (asl::VarArgExpr::iterator i = e.begin(), end = e.end(); i != end; ++i)
      env_.SetObjectArrayElement(args, index++, Visit(*i));
    return CreateCon(cls, args);
  }

  // Converts a binary logical expression to a JaCoP constraint of class cls.
  template <typename ExprType>
  jobject Convert(ExprType e, ClassBase &cls) {
    return cls.NewObject(env_, Visit(e.lhs()), Visit(e.rhs()));
  }

  // Converts an iterated logical expression.
  jobject Convert(asl::IteratedLogicalExpr e, ClassBase &cls, jmethodID &ctor);

  void Impose(jobject constraint) {
    env_.CallVoidMethod(store_, impose_, constraint);
  }

  static void RequireZeroRHS(asl::BinaryExpr e, const std::string &func_name);

  template<typename Term>
  void ConvertExpr(asl::LinearExpr<Term> linear,
      asl::NumericExpr nonlinear, jobject result_var);

  jobject Convert(ClassBase &logop_class, jmethodID &logop_ctor,
                  ClassBase &eq_class, asl::PairwiseExpr e);

 public:
  explicit NLToJaCoPConverter();

  // Converts a logical constraint.
  void ConvertLogicalCon(asl::LogicalExpr e);

  void Convert(const ASLProblem &p);

  jobject store() const { return store_; }
  jobjectArray var_array() const { return var_array_; }
  const std::vector<jobject> &vars() const { return vars_; }
  Class<IntVar> &var_class() { return var_class_; }
  jobject obj() const { return obj_; }

  // The methods below perform conversion of AMPL NL expressions into
  // equivalent JaCoP expressions. JaCoP doesn't support the following
  // expressions/functions:
  // * division other than the integer one
  // * trigonometric & hyperbolic functions
  // * log, log10, exp, sqrt

  jobject VisitAdd(asl::BinaryExpr e);

  jobject VisitSub(asl::BinaryExpr e) {
    return CreateMinus(Visit(e.lhs()), Visit(e.rhs()));
  }

  jobject VisitMul(asl::BinaryExpr e) {
    return CreateCon(mul_class_, Visit(e.lhs()), Visit(e.rhs()));
  }

  jobject VisitMod(asl::BinaryExpr e) {
    return CreateCon(mod_class_, Visit(e.lhs()), Visit(e.rhs()));
  }

  jobject VisitPow(asl::BinaryExpr e) {
    return CreateCon(exp_class_, Visit(e.lhs()), Visit(e.rhs()));
  }

  jobject VisitLess(asl::BinaryExpr e);

  jobject VisitMin(asl::VarArgExpr e) {
    return Convert(e, min_class_);
  }

  jobject VisitMax(asl::VarArgExpr e) {
    return Convert(e, max_class_);
  }

  jobject VisitMinus(asl::UnaryExpr e) {
    return CreateCon(mul_const_class_, Visit(e.arg()), -1);
  }

  jobject VisitAbs(asl::UnaryExpr e) {
    return CreateCon(abs_class_, Visit(e.arg()));
  }

  jobject VisitFloor(asl::UnaryExpr e) {
    // floor does nothing because JaCoP supports only integer expressions.
    return Visit(e.arg());
  }

  jobject VisitCeil(asl::UnaryExpr e) {
    // ceil does nothing because JaCoP supports only integer expressions.
    return Visit(e.arg());
  }

  jobject VisitIf(asl::IfExpr e);

  jobject VisitSum(asl::SumExpr e);

  jobject VisitIntDiv(asl::BinaryExpr e) {
    return CreateCon(div_class_, Visit(e.lhs()), Visit(e.rhs()));
  }

  jobject VisitRound(asl::BinaryExpr e) {
    // round does nothing because JaCoP supports only integer expressions.
    RequireZeroRHS(e, "round");
    return Visit(e.lhs());
  }

  jobject VisitTrunc(asl::BinaryExpr e) {
    // trunc does nothing because JaCoP supports only integer expressions.
    RequireZeroRHS(e, "trunc");
    return Visit(e.lhs());
  }

  jobject VisitCount(asl::CountExpr e);

  jobject VisitNumberOf(asl::NumberOfExpr e);

  jobject VisitPowConstExp(asl::BinaryExpr e) {
    return VisitPow(e);
  }

  jobject VisitPow2(asl::UnaryExpr e) {
    jobject arg = Visit(e.arg());
    return CreateCon(mul_class_, arg, arg);
  }

  jobject VisitPowConstBase(asl::BinaryExpr e) {
    return VisitPow(e);
  }

  jobject VisitNumericConstant(asl::NumericConstant c) {
    return CreateConst(CastToInt(c.value()));
  }

  jobject VisitVariable(asl::Variable v) {
    return vars_[v.index()];
  }

  jobject VisitOr(asl::BinaryLogicalExpr e) {
    return Convert(e, or_class_);
  }

  jobject VisitAnd(asl::BinaryLogicalExpr e) {
    return Convert(e, and_class_);
  }

  jobject VisitLT(asl::RelationalExpr e) {
    return Convert(e, lt_class_);
  }

  jobject VisitLE(asl::RelationalExpr e) {
    return Convert(e, le_class_);
  }

  jobject VisitEQ(asl::RelationalExpr e) {
    return Convert(e, eq_class_);
  }

  jobject VisitGE(asl::RelationalExpr e) {
    return Convert(e, ge_class_);
  }

  jobject VisitGT(asl::RelationalExpr e) {
    return Convert(e, gt_class_);
  }

  jobject VisitNE(asl::RelationalExpr e) {
    return Convert(e, ne_class_);
  }

  jobject VisitNot(asl::NotExpr e) {
    return not_class_.NewObject(env_, Visit(e.arg()));
  }

  jobject VisitForAll(asl::IteratedLogicalExpr e) {
    return Convert(e, and_class_, and_array_ctor_);
  }

  jobject VisitExists(asl::IteratedLogicalExpr e) {
    return Convert(e, or_class_, or_array_ctor_);
  }

  jobject VisitImplication(asl::ImplicationExpr e);

  jobject VisitIff(asl::BinaryLogicalExpr e) {
    return eq_con_class_.NewObject(env_, Visit(e.lhs()), Visit(e.rhs()));
  }

  jobject VisitAllDiff(asl::PairwiseExpr e) {
    return Convert(and_class_, and_array_ctor_, ne_class_, e);
  }

  jobject VisitNotAllDiff(asl::PairwiseExpr e) {
    return Convert(or_class_, or_array_ctor_, eq_class_, e);
  }

  jobject VisitLogicalConstant(asl::LogicalConstant c) {
    if (!one_var_)
      one_var_ = var_class_.NewObject(env_, store_, 1, 1);
    return eq_const_class_.NewObject(env_, one_var_, c.value());
  }
};

// JaCoP solver.
class JaCoPSolver : public ASLSolver {
 private:
  std::vector<std::string> jvm_options_;
  jlong outlev_;
  double output_frequency_;
  steady_clock::time_point next_output_time_;
  unsigned output_count_;
  std::string header_;
  const char *var_select_;
  const char *val_select_;

  // The limits must be jlong to comply with setters' signatures.
  jlong time_limit_;  // Time limit in seconds.
  jlong node_limit_;
  jlong fail_limit_;
  jlong backtrack_limit_;
  jlong decision_limit_;
  jlong solution_limit_;

  int solve_code_;
  std::string status_;

  void SetStatus(int solve_code, const char *status) {
    solve_code_ = solve_code;
    status_ = status;
  }

  Env env_;
  GlobalRef search_;
  jmethodID get_depth_;
  jmethodID get_nodes_;
  jmethodID get_fails_;
  jmethodID value_;

  int DoGetIntOption(const SolverOption &, jlong *option) const {
    return static_cast<int>(*option);
  }

  void DoSetIntOption(const SolverOption &opt, int value, jlong *option) {
    if (value < 0)
      throw InvalidOptionValue(opt, value);
    *option = value;
  }

  void SetBoolOption(const SolverOption &opt, int value, jlong *option) {
    if (value != 0 && value != 1)
      throw InvalidOptionValue(opt, value);
    *option = value;
  }

  std::string GetEnumOption(
      const SolverOption &opt, const char **ptr) const;
  void SetEnumOption(const SolverOption &opt,
      fmt::StringRef value, const char **ptr);

  double GetOutputFrequency(const SolverOption &) const {
    return output_frequency_;
  }
  void SetOutputFrequency(const SolverOption &opt, double value);

  void Output(fmt::StringRef format, const fmt::ArgList &args);
  FMT_VARIADIC(void, Output, fmt::StringRef)

  steady_clock::duration GetOutputInterval() const {
    return steady_clock::duration(
          static_cast<steady_clock::rep>(output_frequency_ *
              steady_clock::period::den / steady_clock::period::num));
  }

  // Prints the solution log entry if the time is right.
  void PrintLogEntry();

  static JNIEXPORT jboolean JNICALL Stop(JNIEnv *, jobject, jlong data);

  // Relays the solution from JaCoP to the solution handler.
  class SolutionRelay {
   private:
    JaCoPSolver &solver_;
    SolutionHandler &sol_handler_;
    ASLProblem &problem_;
    const jobject *vars_;
    jobject obj_var_;
    bool multiple_sol_;
    jlong num_solutions_;
    std::string feasible_sol_message_;
    std::vector<double> solution_;

    bool DoHandleSolution();

   public:
    SolutionRelay(JaCoPSolver &s, SolutionHandler &sh,
        ASLProblem &p, const jobject *vars, jobject obj_var)
    : solver_(s), sol_handler_(sh), problem_(p), vars_(vars), obj_var_(obj_var),
      multiple_sol_(s.need_multiple_solutions()), num_solutions_(0) {
      if (multiple_sol_) {
        feasible_sol_message_ =
            fmt::format("{}: feasible solution", s.long_name());
        solution_.resize(p.num_vars());
      }
    }

    static JNIEXPORT jboolean JNICALL HandleSolution(
        JNIEnv *, jobject, jlong data) {
      return reinterpret_cast<SolutionRelay*>(data)->DoHandleSolution();
    }
  };

 protected:
  void DoSolve(ASLProblem &p, SolutionHandler &sh);

  void HandleUnknownOption(const char *name);

 public:
  JaCoPSolver();
};
}

#endif  // MP_SOLVERS_JACOP_H_
