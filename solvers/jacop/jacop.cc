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

#include "jacop/jacop.h"
#include "mp/os.h"

#include <iterator>

using namespace mp::asl;

namespace {

// Casts value to intptr_t.
#define CAST(value) reinterpret_cast<intptr_t>(value)

const mp::OptionValueInfo VAR_SELECT_VALUES[] = {
  {
    "largestdomain",
    "select the variable which has the largest domain size",
    CAST("LargestDomain")
  },
  {
    "largestmax",
    "select the variable with the largest maximal value in its domain",
    CAST("LargestMax")
  },
  {
    "largestmin",
    "select the variable with the largest minimal value in its domain",
    CAST("LargestMin")
  },
  {
    "maxregret",
    "max regret selector",
    CAST("MaxRegret")
  },
  {
    "mindomainoverdegree",
    "select the variable based on the minimal value of domain size divided "
    "by the number of constraints currently attached to a variable",
    CAST("MinDomainOverDegree")
  },
  {
    "mostconstraineddynamic",
    "select the variable which has the most pending constraints assigned to it",
    CAST("MostConstrainedDynamic")
  },
  {
    "mostconstrainedstatic",
    "select the variable which has the most constraints assigned to it",
    CAST("MostConstrainedStatic")
  },
  {
    "smallestdomain",
    "select the variable which has the smallest domain size",
    CAST("SmallestDomain")
  },
  {
    "smallestmax",
    "select the variable with the smallest maximal value in its domain",
    CAST("SmallestMax")
  },
  {
    "smallestmin",
    "select the variable with the smallest minimal value in its domain",
    CAST("SmallestMin")
  },
  {
    "weighteddegree",
    "select the variable with the highest weight divided by its size; "
    "every time a constraint failure is encountered all variables within "
    "the scope of that constraints have increased weight",
    CAST("WeightedDegree")
  }
};

const mp::OptionValueInfo VAL_SELECT_VALUES[] = {
  {
    "indomainmax",
    "select the maximal value in the domain of the variable",
    CAST("IndomainMax")
  },
  {
    "indomainmedian",
    "select the median value in the domain of the variable and then right "
    "and left values",
    CAST("IndomainMedian")
  },
  {
    "indomainmiddle",
    "select the middle value in the domain of the variable and then right "
    "and left values",
    CAST("IndomainMiddle")
  },
  {
    "indomainmin",
    "select the minimal value in the domain of the variable",
    CAST("IndomainMin")
  },
  {
    "indomainrandom",
    "select the random value in the domain of the variable; can split "
    "domains into multiple intervals",
    CAST("IndomainRandom")
  },
  {
    "indomainsimplerandom",
    "similar to indomainrandom, but faster and does not achieve uniform "
    "probability",
    CAST("IndomainSimpleRandom")
  }
};
}

namespace mp {

jint NLToJaCoPConverter::CastToInt(double value) const {
  jint int_value = static_cast<jint>(value);
  if (int_value != value)
    throw Error("value {} can't be represented as int", value);
  if (int_value < min_int_ || int_value > max_int_)
    throw Error("value {} is out of bounds", value);
  return int_value;
}

jobject NLToJaCoPConverter::Convert(
    IteratedLogicalExpr e, ClassBase &cls, jmethodID &ctor) {
  if (!ctor) {
    cls.Init(env_);
    ctor = env_.GetMethod(cls.get(),
        "<init>", "([Lorg/jacop/constraints/PrimitiveConstraint;)V");
  }
  if (!constraint_class_) {
    constraint_class_ = env_.FindClass(
        "org/jacop/constraints/PrimitiveConstraint");
  }
  int num_args = e.num_args();
  jobjectArray args = env_.NewObjectArray(num_args, constraint_class_, 0);
  for (int i = 0; i < num_args; ++i)
    env_.SetObjectArrayElement(args, i, Visit(e[i]));
  return env_.NewObject(cls.get(), ctor, args);
}

void NLToJaCoPConverter::RequireZeroRHS(
    BinaryExpr e, const std::string &func_name) {
  if (!IsZero(e.rhs()))
    throw MakeUnsupportedError("{} with nonzero second parameter", func_name);
}

template <typename Term>
void NLToJaCoPConverter::ConvertExpr(
    LinearExpr<Term> linear, NumericExpr nonlinear, jobject result_var) {
  jsize num_terms = static_cast<jsize>(
      std::distance(linear.begin(), linear.end()));
  if (nonlinear) {
    NumericConstant n = Cast<NumericConstant>(nonlinear);
    if (n && n.value() == 0)
      nonlinear = NumericExpr();
  }
  if (num_terms != 0) {
    if (nonlinear)
      ++num_terms;
    std::vector<jint> coefs(num_terms);
    jobjectArray vars = CreateVarArray(num_terms);
    int index = 0;
    for (typename LinearExpr<Term>::iterator
        i = linear.begin(), end = linear.end(); i != end; ++i, ++index) {
      coefs[index] = CastToInt(i->coef());
      env_.SetObjectArrayElement(vars, index, vars_[i->var_index()]);
    }
    if (nonlinear) {
      assert(index == num_terms - 1);
      coefs[index] = 1;
      env_.SetObjectArrayElement(vars, index, Visit(nonlinear));
    }
    jintArray coefs_array = env_.NewIntArray(num_terms);
    env_.SetIntArrayRegion(coefs_array, 0, num_terms, &coefs[0]);
    Impose(sum_weight_class_.NewObject(env_, vars, coefs_array, result_var));
  } else if (nonlinear)
    Impose(eq_class_.NewObject(env_, Visit(nonlinear), result_var));
}

void NLToJaCoPConverter::ConvertLogicalCon(LogicalExpr e) {
  if (e.kind() != expr::ALLDIFF) {
    Impose(Visit(e));
    return;
  }
  PairwiseExpr alldiff = Cast<PairwiseExpr>(e);
  int num_args = alldiff.num_args();
  jobjectArray args = CreateVarArray(num_args);
  for (int i = 0; i < num_args; ++i) {
    NumericExpr arg = alldiff[i];
    jobject result_var = 0;
    if (Variable var = Cast<Variable>(arg))
      result_var = vars_[var.index()];
    else
      result_var = Visit(arg);
    env_.SetObjectArrayElement(args, i, result_var);
  }
  Impose(alldiff_class_.NewObject(env_, args));
}

jobject NLToJaCoPConverter::Convert(
    ClassBase &logop_class, jmethodID &logop_ctor,
    ClassBase &eq_class, PairwiseExpr e) {
  if (!logop_ctor) {
    logop_class.Init(env_);
    logop_ctor = env_.GetMethod(logop_class.get(),
        "<init>", "([Lorg/jacop/constraints/PrimitiveConstraint;)V");
  }
  int n = e.num_args();
  std::vector<jobject> args(n);
  int index = 0;
  for (PairwiseExpr::iterator i = e.begin(), end = e.end(); i != end; ++i)
    args[index++] = Visit(*i);
  if (!constraint_class_) {
    constraint_class_ = env_.FindClass(
        "org/jacop/constraints/PrimitiveConstraint");
  }
  jobjectArray and_args = env_.NewObjectArray(
        n * (n - 1) / 2, constraint_class_, 0);
  index = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      env_.SetObjectArrayElement(
            and_args, index++, eq_class.NewObject(env_, args[i], args[j]));
    }
  }
  return env_.NewObject(logop_class.get(), logop_ctor, and_args);
}

NLToJaCoPConverter::NLToJaCoPConverter()
: env_(JVM::env()), store_(), impose_(), var_array_(), obj_(),
  constraint_class_(), or_array_ctor_(), and_array_ctor_(), one_var_() {
  var_class_.Init(env_);
  jclass domain_class = env_.FindClass("org/jacop/core/IntDomain");
  min_int_ = env_.GetStaticIntField(
      domain_class, env_.GetStaticFieldID(domain_class, "MinInt", "I"));
  max_int_ = env_.GetStaticIntField(
      domain_class, env_.GetStaticFieldID(domain_class, "MaxInt", "I"));
}

jobject NLToJaCoPConverter::VisitAdd(BinaryExpr e) {
  NumericExpr lhs = e.lhs(), rhs = e.rhs();
  if (NumericConstant c = Cast<NumericConstant>(lhs))
    return CreateCon(plus_const_class_, Visit(e.rhs()), CastToInt(c.value()));
  if (NumericConstant c = Cast<NumericConstant>(rhs))
    return CreateCon(plus_const_class_, Visit(e.lhs()), CastToInt(c.value()));
  return CreateCon(plus_class_, Visit(lhs), Visit(rhs));
}

jobject NLToJaCoPConverter::VisitLess(BinaryExpr e) {
  jobjectArray args = CreateVarArray(2);
  env_.SetObjectArrayElement(args, 0,
      CreateMinus(Visit(e.lhs()), Visit(e.rhs())));
  env_.SetObjectArrayElement(args, 1, CreateConst(0));
  return CreateCon(max_class_, args);
}

void NLToJaCoPConverter::Convert(const ASLProblem &p) {
  if (p.num_continuous_vars() != 0)
    throw Error("JaCoP doesn't support continuous variables");

  jclass store_class = env_.FindClass("org/jacop/core/Store");
  store_ = env_.NewObject(store_class,
      env_.GetMethod(store_class, "<init>", "()V"));
  impose_ = env_.GetMethod(store_class,
      "impose", "(Lorg/jacop/constraints/Constraint;)V");

  int num_vars = p.num_integer_vars();
  var_array_ = CreateVarArray(num_vars);
  vars_.resize(num_vars);
  for (int j = 0; j < num_vars; ++j) {
    ASLProblem::Variable var = p.var(j);
    double lb = var.lb(), ub = var.ub();
    jobject jvar = var_class_.NewObject(env_, store_,
        lb <= negInfinity ? min_int_ : CastToInt(lb),
        ub >= Infinity ? max_int_ : CastToInt(ub));
    vars_[j] = jvar;
    env_.SetObjectArrayElement(var_array_, j, jvar);
  }

  if (p.num_objs() > 0) {
    ASLProblem::Objective obj = p.obj(0);
    jobject result_var = var_class_.NewObject(env_, store_, min_int_, max_int_);
    ConvertExpr(obj.linear_expr(), obj.nonlinear_expr(), result_var);
    obj_ = obj.type() == obj::MIN ?
        result_var : CreateCon(mul_const_class_, result_var, -1);
  }

  // Convert algebraic constraints.
  for (int i = 0, n = p.num_algebraic_cons(); i < n; ++i) {
    double lb = p.con_lb(i), ub = p.con_ub(i);
    jint int_lb = lb <= negInfinity ? min_int_ : CastToInt(lb);
    jint int_ub = ub >= Infinity ? max_int_ : CastToInt(ub);
    jobject result_var = var_class_.NewObject(env_, store_, int_lb, int_ub);
    ConvertExpr(p.linear_con_expr(i), p.nonlinear_con_expr(i), result_var);
  }

  // Convert logical constraints.
  for (int i = 0, n = p.num_logical_cons(); i < n; ++i)
    ConvertLogicalCon(p.logical_con_expr(i));
}

jobject NLToJaCoPConverter::VisitIf(IfExpr e) {
  jobject result_var = CreateVar();
  Impose(if_else_class_.NewObject(env_, Visit(e.condition()),
      eq_class_.NewObject(env_, result_var, Visit(e.true_expr())),
      eq_class_.NewObject(env_, result_var, Visit(e.false_expr()))));
  return result_var;
}

jobject NLToJaCoPConverter::VisitSum(SumExpr e) {
  jobjectArray args = CreateVarArray(e.num_args());
  int index = 0;
  for (SumExpr::iterator i = e.begin(), end = e.end(); i != end; ++i)
    env_.SetObjectArrayElement(args, index++, Visit(*i));
  return CreateCon(sum_class_, args);
}

jobject NLToJaCoPConverter::VisitCount(CountExpr e) {
  jobjectArray args = CreateVarArray(e.num_args());
  int index = 0;
  for (CountExpr::iterator i = e.begin(), end = e.end(); i != end; ++i) {
    jobject result_var = CreateVar();
    Impose(if_else_class_.NewObject(env_, Visit(*i),
        eq_const_class_.NewObject(env_, result_var, 1),
        eq_const_class_.NewObject(env_, result_var, 0)));
    env_.SetObjectArrayElement(args, index++, result_var);
  }
  return CreateCon(sum_class_, args);
}

jobject NLToJaCoPConverter::VisitNumberOf(NumberOfExpr e) {
  // JaCoP only supports count constraints with constant value.
  NumericConstant num = Cast<NumericConstant>(e[0]);
  if (!num)
    throw MakeUnsupportedError("numberof with variable value");
  jobject result_var = CreateVar();
  int num_args = e.num_args();
  jobjectArray args = CreateVarArray(num_args - 1);
  for (int i = 1; i < num_args; ++i)
    env_.SetObjectArrayElement(args, i - 1, Visit(e[i]));
  Impose(count_class_.NewObject(
      env_, args, result_var, CastToInt(num.value())));
  return result_var;
}

jobject NLToJaCoPConverter::VisitImplication(ImplicationExpr e) {
  jobject condition = Visit(e.condition());
  return if_else_class_.NewObject(env_, condition,
      Visit(e.true_expr()), Visit(e.false_expr()));
}

void JaCoPSolver::SetOutputFrequency(const SolverOption &opt, double value) {
  if (value <= 0)
    throw InvalidOptionValue(opt, value);
  output_frequency_ = value;
}

void JaCoPSolver::HandleUnknownOption(const char *name) {
  if (name[0] == '-') {
    Print("{}\n", name);
    jvm_options_.push_back(name);
  } else {
    Solver::HandleUnknownOption(name);
  }
}

JaCoPSolver::JaCoPSolver()
: ASLSolver("jacop", "jacop " JACOP_VERSION, 20141107, MULTIPLE_SOL),
  outlev_(0), output_frequency_(1), output_count_(0),
  var_select_("SmallestDomain"), val_select_("IndomainMin"),
  time_limit_(-1), node_limit_(-1), fail_limit_(-1),
  backtrack_limit_(-1), decision_limit_(-1),
  solution_limit_(-1), solve_code_(-1),
  get_depth_(), get_nodes_(), get_fails_(), value_() {

  set_version("JaCoP " JACOP_VERSION);

  set_option_header(
      "JaCoP Options for AMPL\n"
      "----------------------\n"
      "\n"
      "To set these options, assign a string specifying their values to "
      "the AMPL option ``jacop_options``. For example::\n"
      "\n"
      "  ampl: option jacop_options 'version nodelimit=30000';\n");

  AddDblOption("outfreq",
      "Output frequency in seconds. The value should be a positive integer.",
      &JaCoPSolver::GetOutputFrequency, &JaCoPSolver::SetOutputFrequency);

  AddIntOption("outlev", "0 or 1 (default 0):  Whether to print solution log.",
      &JaCoPSolver::DoGetIntOption, &JaCoPSolver::SetBoolOption, &outlev_);

  AddStrOption("var_select",
      "Variable selector. Possible values:\n"
      "\n"
      ".. value-table::\n"
      "\n"
      "The default value is ``smallestdomain``.",
      &JaCoPSolver::GetEnumOption, &JaCoPSolver::SetEnumOption,
      &var_select_, VAR_SELECT_VALUES);

  AddStrOption("val_select",
      "Value selector. Possible values:\n"
      "\n"
      ".. value-table::\n"
      "\n"
      "The default value is ``indomainmin``.",
      &JaCoPSolver::GetEnumOption, &JaCoPSolver::SetEnumOption,
      &val_select_, VAL_SELECT_VALUES);

  AddIntOption("timelimit", "Time limit in seconds.",
      &JaCoPSolver::DoGetIntOption, &JaCoPSolver::DoSetIntOption, &time_limit_);
  AddIntOption("nodelimit", "Node limit.",
      &JaCoPSolver::DoGetIntOption, &JaCoPSolver::DoSetIntOption, &node_limit_);
  AddIntOption("faillimit", "Fail (wrong decision) limit.",
      &JaCoPSolver::DoGetIntOption, &JaCoPSolver::DoSetIntOption, &fail_limit_);
  AddIntOption("backtracklimit", "Backtrack limit.",
      &JaCoPSolver::DoGetIntOption, &JaCoPSolver::DoSetIntOption,
      &backtrack_limit_);
  AddIntOption("decisionlimit", "Decision limit.",
      &JaCoPSolver::DoGetIntOption, &JaCoPSolver::DoSetIntOption,
      &decision_limit_);
  AddIntOption("solutionlimit",
      "Limit on the number of feasible solutions found before terminating "
      "a search. Leaving the solution limit unspecified will make the "
      "optimizer search for an optimal solution if there is an objective "
      "function or for a feasible solution otherwise.",
      &JaCoPSolver::DoGetIntOption, &JaCoPSolver::DoSetIntOption,
      &solution_limit_);
}

std::string JaCoPSolver::GetEnumOption(
    const SolverOption &, const char **ptr) const {
  std::string value = *ptr;
  std::transform(value.begin(), value.end(), value.begin(), ::tolower);
  return value;
}

void JaCoPSolver::SetEnumOption(
    const SolverOption &opt, fmt::StringRef value, const char **ptr) {
  for (mp::ValueArrayRef::iterator
      i = opt.values().begin(), e = opt.values().end(); i != e; ++i) {
    if (std::strcmp(value.c_str(), i->value) == 0) {
      *ptr = reinterpret_cast<const char*>(i->data);
      return;
    }
  }
  throw InvalidOptionValue(opt, value);
}

void JaCoPSolver::Output(fmt::StringRef format, const fmt::ArgList &args) {
  if (output_count_ == 0)
    Print("{}", header_);
  output_count_ = (output_count_ + 1) % 20;
  Print(format, args);
}

void JaCoPSolver::PrintLogEntry() {
  if (outlev_ == 0 || steady_clock::now() < next_output_time_)
    return;
  Output("{:10} {:10} {:10}\n",
      env_.CallIntMethodKeepException(search_.get(), get_depth_),
      env_.CallIntMethodKeepException(search_.get(), get_nodes_),
      env_.CallIntMethodKeepException(search_.get(), get_fails_));
  next_output_time_ += GetOutputInterval();
}

bool JaCoPSolver::SolutionRelay::DoHandleSolution() {
  try {
    ++num_solutions_;
    if (solver_.outlev_ != 0 && obj_var_) {
      jint value = solver_.env_.CallIntMethodKeepException(
          obj_var_, solver_.value_);
      solver_.Output("{:46}\n",
                     (problem_.obj(0).type() == obj::MIN ? value : -value));
    }
    if (multiple_sol_) {
      double obj_value = obj_var_ ?
        solver_.env_.CallIntMethod(obj_var_, solver_.value_) : 0;
      for (int j = 0, n = problem_.num_vars(); j < n; ++j)
        solution_[j] = solver_.env_.CallIntMethod(vars_[j], solver_.value_);
      sol_handler_.HandleFeasibleSolution(feasible_sol_message_,
          solution_.empty() ? 0 : solution_.data(), 0, obj_value);
    }
    if (solver_.solution_limit_ != -1 &&
        num_solutions_ >= solver_.solution_limit_) {
      solver_.solve_code_ = 403;
      solver_.status_ = "solution limit";
      return true;
    }
  } catch (const JavaError &) {
    // This indicates that a Java exception has occurred and it will be
    // re-thrown when the control returns to the Java code. Therefore the
    // C++ JavaError exception is ignored here.
  }
  return false;
}

JNIEXPORT jboolean JNICALL JaCoPSolver::Stop(JNIEnv *, jobject, jlong data) {
  try {
    JaCoPSolver* solver = reinterpret_cast<JaCoPSolver*>(data);
    solver->PrintLogEntry();
    if (solver->interrupter()->Stop()) {
      solver->solve_code_ = 600;
      solver->status_ = "interrupted";
      return JNI_TRUE;
    }
  } catch (const JavaError &) {
    // This indicates that a Java exception has occurred and it will be
    // re-thrown when the control returns to the Java code. Therefore the
    // C++ JavaError exception is ignored here.
  }
  return JNI_FALSE;
}

void JaCoPSolver::DoSolve(ASLProblem &p, SolutionHandler &sh) {
  steady_clock::time_point time = steady_clock::now();

  std::vector<const char*> jvm_options(jvm_options_.size() + 2);
  for (size_t i = 0, n = jvm_options_.size(); i != n; ++i)
    jvm_options[i] = jvm_options_[i].c_str();
  std::string exe_dir = GetExecutablePath().remove_filename().string();
  std::string classpath =
      "-Djava.class.path=" + exe_dir + "/jacop-" JACOP_VERSION ".jar"
      AMPL_CLASSPATH_SEP + exe_dir + "/lib/jacop-" JACOP_VERSION ".jar"
      AMPL_CLASSPATH_SEP + exe_dir + "/ampljacop.jar";
  jvm_options[jvm_options_.size()] = classpath.c_str();
  env_ = JVM::env(&jvm_options[0]);

  // Set up an optimization problem in JaCoP.
  NLToJaCoPConverter converter;
  converter.Convert(p);

  Class<DepthFirstSearch> dfs_class;
  search_ = env_.NewGlobalRef(dfs_class.NewObject(env_));
  if (!get_depth_) {
    get_depth_ = env_.GetMethod(dfs_class.get(), "getMaximumDepth", "()I");
    get_nodes_ = env_.GetMethod(dfs_class.get(), "getNodes", "()I");
    get_fails_ = env_.GetMethod(dfs_class.get(), "getWrongDecisions", "()I");
  }
  jmethodID setPrintInfo =
      env_.GetMethod(dfs_class.get(), "setPrintInfo", "(Z)V");
  env_.CallVoidMethod(search_.get(), setPrintInfo, JNI_FALSE);
  jobject var_select = env_.NewObject(
      (std::string("org/jacop/search/") + var_select_).c_str(), "()V");
  jobject val_select = env_.NewObject(
      (std::string("org/jacop/search/") + val_select_).c_str(), "()V");
  jobject select = env_.NewObject("org/jacop/search/SimpleSelect",
      "([Lorg/jacop/core/Var;Lorg/jacop/search/ComparatorVariable;"
      "Lorg/jacop/search/Indomain;)V", converter.var_array(),
      var_select, val_select);
  double obj_val = std::numeric_limits<double>::quiet_NaN();
  bool has_obj = p.num_objs() > 0;

  Class<InterruptingListener> interrupter_class;
  interrupter_class.Init(env_);
  {
    char NAME[] = "stop";
    char SIG[] = "(J)Z";
    JNINativeMethod method = {
        NAME, SIG, reinterpret_cast<void*>(reinterpret_cast<jlong>(Stop))
    };
    env_.RegisterNatives(interrupter_class.get(), &method, 1);
  }
  jobject interrupter =
      interrupter_class.NewObject(env_, reinterpret_cast<jlong>(this));
  env_.CallVoidMethod(search_.get(),
       env_.GetMethod(dfs_class.get(), "setConsistencyListener",
           "(Lorg/jacop/search/ConsistencyListener;)V"), interrupter);

  Class<SolutionListener> solution_listener_class;
  solution_listener_class.Init(env_);
  {
    char NAME[] = "handleSolution";
    char SIG[] = "(J)Z";
    JNINativeMethod method = {
        NAME, SIG,
        reinterpret_cast<void*>(
            reinterpret_cast<jlong>(SolutionRelay::HandleSolution))
    };
    env_.RegisterNatives(solution_listener_class.get(), &method, 1);
  }
  GlobalRef obj_var;  // The variable holding the objective value.
  if (has_obj)
    obj_var = env_.NewGlobalRef(converter.obj());
  jclass var_class = converter.var_class().get();
  value_ = env_.GetMethod(var_class, "value", "()I");
  SolutionRelay sol_relay(*this, sh, p, converter.vars().data(), obj_var.get());
  jobject solution_listener = solution_listener_class.NewObject(
      env_, reinterpret_cast<jlong>(&sol_relay));
  env_.CallVoidMethod(solution_listener, env_.GetMethod(
      solution_listener_class.get(),
      "setSolutionLimit", "(I)V"), std::numeric_limits<jint>::max());
  env_.CallVoidMethod(search_.get(),
      env_.GetMethod(dfs_class.get(), "setSolutionListener",
          "(Lorg/jacop/search/SolutionListener;)V"), solution_listener);

  // Set the limits.
  Class<SimpleTimeOut> timeout_class;
  jobject timeout = timeout_class.NewObject(env_);
  env_.CallVoidMethod(search_.get(),
       env_.GetMethod(dfs_class.get(), "setTimeOutListener",
           "(Lorg/jacop/search/TimeOutListener;)V"), timeout);
  if (time_limit_ != -1) {
    env_.CallVoidMethod(search_.get(),
         env_.GetMethod(dfs_class.get(), "setTimeOut", "(J)V"), time_limit_);
  }
  if (node_limit_ != -1) {
    env_.CallVoidMethod(search_.get(),
         env_.GetMethod(dfs_class.get(), "setNodesOut", "(J)V"), node_limit_);
  }
  if (fail_limit_ != -1) {
    env_.CallVoidMethod(search_.get(), env_.GetMethod(
        dfs_class.get(), "setWrongDecisionsOut", "(J)V"), fail_limit_);
  }
  if (backtrack_limit_ != -1) {
    env_.CallVoidMethod(search_.get(), env_.GetMethod(
        dfs_class.get(), "setBacktracksOut", "(J)V"), backtrack_limit_);
  }
  if (decision_limit_ != -1) {
    env_.CallVoidMethod(search_.get(), env_.GetMethod(
        dfs_class.get(), "setDecisionsOut", "(J)V"), decision_limit_);
  }

  double setup_time = GetTimeAndReset(time);

  // Solve the problem.
  solve_code_ = -1;
  header_ = fmt::format("{:>10} {:>10} {:>10} {:>13}\n",
    "Max Depth", "Nodes", "Fails", (has_obj ? "Best Obj" : ""));
  jboolean found = false;
  output_count_ = 0;
  next_output_time_ = steady_clock::now() + GetOutputInterval();
  try {
    if (has_obj) {
      jmethodID labeling = env_.GetMethod(dfs_class.get(), "labeling",
          "(Lorg/jacop/core/Store;Lorg/jacop/search/SelectChoicePoint;"
          "Lorg/jacop/core/Var;)Z");
      found = env_.CallBooleanMethod(
          search_.get(), labeling, converter.store(), select, obj_var.get());
    } else {
      jmethodID labeling = env_.GetMethod(dfs_class.get(), "labeling",
          "(Lorg/jacop/core/Store;Lorg/jacop/search/SelectChoicePoint;)Z");
      found = env_.CallBooleanMethod(
          search_.get(), labeling, converter.store(), select);
    }
  } catch (const JavaError &e) {
    // Check if exception is of class InterruptSearch which is used to
    // interrupt search.
    bool interrupted = false;
    if (jthrowable throwable = e.exception()) {
      Class<InterruptSearch> interrupt_class;
      interrupt_class.Init(env_);
      if (env_.IsInstanceOf(throwable, interrupt_class.get()))
        interrupted = true;
    }
    if (!interrupted)
      throw;
    if (solve_code_ == 403)
      found = true;
  }

  // Convert the solution status.
  if (env_.GetBooleanField(timeout,
      env_.GetFieldID(timeout_class.get(), "timeOutOccurred", "Z"))) {
    solve_code_ = 400;
    status_ = "limit";
  }
  if (found) {
    if (!has_obj) {
      // If the problem has no objectives and there is a solution, report
      // it as solved even if some limit was reached.
      solve_code_ = 0;
      status_ = "feasible solution";
    } else if (solve_code_ == -1) {
      solve_code_ = 0;
      obj_val = env_.CallIntMethod(obj_var.get(), value_);
      if (p.obj(0).type() == obj::MAX)
        obj_val = -obj_val;
      status_ = "optimal solution";
    }
  } else if (solve_code_ == -1) {
    solve_code_ = 200;
    status_ = "infeasible problem";
  }

  std::vector<double> final_solution;
  if (found) {
    const std::vector<jobject> &vars = converter.vars();
    int num_vars = p.num_vars();
    final_solution.resize(num_vars);
    for (int j = 0; j < num_vars; ++j)
      final_solution[j] = env_.CallIntMethod(vars[j], value_);
  }

  double solution_time = GetTimeAndReset(time);

  fmt::MemoryWriter w;
  w.write("{}: {}\n", long_name(), status_);
  w.write("{} nodes, {} fails",
      env_.CallIntMethod(search_.get(), get_nodes_),
      env_.CallIntMethod(search_.get(), get_fails_));
  if (has_obj && found)
    w.write(", objective {}", FormatObjValue(obj_val));
  sh.HandleSolution(solve_code_, w.c_str(),
      final_solution.empty() ? 0 : final_solution.data(), 0, obj_val);

  double output_time = GetTimeAndReset(time);

  if (timing()) {
    Print("Setup time = {:.6f}s\n"
          "Solution time = {:.6f}s\n"
          "Output time = {:.6f}s\n",
          setup_time, solution_time, output_time);
  }
}

SolverPtr CreateSolver(const char *) { return SolverPtr(new JaCoPSolver()); }
}  // namespace mp
