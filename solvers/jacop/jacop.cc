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

#include "solvers/jacop/jacop.h"
#include "solvers/util/filesystem.h"

#include <iterator>

#define JACOP_VERSION "3.2"

namespace {

const char *const VAR_SELECT[] = {
  "LargestDomain",
  "LargestMax",
  "LargestMin",
  "MaxRegret",
  "MinDomainOverDegree",
  "MostConstrainedDynamic",
  "MostConstrainedStatic",
  "SmallestDomain",
  "SmallestMax",
  "SmallestMin",
  "WeightedDegree",
  0
};

const char *const VAL_SELECT[] = {
  "IndomainMax",
  "IndomainMedian",
  "IndomainMiddle",
  "IndomainMin",
  "IndomainRandom",
  "IndomainSimpleRandom",
  0
};

bool Match(const char *value, const char *s) {
  while (*value && *value == std::tolower(*s))
    ++value, ++s;
  return *value == std::tolower(*s);
}
}

namespace ampl {

jint NLToJaCoPConverter::CastToInt(double value) const {
  jint int_value = static_cast<jint>(value);
  if (int_value != value) {
    throw UnsupportedExprError::CreateFromMessage(
        fmt::Format("value {} can't be represented as int") << value);
  }
  if (int_value < min_int_ || int_value > max_int_)
    throw Error(str(fmt::Format("value {} is out of bounds") << value));
  return int_value;
}

jobject NLToJaCoPConverter::Convert(
    IteratedLogicalExpr e, ClassBase &cls, jmethodID &ctor) {
  if (!ctor) {
    cls.Init(env_);
    ctor = env_.GetMethod(cls.get(),
        "<init>", "([LJaCoP/constraints/PrimitiveConstraint;)V");
  }
  if (!constraint_class_) {
    constraint_class_ = env_.FindClass(
        "JaCoP/constraints/PrimitiveConstraint");
  }
  int num_args = e.num_args();
  jobjectArray args = env_.NewObjectArray(num_args, constraint_class_, 0);
  for (int i = 0; i < num_args; ++i)
    env_.SetObjectArrayElement(args, i, Visit(e[i]));
  return env_.NewObject(cls.get(), ctor, args);
}

void NLToJaCoPConverter::RequireZeroRHS(
    BinaryExpr e, const std::string &func_name) {
  NumericConstant num = Cast<NumericConstant>(e.rhs());
  if (!num || num.value() != 0) {
    throw UnsupportedExprError::CreateFromExprString(
        func_name + " with nonzero second parameter");
  }
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
  AllDiffExpr alldiff = Cast<AllDiffExpr>(e);
  if (!alldiff) {
    Impose(Visit(e));
    return;
  }
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

NLToJaCoPConverter::NLToJaCoPConverter()
: env_(JVM::env()), store_(), impose_(), var_array_(), obj_(),
  constraint_class_(), or_array_ctor_(), and_array_ctor_(), one_var_() {
  var_class_.Init(env_);
  jclass domain_class = env_.FindClass("JaCoP/core/IntDomain");
  min_int_ = env_.GetStaticIntField(
      domain_class, env_.GetStaticFieldID(domain_class, "MinInt", "I"));
  max_int_ = env_.GetStaticIntField(
      domain_class, env_.GetStaticFieldID(domain_class, "MaxInt", "I"));
}

jobject NLToJaCoPConverter::VisitPlus(BinaryExpr e) {
  NumericExpr lhs = e.lhs(), rhs = e.rhs();
  if (NumericConstant c = Cast<NumericConstant>(lhs))
    return CreateCon(plus_const_class_, Visit(e.rhs()), CastToInt(c.value()));
  if (NumericConstant c = Cast<NumericConstant>(rhs))
    return CreateCon(plus_const_class_, Visit(e.lhs()), CastToInt(c.value()));
  return CreateCon(plus_class_, Visit(lhs), Visit(rhs));
}

jobject NLToJaCoPConverter::VisitNumericLess(BinaryExpr e) {
  jobjectArray args = CreateVarArray(2);
  env_.SetObjectArrayElement(args, 0,
      CreateMinus(Visit(e.lhs()), Visit(e.rhs())));
  env_.SetObjectArrayElement(args, 1, CreateConst(0));
  return CreateCon(max_class_, args);
}

void NLToJaCoPConverter::Convert(const Problem &p) {
  if (p.num_continuous_vars() != 0)
    throw Error("JaCoP doesn't support continuous variables");

  jclass store_class = env_.FindClass("JaCoP/core/Store");
  store_ = env_.NewObject(store_class,
      env_.GetMethod(store_class, "<init>", "()V"));
  impose_ = env_.GetMethod(store_class,
      "impose", "(LJaCoP/constraints/Constraint;)V");

  int num_vars = p.num_integer_vars();
  var_array_ = CreateVarArray(num_vars);
  vars_.resize(num_vars);
  for (int j = 0; j < num_vars; ++j) {
    double lb = p.var_lb(j), ub = p.var_ub(j);
    jobject var = var_class_.NewObject(env_, store_,
        lb <= negInfinity ? min_int_ : CastToInt(lb),
        ub >= Infinity ? max_int_ : CastToInt(ub));
    vars_[j] = var;
    env_.SetObjectArrayElement(var_array_, j, var);
  }

  if (p.num_objs() != 0) {
    jobject result_var = var_class_.NewObject(env_, store_, min_int_, max_int_);
    ConvertExpr(p.linear_obj_expr(0), p.nonlinear_obj_expr(0), result_var);
    obj_ = p.obj_type(0) == MIN ?
        result_var : CreateCon(mul_const_class_, result_var, -1);
  }

  // Convert constraints.
  for (int i = 0, n = p.num_cons(); i < n; ++i) {
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
  NumericConstant num = Cast<NumericConstant>(e.value());
  if (!num) {
    throw UnsupportedExprError::CreateFromExprString(
        "numberof with variable value");
  }
  jobject result_var = CreateVar();
  int num_args = e.num_args();
  jobjectArray args = CreateVarArray(num_args);
  for (int i = 0; i < num_args; ++i)
    env_.SetObjectArrayElement(args, i, Visit(e[i]));
  Impose(count_class_.NewObject(
      env_, args, result_var, CastToInt(num.value())));
  return result_var;
}

jobject NLToJaCoPConverter::VisitImplication(ImplicationExpr e) {
  jobject condition = Visit(e.condition());
  return if_else_class_.NewObject(env_, condition,
      Visit(e.true_expr()), Visit(e.false_expr()));
}

void JaCoPSolver::SetOutputFrequency(const char *name, double value) {
  if (value <= 0)
    throw InvalidOptionValue(name, value);
  output_frequency_ = value;
}

void JaCoPSolver::HandleUnknownOption(const char *name) {
  if (name[0] == '-') {
    Print("{}\n") << name;
    jvm_options_.push_back(name);
  } else {
    BasicSolver::HandleUnknownOption(name);
  }
}

JaCoPSolver::JaCoPSolver()
: Solver<JaCoPSolver>("jacop", "jacop " JACOP_VERSION, 20130820),
  outlev_(0), output_frequency_(1), output_count_(0),
  var_select_("SmallestDomain"), val_select_("IndomainMin"),
  time_limit_(-1), node_limit_(-1), fail_limit_(-1),
  backtrack_limit_(-1), decision_limit_(-1),
  get_depth_(), get_nodes_(), get_fails_(), obj_type_(MIN), value_() {

  set_version("JaCoP " JACOP_VERSION);

  AddDblOption("outfreq",
      "Output frequency in seconds.  The value should be a positive integer.",
      &JaCoPSolver::GetOutputFrequency, &JaCoPSolver::SetOutputFrequency);

  AddIntOption("outlev", "0 or 1 (default 0):  Whether to print solution log.",
      &JaCoPSolver::DoGetIntOption, &JaCoPSolver::SetBoolOption, &outlev_);

  AddStrOption("var_select",
      "Variable selector.  Possible values:\n"
      "      largestdomain          - select the variable which has the\n"
      "                               largest domain size\n"
      "      largestmax             - select the variable with the largest\n"
      "                               maximal value in its domain\n"
      "      largestmin             - select the variable with the largest\n"
      "                               minimal value in its domain\n"
      "      maxregret              - max regret selector\n"
      "      mindomainoverdegree    - select the variable based on the\n"
      "                               minimal value of domain size divided\n"
      "                               by the number of constraints currently\n"
      "                               attached to a variable\n"
      "      mostconstraineddynamic - select the variable which has the\n"
      "                               most pending constraints assign to it\n"
      "      mostconstrainedstatic  - select the variable which has the\n"
      "                               most constraints assign to it\n"
      "      smallestdomain         - select the variable which has the\n"
      "                               smallest domain size (default)\n"
      "      smallestmax            - select the variable with the smallest\n"
      "                               maximal value in its domain\n"
      "      smallestmin            - select the variable with the smallest\n"
      "                               minimal value in its domain\n"
      "      weighteddegree         - select the variable with the highest\n"
      "                               weight divided by its size; every time\n"
      "                               a constraint failure is encountered\n"
      "                               all variables within the scope of that\n"
      "                               constraints have increased weight. \n",
      &JaCoPSolver::GetEnumOption, &JaCoPSolver::SetEnumOption,
      OptionInfo(VAR_SELECT, var_select_));

  AddStrOption("val_select",
      "Value selector.  Possible values:\n"
      "      indomainmax          - select the maximal value in the domain\n"
      "                             of the variable\n"
      "      indomainmedian       - select the median value in the domain\n"
      "                             of the variable and then right and left\n"
      "                             values\n"
      "      indomainmiddle       - select the middle value in the domain\n"
      "                             of the variable and then right and left\n"
      "                             values\n"
      "      indomainmin          - select the minimal value in the domain\n"
      "                             of the variable (default)\n"
      "      indomainrandom       - select the random value in the domain\n"
      "                             of the variable; can split domains into\n"
      "                             multiple intervals\n"
      "      indomainsimplerandom - similar to indomainrandom, but faster\n"
      "                             and does not achieve uniform probability\n",
      &JaCoPSolver::GetEnumOption, &JaCoPSolver::SetEnumOption,
      OptionInfo(VAL_SELECT, val_select_));

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
}

std::string JaCoPSolver::GetEnumOption(
    const char *, const OptionInfo &info) const {
  std::string value = info.value;
  std::transform(value.begin(), value.end(), value.begin(), ::tolower);
  return value;
}

void JaCoPSolver::SetEnumOption(
    const char *name, const char *value, const OptionInfo &info) {
  for (const char *const *v = info.values; *v; ++v) {
    if (Match(value, *v)) {
      info.value = *v;
      return;
    }
  }
  throw InvalidOptionValue(name, value);
}

fmt::Formatter<BasicSolver::Printer>
    JaCoPSolver::Output(fmt::StringRef format) {
  if (output_count_ == 0)
    Print("{}") << header_;
  output_count_ = (output_count_ + 1) % 20;
  return fmt::Formatter<Printer>(format, MakePrinter());
}

void JaCoPSolver::PrintLogEntry() {
  if (outlev_ == 0 || steady_clock::now() < next_output_time_)
    return;
  Output("{:10} {:10} {:10}\n")
      << env_.CallIntMethodKeepException(search_.get(), get_depth_)
      << env_.CallIntMethodKeepException(search_.get(), get_nodes_)
      << env_.CallIntMethodKeepException(search_.get(), get_fails_);
  next_output_time_ += GetOutputInterval();
}

void JaCoPSolver::PrintObjValue() {
  try {
    if (outlev_ == 0)
      return;
    jint value = env_.CallIntMethodKeepException(obj_var_.get(), value_);
    Output("{:46}\n") << (obj_type_ == MIN ? value : -value);
  } catch (const JavaError &) {
    // This indicates that a Java exception has occurred and it will be
    // re-thrown when the control returns to the Java code. Therefore the
    // C++ JavaError exception is ignored here.
  }
}

JNIEXPORT jboolean JNICALL JaCoPSolver::Stop(JNIEnv *, jobject, jlong data) {
  try {
    reinterpret_cast<JaCoPSolver*>(data)->PrintLogEntry();
    return SignalHandler::stop();
  } catch (const JavaError &) {
    // This indicates that a Java exception has occurred and it will be
    // re-thrown when the control returns to the Java code. Therefore the
    // C++ JavaError exception is ignored here.
  }
  return JNI_FALSE;
}

std::string JaCoPSolver::GetOptionHeader() {
  return
      "JaCoP Directives for AMPL\n"
      "--------------------------\n"
      "\n"
      "To set these directives, assign a string specifying their values to "
      "the AMPL option jacop_options.  For example:\n"
      "\n"
      "  ampl: option jacop_options 'version nodelimit=30000';\n";
}

void JaCoPSolver::Solve(Problem &p) {
  steady_clock::time_point time = steady_clock::now();

  std::vector<const char*> jvm_options(jvm_options_.size() + 2);
  for (size_t i = 0, n = jvm_options_.size(); i != n; ++i)
    jvm_options[i] = jvm_options_[i].c_str();
  std::string exe_dir = GetExecutablePath().remove_filename().string();
  std::string classpath =
      "-Djava.class.path=" + exe_dir + "/JaCoP-" JACOP_VERSION ".jar"
      AMPL_CLASSPATH_SEP + exe_dir + "/lib/JaCoP-3.2.jar"
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
      (std::string("JaCoP/search/") + var_select_).c_str(), "()V");
  jobject val_select = env_.NewObject(
      (std::string("JaCoP/search/") + val_select_).c_str(), "()V");
  jobject select = env_.NewObject("JaCoP/search/SimpleSelect",
      "([LJaCoP/core/Var;LJaCoP/search/ComparatorVariable;"
      "LJaCoP/search/Indomain;)V", converter.var_array(),
      var_select, val_select);
  double obj_val = std::numeric_limits<double>::quiet_NaN();
  bool has_obj = p.num_objs() != 0;

  SignalHandler signal_handler(*this);
  Class<Interrupter> interrupter_class;
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
           "(LJaCoP/search/ConsistencyListener;)V"), interrupter);

  if (has_obj) {
    Class<SolutionListener> solution_listener_class;
    solution_listener_class.Init(env_);
    {
      char NAME[] = "handleSolution";
      char SIG[] = "(J)V";
      JNINativeMethod method = {
          NAME, SIG,
          reinterpret_cast<void*>(reinterpret_cast<jlong>(DoHandleSolution))
      };
      env_.RegisterNatives(solution_listener_class.get(), &method, 1);
    }
    jobject solution_listener =
        solution_listener_class.NewObject(env_, reinterpret_cast<jlong>(this));
    env_.CallVoidMethod(search_.get(),
         env_.GetMethod(dfs_class.get(), "setSolutionListener",
             "(LJaCoP/search/SolutionListener;)V"), solution_listener);
    obj_var_ = env_.NewGlobalRef(converter.obj());
    obj_type_ = p.obj_type(0);
  }
  jclass var_class = converter.var_class().get();
  value_ = env_.GetMethod(var_class, "value", "()I");

  // Set the limits.
  Class<SimpleTimeOut> timeout_class;
  jobject timeout = timeout_class.NewObject(env_);
  env_.CallVoidMethod(search_.get(),
       env_.GetMethod(dfs_class.get(), "setTimeOutListener",
           "(LJaCoP/search/TimeOutListener;)V"), timeout);
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
  header_ = str(fmt::Format("{:>10} {:>10} {:>10} {:>13}\n")
    << "Max Depth" << "Nodes" << "Fails" << (has_obj ? "Best Obj" : ""));
  jboolean found = false;
  bool interrupted = false;
  output_count_ = 0;
  next_output_time_ = steady_clock::now() + GetOutputInterval();
  try {
    if (has_obj) {
      jmethodID labeling = env_.GetMethod(dfs_class.get(), "labeling",
          "(LJaCoP/core/Store;LJaCoP/search/SelectChoicePoint;"
          "LJaCoP/core/IntVar;)Z");
      found = env_.CallBooleanMethod(
          search_.get(), labeling, converter.store(), select, obj_var_.get());
    } else {
      jmethodID labeling = env_.GetMethod(dfs_class.get(), "labeling",
          "(LJaCoP/core/Store;LJaCoP/search/SelectChoicePoint;)Z");
      found = env_.CallBooleanMethod(
          search_.get(), labeling, converter.store(), select);
    }
  } catch (const JavaError &e) {
    // Check if exception is of class InterruptSearch which is used to
    // interrupt search.
    if (jthrowable throwable = e.exception()) {
      Class<InterruptSearch> interrupt_class;
      interrupt_class.Init(env_);
      if (env_.IsInstanceOf(throwable, interrupt_class.get()))
        interrupted = true;
    }
    if (!interrupted)
      throw;
  }

  // Convert solution status.
  const char *status = 0;
  int solve_code = 0;
  if (interrupted || env_.GetBooleanField(timeout,
      env_.GetFieldID(timeout_class.get(), "timeOutOccurred", "Z"))) {
    solve_code = 600;
    status = "interrupted";
  } else if (!found) {
    solve_code = 200;
    status = "infeasible problem";
  } else if (has_obj) {
    solve_code = 0;
    status = "optimal solution";
  } else {
    solve_code = 100;
    status = "feasible solution";
  }
  p.set_solve_code(solve_code);

  std::vector<double> final_solution;
  if (found) {
    if (has_obj) {
      obj_val = env_.CallIntMethod(obj_var_.get(), value_);
      if (p.obj_type(0) == MAX)
        obj_val = -obj_val;
    }
    const std::vector<jobject> &vars = converter.vars();
    int num_vars = p.num_vars();
    final_solution.resize(num_vars);
    for (int j = 0; j < num_vars; ++j)
      final_solution[j] = env_.CallIntMethod(vars[j], value_);
  }

  double solution_time = GetTimeAndReset(time);

  fmt::Writer w;
  w.Format("{}: {}\n") << long_name() << status;
  w.Format("{} nodes, {} fails")
      << env_.CallIntMethod(search_.get(), get_nodes_)
      << env_.CallIntMethod(search_.get(), get_fails_);
  if (has_obj && found)
    w.Format(", objective {}") << ObjPrec(obj_val);
  BasicSolver::DoHandleSolution(p, w.c_str(),
      final_solution.empty() ? 0 : &final_solution[0], 0, obj_val);

  double output_time = GetTimeAndReset(time);

  if (timing()) {
    Print("Setup time = {:.6f}s\n"
          "Solution time = {:.6f}s\n"
          "Output time = {:.6f}s\n")
            << setup_time << solution_time << output_time;
  }
}
}
