/*
 LocalSolver tests

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

#include "localsolver/localsolver.h"
#include "feature.h"

#include <gmock/gmock.h>

using testing::ElementsAre;

typedef mp::LocalSolver Solver;
unsigned FEATURES = localsolver::LSVersion::getMajorVersionNumber() < 5 ?
      ~feature::PLTERM : feature::ALL;

#define MP_THREAD_INTERRUPT 1

// Demo version of LocalSolver cannot handle a TSP problem with n > 9.
#define MP_TSP_SIZE 9

#include "nl-solver-test.h"

using localsolver::LSParam;
using localsolver::LSPhase;

int ls_version() { return localsolver::LSVersion::getMajorVersionNumber(); }

// Creates a test LocalSolver model.
void CreateTestModel(localsolver::LocalSolver &s) {
  s.getParam().setVerbosity(0);
  localsolver::LSModel model = s.getModel();
  model.addObjective(model.createConstant(localsolver::lsint(0)),
                     localsolver::OD_Minimize);
  model.close();
  s.createPhase();
}

// Basic LocalSolver option
template <typename T>
struct BasicOption {
  const char *name;
  T default_value;
  T lb;
  T ub;
  bool check_default;

  BasicOption(const char *name, T default_value,
              T lb, T ub, bool check_default = true)
    : name(name), default_value(default_value),
      lb(lb), ub(ub), check_default(check_default) {}
  virtual ~BasicOption() {}

  typedef localsolver::LocalSolver Solver;

  virtual T get(Solver &s) const = 0;
  virtual void set(Solver &s, T value) const = 0;
};

template <typename T>
T ConvertOptionValue(T value) { return value; }

std::string ConvertOptionValue(const fmt::StringRef &value) {
  return value.to_string();
}

int ConvertOptionValue(fmt::LongLong value) {
  return static_cast<int>(value);
}

// LocalSolver option stored in LSParam
template <typename T, typename ParamT = T>
struct Option : public BasicOption<T> {
  typedef ParamT (LSParam::*Getter)() const;
  Getter get_;
  typedef void (LSParam::*Setter)(ParamT);
  Setter set_;

  explicit Option(Getter get)
    : BasicOption<T>(0, T(), T(), T()), get_(get), set_(0) {}

  Option(const char *name, T default_value, T lb, T ub, Getter get, Setter set)
    : BasicOption<T>(name, default_value, lb, ub), get_(get), set_(set) {}

  typedef localsolver::LocalSolver Solver;

  T get(Solver &s) const { return (s.getParam().*get_)(); }
  void set(Solver &s, T value) const {
    (s.getParam().*set_)(ConvertOptionValue(value));
  }
};

typedef Option<fmt::LongLong, int> IntOption;

template <typename T = fmt::LongLong>
struct TestLocalSolver : mp::LocalSolver {
  const BasicOption<T> *opt;
  T value;
  void DoSolve(localsolver::LocalSolver &s) {
    value = opt->get(s);
  }
  TestLocalSolver() : opt(), value() {}
};

// Creates and solves a test problem.
template <typename T>
T GetOption(TestLocalSolver<T> &solver, const BasicOption<T> &opt) {
  mp::LSProblemBuilder pb(solver);
  MakeAllDiffProblem(pb);
  TestSolutionHandler sh;
  solver.opt = &opt;
  solver.Solve(pb, sh);
  return solver.value;
}

class OptionTest :
    public testing::TestWithParam<BasicOption<fmt::LongLong>*> {};

// Test that the default value agrees with LocalSolver.
TEST_P(OptionTest, DefaultValue) {
  auto &opt = *GetParam();
  if (!opt.check_default)
    return;
  localsolver::LocalSolver ls;
  CreateTestModel(ls);
  EXPECT_EQ(opt.get(ls), mp::LocalSolver().GetIntOption(opt.name));
}

// Test that the range agrees with LocalSolver.
TEST_P(OptionTest, Range) {
  auto &opt = *GetParam();
  localsolver::LocalSolver ls;
  CreateTestModel(ls);
  mp::LocalSolver solver;
  opt.set(ls, opt.lb);
  solver.SetIntOption(opt.name, opt.lb);
  opt.set(ls, opt.ub);
  solver.SetIntOption(opt.name, opt.ub);
  EXPECT_THROW(opt.set(ls, opt.lb - 1), localsolver::LSException);
  EXPECT_THROW(solver.SetIntOption(opt.name, opt.lb - 1),
               mp::InvalidOptionValue);
  if (opt.ub != INT_MAX) {
    EXPECT_THROW(opt.set(ls, opt.ub + 1), localsolver::LSException);
    EXPECT_THROW(solver.SetIntOption(opt.name, opt.ub + 1),
                 mp::InvalidOptionValue);
  }
}

// Test that the option value is passed to LocalSolver.
TEST_P(OptionTest, PassValue) {
  auto &opt = *GetParam();
  TestLocalSolver<> solver;
  enum {TEST_VALUE = 7};
  solver.SetIntOption(opt.name, TEST_VALUE);
  EXPECT_EQ(TEST_VALUE, solver.GetIntOption(opt.name));
  EXPECT_EQ(TEST_VALUE, GetOption(solver, opt));
  EXPECT_THROW(solver.SetStrOption(opt.name, "oops"), mp::OptionError);
}

IntOption seed("seed", 0, 0, INT_MAX, &LSParam::getSeed, &LSParam::setSeed);
IntOption threads(
    "threads", 2, 1, 1024, &LSParam::getNbThreads, &LSParam::setNbThreads);
IntOption annealing_level(
    "annealing_level", 1, 0, 9,
    &LSParam::getAnnealingLevel, &LSParam::setAnnealingLevel);
IntOption time_between_displays(
    "time_between_displays", 1, 1, 65535,
    &LSParam::getTimeBetweenDisplays, &LSParam::setTimeBetweenDisplays);
INSTANTIATE_TEST_CASE_P(, OptionTest, testing::Values(
                          &seed, &annealing_level, &time_between_displays));

TEST(LocalSolverTest, VerbosityOption) {
  TestLocalSolver<> solver;

  // Test default value.
  EXPECT_EQ("quiet", solver.GetStrOption("verbosity"));
  mp::LSProblemBuilder builder(solver);
  EXPECT_EQ(0, builder.solver().getParam().getVerbosity());

  // Test option values.
  auto values = solver.FindOption("verbosity")->values();
  EXPECT_EQ(4, values.size());
  struct ValueInfo {
    const char *value;
    int data;
  };
  int detailed_verbosity = ls_version() >= 5 ? 2 : 10;
  ValueInfo value_info[] = {
    {"terse", -1}, {"quiet", 0}, {"normal", 1}, {"detailed", detailed_verbosity}
  };
  int index = 0;
  for (auto i = values.begin(), end = values.end(); i != end; ++i) {
    const auto &info = value_info[index++];
    EXPECT_STREQ(info.value, i->value);
    EXPECT_EQ(info.data, i->data);
    solver.SetStrOption("verbosity", info.value);
    EXPECT_EQ(info.value, solver.GetStrOption("verbosity"));
    solver.SetStrOption("verbosity", fmt::format("{}", info.data));
    EXPECT_EQ(info.value, solver.GetStrOption("verbosity"));
  }

  // Test that value is passed to LocalSolver.
  solver.SetStrOption("verbosity", "detailed");
  EXPECT_EQ(detailed_verbosity,
            GetOption(solver, IntOption(&LSParam::getVerbosity)));

  // Terse verbosity triggers custom output and the LocalSolver output
  // is disabled (verbosity = 0).
  solver.SetStrOption("verbosity", "terse");
  EXPECT_EQ(0, GetOption(solver, IntOption(&LSParam::getVerbosity)));

  EXPECT_THROW(solver.SetDblOption("verbosity", 1.2), mp::OptionError);
  EXPECT_THROW(solver.SetStrOption("verbosity", "oops"),
               mp::InvalidOptionValue);
}

TEST(LocalSolverTest, LogFileOption) {
  TestLocalSolver<std::string> solver;
  EXPECT_EQ("", solver.GetStrOption("logfile"));
  solver.SetStrOption("logfile", "testlog");
  EXPECT_EQ("testlog",
            GetOption(solver, Option<std::string>(&LSParam::getLogFile)));
  EXPECT_THROW(solver.SetIntOption("logfile", 1), mp::OptionError);
}

TEST(LocalSolverTest, EnvFileOption) {
  mp::LocalSolver solver;
  EXPECT_EQ("", solver.GetStrOption("envfile"));
  const char *envfile = "test.lsm";
  solver.SetStrOption("envfile", envfile);
  EXPECT_EQ(envfile, solver.GetStrOption("envfile"));
  EXPECT_THROW(solver.SetIntOption("envfile", 1), mp::OptionError);

  mp::LSProblemBuilder pb(solver);
  pb.AddVar(0, 1, var::INTEGER);
  pb.AddObj(obj::MIN, pb.MakeVariable(0));
  TestSolutionHandler sh;
  std::remove(envfile);
  solver.Solve(pb, sh);
  // Check if file exists by opening it.
  fmt::File(envfile, fmt::File::RDONLY);
}

// LocalSolver option stored in LSPhase
template <typename T>
struct PhaseOption : public BasicOption<fmt::LongLong> {
  typedef T (LSPhase::*Getter)() const;
  Getter get_;
  typedef void (LSPhase::*Setter)(T);
  Setter set_;

  PhaseOption(const char *name, fmt::LongLong default_value,
              int lb, int ub, Getter get, Setter set, bool check_default = true)
    : BasicOption<fmt::LongLong>(name, default_value, lb, ub, check_default),
      get_(get), set_(set) {}

  typedef localsolver::LocalSolver Solver;

  fmt::LongLong get(Solver &s) const { return (s.getPhase(0).*get_)(); }
  void set(Solver &s, fmt::LongLong value) const {
    (s.getPhase(0).*set_)(static_cast<T>(value));
  }
};

PhaseOption<int> timelimit(
    "timelimit", 10, ls_version() >= 5 ? 0 : 1, INT_MAX,
    &LSPhase::getTimeLimit, &LSPhase::setTimeLimit, false);
const fmt::LongLong long_long_max = std::numeric_limits<fmt::LongLong>::max();
PhaseOption<fmt::LongLong> iterlimit(
    "iterlimit", long_long_max, ls_version() >= 5 ? 0 : 1, INT_MAX,
    &LSPhase::getIterationLimit, &LSPhase::setIterationLimit);
INSTANTIATE_TEST_CASE_P(
    Limit, OptionTest, testing::Values(&timelimit, &iterlimit));

template <typename SuffixInfo>
class SuffixTest : public ::testing::Test{};

struct IntBoundInfo {
  typedef localsolver::lsint Type;

  static var::Type type() { return var::INTEGER; }
  static int value() { return 42; }
  static Type get(const localsolver::LSParam &p) {
    return p.getIntObjectiveThreshold(0);
  }
  static mp::LSProblemBuilder::IntSuffixHandler AddSuffix(
      mp::LSProblemBuilder &pb) {
    return pb.AddIntSuffix("bound", mp::suf::OBJ, 1);
  }
};

struct DblBoundInfo {
  typedef double Type;

  static var::Type type() { return var::CONTINUOUS; }
  static Type value() { return 4.2; }
  static Type get(const localsolver::LSParam &p) {
    return p.getDoubleObjectiveThreshold(0);
  }
  static mp::LSProblemBuilder::DblSuffixHandler AddSuffix(
      mp::LSProblemBuilder &pb) {
    return pb.AddDblSuffix("bound", mp::suf::OBJ, 1);
  }
};

typedef ::testing::Types<IntBoundInfo, DblBoundInfo> SuffixTestArgs;
TYPED_TEST_CASE(SuffixTest, SuffixTestArgs);

TYPED_TEST(SuffixTest, ObjBound) {
  struct TestSolver : mp::LocalSolver {
    typename TypeParam::Type bound;
    TestSolver() : bound(0) {}
    void DoSolve(localsolver::LocalSolver &s) {
      bound = TypeParam::get(s.getParam());
    }
  };
  TestSolver solver;
  mp::LSProblemBuilder pb(solver);
  pb.AddVar(0, 1, TypeParam::type());
  pb.AddObj(obj::MIN, pb.MakeVariable(0));
  TypeParam::AddSuffix(pb).SetValue(0, TypeParam::value());
  TestSolutionHandler sh;
  solver.Solve(pb.problem(), sh);
  EXPECT_EQ(TypeParam::value(), solver.bound);
}

std::vector<double> LSArrayToVector(localsolver::LSExpression array) {
  std::vector<double> result;
  for (int i = 0, n = array.getNbOperands(); i < n; ++i)
    result.push_back(array.getOperand(i).getDoubleValue());
  return result;
}

class PLTermTest : public ::testing::Test {
 protected:
  mp::LocalSolver solver;
  mp::LSProblemBuilder pb;

  PLTermTest() : pb(solver) {}

  localsolver::LSExpression AddPLTerm(
      mp::LSProblemBuilder &pb, double lb, double ub, double bp = 0) {
    pb.AddVar(lb, ub, var::CONTINUOUS);
    auto pl_builder = pb.BeginPLTerm(1);
    pl_builder.AddSlope(-1);
    pl_builder.AddBreakpoint(bp);
    pl_builder.AddSlope(1);
    return pb.EndPLTerm(pl_builder, pb.MakeVariable(0));
  }

  localsolver::LSExpression AddPLTerm(double lb, double ub, double bp = 0) {
    return AddPLTerm(pb, lb, ub, bp);
  }
};

TEST_F(PLTermTest, Bounds) {
  if ((FEATURES & feature::PLTERM) == 0) return;
  auto plterm = AddPLTerm(-1.1, 22.2);
  ASSERT_THAT(LSArrayToVector(plterm.getOperand(0)),
              ElementsAre(-1.1, 0, 22.2));
  ASSERT_THAT(LSArrayToVector(plterm.getOperand(1)),
              ElementsAre(1.1, 0, 22.2));
}

TEST_F(PLTermTest, BoundsBelowBreakpoints) {
  if ((FEATURES & feature::PLTERM) == 0) return;
  auto plterm = AddPLTerm(-20, -10);
  ASSERT_THAT(LSArrayToVector(plterm.getOperand(0)), ElementsAre(-20, -10));
  ASSERT_THAT(LSArrayToVector(plterm.getOperand(1)), ElementsAre(20, 10));
}

TEST_F(PLTermTest, BoundsAboveBreakpoints) {
  if ((FEATURES & feature::PLTERM) == 0) return;
  auto plterm = AddPLTerm(10, 20);
  ASSERT_THAT(LSArrayToVector(plterm.getOperand(0)), ElementsAre(10, 20));
  ASSERT_THAT(LSArrayToVector(plterm.getOperand(1)), ElementsAre(10, 20));
}

TEST_F(PLTermTest, InfiniteLB) {
  if ((FEATURES & feature::PLTERM) == 0) return;
  auto plterm = AddPLTerm(-std::numeric_limits<double>::infinity(), 10);
  ASSERT_THAT(LSArrayToVector(plterm.getOperand(0)), ElementsAre(-1e6, 0, 10));
  ASSERT_THAT(LSArrayToVector(plterm.getOperand(1)), ElementsAre(1e6, 0, 10));
}

TEST_F(PLTermTest, InfiniteUB) {
  if ((FEATURES & feature::PLTERM) == 0) return;
  auto plterm = AddPLTerm(-10, std::numeric_limits<double>::infinity());
  ASSERT_THAT(LSArrayToVector(plterm.getOperand(0)), ElementsAre(-10, 0, 1e6));
  ASSERT_THAT(LSArrayToVector(plterm.getOperand(1)), ElementsAre(10, 0, 1e6));
}

TEST_F(PLTermTest, BoundsFromBreakpoint) {
  if ((FEATURES & feature::PLTERM) == 0) return;
  double bp = -1234321, inf = std::numeric_limits<double>::infinity();
  auto plterm = AddPLTerm(-inf, inf, bp);
  ASSERT_THAT(LSArrayToVector(plterm.getOperand(0)),
              ElementsAre(2 * bp, bp, -2 * bp));
  ASSERT_THAT(LSArrayToVector(plterm.getOperand(1)),
              ElementsAre(0, bp, -2 * bp));
}

TEST_F(PLTermTest, BoundsFromLargeBound) {
  if ((FEATURES & feature::PLTERM) == 0) return;
  double bound = -1234321;
  auto plterm = AddPLTerm(bound, std::numeric_limits<double>::infinity(), 0);
  ASSERT_THAT(LSArrayToVector(plterm.getOperand(0)),
              ElementsAre(bound, 0, -bound));
  ASSERT_THAT(LSArrayToVector(plterm.getOperand(1)),
              ElementsAre(-bound, 0, -bound));
}

TEST_F(PLTermTest, PLBigMOption) {
  if ((FEATURES & feature::PLTERM) == 0) return;
  if ((FEATURES & feature::PLTERM) == 0) return;
  EXPECT_EQ(1e6, solver.GetDblOption("pl_bigm"));
  double pl_bigm = 100;
  solver.SetDblOption("pl_bigm", pl_bigm);
  EXPECT_EQ(pl_bigm, solver.GetDblOption("pl_bigm"));
  EXPECT_EQ(pl_bigm, solver.pl_bigm());
  mp::LSProblemBuilder pb(solver);
  auto plterm = AddPLTerm(pb, -10, std::numeric_limits<double>::infinity());
  ASSERT_THAT(LSArrayToVector(plterm.getOperand(0)),
              ElementsAre(-10, 0, pl_bigm));
  ASSERT_THAT(LSArrayToVector(plterm.getOperand(1)),
              ElementsAre(10, 0, pl_bigm));
}
