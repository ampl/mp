// Tests of mp::NLSolver's temporary files and solver command line.

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef _WIN32
#include <unistd.h>    // geteuid
#endif

#include "gtest/gtest.h"
#include "mp/nl-solver.h"

namespace fs = std::filesystem;

namespace {

/// What NLUtils::myexit() raises here instead of ending the test.
struct ExitCalled : std::runtime_error {
  using std::runtime_error::runtime_error;
};

class ThrowingUtils : public mp::NLUtils {
public:
  void myexit(const std::string& msg) override { throw ExitCalled(msg); }
};

#ifdef _WIN32
const char* const kTempVars[] = {"TMP", "TEMP"};
#else
const char* const kTempVars[] = {"TMPDIR"};
#endif

void SetEnv(const char* name, const char* value) {
#ifdef _WIN32
  _putenv_s(name, value ? value : "");
#else
  if (value)
    setenv(name, value, 1);
  else
    unsetenv(name);
#endif
}

/// Points the temporary directory at \a dir, restores it on exit.
class TempDirOverride {
public:
  explicit TempDirOverride(const fs::path& dir) {
    for (int i = 0; i < kNumVars; ++i) {
      auto old = std::getenv(kTempVars[i]);
      had_[i] = old != nullptr;
      old_[i] = old ? old : "";
      SetEnv(kTempVars[i], dir.string().c_str());
    }
  }
  ~TempDirOverride() {
    for (int i = 0; i < kNumVars; ++i)
      SetEnv(kTempVars[i], had_[i] ? old_[i].c_str() : nullptr);
  }

private:
  static constexpr int kNumVars = int(std::size(kTempVars));
  bool had_[kNumVars] {};
  std::string old_[kNumVars];
};

/// A fresh directory under the system temporary directory,
/// removed on exit.
class ScratchDir {
public:
  explicit ScratchDir(const std::string& name)
    : path_(fs::temp_directory_path() / name) {
    fs::remove_all(path_);
    fs::create_directories(path_);
  }
  ~ScratchDir() {
    std::error_code ec;
    fs::remove_all(path_, ec);
  }
  const fs::path& path() const { return path_; }

private:
  fs::path path_;
};

/// minimize x, 0 <= x <= 1.
class OneVarModel {
public:
  mp::NLModel Get() {
    mp::NLModel mdl("one variable");
    mdl.SetCols({1, lb_, ub_, type_});
    mdl.SetLinearObjective(NLW2_ObjSenseMinimize, 0.0, c_);
    return mdl;
  }

private:
  double lb_[1] {0.0};
  double ub_[1] {1.0};
  int type_[1] {0};
  double c_[1] {1.0};
};

std::string Solver() { return mp::QuoteCommandArg(NLSOL_FAKE_SOLVER); }

std::string ReadFile(const fs::path& path) {
  std::ifstream in(path, std::ios::binary);
  return {std::istreambuf_iterator<char>(in), {}};
}

/// Directory name with characters that a shell would interpret.
std::string ShellCharsName() {
#ifdef _WIN32
  return "nlw2 test dir & (x) ^y";
#else
  return "nlw2 test dir it's $HOME `x` \"q\"";
#endif
}

}  // namespace

TEST(NLSolverQuoteTest, QuotesArgument) {
#ifdef _WIN32
  EXPECT_EQ("\"C:\\a b\\c\"", mp::QuoteCommandArg("C:\\a b\\c"));
  EXPECT_EQ("\"C:\\a b\\\\\"", mp::QuoteCommandArg("C:\\a b\\"));
  EXPECT_EQ("\"\"", mp::QuoteCommandArg(""));
  EXPECT_EQ("", mp::QuoteCommandArg("a\"b"));
#else
  EXPECT_EQ("'/a b/c'", mp::QuoteCommandArg("/a b/c"));
  EXPECT_EQ("'it'\\''s'", mp::QuoteCommandArg("it's"));
  EXPECT_EQ("'$HOME `x` \"q\"'", mp::QuoteCommandArg("$HOME `x` \"q\""));
  EXPECT_EQ("''", mp::QuoteCommandArg(""));
#endif
}

TEST(NLSolverTest, CustomStubWithShellCharacters) {
  ScratchDir dir(ShellCharsName());
  auto stub = (dir.path() / "model").string();
  std::ofstream(stub + ".nl") << "not read\n";
  ThrowingUtils utils;
  mp::NLSolver nlsol(&utils);
  nlsol.SetFileStub(stub);
  ASSERT_TRUE(nlsol.Solve(Solver(), "")) << nlsol.GetErrorMessage();
  EXPECT_EQ(stub, ReadFile(stub + ".ok"));
}

TEST(NLSolverTest, AutoStubInTempDirWithShellCharacters) {
  ScratchDir dir(ShellCharsName());
  TempDirOverride tmp(dir.path());
  ThrowingUtils utils;
  OneVarModel model;
  mp::NLSolver nlsol(&utils);
  ASSERT_TRUE(nlsol.LoadModel(model.Get())) << nlsol.GetErrorMessage();
  const auto stub = nlsol.GetFileStub();
  EXPECT_EQ(0u, stub.find(dir.path().string()));
  ASSERT_TRUE(nlsol.Solve(Solver(), "")) << nlsol.GetErrorMessage();
  EXPECT_EQ(stub, ReadFile(stub + ".ok"));
}

/// A fresh scratch stub with a model file, for the solver call tests.
struct StubFixture {
  ScratchDir dir {"nlw2 call test"};
  std::string stub {(dir.path() / "model").string()};
  StubFixture() { std::ofstream(stub + ".nl") << "not read\n"; }
};

TEST(NLSolverTest, SolutionOfPreviousSolveIsNotReused) {
  StubFixture f;
  std::ofstream(f.stub + ".sol") << "stale result\n";
  ThrowingUtils utils;
  mp::NLSolver nlsol(&utils);
  nlsol.SetFileStub(f.stub);
  // The solver exits with 0 but writes no solution
  EXPECT_FALSE(nlsol.Solve(Solver(), "nosol"));
  EXPECT_NE(std::string::npos,
            std::string(nlsol.GetErrorMessage()).find("wrote no solution"));
  EXPECT_FALSE(fs::exists(f.stub + ".sol"));
}

TEST(NLSolverTest, CannotRemovePreviousSolution) {
  StubFixture f;
  const fs::path sol(f.stub + ".sol");
  // A nonempty directory cannot be removed as a file, even by an
  // administrator or by a test run with write access to the parent.
  fs::create_directory(sol);
  std::ofstream(sol / "entry") << "keep\n";
  ThrowingUtils utils;
  mp::NLSolver nlsol(&utils);
  nlsol.SetFileStub(f.stub);
  EXPECT_FALSE(nlsol.Solve(Solver(), "nosol_exit7"));
  EXPECT_NE(std::string::npos,
            std::string(nlsol.GetErrorMessage()).find("remove"));
  EXPECT_FALSE(fs::exists(f.stub + ".ok"));
  EXPECT_TRUE(fs::exists(sol / "entry"));
}

TEST(NLSolverTest, FailedCallReportsExitCode) {
  StubFixture f;
  ThrowingUtils utils;
  mp::NLSolver nlsol(&utils);
  nlsol.SetFileStub(f.stub);
  EXPECT_FALSE(nlsol.Solve(Solver(), "nosol_exit7"));
  // The exit code itself, not the raw wait status (1792 for 7)
  EXPECT_NE(std::string::npos,
            std::string(nlsol.GetErrorMessage()).find("(code 7)"));
}

TEST(NLSolverTest, SolutionWrittenBeforeNonzeroExitIsKept) {
  StubFixture f;
  ThrowingUtils utils;
  mp::NLSolver nlsol(&utils);
  nlsol.SetFileStub(f.stub);
  EXPECT_TRUE(nlsol.Solve(Solver(), "exit7"));
  EXPECT_TRUE(fs::exists(f.stub + ".sol"));
  EXPECT_NE(std::string::npos,
            std::string(nlsol.GetErrorMessage()).find("code 7"));
}

TEST(NLSolverTest, NoteOfPreviousSolveIsCleared) {
  StubFixture f;
  ThrowingUtils utils;
  mp::NLSolver nlsol(&utils);
  nlsol.SetFileStub(f.stub);
  ASSERT_TRUE(nlsol.Solve(Solver(), "exit7"));
  ASSERT_NE(std::string::npos,
            std::string(nlsol.GetErrorMessage()).find("code 7"));
  // A clean solve after it must not report the old exit code
  ASSERT_TRUE(nlsol.Solve(Solver(), ""));
  EXPECT_STREQ("", nlsol.GetErrorMessage());
}

TEST(NLSolverTest, SolverOptionsAreNotInterpretedByShell) {
  StubFixture f;
  ThrowingUtils utils;
  mp::NLSolver nlsol(&utils);
  nlsol.SetFileStub(f.stub);
  const auto injected = (f.dir.path() / "INJECTED").string();
  // Each option is one argument; quotes group words, nothing else is special
  ASSERT_TRUE(nlsol.Solve(
      Solver(), "outlev=0;touch  \"" + injected + "\"  \"log=a b.txt\"  $HOME"))
      << nlsol.GetErrorMessage();
  EXPECT_FALSE(fs::exists(injected));
  EXPECT_EQ("outlev=0;touch\n" + injected + "\nlog=a b.txt\n$HOME\n",
            ReadFile(f.stub + ".args"));
}

TEST(NLSolverTest, UnclosedQuoteInSolverOptionsIsAnError) {
  StubFixture f;
  ThrowingUtils utils;
  mp::NLSolver nlsol(&utils);
  nlsol.SetFileStub(f.stub);
  // Everything after the quote would become one argument, and timelim
  // would be lost without a message
  EXPECT_FALSE(nlsol.Solve(Solver(), "logfile=\"run.log timelim=10"));
  EXPECT_NE(std::string::npos,
            std::string(nlsol.GetErrorMessage()).find("unclosed quote"))
      << nlsol.GetErrorMessage();
  EXPECT_FALSE(fs::exists(f.stub + ".args"));   // the solver was not run
}

TEST(NLSolverTest, MissingTempDirGivesMessage) {
  ScratchDir parent("nlw2 test parent");
  TempDirOverride tmp(parent.path() / "missing");
  ThrowingUtils utils;
  OneVarModel model;
  mp::NLSolver nlsol(&utils);
  try {
    nlsol.LoadModel(model.Get());
    FAIL() << "the missing temp directory was not reported";
  } catch (const ExitCalled& e) {
    EXPECT_NE(nullptr, std::strstr(e.what(), "temp dir")) << e.what();
  }
}

namespace {

/// An error handler that returns, as a C API handler may:
/// it records the messages.
class ReturningUtils : public mp::NLUtils {
public:
  void myexit(const std::string& msg) override { msgs_.push_back(msg); }
  const std::vector<std::string>& Messages() const { return msgs_; }

private:
  std::vector<std::string> msgs_;
};

/// Makes \a dir the current directory, restores it on exit.
class CurrentDirOverride {
public:
  explicit CurrentDirOverride(const fs::path& dir)
    : old_(fs::current_path()) { fs::current_path(dir); }
  ~CurrentDirOverride() {
    std::error_code ec;
    fs::current_path(old_, ec);
  }

private:
  fs::path old_;
};

/// With a returning myexit(), a failed temporary directory must not
/// make NLSolver write the NL file elsewhere (before, a relative
/// directory and file appeared in the current directory).
void ExpectNothingWrittenWithTempDir(const fs::path& tmpdir,
                                     const fs::path& cwd) {
  TempDirOverride tmp(tmpdir);
  CurrentDirOverride here(cwd);
  ReturningUtils utils;
  OneVarModel model;
  {
    mp::NLSolver nlsol(&utils);
    EXPECT_FALSE(nlsol.LoadModel(model.Get()));
    EXPECT_NE(nullptr, std::strstr(nlsol.GetErrorMessage(), "temporary"))
        << nlsol.GetErrorMessage();
    EXPECT_FALSE(nlsol.Solve(Solver(), ""));
  }
  ASSERT_EQ(1u, utils.Messages().size());
  EXPECT_NE(std::string::npos, utils.Messages()[0].find("temp dir"))
      << utils.Messages()[0];
  EXPECT_TRUE(fs::is_empty(cwd)) << fs::directory_iterator(cwd)->path();
}

}  // namespace

TEST(NLSolverTest, MissingTempDirWithReturningHandlerWritesNothing) {
  ScratchDir parent("nlw2 test returning");
  const auto cwd = parent.path() / "cwd";
  fs::create_directories(cwd);
  ExpectNothingWrittenWithTempDir(parent.path() / "missing", cwd);
}

#ifndef _WIN32
// mkdtemp() fails in a read-only temporary directory
// (not as root, who may write there anyway)
TEST(NLSolverTest, UnwritableTempDirWithReturningHandlerWritesNothing) {
  if (0 == geteuid())
    GTEST_SKIP() << "root can write to a read-only directory";
  ScratchDir parent("nlw2 test readonly");
  const auto cwd = parent.path() / "cwd", ro = parent.path() / "ro";
  fs::create_directories(cwd);
  fs::create_directories(ro);
  fs::permissions(ro, fs::perms::owner_read | fs::perms::owner_exec);
  ExpectNothingWrittenWithTempDir(ro, cwd);
  fs::permissions(ro, fs::perms::owner_all);   // so that it can be removed
}
#endif

#ifdef _WIN32
// Windows refuses to create a directory whose path exceeds 248
// characters: here the temporary directory itself fits, but
// NLSolver's subdirectory in it does not.
TEST(NLSolverTest, TooLongTempDirGivesMessage) {
  ScratchDir parent("nlw2 test long");
  auto dir = parent.path();
  while (dir.string().size() < 230)
    dir /= std::string(std::min<size_t>(40, 230 - dir.string().size()), 'd');
  fs::create_directories(dir);
  TempDirOverride tmp(dir);
  ThrowingUtils utils;
  OneVarModel model;
  mp::NLSolver nlsol(&utils);
  try {
    nlsol.LoadModel(model.Get());
    FAIL() << "the too long temp directory was not reported";
  } catch (const ExitCalled& e) {
    EXPECT_NE(nullptr, std::strstr(e.what(), "Could not create temp dir"))
        << e.what();
  }
}
#endif
