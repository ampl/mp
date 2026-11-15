// Tests of mp::NLSolver's temporary files and solver command line.

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>

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
