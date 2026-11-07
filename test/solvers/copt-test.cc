#include <cfloat>
#include <climits>
#include <string>

#include "copt/coptbackend.h"
#include "../gtest-extra.h"

namespace {

class CoptOptionsTest : public ::testing::Test {
 protected:
  mp::CoptBackend backend_;

  void SetUp() override {
    backend_.GetCallbacks() = {};  // Match the driver's AMPLS initialization.
    backend_.InitMetaInfoAndOptions();
  }

  bool HasValue(const char* option, const char* value) {
    for (const auto& item : backend_.GetOption(option)->values())
      if (std::string(item.value) == value)
        return true;
    return false;
  }
};

// Registering and inspecting options must not open a licensed environment.
TEST_F(CoptOptionsTest, VersionedAlgorithmChoices) {
  const bool has_v9 = COPT_VERSION_MAJOR >= 9;
  EXPECT_EQ(has_v9, HasValue("alg:method", "7"));
  EXPECT_EQ(has_v9, HasValue("lp:concurrentmode", "4"));
  EXPECT_EQ(has_v9, HasValue("lp:barorder", "3"));
  EXPECT_EQ(has_v9, HasValue("lp:barorder", "4"));
  EXPECT_EQ(has_v9, HasValue("lp:barorder", "5"));
  EXPECT_EQ(has_v9, HasValue("lp:barorder", "6"));
  EXPECT_EQ(has_v9, HasValue("mip:startmode", "3"));
}

TEST_F(CoptOptionsTest, LegacyChoicesRemainAvailable) {
  for (const char* option : {"alg:method", "lp:concurrentmode",
                             "lp:barorder", "mip:startmode"}) {
    EXPECT_TRUE(HasValue(option, "-1"));
    EXPECT_TRUE(HasValue(option, "1"));
    EXPECT_TRUE(HasValue(option, "2"));
  }
  EXPECT_TRUE(HasValue("alg:method", "6"));
  EXPECT_TRUE(HasValue("lp:concurrentmode", "0"));
}

TEST_F(CoptOptionsTest, ConcurrentMipAvailability) {
#if COPT_VERSION_MAJOR >= 9
  EXPECT_EQ(backend_.GetOption("mip:concurrent"),
            backend_.GetOption("concurrentmip"));
#else
  EXPECT_THROW(backend_.GetOption("mip:concurrent"), mp::OptionError);
  EXPECT_THROW(backend_.GetOption("concurrentmip"), mp::OptionError);
#endif
}

}  // namespace
