/*
 SafeInt tests.

 Copyright (c) 2012, Victor Zverovich
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice, this
    list of conditions and the following disclaimer.
 2. Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <gtest/gtest.h>

#include <climits>
#include "solvers/util/safeint.h"

TEST(SafeIntTest, Ctor) {
  SafeInt<int> n(42);
  EXPECT_EQ(42, n.value());
  EXPECT_THROW(SafeInt<int>(INT_MAX + 1u), OverflowError);
  EXPECT_EQ(42, SafeInt<int>(42u).value());
}

TEST(SafeIntTest, Add) {
  EXPECT_EQ(42, (SafeInt<int>(40) + SafeInt<int>(2)).value());
  EXPECT_EQ(42, (SafeInt<int>(40) + 2).value());
  EXPECT_EQ(42, (40 + SafeInt<int>(2)).value());
  EXPECT_EQ(INT_MAX, (SafeInt<int>(INT_MAX - 1) + 1).value());
}

TEST(SafeIntTest, AddOverflow) {
  EXPECT_THROW(SafeInt<int>(1) + INT_MAX, OverflowError);
  EXPECT_THROW(SafeInt<int>(-1) + INT_MIN, OverflowError);
  EXPECT_THROW(SafeInt<int>(INT_MAX) + 1, OverflowError);
  EXPECT_THROW(SafeInt<int>(INT_MAX / 2) + (INT_MAX / 2 + 2), OverflowError);
  EXPECT_THROW(INT_MAX + SafeInt<int>(INT_MAX), OverflowError);
  EXPECT_THROW(SafeInt<int>(INT_MIN) + -1, OverflowError);
  EXPECT_THROW(SafeInt<int>(INT_MAX) + INT_MAX, OverflowError);
}

TEST(SafeIntTest, Sub) {
  EXPECT_EQ(42, (SafeInt<int>(44) - SafeInt<int>(2)).value());
  EXPECT_EQ(42, (SafeInt<int>(44) - 2).value());
  EXPECT_EQ(42, (44 - SafeInt<int>(2)).value());
  EXPECT_EQ(INT_MIN, (SafeInt<int>(INT_MIN + 1) - 1).value());
  EXPECT_EQ(INT_MAX, (SafeInt<int>(-1) - INT_MIN).value());
}

TEST(SafeIntTest, SubOverflow) {
  EXPECT_THROW(SafeInt<int>(0) - INT_MIN, OverflowError);
  EXPECT_THROW(SafeInt<int>(INT_MIN) - 1, OverflowError);
  EXPECT_THROW(SafeInt<int>(-2) - INT_MAX, OverflowError);
  EXPECT_THROW(SafeInt<int>(100) - (INT_MIN + 10), OverflowError);
}

TEST(SafeIntTest, Mul) {
  EXPECT_EQ(42, (SafeInt<int>(6) * SafeInt<int>(7)).value());
  EXPECT_EQ(42, (SafeInt<int>(6) * 7).value());
  EXPECT_EQ(42, (6 * SafeInt<int>(7)).value());
  EXPECT_EQ(INT_MAX & ~1, (SafeInt<int>(INT_MAX / 2) * 2).value());
}

TEST(SafeIntTest, MulOverflow) {
  EXPECT_THROW(SafeInt<int>(INT_MAX / 2) * 3, OverflowError);
  EXPECT_THROW(SafeInt<int>(101) * (INT_MAX / 100), OverflowError);
}
