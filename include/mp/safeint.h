/*
 An integer with overflow check.

 Copyright (c) 2013, Victor Zverovich
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

#ifndef MP_SAFEINT_H_
#define MP_SAFEINT_H_

#include <exception>
#include <limits>

namespace mp {

class OverflowError : public std::exception {
 public:
  const char *what() const throw() { return "integer overflow"; }
};

template <typename T>
struct MakeUnsigned { typedef T Type; };

template <>
struct MakeUnsigned<short> { typedef unsigned short Type; };

template <>
struct MakeUnsigned<int> { typedef unsigned Type; };

template <>
struct MakeUnsigned<long> { typedef unsigned long Type; };

// Safe std::abs replacement. Unlike std::abs, SafeAbs doesn't result
// in undefined behavior for std::numeric_limits<T>::min().
template <typename T>
inline typename MakeUnsigned<T>::Type SafeAbs(T value) {
  typename MakeUnsigned<T>::Type result = value;
  if (value < 0)
    result = 0 - result;
  return result;
}

template <typename T>
class SafeInt {
 private:
  T value_;

 public:
  SafeInt(T value) : value_(value) {}

  template <typename U>
  explicit SafeInt(U value) : value_(value) {
    if ((std::numeric_limits<U>::is_signed &&
         value < std::numeric_limits<T>::min()) ||
        value > std::numeric_limits<T>::max()) {
      throw OverflowError();
    }
  }
  
  T value() const { return value_; }
};

template <typename T>
inline SafeInt<T> operator+(SafeInt<T> a, SafeInt<T> b) {
  T a_value = a.value(), b_value = b.value();
  if (a_value >= 0) {
    if (b_value > std::numeric_limits<T>::max() - a_value)
      throw OverflowError();
  } else if (b_value < std::numeric_limits<T>::min() - a_value)
    throw OverflowError();
  return a_value + b_value;
}

template <typename T1, typename T2>
inline SafeInt<T1> operator+(SafeInt<T1> a, T2 b) { return a + SafeInt<T1>(b); }

template <typename T1, typename T2>
inline SafeInt<T2> operator+(T1 a, SafeInt<T2> b) { return SafeInt<T2>(a) + b; }

template <typename T>
inline SafeInt<T> operator-(SafeInt<T> a, SafeInt<T> b) {
  T a_value = a.value(), b_value = b.value();
  if (a_value >= 0) {
    if (b_value < a_value - std::numeric_limits<T>::max())
      throw OverflowError();
  } else if (b_value > a_value - std::numeric_limits<T>::min())
    throw OverflowError();
  return a_value - b_value;
}

template <typename T1, typename T2>
inline SafeInt<T1> operator-(SafeInt<T1> a, T2 b) { return a - SafeInt<T1>(b); }

template <typename T1, typename T2>
inline SafeInt<T2> operator-(T1 a, SafeInt<T2> b) { return SafeInt<T2>(a) - b; }

template <typename T>
inline SafeInt<T> operator*(SafeInt<T> a, SafeInt<T> b) {
  T a_value = a.value(), b_value = b.value();
  if (b_value != 0 &&
      SafeAbs(a_value) > std::numeric_limits<T>::max() / SafeAbs(b_value)) {
    throw OverflowError();
  }
  return a_value * b_value;
}

template <typename T1, typename T2>
inline SafeInt<T1> operator*(SafeInt<T1> a, T2 b) { return a * SafeInt<T1>(b); }

template <typename T1, typename T2>
inline SafeInt<T2> operator*(T1 a, SafeInt<T2> b) { return SafeInt<T2>(a) * b; }

}  // namespace mp

#endif  // MP_SAFEINT_H_
