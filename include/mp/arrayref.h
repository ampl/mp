/*
 Array reference class.

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

#ifndef MP_ARRAYREF_H_
#define MP_ARRAYREF_H_

#include <cstddef>  // for std::size_t

namespace mp {

// A reference to an immutable array.
template <typename T>
class ArrayRef {
 private:
  const T *data_;
  std::size_t size_;

 public:
  ArrayRef(const T *data, std::size_t size) : data_(data), size_(size) {}

  template <typename U>
  ArrayRef(ArrayRef<U> other) : data_(other.data()), size_(other.size()) {}

  template <typename Vector>
  ArrayRef(const Vector &other) : data_(other.data()), size_(other.size()) {}

  template <std::size_t SIZE>
  ArrayRef(const T (&data)[SIZE]) : data_(data), size_(SIZE) {}

  const T *data() const { return data_; }
  std::size_t size() const { return size_; }

  const T &operator[](std::size_t i) const { return data_[i]; }
};

template <typename T>
ArrayRef<T> MakeArrayRef(const T *data, std::size_t size) {
  return ArrayRef<T>(data, size);
}
}  // namespace mp

#endif  // MP_ARRAYREF_H_
