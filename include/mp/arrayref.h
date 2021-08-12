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

/// A reference to an immutable array which can be
/// stored inside if ArrayRef<> is constructed
/// from an rvalue std::vector.
template <typename T>
class ArrayRef {
  std::vector<T> save_;
  const T *data_ = nullptr;
  std::size_t size_ = 0;

 public:
  ArrayRef() { }

  ArrayRef(const T *data, std::size_t size) : data_(data), size_(size) {}

  template <std::size_t SIZE>
  ArrayRef(const T (&data)[SIZE]) : data_(data), size_(SIZE) {}

  /// Rvalue ArrayRef, take over vector if any
  template <typename U>
  ArrayRef(ArrayRef<U>&& other) { init_from_rvalue(std::move(other)); }

  /// Lvalue ArrayRef, pure reference
  template <typename U>
  ArrayRef(const ArrayRef<U>& other) : data_(other.data()), size_(other.size()) {}

  /// Rvalue std::vector, take over
  template <typename TT>
  ArrayRef(std::vector<TT> &&other) :
    save_(std::move(other)), data_(save_.data()), size_(save_.size()) {}

  /// Lvalue Vector, pure reference
  template <typename Vector>
  ArrayRef(const Vector &other) : data_(other.data()), size_(other.size()) {}

  /// = Rvalue, take over vector if any
  template <class U>
  ArrayRef<T>& operator=(ArrayRef<U>&& other) {
    init_from_rvalue(std::move(other));
    return *this;
  }

  /// = Lvalue, pure reference
  template <class U>
  ArrayRef<T>& operator=(const ArrayRef<U>& other) {
    data_ = other.data();
    size_ = other.size();
    return *this;
  }

  operator bool() const { return !empty(); }
  bool empty() const { return 0==size(); }

  const T *data() const { return data_; }
  std::size_t size() const { return size_; }

  const T &operator[](std::size_t i) const { return data_[i]; }

protected:
  template <class AR>
  void init_from_rvalue(AR&& other) {
    if (other.save_.size()) {
      save_ = std::move(other.save_);
      data_ = save_.data();
      size_ = save_.size();
    } else {
      data_ = other.data();
      size_ = other.size();
    }
  }
};

template <typename T>
ArrayRef<T> MakeArrayRef(const T *data, std::size_t size) {
  return ArrayRef<T>(data, size);
}
}  // namespace mp

#endif  // MP_ARRAYREF_H_
