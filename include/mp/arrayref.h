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

#include <vector>
#include <cstddef>  // for std::size_t

namespace mp {

/// A reference to a mutable array which can be
/// stored inside if MutArrayRef<> is constructed
/// from an rvalue `~std::vector`.
///
/// @see `~mp::ArrayRef`
template <typename T>
class MutArrayRef {
public:
  /// Typedef internally stored type.
  /// std::vector<const T> is illegal?
  using TS = typename std::remove_const<T>::type;

  /// Construct
  MutArrayRef() { }

  /// From pointer + size
  MutArrayRef(T *data, std::size_t size) noexcept :
    data_(data), size_(size) {}

  /// From C-array
  template <std::size_t SIZE>
  MutArrayRef(T (&data)[SIZE]) noexcept :
    data_(data), size_(SIZE) {}

  /// Rvalue ArrayRef, take over stored vector if any
  MutArrayRef(MutArrayRef&& other) noexcept
  { init_from_rvalue(std::move(other)); }

  /// Lvalue ArrayRef, pure reference
  MutArrayRef(const MutArrayRef& other) noexcept :
    data_(other.data()), size_(other.size()) {}

  /// Rvalue std::vector, take over
  MutArrayRef(std::vector<TS> &&other) noexcept :
    save_(std::move(other)),
    data_(save_.data()), size_(save_.size()) {}

  /// Lvalue std::vector, pure reference
  MutArrayRef(std::vector<TS> &other) noexcept :
    data_(other.data()), size_(other.size()) {}

  /// = Rvalue, take over vector if any
  MutArrayRef& operator=(MutArrayRef&& other) noexcept {
    init_from_rvalue(std::move(other));
    return *this;
  }

  /// = Lvalue, pure reference
  MutArrayRef& operator=(const MutArrayRef& other) {
    data_ = other.data();
    size_ = other.size();
    return *this;
  }

  // operator bool()
  operator bool() const { return !empty(); }
  // bool empty()
  bool empty() const { return 0==size(); }

  operator std::vector<TS>() & { return {begin(), end()}; }
  operator std::vector<TS>() && { return move_or_copy(); }

  const T *data() const { return data_; }
  T *data() { return data_; }
  std::size_t size() const { return size_; }

  const T* begin() const { return data(); }
  const T* end() const { return data()+size(); }
  T* begin() { return data(); }
  T* end() { return data()+size(); }

  const T &operator[](std::size_t i) const { return data_[i]; }

protected:
  std::vector<TS> move_or_copy() {
    if (save_.size()) {
      data_ = nullptr;
      size_ = 0;
      return std::move(save_);
    }
    return {begin(), end()};
  }

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

private:
  std::vector<TS> save_;
  T *data_ = nullptr;
  std::size_t size_ = 0;
};


/// A reference to an immutable array which can be
/// stored inside if ArrayRef<> is constructed
/// from an rvalue `~std::vector`.
///
/// @see `~mp::MutArrayRef`
///
/// @note `~mp::ArrayRef<T>` inherits from `~mp::MutArrayRef<T>`,
///   so easy to remove immutability.
template <typename T>
class ArrayRef : public MutArrayRef<const T> {
public:
  /// Typedef base class
  using Base = MutArrayRef<const T>;

  /// Construct
  ArrayRef() : Base() { }

  /// From pointer + size
  ArrayRef(const T *data, std::size_t size) noexcept :
      Base(data, size) {}

  /// From C-array
  template <std::size_t SIZE>
  ArrayRef(const T (&data)[SIZE]) noexcept :
      Base(data) {}

  /// Rvalue ArrayRef, take over stored vector if any
  ArrayRef(ArrayRef&& other) noexcept :
      Base(other) { }

  /// Lvalue ArrayRef, pure reference
  ArrayRef(const ArrayRef& other) noexcept :
      Base(other) {}

  /// Rvalue std::vector, take over
  ArrayRef(std::vector<T> &&other) noexcept :
      Base(std::move(other)) {}

  /// Lvalue std::vector, pure reference
  ArrayRef(const std::vector<T> &other) noexcept :
      Base(other.data(), other.size()) {}

  /// = Rvalue, take over vector if any
  ArrayRef& operator=(ArrayRef&& other) noexcept {
    Base::operator=(std::move(other));
    return *this;
  }

  /// = Lvalue, pure reference
  ArrayRef& operator=(const ArrayRef& other) {
    Base::operator=(other);
    return *this;
  }

  // reuse stuff
  const T *data() const { return Base::data(); }
  const T* begin() const { return Base::begin(); }
  const T* end() const { return Base::end(); }
  using Base::operator std::vector<T>;
};

template <typename T> inline
ArrayRef<T> MakeArrayRef(const T *data, std::size_t size) {
  return ArrayRef<T>(data, size);
}

/// std::vector::data() might not return nullptr when empty
template <class Vec> inline
auto data_or_null(const Vec& v) -> decltype(v.data()) {
  return v.empty() ? nullptr : v.data();
}

}  // namespace mp

#endif  // MP_ARRAYREF_H_
