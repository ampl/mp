/*
 AMPL suffix support
 Suffixes are values associated with model components.
 See http://www.ampl.com/NEW/suffixes.html

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
 whatsoever resulting from loss of use	, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Victor Zverovich
 */

#ifndef MP_SUFFIX_H_
#define MP_SUFFIX_H_

#include <cstddef>     // for std::size_t
#include <algorithm>   // for std::fill_n
#include <iterator>
#include <set>
#include <string>      // for std::char_traits

#include "mp/common.h"
#include "mp/error.h"  // for MP_ASSERT
#include "mp/format.h"
#include "mp/arrayref.h"

namespace mp {

template <typename T>
class BasicSuffix;

template <typename T>
class BasicMutSuffix;

template <typename Alloc>
class BasicSuffixSet;

/// Suffix table is a newline-separated string
using SuffixTable = std::string;

namespace internal {

class SuffixBase {
 protected:
  struct Impl {
    // Name is stored as a StringRef rather than std::string to avoid
    // dynamic memory allocation when using set::find.
    fmt::StringRef name;
    mutable int kind;
    int num_values;
    union {
      void *values;
      int *int_values;
      double *dbl_values;
    };
    SuffixTable table;

    explicit Impl(fmt::StringRef name, int kind = 0, int num_values = 0,
                  const SuffixTable& tab = {})
      : name(name), kind(kind), num_values(num_values), int_values(0),
    table(tab) {}

    int kind_full() const { return kind; }
    int kind_pure() const { return kind & suf::KIND_MASK; }
  };

  template <typename SuffixType>
  explicit SuffixBase(SuffixType s) : impl_(s.impl()) {}

  void get_value(int index, int &value) const {
    value = impl_->int_values[index];
  }
  void set_value(int index, int value) { impl_->int_values[index] = value; }

  void get_value(int index, double &value) const {
    value = impl_->dbl_values[index];
  }
  void set_value(int index, double value) { impl_->dbl_values[index] = value; }

  const Impl *impl() const { return impl_; }

  explicit SuffixBase(const Impl *impl) : impl_(impl) {}

  // Safe bool type.
  typedef void (SuffixBase::*SafeBool)() const;

 private:
  const Impl *impl_;

  // A member function representing the true value of SafeBool.
  void True() const {}

 public:
  // Constructs a Suffix object representing a null reference to a
  // suffix. The only operation permitted for such object is copying,
  // assignment and check whether it is null using operator SafeBool.
  SuffixBase() : impl_() {}

  // Returns the suffix name.
  const char *name() const { return impl_->name.data(); }

  // Returns the suffix pure kind (var/con/prob/obj).
  int kind_pure() const { return impl_->kind_pure(); }

  // Returns the suffix kind.
  int kind() const { return impl_->kind_full(); }

  /// Or's the kind with a given int argument
  void or_kind(int flg) { impl_->kind |= flg; }

  int num_values() const { return impl_->num_values; }

  template <class T>
  ArrayRef<T> get_values() const { assert(0); return {}; }

  const SuffixTable& table() const { return impl_->table; }

  // Returns a value convertible to bool that can be used in conditions but not
  // in comparisons and evaluates to "true" if this suffix is not null
  // and "false" otherwise.
  // Example:
  //   if (s) {
  //     // Do something if s is not null.
  //   }
  operator SafeBool() const { return impl_ != 0 ? &SuffixBase::True : 0; }
};

template <> inline
ArrayRef<int> SuffixBase::get_values<int>() const {
  assert(0 == (suf::FLOAT & kind()));
  return { impl_->int_values, (std::size_t)num_values() };
}

template <> inline
ArrayRef<double> SuffixBase::get_values<double>() const {
  assert(0 != (suf::FLOAT & kind()));
  return { impl_->dbl_values, (std::size_t)num_values() };
}


}  // namespace internal

// A suffix.
// Suffixes are data that can be attached to variables, objectives,
// constraints and problems.
class Suffix : private internal::SuffixBase {
 private:
  // SuffixBase is a friend because it needs access to SuffixBase::impl_ via
  // a private base class.
  friend class internal::SuffixBase;
  friend class MutSuffix;

  template <typename Alloc>
  friend class BasicSuffixSet;

  explicit Suffix(const Impl *impl) : SuffixBase(impl) {}

 public:
  Suffix() {}

  template <typename T>
  Suffix(BasicSuffix<T> other) : SuffixBase(other) {}

  using SuffixBase::name;
  using SuffixBase::kind;
  using SuffixBase::num_values;
  using SuffixBase::table;
  using SuffixBase::impl;
  using SuffixBase::operator SafeBool;

  // Iterates over nonzero suffix values and sends them to the visitor.
  template <typename Visitor>
  void VisitValues(Visitor &v) const;
};

class MutSuffix : public Suffix {
 private:
  template <typename Alloc>
  friend class BasicSuffixSet;

  explicit MutSuffix(const Impl *impl) : Suffix(impl) {}

 public:
  MutSuffix() {}

  template <typename T>
  MutSuffix(BasicMutSuffix<T> other) : Suffix(other) {}
};

// "inline" is used here instead of the definition to suppress bogus C4396
// warnings in MSVC.
template <typename SuffixType>
inline SuffixType Cast(Suffix s);

template <typename SuffixType>
inline SuffixType Cast(MutSuffix s);

template <typename T>
class BasicSuffix : private internal::SuffixBase {
 private:
  // SuffixBase is a friend because it needs access to SuffixBase::impl_ via
  // a private base class.
  friend class internal::SuffixBase;
  friend class BasicMutSuffix<T>;

  template <typename Alloc>
  friend class BasicSuffixSet;

  friend BasicSuffix Cast<BasicSuffix>(Suffix s);

  explicit BasicSuffix(const Impl *impl) : SuffixBase(impl) {}
  explicit BasicSuffix(Suffix other) : SuffixBase(other) {}

 public:
  typedef T Type;

  BasicSuffix() {}

  using SuffixBase::name;
  using SuffixBase::kind;
  using SuffixBase::or_kind;
  using SuffixBase::num_values;
  using SuffixBase::operator SafeBool;

  ArrayRef<T> get_values() const {
    return SuffixBase::get_values<T>();
  }

  T value(int index) const {
    MP_ASSERT(index >= 0 && index < impl()->num_values, "index out of bounds");
    T result = T();
    get_value(index, result);
    return result;
  }

  template <typename Visitor>
  void VisitValues(Visitor &v) const {
    for (int i = 0, n = num_values(); i < n; ++i) {
      if (T value = this->value(i))
        v.Visit(i, value);
    }
  }
};

// A mutable suffix.
template <typename T>
class BasicMutSuffix : public BasicSuffix<T> {
 private:
  template <typename Alloc>
  friend class BasicSuffixSet;

  friend BasicMutSuffix Cast<BasicMutSuffix>(MutSuffix s);

  explicit BasicMutSuffix(const typename BasicSuffix<T>::Impl *impl)
    : BasicSuffix<T>(impl) {}
  explicit BasicMutSuffix(MutSuffix other) : BasicSuffix<T>(other) {}

 public:
  BasicMutSuffix() {}

  using BasicSuffix<T>::get_values;

  void set_value(int index, T value) {
    MP_ASSERT(index >= 0 &&
              index < this->impl()->num_values, "index out of bounds");
    BasicSuffix<T>::set_value(index, value);
  }
};

typedef BasicSuffix<int> IntSuffix;
typedef BasicSuffix<double> DoubleSuffix;

typedef BasicMutSuffix<int> MutIntSuffix;
typedef BasicMutSuffix<double> MutDoubleSuffix;

namespace internal {

template <typename T>
struct SuffixInfo {
  enum { KIND = 0 };
};

template <>
struct SuffixInfo<double> {
  enum { KIND = suf::FLOAT };
};

// Returns true if s is of type SuffixType.
template <typename SuffixType>
inline bool Is(Suffix s) {
  return (s.kind() & suf::FLOAT) == SuffixInfo<typename SuffixType::Type>::KIND;
}
}

// Casts a suffix to type SuffixType which must be a valid suffix type.
// Returns a null suffix if s is not convertible to SuffixType.
template <typename SuffixType>
SuffixType Cast(Suffix s) {
  return internal::Is<SuffixType>(s) ? SuffixType(s) : SuffixType();
}
template <typename SuffixType>
SuffixType Cast(MutSuffix s) {
  return internal::Is<SuffixType>(s) ? SuffixType(s) : SuffixType();
}

template <typename Visitor>
inline void Suffix::VisitValues(Visitor &v) const {
  if (IntSuffix int_suffix = Cast<IntSuffix>(*this))
    int_suffix.VisitValues(v);
  else
    Cast<DoubleSuffix>(*this).VisitValues(v);
}

// A set of suffixes.
template <typename Alloc>
class BasicSuffixSet : private Alloc {
 private:
  typedef Suffix::Impl SuffixImpl;

  struct SuffixNameLess {
    bool operator()(const SuffixImpl &lhs, const SuffixImpl &rhs) const {
      std::size_t lhs_size = lhs.name.size(), rhs_size = rhs.name.size();
      if (lhs_size != rhs_size)
        return lhs_size < rhs_size;
      auto cmp = std::char_traits<char>::compare(
            lhs.name.data(), rhs.name.data(), lhs_size);
      return cmp<0;
      /// Don't compare kind_pure() because BasicProb
      /// has own suffix set for each kind
    }
  };

  typedef std::set<SuffixImpl, SuffixNameLess> Set;
  Set set_;

  FMT_DISALLOW_COPY_AND_ASSIGN(BasicSuffixSet);

  SuffixImpl *DoAdd(fmt::StringRef name, int kind, int num_values,
                    const SuffixTable& table);

  template <typename T>
  T *Allocate(std::size_t size) {
    return typename Alloc::template rebind<T>::other(*this).allocate(size);
  }

  template <typename T>
  void Deallocate(T *values, std::size_t size) {
    typename Alloc::template rebind<T>::other(*this).deallocate(values, size);
  }

 public:
  explicit BasicSuffixSet(Alloc alloc = Alloc()) : Alloc(alloc) {}
  ~BasicSuffixSet();

  // Adds a suffix throwing Error if another suffix with the same name is
  // in the set.
  template <typename T>
  BasicMutSuffix<T> Add(fmt::StringRef name, int kind, int num_values,
                        const SuffixTable& table = {}) {
    MP_ASSERT((kind & suf::FLOAT) == 0 ||
              (kind & suf::FLOAT) == internal::SuffixInfo<T>::KIND,
              "invalid suffix kind");
    SuffixImpl *impl = DoAdd(
          name, kind | internal::SuffixInfo<T>::KIND, num_values, table);
    if (num_values != 0) {
      T *values = Allocate<T>(num_values);
      std::fill_n(fmt::internal::make_ptr(values, num_values), num_values, 0);
      impl->values = values;
    }
    return BasicMutSuffix<T>(impl);
  }

  // Finds a suffix with the specified name.
  Suffix Find(fmt::StringRef name) const {
    typename Set::iterator i = set_.find(SuffixImpl(name));
    return Suffix(i != set_.end() ? &*i : 0);
  }
  MutSuffix Find(fmt::StringRef name) {
    typename Set::iterator i = set_.find(SuffixImpl(name));
    return MutSuffix(i != set_.end() ? &*i : 0);
  }

  // Finds a suffix with the specified name and value type.
  template <typename T>
  BasicSuffix<T> Find(fmt::StringRef name) const {
    Suffix s = Find(name);
    return s ? Cast< BasicSuffix<T> >(s) : BasicSuffix<T>();
  }

  // A suffix iterator.
  class iterator : public std::iterator<std::forward_iterator_tag, Suffix> {
   private:
    typename Set::const_iterator it_;

    // A suffix proxy used for implementing operator->.
    class Proxy {
     private:
      Suffix suffix_;

     public:
      explicit Proxy(const SuffixImpl *impl) : suffix_(impl) {}

      const Suffix *operator->() const { return &suffix_; }
    };

   public:
    iterator(typename Set::const_iterator it) : it_(it) {}

    Suffix operator*() const { return Suffix(&*it_); }

    Proxy operator->() const { return Proxy(&*it_); }

    iterator &operator++() {
      ++it_;
      return *this;
    }

    iterator operator++(int ) {
      iterator it(*this);
      ++it_;
      return it;
    }

    bool operator==(iterator other) const { return it_ == other.it_; }
    bool operator!=(iterator other) const { return it_ != other.it_; }
  };

  iterator begin() const { return iterator(set_.begin()); }
  iterator end() const { return iterator(set_.end()); }
};

template <typename Alloc>
BasicSuffixSet<Alloc>::~BasicSuffixSet() {
  // Deallocate names and values.
  for (typename Set::iterator i = set_.begin(), e = set_.end(); i != e; ++i) {
    Deallocate(const_cast<char*>(i->name.data()), i->name.size());
    if ((i->kind & suf::FLOAT) != 0)
      Deallocate(i->dbl_values, i->num_values);
    else
      Deallocate(i->int_values, i->num_values);
  }
}

template <typename Alloc>
typename BasicSuffixSet<Alloc>::SuffixImpl *BasicSuffixSet<Alloc>::DoAdd(
    fmt::StringRef name, int kind, int num_values,
    const SuffixTable& table) {
  std::pair<typename Set::iterator, bool> result =
      set_.insert(Suffix::Impl(name, kind));
  if (!result.second)
    throw Error("duplicate suffix '{}'", name);
  Suffix::Impl *impl = const_cast<SuffixImpl*>(&*result.first);
  // Set name to empty string so that it is not deleted if new throws.
  std::size_t size = name.size();
  impl->name = fmt::StringRef(0, 0);
  char *name_copy = Allocate<char>(size + 1);
  const char *s = name.data();
  std::copy(s, s + size, fmt::internal::make_ptr(name_copy, size));
  name_copy[size] = 0;
  impl->name = name_copy;
  impl->num_values = num_values;
  impl->table = table;
  return impl;
}

typedef BasicSuffixSet< std::allocator<char> > SuffixSet;

class SuffixManager {
 private:
  mp::SuffixSet suffixes_[internal::NUM_SUFFIX_KINDS];

  static void Check(suf::Kind kind) {
    // Assign to an int variable to avoid warning about comparing enums.
    int num_kinds = internal::NUM_SUFFIX_KINDS;
    internal::Unused(kind, num_kinds);
    MP_ASSERT(kind >= 0 && kind < num_kinds, "invalid suffix kind");
  }

 public:
  virtual ~SuffixManager() {}

  typedef MutSuffix Suffix;
  typedef MutIntSuffix IntSuffix;
  typedef mp::SuffixSet SuffixSet;

  // Returns a set of suffixes. TODO hide
  SuffixSet &suffixes(suf::Kind kind) {
    Check(kind);
    return suffixes_[kind];
  }
  const SuffixSet &suffixes(suf::Kind kind) const {
    Check(kind);
    return suffixes_[kind];
  }
};

/// High-level suffix description
template <class T>
class SuffixDef {
  fmt::StringRef name_;
  int kind_;
  SuffixTable tab_;
public:
  SuffixDef(fmt::StringRef nm, int ki, const SuffixTable& st={}) :
    name_(nm), kind_(ki), tab_(st) { }

  using value_type = T;
  fmt::StringRef name() const { return name_; }
  int kind() const { return kind_; }
  const SuffixTable& table() const { return tab_; }
};

}  // namespace mp

#endif  // MP_SUFFIX_H_
