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

#include <cstddef>  // for std::size_t
#include <cstring>
#include <set>
#include <string>
#include <vector>

#include "mp/format.h"

namespace mp {

template <typename SuffixPtr>
class SuffixData;

// A suffix.
// Suffixes are arbitrary metadata that can be attached to variables,
// objectives, constraints and problems.
class Suffix {
 private:
  std::string name_;
  int kind_;
  int *values_;
  std::size_t size_;

  void set_data(int *values, std::size_t size) {
    size_ = size;
    values_ = values;
  }

  friend class SuffixData<Suffix*>;

 public:
  Suffix(fmt::StringRef name, int kind)
    : name_(name.c_str(), name.size()), kind_(kind), values_(0), size_(0) {}

  // Returns the suffix name.
  const char *name() const { return name_.c_str(); }

  // Returns the suffix kind.
  int kind() const { return kind_; }

  int value(std::size_t index) const {
    assert(index <= size_);
    return values_[index];
  }

  // Iterates over nonzero suffix values and sends them to the visitor.
  template <typename Visitor>
  void VisitValues(Visitor &visitor) const {
    for (std::size_t i = 0; i < size_; ++i) {
      int value = values_[i];
      if (value != 0)
        visitor.Visit(i, value);
    }
  }
};

// Suffix data.
// For compatibility with ASL, suffix data is stored separately.
// This class is used to control the lifetime of the data and
// automatically detach it from the suffix once it is destroyed.
template <typename SuffixPtr>
class SuffixData {
 private:
  SuffixPtr suffix_;
  std::vector<int> values_;

  FMT_DISALLOW_COPY_AND_ASSIGN(SuffixData);

 public:
  SuffixData() : suffix_() {}
  ~SuffixData() {
    if (suffix_)
      suffix_->set_data(0, 0);
  }

  // Attaches data to a suffix.
  void Attach(SuffixPtr suffix, std::size_t num_values) {
    suffix_ = suffix;
    values_.resize(num_values);
    suffix->set_data(&values_[0], values_.size());
  }

  int value(std::size_t index) const {
    assert(index <= values_.size());
    return values_[index];
  }

  void set_value(std::size_t index, int value) {
    assert(index <= values_.size());
    values_[index] = value;
  }
};

// A set of suffixes.
class SuffixSet {
 private:
  struct SuffixNameLess {
    bool operator()(const Suffix &lhs, const Suffix &rhs) const {
      return std::strcmp(lhs.name(), rhs.name()) < 0;
    }
  };

  typedef std::set<Suffix, SuffixNameLess> Set;
  Set set_;

  friend class Problem;

 public:
  // Finds a suffix with specified name.
  Suffix *Find(fmt::StringRef name) {
    Set::iterator i = set_.find(Suffix(name, 0));
    return const_cast<Suffix*>(i != set_.end() ? &*i : 0);
  }
  const Suffix *Find(fmt::StringRef name) const {
    Set::const_iterator i = set_.find(Suffix(name, 0));
    return i != set_.end() ? &*i : 0;
  }

  typedef Set::const_iterator iterator;

  iterator begin() const { return set_.begin(); }
  iterator end() const { return set_.end(); }
};
}  // namespace mp

#endif  // MP_SUFFIX_H_
