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

#include "mp/common.h"
#include "mp/format.h"

namespace mp {

template <typename Alloc>
class BasicProblem;

// A suffix.
// Suffixes are arbitrary metadata that can be attached to variables,
// objectives, constraints and problems.
class Suffix {
 private:
  std::string name_;
  int kind_;
  int *values_;
  int size_;

  template <typename Alloc>
  friend class BasicProblem;

  void InitValues(int size) {
    assert(!values_);
    values_ = new int[size];
    size_ = size;
  }

 public:
  Suffix(fmt::StringRef name, int kind)
    : name_(name.c_str(), name.size()), kind_(kind), values_(0), size_(0) {}
  ~Suffix() {
    delete [] values_;
  }

  // Returns the suffix name.
  const char *name() const { return name_.c_str(); }

  // Returns the suffix kind.
  int kind() const { return kind_; }

  int value(int index) const {
    assert(index < size_);
    return values_[index];
  }

  void set_value(int index, int value) {
    assert(index < size_);
    values_[index] = value;
  }

  void set_value(int index, double value) {
    assert(index < size_);
    // TODO: set double value
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

// A set of suffixes.
class SuffixSet {
 private:
  struct SuffixNameLess {
    bool operator()(const Suffix &lhs, const Suffix &rhs) const {
      return std::strcmp(lhs.name(), rhs.name()) < 0;
    }
  };

  // Suffixes are stored in a set which is not efficient, but it doesn't
  // matter because the number of suffixes is normally small.
  typedef std::set<Suffix, SuffixNameLess> Set;
  Set set_;

  template <typename Alloc>
  friend class BasicProblem;

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

class SuffixManager {
 private:
  mp::SuffixSet suffixes_[suf::NUM_KINDS];

 public:
  virtual ~SuffixManager() {}

  typedef mp::Suffix *SuffixPtr;
  typedef mp::SuffixSet SuffixSet;

  // Returns a set of suffixes.
  SuffixSet &suffixes(int kind) {
    assert(kind < suf::NUM_KINDS);
    return suffixes_[kind];
  }
};

}  // namespace mp

#endif  // MP_SUFFIX_H_
