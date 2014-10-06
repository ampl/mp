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

namespace mp {

template <typename SuffixPtr>
class SuffixData;

class Suffix {
 private:
  std::string name_;
  int kind_;
  int *values_;
  std::size_t size_;

  friend class SuffixData<Suffix*>;

  void set_data(int *values, std::size_t size) {
    size_ = size;
    values_ = values;
  }

 public:
  explicit Suffix(const std::string &name, int kind = 0)
    : name_(name), kind_(kind), values_(0), size_(0) {}

  const char *name() const { return name_.c_str(); }

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

class SuffixSet {
 private:
  int kind_;

  struct SuffixLess {
    bool operator()(const Suffix &lhs, const Suffix &rhs) const {
      return std::strcmp(lhs.name(), rhs.name()) < 0;
    }
  };

  typedef std::set<Suffix, SuffixLess> Set;
  Set set_;

 public:
  explicit SuffixSet(int kind = 0) : kind_(kind) {
    assert(0 <= kind && kind <= suf::NUM_KINDS);
  }

  // Finds a suffix with specified name.
  Suffix *Find(const char *name) {
    Set::iterator i = set_.find(Suffix(name));
    return const_cast<Suffix*>(i != set_.end() ? &*i : 0);
  }
  const Suffix *Find(const char *name) const {
    Set::const_iterator i = set_.find(Suffix(name));
    return i != set_.end() ? &*i : 0;
  }

  Suffix &Add(const char *name) {
    return const_cast<Suffix&>(*set_.insert(Suffix(name, kind_)).first);
  }

  typedef Set::const_iterator iterator;

  iterator begin() const { return set_.begin(); }
  iterator end() const { return set_.end(); }
};

class SuffixManager {
 private:
  SuffixSet suffixes_[suf::NUM_KINDS];

 public:
  SuffixManager() {
    for (int kind = 0; kind < suf::NUM_KINDS; ++kind)
      suffixes_[kind] = SuffixSet(kind);
  }

  SuffixSet &get(int kind) {
    assert(kind < suf::NUM_KINDS);
    return suffixes_[kind];
  }
};

// TODO: test suffixes

}  // namespace mp

#endif  // MP_SUFFIX_H_
