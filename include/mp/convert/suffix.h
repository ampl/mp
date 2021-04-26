#ifndef SUFFIX_H
#define SUFFIX_H

#include "mp/format.h"

/// High-level suffix description
template <class T>
class SuffixDef {
  fmt::StringRef name_;
  int kind_;
public:
  SuffixDef(fmt::StringRef nm, int ki) : name_(nm), kind_(ki) { }

  using value_type = T;
  fmt::StringRef name() const { return name_; }
  int kind() const { return kind_; }
};

#endif // SUFFIX_H
