#ifndef SUFFIX_H
#define SUFFIX_H

#include <string>

#include "mp/format.h"
#include "mp/suffix.h"

namespace mp {

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

} // namespace mp

#endif // SUFFIX_H
