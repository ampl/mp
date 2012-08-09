// Function adapters, binders, numerical differentiator and other
// function-related stuff.

#include "tests/functional.h"

std::ostream &fun::operator<<(std::ostream &os, const Tuple &t) {
  os << "(";
  if (unsigned size = t.size()) {
    os << t[0];
    for (size_t i = 1; i < size; ++i)
      os << ", " << t[i];
  }
  os << ")";
  return os;
}
