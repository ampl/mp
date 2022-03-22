#ifndef UTILSHASH_H
#define UTILSHASH_H

#include <functional>


namespace mp {
namespace internal {
/// HashCombine template, 32-bit
template <class T>
inline std::size_t HashCombine(std::size_t seed, const T &v) {
  return seed ^ (std::hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2));
}
}
}


#endif // UTILSHASH_H
