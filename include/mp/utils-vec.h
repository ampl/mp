#ifndef UTILSVEC_H
#define UTILSVEC_H

#include <cmath>

// #include "../../thirdparty/martinus/svector/svector.h"
#include "../../thirdparty/gharveymn/small_vector/small_vector.hpp"


namespace mp {

/// Typedef small vector
// template <class T, size_t N>
// using SmallVec = ankerl::svector<T, N>;

/// Typedef small vector
template <class T, size_t N>
using SmallVec = gch::small_vector<T, N>;

/// Typedef small vector, default size of 64 bytes
template <class T>
using SmallVecDefSz = gch::small_vector<T>;

/// Grow vector capacity by a factor if needed;
/// @note better preallocate, or call in the reverse order of indexes.
template <class Vec>
void ResizePlus(Vec& vec, typename Vec::size_type i) {
  if (vec.size()<=i) {
    if (vec.capacity()<=i)
      vec.reserve(
          typename Vec::size_type(std::ceil((i+1)*1.3)));
    vec.resize(i+1);
  }
}

/// Grow vector capacity by a factor if needed;
/// set vec[i] = v.
/// @note better preallocate, or call in the reverse order of indexes.
template <class Vec>
void AutoExpand(
    Vec& vec, typename Vec::size_type i, typename Vec::value_type v) {
  ResizePlus(vec, i);
  vec[i] = std::move(v);
}

}  // namespace mp


#endif // UTILSVEC_H
