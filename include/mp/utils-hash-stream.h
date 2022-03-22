#ifndef UTILSHASHSTREAM_H
#define UTILSHASHSTREAM_H

#include "mp/utils-hash.h"

/// Wrappers implementing a HashStreamer interface.
/// HashStreamer can
///    1. Compute hash for a given array in 1 go;
///    2. "Accumulate" new values / arrays
///       to an intermediate "hash state".
///
/// Various implementations are possible,
/// one is typedef'd as HashStreamer.

namespace mp {

/// HashStreamerCombine
///
/// Implements the HashStreamer interface using HashCombine()
class HashStreamerCombine {
public:
  /// HashArray, one-off
  template <class T>
  static std::size_t HashArray(std::size_t seed, const T &v) {
    return HashRange(seed, std::begin(v), std::end(v));
  }

  /// HashRange, one-off
  template <class IT>
  static std::size_t HashRange(std::size_t seed, IT beg, const IT end) {
    for ( ; beg!=end; ++beg)
      seed = mp::internal::HashCombine(seed, *beg);
    return seed;
  }

  /// Accumulation (stateful combination)

  /// Constructor
  HashStreamerCombine(std::size_t sd=0) : seed_(sd), seed00_(sd) { }

  /// Add value
  template <class T>
  void Add(const T& v) { seed_ = mp::internal::HashCombine(seed_, v); }

  /// Add array
  template <class T>
  void AddArray(const T& v) { seed_ = HashArray(seed_, v); }

  /// Add range
  template <class IT>
  void AddRange(IT beg, const IT end) { seed_ = HashRange(seed_, beg, end); }

  /// Finalize hash value
  /// @return the value
  std::size_t FinalizeHashValue() { return seed_; }

  /// Reset state
  void Reset() { seed_ = seed00_; }

  /// Finalize & reset
  /// @return hash value
  std::size_t FinalizeAndReset() { auto s = seed_; Reset(); return s; }


private:
  std::size_t seed_, seed00_;
};


/// Selecting one particular implementation as HashStreamer
using HashStreamer = HashStreamerCombine;

} // namespace mp

#endif // UTILSHASHSTREAM_H
