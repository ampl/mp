#ifndef BUCKETACCUM_H
#define BUCKETACCUM_H

#include <array>
#include <climits>
#include <type_traits>
#include <utility>

namespace mp {

/// Accumulate single body type
/// (linear/quadratic).
///
/// It performs a MergeSort-like accumulation.
template <class Body>
class BucketAccum1Type {
public:
  /// Construct
  BucketAccum1Type(size_t n_terms = 2);
  /// Destruct
  ~BucketAccum1Type();
  /// Insert a Body
  void Add(Body b);
  /// Extract the sum
  /// @note should be done exactly once
  Body ExtractSum();
protected:
  size_t FindBucket(const Body& ) const;
  bool HasBucket(size_t i) const;
  /// If the contents are non-empty
  /// and in the wrong place
  bool WrongBucket(size_t i) const;
  void Merge2Bucket(size_t i, const Body& terms);
  void Assign2Bucket(size_t i, Body terms);

  const Body& GetBucket(size_t i) const { return buckets_[i]; }
  /// Move out (and delete the contents)
  Body ExtractBucket(size_t i);
private:
  static constexpr double bucket_factor_ {2.0};  // ignored
  static constexpr size_t sz_max_ = size_t(-1);  // SIZE_MAX
  static constexpr size_t n_buckets_max_
      = sizeof(size_t) * CHAR_BIT;

  std::array<Body, n_buckets_max_> buckets_;
};


/// Accumulate EExpr's,
/// or anything with its interface
template <class EExprLike>
class BucketAccumulator {
public:
  /// Alg expr body type
  using AlgExprBodyType = std::decay_t<
      decltype(std::declval<EExprLike>().GetBody()) >;
  /// Linear expr type
  using LinTermsType = std::decay_t<
      decltype(std::declval<AlgExprBodyType>().GetLinTerms()) >;
  /// Quad expr type
  using QuadTermsType = std::decay_t<
      decltype(std::declval<AlgExprBodyType>().GetQPTerms()) >;

  /// Construct
  BucketAccumulator(size_t n_terms = 2)
      : ba_l_(n_terms), ba_q_(n_terms) { }
  /// Insert EExpr
  void Add(EExprLike ee);
  /// Extract the sum
  /// @note should happen exactly 1x.
  EExprLike ExtractSum() {
    return EExprLike
        {ba_l_.ExtractSum(), ba_q_.ExtractSum(), ct_};
  }
private:
  BucketAccum1Type<LinTermsType> ba_l_;
  BucketAccum1Type<QuadTermsType> ba_q_;
  double ct_ {0.0};
};

}  // namespace mp

#endif // BUCKETACCUM_H
