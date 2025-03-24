#ifndef BUCKETACCUM_HPP
#define BUCKETACCUM_HPP

#include <cmath>

#include "mp/flat/bucketaccum.h"

namespace mp {

template <class Body>
void BucketAccum1Type<Body>::Add(Body terms) {
  if (terms.empty())
    return;

  // assert(terms.is_sorted());   // We only add sorted summands
  if (!terms.is_sorted())
    terms.sort_terms();           // now do it if necessary @todo

  auto i_bucket = FindBucket(terms);
  if (HasBucket(i_bucket)) {
    Merge2Bucket(i_bucket, terms);
    if (WrongBucket(i_bucket))
      Add(ExtractBucket(i_bucket));     // re-insert
  } else {
    Assign2Bucket(i_bucket, std::move(terms));
  }
}

template <class Body>
size_t BucketAccum1Type<Body>::FindBucket(
    const Body& terms) const {
  assert(terms.size());
  return size_t(std::floor(
      std::log2((long double)terms.size())));
}

template <class Body>
bool BucketAccum1Type<Body>::HasBucket(
    size_t i_bucket) const
{ return buckets_[i_bucket].size(); }

template <class Body>
bool BucketAccum1Type<Body>::WrongBucket(
    size_t i_bucket) const {
  return HasBucket(i_bucket)
      && i_bucket!=FindBucket(GetBucket(i_bucket));
}

template <class Terms>
Terms MergeSorted(const Terms& t1, const Terms& t2);

template <class Body>
void BucketAccum1Type<Body>::Merge2Bucket(
    size_t i_bucket, const Body& terms) {
  buckets_[i_bucket]
      = MergeSorted(GetBucket(i_bucket), terms);
}

template <class Body>
void BucketAccum1Type<Body>::Assign2Bucket(
    size_t i_bucket, Body terms) {
  buckets_[i_bucket] = std::move(terms);
}

template <class Body>
Body BucketAccum1Type<Body>::ExtractBucket(size_t i) {
  Body result = std::move(buckets_.at(i));
  buckets_[i].clear();            // std::move does not
  return result;
}

template <class Body>
Body BucketAccum1Type<Body>::ExtractSum() {
  Body result;

  size_t i_bucket = 0;
  for ( ; i_bucket<buckets_.size(); ++i_bucket)
    if (HasBucket(i_bucket)) {
      result = ExtractBucket(i_bucket);
      break;
    }

  while ( ++i_bucket<buckets_.size())
    if (HasBucket(i_bucket)) {
      result
          = result.size()
          ? MergeSorted(ExtractBucket(i_bucket), result)
          : ExtractBucket(i_bucket);
    }

  return result;
}

template <class EExprLike>
void BucketAccumulator<EExprLike>::Add(EExprLike ee) {
  ct_ += ee.constant_term();
  ba_l_.Add(std::move(ee.GetBody().GetLinTerms()));
  ba_q_.Add(std::move(ee.GetBody().GetQPTerms()));
}


template <class Body>
BucketAccum1Type<Body>::BucketAccum1Type(size_t ) { }

template <class Body>
BucketAccum1Type<Body>::~BucketAccum1Type() {
#ifndef NDEBUG
  for (const auto& b: buckets_)
    assert(b.empty());
#endif
}

}  // namespace mp

#endif // BUCKETACCUM_HPP
