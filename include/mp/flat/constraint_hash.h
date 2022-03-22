#ifndef CONSTRAINT_HASH_H
#define CONSTRAINT_HASH_H

/// Specialize std::hash<> for some standard constraints

#include <functional>
#include <array>
#include <vector>

#include "mp/utils-hash-stream.h"
#include "mp/flat/constraint_base.h"
#include "mp/flat/expr_affine.h"

namespace std {

/// Specialize std::hash<> for CustomFunctionalConstraint<>
///
/// Assumes that std::hash<> is specialized
/// for Arguments and Parameters
template <class Args, class Params, class NumOrLogic, class Id>
struct hash<
    std::reference_wrapper< const
      mp::CustomFunctionalConstraint<Args, Params, NumOrLogic, Id> > >
{
  size_t operator()( std::reference_wrapper<
        const mp::CustomFunctionalConstraint
          <Args, Params, NumOrLogic, Id> > x) const
  {
    mp::HashStreamer hs;
    hs.Add(std::hash<Args>{}(x.get().GetArguments()));
    hs.Add(std::hash<Params>{}(x.get().GetParameters()));
    return hs.FinalizeHashValue();
  }
};


/// Specialize std::hash<> for std::array<>
///
/// Might assume std::hash<> specialized for elements
template <class Element, std::size_t N>
struct hash< std::array<Element, N> >
{
  size_t operator()(
      const std::array<Element, N>& x) const
  {
    return mp::HashStreamer::HashArray(0, x);
  }
};


/// Specialize std::hash<> for std::vector<>
///
/// Might assume std::hash<> specialized for elements
template <class Element, class Allocator>
struct hash< std::vector<Element, Allocator> >
{
  size_t operator()(
      const std::vector<Element, Allocator>& x) const
  {
    return mp::HashStreamer::HashArray(0, x);
  }
};


/// Specialize std::hash<> for mp::LinTerms
template <>
struct hash< mp::LinTerms >
{
  size_t operator()(
      const mp::LinTerms& x) const
  {
    mp::HashStreamer hs;
    hs.Add(x.vars());
    hs.Add(x.coefs());
    return hs.FinalizeHashValue();
  }
};


/// Specialize std::hash<> for mp::AffExp
template <>
struct hash< mp::AffExp >
{
  size_t operator()(
      const mp::AffExp& x) const
  {
    mp::HashStreamer hs;
    hs.Add(std::hash<mp::LinTerms>{}(x.get_lin_exp()));
    hs.Add(x.constant_term());
    return hs.FinalizeHashValue();
  }
};

}  // namespace std

#endif // CONSTRAINT_HASH_H
