#ifndef CONSTRAINT_HASH_H
#define CONSTRAINT_HASH_H

/// Specialize std::hash<> for some standard constraints

#include <functional>
#include <array>
#include <vector>

#include "mp/utils-hash-stream.h"
#include "mp/flat/constr_std.h"
#include "mp/flat/expr_affine.h"

namespace std {

/// Specialize std::hash<> for various expressions and constraints.
/// Remember we also need operator== for [std::reference_wrapper<> of]
/// the compared types. Operator=='s are class memebrs, or
/// defined in ``constraint_keeper.h``.

/// Partially specialize std::hash<> for CustomFunctionalConstraint<>
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


/// Partially specialize std::hash<> for ConditionalConstraint<>
///
/// Assumes that std::hash<> is specialized
/// for \a Con
template <class Con>
struct hash<
    std::reference_wrapper< const
      mp::ConditionalConstraint<Con> > >
{
  size_t operator()( std::reference_wrapper<
        const mp::ConditionalConstraint<Con> > x) const
  {
    return std::hash<Con>{}(x.get().GetConstraint());
  }
};


/// ... which happens here:
/// Partially specialize std::hash<> for AlgebraicConstraint<>
///
/// Assumes that std::hash<> is specialized
/// for \a Body and \a RangeOrRhs
template <class Body, class RangeOrRhs>
struct hash<
      mp::AlgebraicConstraint<Body, RangeOrRhs> >
{
  size_t operator()(
        const mp::AlgebraicConstraint<Body, RangeOrRhs>& x) const
  {
    mp::HashStreamer hs;
    hs.Add(std::hash<Body>{}(x.GetBody()));
    hs.Add(std::hash<RangeOrRhs>{}(x.GetRhsOrRange()));
    return hs.FinalizeHashValue();
  }
};


/// Partially specialize std::hash<> for mp::AlgConRhs<>
template <int kind>
struct hash< mp::AlgConRhs<kind> >
{
  size_t operator()(
      mp::AlgConRhs<kind> x) const
  {
    return std::hash<double>{}(x.rhs());
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


/// Partially specialize std::hash<> for mp::AlgebraicExpr
template <class Body>
struct hash< mp::AlgebraicExpression<Body> >
{
  size_t operator()(
      const mp::AlgebraicExpression<Body>& x) const
  {
    mp::HashStreamer hs;
    hs.Add(std::hash<Body>{}(x.GetBody()));
    hs.Add(x.constant_term());
    return hs.FinalizeHashValue();
  }
};


/// Specialize std::hash<> for mp::QuadTerms
template <>
struct hash< mp::QuadTerms >
{
  size_t operator()(
      const mp::QuadTerms& qt) const
  {
    mp::HashStreamer hs;
    hs.Add(qt.vars1());
    hs.Add(qt.vars2());
    hs.Add(qt.coefs());
    return hs.FinalizeHashValue();
  }
};


/// Specialize std::hash<> for mp::QuadAndLinTerms
template <>
struct hash< mp::QuadAndLinTerms >
{
  size_t operator()(
      const mp::QuadAndLinTerms& qe) const
  {
    mp::HashStreamer hs;
    hs.Add(std::hash<mp::LinTerms>{}(qe.GetLinTerms()));
    hs.Add(std::hash<mp::QuadTerms>{}(qe.GetQPTerms()));
    return hs.FinalizeHashValue();
  }
};


/// Specialize std::hash<> for mp::LinearFunctionalCon
template <>
struct hash<
    std::reference_wrapper< const mp::LinearFunctionalConstraint > >
{
  size_t operator()(
      std::reference_wrapper<
        const mp::LinearFunctionalConstraint> lfc) const
  {
    return std::hash<mp::AffineExpr>{}(lfc.get().GetAffineExpr());
  }
};


/// Specialize std::hash<> for mp::QuadraticFunctionalCon
template <>
struct hash<
    std::reference_wrapper< const mp::QuadraticFunctionalConstraint > >
{
  size_t operator()(
      std::reference_wrapper<
        const mp::QuadraticFunctionalConstraint > qfc) const
  {
    return std::hash<mp::QuadraticExpr>{}(qfc.get().GetQuadExpr());
  }
};

}  // namespace std

#endif // CONSTRAINT_HASH_H
