#ifndef CONSTRAINTS_ALGEBRAIC_H
#define CONSTRAINTS_ALGEBRAIC_H

/**
  * Static algebraic constraints
  */

#include <string>
#include <cmath>

#include "mp/flat/constr_base.h"
#include "mp/flat/expr_quadratic.h"


namespace mp {

////////////////////////////////////////////////////////////////////////
/// Generic algebraic constraint
/// @param Body: linear or linear and higher-order terms
/// @param RhsOrRange: rhs or range
template <class Body, class RhsOrRange>
class AlgebraicConstraint :
    public BasicConstraint, public Body, public RhsOrRange {
public:
  /// Constraint type name
  static const std::string& GetTypeName() {
    static std::string name {
      std::string("AlgebraicConstraint< ") +
      Body::GetTypeName() + ", " +
      RhsOrRange::GetTypeName() + " >" };
    return name;
  }

  /// Is logical?
  static bool IsLogical() { return false; }

  /// BodyType
  using BodyType = Body;

  /// RhsOrRangeType
  using RhsOrRangeType = RhsOrRange;

  /// Constructor.
  /// By default (\a fSort = true), it sorts terms.
  /// Pass \a fSort = false to skip if you populate the terms
  /// but do sorting later.
  /// @param le: linear / linear + higher-order terms
  /// @param rr: rhs or range
  AlgebraicConstraint(Body le, RhsOrRange rr, bool fSort=true)
    : Body(std::move(le)), RhsOrRange(std::move(rr))
  { if (fSort) sort_terms(); }

  /// Body: linear or linear + higher-order terms, const
  const Body& GetBody() const { return (const Body&)(*this); }

  /// Body: linear or linear + higher-order terms
  Body& GetBody() { return (Body&)(*this); }

  /// Range or RHS. Used for hash<>
  RhsOrRange GetRhsOrRange() const { return (RhsOrRange)(*this); }

  /// Synonym, For PropagateResult()
  const Body& GetArguments() const { return GetBody(); }

  /// If no variable terms in the body
  bool empty() const { return Body::empty(); }

  /// Compute lower slack
  double ComputeLowerSlack(ArrayRef<double> x) const {
    return Body::ComputeValue(x) - RhsOrRange::lb();
  }

  /// Compute violation.
  /// Negative if holds with slack.
  template <class VarInfo>
  double
  ComputeViolation(const VarInfo& x) const {
    auto bd = Body::ComputeValue(x);
    return std::max(      // same for strict cmp?
          RhsOrRange::lb() - bd, bd - RhsOrRange::ub());
  }

  /// Sorting and merging terms, some solvers require
  void sort_terms() { Body::sort_terms(); }

  /// Is Normalized?
  bool is_normalized() {
    sort_terms();
    return GetBody().is_normalized();
  }

  /// Negate all terms
  void negate() { Body::negate(); RhsOrRange::negate(); }

  /// Testing API
  bool operator==(const AlgebraicConstraint& lc) const {
    return Body::equals(lc) && RhsOrRange::equals(lc);
  }
};


/// Algebraic constraint range (template parameter)
class AlgConRange {
public:
  /// Class name
  static constexpr const char* GetTypeName() { return "Range"; }
  /// Constructor
  AlgConRange(double l, double u) : lb_(l), ub_(u) { }
  /// range lb()
  double lb() const { return lb_; }
  /// range ub()
  double ub() const { return ub_; }
  /// negate
  void negate() { auto tmp=ub_; ub_=-lb_; lb_=-tmp; }
  /// operator==
  bool equals(const AlgConRange& r) const
  { return lb()==r.lb() && ub()==r.ub(); }
private:
  double lb_, ub_;
};


/// Algebraic constraint right-hand side (template parameter).
/// Kind: -2/-1/0/1/2 for < / <= / == / >= / >
template <int kind_>
class AlgConRhs {
  static constexpr const char* kind_str_[] =
  { "LT", "LE", "EQ", "GE", "GT" };
public:
  /// name
  static std::string GetTypeName()
  { return std::string("Rhs") + kind_str_[kind_+2]; }
  /// Constructor
  AlgConRhs(double r) : rhs_(r) { }
  /// Kind
  static constexpr int kind() { return kind_; }
  /// rhs()
  double rhs() const { return rhs_; }
  /// lb(): this is a specialization of the range constraint
  double lb() const {
    return kind_<0 ? -INFINITY : rhs();
  }
  /// ub(): this is a specialization of the range constraint
  double ub() const {
    return kind_>0 ? INFINITY : rhs();
  }
  /// Set rhs
  void set_rhs(double v) { rhs_ = v; }
  /// Add to RHS
  void add_to_rhs(double v) { rhs_ += v; }
  /// negate
  void negate() { rhs_ = -rhs_; }
  /// operator==
  bool equals(const AlgConRhs& r) const
  { return rhs()==r.rhs(); }
private:
  double rhs_;
};


////////////////////////////////////////////////////////////////////////
/// Linear range constraint
using LinConRange = AlgebraicConstraint<LinTerms, AlgConRange>;

/// Convenience typedef
template <int sens>
using LinConRhs = AlgebraicConstraint< LinTerms, AlgConRhs<sens> >;
/// Linear constraint c'x <  d
using LinConLT = LinConRhs<-2>;
/// Linear constraint c'x <= d
using LinConLE = LinConRhs<-1>;
/// Linear constraint c'x == d
using LinConEQ = LinConRhs< 0>;
/// Linear constraint c'x >= d
using LinConGE = LinConRhs< 1>;
/// Linear constraint c'x >  d
using LinConGT = LinConRhs< 2>;


////////////////////////////////////////////////////////////////////////
/// Quadratic range constraint
using QuadConRange =
    AlgebraicConstraint<QuadAndLinTerms, AlgConRange>;

/// Convenience typedef
template <int sens>
using QuadConRhs =
    AlgebraicConstraint< QuadAndLinTerms, AlgConRhs<sens> >;
/// Quadratic constraint c'x+x'Qx <  d
using QuadConLT = QuadConRhs<-2>;
/// Quadratic constraint c'x+x'Qx <= d
using QuadConLE = QuadConRhs<-1>;
/// Quadratic constraint c'x+x`Qx == d
using QuadConEQ = QuadConRhs< 0>;
/// Quadratic constraint c'x+x`Qx >= d
using QuadConGE = QuadConRhs< 1>;
/// Quadratic constraint c'x+x`Qx >  d
using QuadConGT = QuadConRhs< 2>;

} // namespace mp

#endif // CONSTRAINTS_ALGEBRAIC_H
