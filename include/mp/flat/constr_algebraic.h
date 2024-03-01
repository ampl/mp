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
  /// For checking solver values,
  /// report violation amount
  /// (negative if holds with slack.)
  /// In logical mode, report 1/0
  /// (what's the violation amount
  ///  for "x>0" when x==0?)
  template <class VarInfo>
  Violation
  ComputeViolation(const VarInfo& x, bool logical=false) const {
    double bd = Body::ComputeValue(x);
    if (!logical) {
      if (RhsOrRange::lb() > bd)
        return {RhsOrRange::lb() - bd, RhsOrRange::lb()};
      if (bd > RhsOrRange::ub())
        return {bd - RhsOrRange::ub(), RhsOrRange::ub()};
      return
      {std::max( // negative. Same for strict cmp?
                 RhsOrRange::lb() - bd, bd - RhsOrRange::ub()),
            0.0};
    }
    return {double(!RhsOrRange::is_valid(bd)), 1.0};
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
  /// kind placeholder
  int kind() const { return -100; }
  /// range lb()
  double lb() const { return lb_; }
  /// range ub()
  double ub() const { return ub_; }
  /// negate
  void negate() { auto tmp=ub_; ub_=-lb_; lb_=-tmp; }
  /// operator==
  bool equals(const AlgConRange& r) const
  { return lb()==r.lb() && ub()==r.ub(); }
  /// validity of the body value
  bool is_valid(double bv) const
  { return bv >= lb() && bv <= ub(); }
private:
  double lb_, ub_;
};


/// Algebraic constraint right-hand side (template parameter).
/// Kind: -2/-1/0/1/2 for < / <= / == / >= / >
template <int kind_>
class AlgConRhs {
  static constexpr const char* kind_str_[] =
  { "LT", "LE", "EQ", "GE", "GT" };
  static constexpr const char* kind_cmp_[] =
  { "<", "<=", "==", ">=", ">" };
public:
  /// name
  static std::string GetTypeName()
  { return std::string("Rhs") + kind_str_[kind_+2]; }
  /// Comparison name as string
  static const char* GetCmpName()
  { return kind_str_[kind_+2]; }
  /// Comparison operator as string
  static const char* GetCmpStr()
  { return kind_cmp_[kind_+2]; }
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
  /// validity of the body value
  bool is_valid(double bv) const {
    switch (kind_) {
    case -2: return bv <  rhs();
    case -1: return bv <= rhs();
    case  0: return bv == rhs();
    case  1: return bv >= rhs();
    case  2: return bv >  rhs();
    default: MP_RAISE("wrong comparison kind");
    }
  }
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


/// Write an algebraic constraint
template <class Body, class RhsOrRange>
inline void WriteJSON(JSONW jw,
                      const AlgebraicConstraint<Body, RhsOrRange>& algc) {
  WriteJSON(jw["body"], algc.GetBody());
  WriteJSON(jw["rhs_or_range"], algc.GetRhsOrRange());
}

/// Write alg con range
inline void WriteJSON(JSONW jw,
                      const AlgConRange& acr) {
  jw << acr.lb() << acr.ub();
}

/// Write alg con rhs
template <int kind>
inline void WriteJSON(JSONW jw,
                      const AlgConRhs<kind>& acrhs) {
  jw << acrhs.GetCmpName() << acrhs.rhs();
}

} // namespace mp

#endif // CONSTRAINTS_ALGEBRAIC_H
