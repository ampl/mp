#ifndef CONTEXT_H
#define CONTEXT_H

#include "mp/error.h"

namespace mp {

/// Expression context.
///
/// Considering the relation (result) <-> (Expression):
/// - CTX_ROOT: expression is true (even if no result variable)
/// - CTX_POS: expression is implied by the boolean result
/// - CTX_NEG: expression's negation is implied by the neg result
/// - CTX_MIX: expression is equivalent to the result variable
class Context {
public:
  /// Possible values
  enum CtxVal { CTX_NONE, CTX_ROOT, CTX_POS, CTX_NEG, CTX_MIX };

  /// Construct
  Context(CtxVal v=CTX_NONE) noexcept : value_(v) { }

  /// Is CTX_NONE?
  bool IsNone() const { return CTX_NONE==value_; }

  /// Is a subset of or equal to \a other?
  bool IsSubsetOf(Context other) const {
    return
        IsProperSubsetOf(other) || *this==other;
  }

  /// Is a proper subset of \a other?
  bool IsProperSubsetOf(Context other) const {
    return
        (IsNone() && !other.IsNone())
        || (!IsMixed() && other.IsMixed());
  }

  /// Has CTX_ROOT?
  bool HasRoot() const { return CTX_ROOT==value_; }

  /// Has CTX_POS?
  bool HasPositive() const
  { return CTX_ROOT==value_ || CTX_POS==value_ || CTX_MIX==value_; }

  /// Has CTX_NEG?
  bool HasNegative() const { return CTX_NEG==value_ || CTX_MIX==value_; }

  /// Is CTX_ROOT?
  bool IsRoot() const { return CTX_ROOT==value_; }

  /// Is CTX_POS?
  bool IsPositive() const { return CTX_POS==value_; }

  /// Is CTX_NEG?
  bool IsNegative() const { return CTX_NEG==value_; }

  /// Has CTX_MIX?
  bool IsMixed() const { return CTX_MIX==value_; }

  /// Get value
  CtxVal GetValue() const { return value_; }

  /// Equal
  bool operator==(Context ctx) const { return value_==ctx.value_; }

  /// Positivize
  Context operator+() {
    switch (value_) {
    case CTX_NONE:
      return CTX_POS;
    case CTX_ROOT:
      return CTX_POS;
    default:
      return value_;
    }
  }

  /// Negate
  Context operator-() {
    switch (value_) {
    case CTX_NONE:
    case CTX_ROOT:
    case CTX_POS:
      return CTX_NEG;
    case CTX_NEG:
      return CTX_POS;
    default:
      return value_;
    }
  }

  /// Add (merge) new context
  Context& Add(Context ct) {
    switch (value_) {          // Consider already stored value
    case CTX_NONE:
      value_ = ct.value_;
      break;
    case CTX_ROOT:
      MP_ASSERT(0, "Should not extend root context");
      [[fallthrough]];           // in release
    case CTX_POS:
      if (ct.HasNegative())
        value_ = CTX_MIX;
      break;
    case CTX_NEG:
      if (ct.HasPositive())
        value_ = CTX_MIX;
      break;
    default:
      break;
    }
    return *this;
  }


private:
  CtxVal value_=CTX_NONE;
};

} // namespace mp

#endif // CONTEXT_H
