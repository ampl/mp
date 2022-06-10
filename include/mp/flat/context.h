#ifndef CONTEXT_H
#define CONTEXT_H

namespace mp {

/// Expression context
///
/// <result> <-> <Expression>
/// CTX_POS: expression is implied by the boolean result
/// CTX_NEG: expression's negation is implied by the neg result
/// CTX_MIX: expression is equivalent to the result variable
class Context {
public:
  /// Possible values
  enum CtxVal { CTX_NONE, CTX_POS, CTX_NEG, CTX_MIX };

  /// Construct
  Context(CtxVal v=CTX_NONE) noexcept : value_(v) { }

  /// Is CTX_NONE?
  bool IsNone() const { return CTX_NONE==value_; }

  /// Has CTX_POS?
  bool HasPositive() const { return CTX_POS==value_ || CTX_MIX==value_; }

  /// Has CTX_NEG?
  bool HasNegative() const { return CTX_NEG==value_ || CTX_MIX==value_; }

  /// Has CTX_MIX?
  bool IsMixed() const { return CTX_MIX==value_; }

  /// Get value
  CtxVal GetValue() const { return value_; }

  /// Positivize
  Context operator+() {
    switch (value_) {
    case CTX_NONE:
      return CTX_POS;
    default:
      return value_;
    }
  }

  /// Negate
  Context operator-() {
    switch (value_) {
    case CTX_NONE:
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
