#ifndef CONTEXT_H
#define CONTEXT_H

namespace mp {

class Context {
  enum CtxVal { CTX_NONE, CTX_POS, CTX_NEG, CTX_MIX };
  CtxVal value_ = CTX_NONE;
  Context(CtxVal v) : value_(v) { }
public:
  Context() { }
  bool IsNone() const { return CTX_NONE==value_; }
  bool HasPositive() const { return CTX_POS==value_ || CTX_MIX==value_; }
  bool HasNegative() const { return CTX_NEG==value_ || CTX_MIX==value_; }
  bool IsMixed() const { CTX_MIX==value_; }

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
  /// Add
  Context& Add(Context ct) {
    switch (value_) {
    case CTX_NONE:
      value_ = ct.value_;
      break;
    default:
      throw std::logic_error("Adding context to non-void not implemented");
    }
    return *this;
  }
};

} // namespace mp

#endif // CONTEXT_H
