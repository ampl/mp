#ifndef CONTEXT_H
#define CONTEXT_H

#include <stdexcept>

namespace mp {

class Context {
public:
  enum CtxVal { CTX_NONE, CTX_POS, CTX_NEG, CTX_MIX };
private:
  CtxVal value_ = CTX_NONE;
public:
  Context() { }
  Context(CtxVal v) : value_(v) { }
  bool IsNone() const { return CTX_NONE==value_; }
  bool HasPositive() const { return CTX_POS==value_ || CTX_MIX==value_; }
  bool HasNegative() const { return CTX_NEG==value_ || CTX_MIX==value_; }
  bool IsMixed() const { return CTX_MIX==value_; }

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
      throw std::logic_error("Adding context to non-empty context not implemented");
    }
    return *this;
  }
};

} // namespace mp

#endif // CONTEXT_H
