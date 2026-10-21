#ifndef VAR_BND_IMPL_H
#define VAR_BND_IMPL_H

/**
 * Comparisons decided by variable-bound constraints
 * (option cvt:pre:vbnd).
 *
 * A pair of variable-bound constraints  x >= L*b,  x <= U*b,
 * with b binary and L > 0, gives  b=0 ==> x <= 0,  b=1 ==> x >= L.
 * This is AMPL's encoding of  x in {0} union interval[L, U]
 * (rows "x@IUlb", "x@IUub" with binary "x@b"), and a common
 * hand-written semicontinuous / on-off formulation.
 * A comparison x (cmp) c with threshold c in the gap of these
 * implications is equivalent to b or to !b, e.g.,
 *   x >= c, 0 < c <= L  <==>  b == 1;
 *   x <= c, 0 <= c < L  <==>  b == 0.
 * Reusing b avoids a new reified comparison with its own binary.
 * Savelsbergh, "Preprocessing and probing techniques for mixed integer
 * programming problems", ORSA J. Computing 6(4), 1994, derives such
 * implications by probing on b; here they are read off the rows.
 */

#include <algorithm>
#include <unordered_map>

namespace mp {

/// Implications of variable-bound constraints on a variable x
/// controlled by a binary b
class VarBinaryBounds {
public:
  /// Literal: (var_ == val_); var_ < 0 if none
  struct Literal {
    int var_ {-1};
    int val_ {0};
  };

  /// Note row  ax*x + ab*b (in) [lo, up]  with binary b.
  /// Records x >= L*b (L > 0) and x <= U*b (U > 0) implications.
  void NoteRow(int x, double ax, int b, double ab, double lo, double up) {
    if (0.0 == ax || 0.0 == ab)
      return;
    const double k = -ab / ax;                  // ax*(x - k*b) (in) [lo, up]
    if (0.0 == lo)                              // ax*(x - k*b) >= 0
      NoteSide(x, b, k, ax > 0.0);
    if (0.0 == up)                              // ax*(x - k*b) <= 0
      NoteSide(x, b, k, ax < 0.0);
  }

  /// The literal deciding  a*x (kind) rhs,  kind -2..2 for <,<=,==,>=,>.
  /// @param x_lb: lower bound of x (for ==)
  Literal Decide(int x, double a, int kind, double rhs, double x_lb) const {
    auto it = map_.find(x);
    if (map_.end() == it || !it->second.IsComplete() || 0.0 == a)
      return {};
    const double c = rhs / a;
    return DecideNormalized(it->second, a > 0.0 ? kind : -kind, c, x_lb);
  }

  /// Any implications recorded?
  bool empty() const { return map_.empty(); }


private:
  /// Implications for one variable
  struct Entry {
    int b_ {-1};              ///< the binary
    double L_ {0.0};          ///< b=1 ==> x >= L_ (best known)
    bool f_vub_ {false};      ///< b=0 ==> x <= 0
    /// Both implications known, L > 0?
    bool IsComplete() const { return b_ >= 0 && f_vub_ && L_ > 0.0; }
  };

  /// Note x >= k*b (if \a f_lower) or x <= k*b
  void NoteSide(int x, int b, double k, bool f_lower) {
    if (k <= 0.0)
      return;
    auto& e = map_[x];
    if (e.b_ >= 0 && e.b_ != b)               // keep the first binary
      return;
    e.b_ = b;
    if (f_lower)
      e.L_ = std::max(e.L_, k);
    else
      e.f_vub_ = true;
  }

  /// Decide x (kind) c for implications \a e
  static Literal DecideNormalized(const Entry& e, int kind, double c,
                                  double x_lb) {
    const double L = e.L_;
    switch (kind) {
    case 1:  return (0.0 < c && c <= L) ? Literal{e.b_, 1} : Literal{};
    case 2:  return (0.0 <= c && c < L) ? Literal{e.b_, 1} : Literal{};
    case -1: return (0.0 <= c && c < L) ? Literal{e.b_, 0} : Literal{};
    case -2: return (0.0 < c && c <= L) ? Literal{e.b_, 0} : Literal{};
    default: return (0.0 == c && x_lb >= 0.0) ? Literal{e.b_, 0} : Literal{};
    }
  }

  std::unordered_map<int, Entry> map_;
};

} // namespace mp

#endif // VAR_BND_IMPL_H
