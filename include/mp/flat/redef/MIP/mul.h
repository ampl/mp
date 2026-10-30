#ifndef MUL_H
#define MUL_H

#include <cassert>
#include <cmath>
#include <utility>

#include "mp/error.h"
#include "mp/flat/redef/redef_base.h"
#include "mp/flat/constr_std.h"

namespace mp {

/// r = (functional constraint \a con), with mixed context:
/// converters creating expressions must give them a context
/// (asserted when explicified for expression-based solvers)
template <class MC, class FuncCon>
int AssignResultVarMixedCtx(MC& mc, FuncCon&& con) {
  const int r = mc.AssignResultVar2Args(std::forward<FuncCon>(con));
  mc.PropagateResultOfInitExpr(r, Context::CTX_MIX);
  return r;
}

/// Number of binary digits to represent 0..range
inline int NumBinaryDigits(double range) {
  int k=0;
  for ( ; range >= 1.0; range = std::floor(range / 2.0))
    ++k;
  return k;
}

/// The factor of x*y to expand in binary digits for an exact MIP
/// reformulation: an integer one with finite bounds and the smaller
/// range, needing at most cvt:mip:intprodbits binaries,
/// while the other factor has finite bounds (for the big-Ms).
/// @return -1 if none.
template <class MC>
int FactorToExpand(MC& mc, int x, int y) {
  int result = -1;
  double best = INFINITY;
  for (int v: {x, y}) {
    const int w = (v==x) ? y : x;
    if (!mc.is_var_integer(v)
        || mc.lb(w) <= mc.PracticallyMinusInf()
        || mc.ub(w) >= mc.PracticallyInf()
        || mc.lb(v) <= mc.PracticallyMinusInf()
        || mc.ub(v) >= mc.PracticallyInf())
      continue;
    const double range = mc.ub(v) - mc.lb(v);
    if (range < best && NumBinaryDigits(range) <= mc.IntProdBits()) {
      best = range;
      result = v;
    }
  }
  return result;
}

/// c*x*y, y integer with finite domain [ly, uy]: binary expansion
/// y = ly + sum_j 2^j b_j, so that
/// x*y = ly*x + sum_j 2^j (b_j * x), each b_j*x an IfThen (exact).
/// See Gupte, Ahmed, Cheon, Dey (2013), SIAM J. Optim. 23(2).
template <class MC>
LinTerms ProductByBinaryExpansion(MC& mc, double c, int x, int y) {
  const double ly = mc.lb(y), range = mc.ub(y) - ly;
  assert(mc.is_var_integer(y));                 // see FactorToExpand()
  assert(ly > mc.PracticallyMinusInf() && mc.ub(y) < mc.PracticallyInf());
  assert(mc.lb(x) > mc.PracticallyMinusInf()    // big-Ms of the IfThens
         && mc.ub(x) < mc.PracticallyInf());
  assert(range >= 0.0 && NumBinaryDigits(range) <= 62);
  LinTerms link { {1.0}, {y} };                 // y - sum 2^j b_j = ly
  LinTerms result;
  if (ly)
    result.add_term(c*ly, x);
  const int zero = int( mc.MakeFixedVar(0.0) );
  for (double w=1.0; w <= range; w *= 2.0) {
    const int b = int( mc.AddVar(0.0, 1.0, var::INTEGER) );
    link.add_term(-w, b);
    const int z = AssignResultVarMixedCtx(mc,  // z = b ? x : 0
          IfThenConstraint{ {b, x, zero} });
    result.add_term(c*w, z);
  }
  assert((int)link.size() == 1 + NumBinaryDigits(range)); // digits cover it
  mc.AddConstraint( LinConEQ( std::move(link), {ly} ) );
  return result;
}

/// c*x*y -> 0.5*c((x+y)^2-x^2-y^2):
/// the squares are piecewise-linearly approximated by MIP solvers
template <class MC>
LinTerms ProductBySquares(MC& mc, double c, int x, int y) {
  LinTerms lt;
  auto x_plus_y = mc.AssignResultVar2Args(
        LinearFunctionalConstraint{{  // = x+y
                                      { {1.0, 1.0}, {x, y} }, 0.0 }} );
  for (auto v_coef: { std::make_pair(x_plus_y, 0.5*c),
                      std::make_pair(x, -0.5*c),
                      std::make_pair(y, -0.5*c) }) {
    auto v_pow2 = mc.AssignResultVar2Args(
          PowConstExpConstraint{  // no context as of now
                          VarArray1{v_coef.first}, {2.0} } );
    lt.add_term(v_coef.second, v_pow2);
  }
  return lt;
}


/// Convert all quadratic terms in the body of a QC.
/// This is the most universal way,
/// in case some converter adds QC.
/// An alternative would be, when flattening, to replace
/// every multiplication by
/// a special new constraint, and convert only these.
template <class ModelConverter, class QuadCon, class LinCon>
class QCConverter_MIP :
    public BasicItemConverter<ModelConverter> {
public:
  /// Base class
  using Base = BasicItemConverter<ModelConverter>;
  /// Constructor
  QCConverter_MIP(ModelConverter& mc) : Base(mc) { }
  /// Converted item type
  using ItemType = QuadCon;

  /// Check whether the constraint
  /// needs to be converted despite being accepted by ModelAPI.
  bool IfNeedsConversion(const ItemType& , int ) {
		return !GetMC().IfPassQuadCon();
  }

  /// Keep for Convert2NL()? Assume QuadExpr accepted.
  /// Otherwise, not accepted QCs are linearized,
  /// also when the user lowered acc:quad...
  bool IfKeepUnconverted(const ItemType& , int ) {
    return GetMC().IfPassQuadConAsNL();
  }

  /// Conversion
  Context Convert(const ItemType& qc, int ) {
		LinearizeQPTerms(qc);
    return Context::CTX_MIX;
  }


protected:
  using Base::GetMC;

  void LinearizeQPTerms(const ItemType& qc) {
    const auto& body = qc.GetBody();
    // Copy lin terms
    auto lin_terms = body.GetLinTerms();
    // Convert quadratic terms
    const auto& qp_terms = body.GetQPTerms();
    for (int i=0; i<(int)qp_terms.size(); ++i) {
      auto c = qp_terms.coef(i);
      auto x = qp_terms.var1(i);
      auto y = qp_terms.var2(i);
      lin_terms.add( LinearizeQPTerm(c, x, y) );
    }
    // Sort linear body.
    lin_terms.sort_terms();
    // Not _AS_ROOT: see #248
    GetMC().AddConstraint( LinCon{
                             lin_terms, qc.GetRhsOrRange() } );
  }

  LinTerms LinearizeQPTerm(double c, int x, int y) {
    if (GetMC().is_binary_var(x) ||
        GetMC().is_binary_var(y))
      return LinearizeProductWithBinaryVar(c, x, y);
    if (GetMC().UserAcceptsAndRecommends((MulConstraint*)nullptr)
        || FactorToExpand(GetMC(), x, y) >= 0)
      return ProductAsMul(c, x, y);              // see MulConverter_MIP
    return ProductBySquares(GetMC(), c, x, y);
  }

  /// c*x*y -> c*r, r = x*y, for solvers taking products natively
  /// (acc:mul), or to be expanded exactly by MulConverter_MIP.
  /// Equal products share r.
  LinTerms ProductAsMul(double c, int x, int y) {
    if (y < x)
      std::swap(x, y);        // so that x*y and y*x are the same
    auto r = AssignResultVarMixedCtx(GetMC(),
          MulConstraint{ {x, y} } );
    return { {c}, {r} };
  }

  LinTerms LinearizeProductWithBinaryVar(double c, int x, int y) {
    LinTerms lt;
    MP_ASSERT_ALWAYS(GetMC().is_binary_var(x) ||
                     GetMC().is_binary_var(y),
                     "Can only convert product with a binary variable");
    bool is_x_bin = GetMC().is_binary_var(x);
    auto i_bin = is_x_bin ? x : y;       // index of the (chosen) binary var
    auto i_other = is_x_bin ? y : x;                  // the other var index
    auto x_ifthen = GetMC().AssignResultVar2Args(
          IfThenConstraint {{      // no context as of now
                                   i_bin,
                                   i_other,
                                   int( GetMC().MakeFixedVar(0.0) )
                            }} );
    lt.add_term(c, x_ifthen);
    return lt;
  }
};


/// Typedef MulCvt<QuadConLE>
template <class MC>
using MulCvtLE_MIP = QCConverter_MIP<MC, QuadConLE, LinConLE>;

/// Typedef MulCvt<QuadConEQ>
template <class MC>
using MulCvtEQ_MIP = QCConverter_MIP<MC, QuadConEQ, LinConEQ>;

/// Typedef MulCvt<QuadConGE>
template <class MC>
using MulCvtGE_MIP = QCConverter_MIP<MC, QuadConGE, LinConGE>;


/// Converts r = x*y (MulConstraint) for MIP:
/// exactly by binary expansion of an integer factor
/// (cvt:mip:intprodbits), otherwise approximately via squares.
template <class ModelConverter>
class MulConverter_MIP :
    public BasicFuncConstrCvt<
      MulConverter_MIP<ModelConverter>, ModelConverter> {
public:
  /// Base class
  using Base = BasicFuncConstrCvt<
    MulConverter_MIP<ModelConverter>, ModelConverter>;
  /// Constructor
  MulConverter_MIP(ModelConverter& mc) : Base(mc) { }
  /// Converted item type
  using ItemType = MulConstraint;

  /// Convert in any context
  Context Convert(const ItemType& mc, int ) {
    const int x = mc.GetArguments()[0], y = mc.GetArguments()[1];
    LinTerms lt;
    if (GetMC().is_binary_var(x) || GetMC().is_binary_var(y)) {
      const bool fx = GetMC().is_binary_var(x);  // b*w = b ? w : 0
      lt.add_term(1.0, AssignResultVarMixedCtx(GetMC(),
                    IfThenConstraint{ {fx ? x : y, fx ? y : x,
                                       int(GetMC().MakeFixedVar(0.0))} }));
    } else {
      const int e = FactorToExpand(GetMC(), x, y);
      lt = e >= 0
          ? ProductByBinaryExpansion(GetMC(), 1.0, e==x ? y : x, e)
          : ProductBySquares(GetMC(), 1.0, x, y);
    }
    lt.add_term(-1.0, mc.GetResultVar());        // r = lt
    GetMC().AddConstraint( LinConEQ( std::move(lt), {0.0} ) );
    return Context::CTX_MIX;
  }

  /// Reuse the stored ModelConverter
  using Base::GetMC;
};


} // namespace mp

#endif // MUL_H
