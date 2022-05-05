#ifndef INDICATOR_QUAD_H
#define INDICATOR_QUAD_H

/**
  * Convert quadratic indicators for MIP
  */

#include "mp/flat/redef/redef_base.h"
#include "mp/flat/constraints_std.h"

namespace mp {

/// Convert IndicatorQuad(LE/EQ)
/// b==val ==> c'x + x'Qx (<=/==) d.
/// The quadratic part is moved into a separate constraint,
/// thus producing a linear indicator of the respective type.
/// TODO make sure the quadratic part is an inequality for LE
template <class ModelConverter, int sens>
class IndicatorQuadConverter_MIP :
    public BasicItemConverter<ModelConverter> {
public:
  /// Base class
  using Base = BasicItemConverter<ModelConverter>;
  /// Constructor
  IndicatorQuadConverter_MIP(ModelConverter& mc) : Base(mc) { }
  /// Underlying algebraic constraint
  using QuadCon = AlgebraicConstraint< QuadAndLinTerms, AlgConRhs<sens> >;
  /// Resulting implied linear constraint
  using LinCon = AlgebraicConstraint< LinTerms, AlgConRhs<sens> >;
  /// Converted item type
  using ItemType = IndicatorConstraint<QuadCon>;
  /// Resulting linear indicator
  using IndicatorLin = IndicatorConstraint<LinCon>;

  /// Conversion
  void Convert(const ItemType& indc, int ) {
    auto binvar=indc.get_binary_var();
    auto body = indc.get_constraint().GetBody();
    assert(!body.GetQPTerms().empty());
    auto auxvar = GetMC().AssignResultVar2Args(
          QuadraticFunctionalConstraint{
            { { {}, std::move(body.GetQPTerms()) }, 0.0} } );
    auto lt = body.GetLinTerms();
    lt.add_term(1.0, auxvar);
    GetMC().AddConstraint( IndicatorLin{binvar, indc.get_binary_value(),
                                        {lt, indc.get_constraint().rhs()}} );
  }

protected:
  using Base::GetMC;

};


/// Typedef IndicatorQuadLEConverter_MIP
template <class MC>
using IndicatorQuadLEConverter_MIP = IndicatorQuadConverter_MIP<MC, -1>;

/// Typedef IndicatorQuadEQConverter_MIP
template <class MC>
using IndicatorQuadEQConverter_MIP = IndicatorQuadConverter_MIP<MC, 0>;

} // namespace mp

#endif // INDICATOR_QUAD_H
