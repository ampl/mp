#ifndef INDICATOR_QUAD_H
#define INDICATOR_QUAD_H

/**
  * Convert quadratic indicators for MIP
  */

#include "mp/flat/redef/redef_base.h"
#include "mp/flat/constraints_std.h"

namespace mp {

/// Convert IndicatorQuad(LE/EQ/GE)
/// b==val ==> c'x + x'Qx (<=/==/>=) d.
/// The quadratic part is moved into a separate constraint,
/// thus producing a linear indicator of the respective type.
/// TODO make sure the quadratic part is an inequality for LE/GE
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

  /// Conversion.
  /// Substitute constraint body by a new variable.
  void Convert(const ItemType& indc, int ) {
    auto binvar=indc.get_binary_var();
    const auto& body = indc.get_constraint().GetBody();
    assert(body.is_quadratic());
    auto auxvar = GetMC().AssignResultVar2Args(  // auxvar = body + 0.0
          QuadraticFunctionalConstraint{ {body, 0.0} } );
    GetMC().AddConstraint( IndicatorLin{binvar, indc.get_binary_value(),
                                        LinCon{ { {1.0}, {auxvar} },
                                          indc.get_constraint().rhs() }} );
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

/// Typedef IndicatorQuadGEConverter_MIP
template <class MC>
using IndicatorQuadGEConverter_MIP = IndicatorQuadConverter_MIP<MC, 1>;

} // namespace mp

#endif // INDICATOR_QUAD_H
