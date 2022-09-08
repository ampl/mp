#ifndef INDICATOR_GE_H
#define INDICATOR_GE_H

#include "mp/flat/redef/redef_base.h"
#include "mp/flat/constr_std.h"

namespace mp {

/// Convert IndicatorLinGE
/// b==val ==> c'x>=d
template <class ModelConverter>
class IndicatorLinGEConverter_MIP :
    public BasicItemConverter<ModelConverter> {
public:
  /// Base class
  using Base = BasicItemConverter<ModelConverter>;
  /// Constructor
  IndicatorLinGEConverter_MIP(ModelConverter& mc) : Base(mc) { }
  /// Converted item type
  using ItemType = IndicatorConstraintLinGE;

  /// Conversion
  void Convert(const ItemType& indc, int ) {
    auto binvar=indc.get_binary_var();
    auto bnds = GetMC().ComputeBoundsAndType(
          indc.get_constraint().GetBody());
    ConvertImplicationGE(binvar, indc.get_binary_value(),
                         bnds.lb(), indc.get_constraint());
  }

protected:
  using Base::GetMC;

  /// Linearize (b==val ==> c'x>=d) via big-M
  void ConvertImplicationGE(int b, int val,
                   double body_lb, LinConGE con) {
    if (body_lb <= GetMC().PracticallyMinusInfty())
      throw ConstraintConversionFailure( "IndicatorInfBound",
          "The redefinition of a (possibly auxiliary) indicator constraint"
          " 'bin_var==value ==> c'x>=d' failed"
          " so it will be passed to the solver natively if supported."
          " Provide tight bounds on variables entering logical expressions,"
          " or, if available, set acc:indge=2 for native handling");
    if (0==val)                                // left condition is b==0
      con.GetBody().add_term(-body_lb+con.rhs(), b);
    else {
      con.GetBody().add_term(body_lb-con.rhs(), b);
      con.set_rhs(body_lb);
    }
    GetMC().AddConstraint(con);                // Big-M constraint
  }

};

} // namespace mp

#endif // INDICATOR_GE_H
