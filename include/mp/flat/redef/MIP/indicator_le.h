#ifndef INDICATOR_LE_H
#define INDICATOR_LE_H

#include "mp/flat/redef/redef_base.h"
#include "mp/flat/constr_std.h"

namespace mp {

/// Convert IndicatorLinLE
/// b==val ==> c'x<=d
template <class ModelConverter>
class IndicatorLinLEConverter_MIP :
    public BasicItemConverter<ModelConverter> {
public:
  /// Base class
  using Base = BasicItemConverter<ModelConverter>;
  /// Constructor
  IndicatorLinLEConverter_MIP(ModelConverter& mc) : Base(mc) { }
  /// Converted item type
  using ItemType = IndicatorConstraintLinLE;

  /// Conversion
  void Convert(const ItemType& indc, int ) {
    auto binvar=indc.get_binary_var();
    auto bnds = GetMC().ComputeBoundsAndType(
          indc.get_constraint().GetBody());
    ConvertImplicationLE(binvar, indc.get_binary_value(),
                         bnds.ub(), indc.get_constraint());
  }

protected:
  using Base::GetMC;

  /// Linearize (b==val ==> c'x<=d) via big-M
  void ConvertImplicationLE(int b, int val,
                   double body_ub, LinConLE con) {
    /// TODO fail if lb>0 +report .iis if requested
    /// TODO skip if ub<0
    if (body_ub >= GetMC().PracticallyInf()) {
      if ((body_ub = GetMC().bigMDefault())<=0.0)
        throw ConstraintConversionFailure(
            "IndicatorInfBound",
            "Set bounds on variables\n"
            "participating in logical expressions,\n"
            "or use option cvt:bigM (with caution).\n"
            "See more: mp.ampl.com/modeling-numerics.html");
    }
    if (body_ub != con.rhs()) {
      if (0==val)                                // left condition is b==0
        con.GetBody().add_term(-body_ub+con.rhs(), b);
      else {
        con.GetBody().add_term(body_ub-con.rhs(), b);
        con.set_rhs(body_ub);
      }
      GetMC().AddConstraint(con);                // Big-M constraint
    }
  }

};

} // namespace mp

#endif // INDICATOR_LE_H
