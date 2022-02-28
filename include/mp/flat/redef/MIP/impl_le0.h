#ifndef IMPL_LE0_H
#define IMPL_LE0_H

#include "mp/flat/redef/redef_base.h"
#include "mp/flat/constraints_std.h"
#include "mp/flat/preprocess.h"

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
    auto ae = indc.to_lhs_expr();
    auto bnds = GetMC().ComputeBoundsAndType(ae);
    ConvertImplicationLE(binvar, indc.get_binary_value(),
                         bnds, std::move(ae));
  }

protected:
  using Base::GetMC;

  /// Linearize (b==val ==> ae<=0) via big-M
  void ConvertImplicationLE(int b, int val,
                   const PreprocessInfoStd& bnds, AffExp ae) {
    /// TODO fail if lb>0 +report .iis if requested
    /// TODO skip if ub<0
    if (bnds.ub() >= GetMC().PracticallyInfty())
      throw ConstraintConversionFailure( "IndicatorInfBound",
          "The redefinition of a (possibly auxiliary) indicator constraint failed"
          " so it had to be passed to the solver."
          " Provide tight bounds on variables entering logical expressions, "
          "or set acc:ind_le=2");
    if (val)                                     // left condition is b==1
      ae += {{bnds.ub(), b}, -bnds.ub()};
    else
      ae += {{-bnds.ub(), b}, 0.0};
    GetMC().AddConstraint(LinConLE(              // Big-M constraint
        (LinTerms&&)ae,               // skip the constant
        -ae.constant_term() ));       // and use it here
  }

};

} // namespace mp

#endif // IMPL_LE0_H
