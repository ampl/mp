#ifndef DIV_H
#define DIV_H

#include "mp/flat/redef/redef_base.h"
#include "mp/flat/constr_std.h"

namespace mp {

/// Converts DIV (float division) for MIP
template <class ModelConverter>
class DivConverter_MIP :
    public BasicFuncConstrCvt<
      DivConverter_MIP<ModelConverter>, ModelConverter> {
public:
  /// Base class
  using Base = BasicFuncConstrCvt<
    DivConverter_MIP<ModelConverter>, ModelConverter>;
  /// Constructor
  DivConverter_MIP(ModelConverter& mc) : Base(mc) { }
  /// Converted item type
  using ItemType = DivConstraint;

  /// Check whether the constraint
  /// needs to be converted despite being accepted by ModelAPI.
  bool IfNeedsConversion(const ItemType& dc, int ) {
    return GetMC().is_fixed(dc.GetResultVar())
        || GetMC().is_fixed(dc.GetArguments()[1]);
  }

  /// Convert in both contexts (full reification)
  Context Convert(const ItemType& dc, int ) {
    if (GetMC().is_fixed(dc.GetResultVar()))
      ConvertWithConstResult(dc);
    else if (GetMC().is_fixed(dc.GetArguments()[1]))
      ConvertWithConstDivisor(dc);
    else
      ConvertWithNonConstDivisor(dc);
    return Context::CTX_MIX;
  }


protected:
  /// Convert with constant result
  void ConvertWithConstResult(const ItemType& dc) {
    auto r = dc.GetResultVar();
    const auto& args = dc.GetArguments();
    /// r = arg0 / arg1,    arg0 = arg1 * r;
    auto lt = LinTerms{
        { GetMC().fixed_value(r), -1.0 },
        { args[1], args[0] }
    };
    GetMC().AddConstraint (LinConEQ( std::move(lt), { 0.0 } ) );
    EnsureDivisorNot0(r, args);
  }


  /// Convert with constant divisor
  void ConvertWithConstDivisor(const ItemType& dc) {
    auto r = dc.GetResultVar();
    const auto& args = dc.GetArguments();
    /// r = arg0 / arg1,    arg0 = arg1 * r;
    auto lt = LinTerms{
      { GetMC().fixed_value(args[1]), -1.0 },
      { r, args[0] }
    };
    GetMC().AddConstraint (LinConEQ( std::move(lt), { 0.0 } ) );
  }

  /// Convert with non-constant divisor
  void ConvertWithNonConstDivisor(const ItemType& dc) {
    auto r = dc.GetResultVar();
    const auto& args = dc.GetArguments();
    /// r = arg0 / arg1,    arg0 = arg1 * r;
    auto lt = LinTerms{ {-1.0}, {args[0]} };
    auto qt = QuadTerms{ { 1.0 }, {args[1]}, {r} };
    auto qlt = QuadAndLinTerms{ std::move(lt), std::move(qt) };
    GetMC().AddConstraint(QuadConRange( std::move(qlt), { 0.0, 0.0 } ));
    EnsureDivisorNot0(r, args);
  }

  /// See if need to constrain divisor!=0
  /// @todo consider resc/ dividend bounds - in Prepro
  void EnsureDivisorNot0(int res, VarArray2 args) {
    auto lb1 = GetMC().lb(args[1]);
    auto ub1 = GetMC().ub(args[1]);
    if (lb1 <= 0.0 && ub1 >= 0.0) {
      /// Creating "args[1] != 0" via "not (args[1]==0)"...
      auto arg1is0 = GetMC().AssignResultVar2Args(// arg1is0 = (arg1==0)
          CondLinConEQ( { {{1.0}, {args[1]}}, 0.0 } ) );
      auto r_not = GetMC().AssignResultVar2Args(  // r_not = (arg1!=0)
          NotConstraint( {arg1is0} ));
      GetMC().FixAsTrue(r_not);                   // TODO Simplify these API
    }
  }

  /// Reuse the stored ModelConverter
  using Base::GetMC;
};

} // namespace mp

#endif // DIV_H
