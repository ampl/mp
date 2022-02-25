#ifndef EQ0_H
#define EQ0_H

#include "mp/flat/redef/redef_base.h"
#include "mp/flat/constraints_std.h"

namespace mp {

/// Converts EQ0Constraint for MIP
template <class ModelConverter>
class EQ0Converter_MIP :
    public BasicFuncConstrCvt<
      EQ0Converter_MIP<ModelConverter>, ModelConverter> {
public:
  /// Base class
  using Base = BasicFuncConstrCvt<
    EQ0Converter_MIP<ModelConverter>, ModelConverter>;
  /// Constructor
  EQ0Converter_MIP(ModelConverter& mc) : Base(mc) { }
  /// Converted item type
  using ItemType = EQ0Constraint;

  /// Reimplement Convert(): add filter
  void Convert(const ItemType& eq0c, int i) {
    assert(!eq0c.GetContext().IsNone());
    const auto& args = eq0c.GetArguments();
    if (1<args.size() ||
        !GetMC().IfUseEqualityEncodingForVar(
          args.var(0))) {
      Base::Convert(eq0c, i);
    } // else, using unary encoding whose flags are,
  }   // in the fixed case, fixed by PropagateResult()

  /// Convert in positive context
  /// resvar==1 => c'x==d
  void ConvertCtxPos(const ItemType& eq0c, int ) {
    const auto& ae = eq0c.GetArguments();
    const auto res = eq0c.GetResultVar();
    if (ae.is_constant()) {         // TODO consider resvar+context
      if (std::fabs(ae.constant_term()) != 0.0)
        GetMC().NarrowVarBounds(res, 0.0, 0.0);
    } else {
      if (GetMC().is_fixed(res)) {
        if (GetMC().fixed_value(res)) {  // fixed to 1
          GetMC().AddConstraint( ExtractConstraint(eq0c) );
        } // else, skip
      }
      GetMC().AddConstraint(IndicatorConstraintLinEQ(
                              res, 1,
                              ExtractConstraint(eq0c)));
    }
  }

  /// Convert in negative context
  /// resvar==0 ==> c'x!=d
  void ConvertCtxNeg(const ItemType& eq0c, int ) {
    const auto& ae = eq0c.GetArguments();
    const auto res = eq0c.GetResultVar();
    if (ae.is_constant()) {
      if (std::fabs(ae.constant_term()) != 0.0)
        GetMC().NarrowVarBounds(res, 1.0, 1.0);
      // TODO use resvar + context
    } else if ( !GetMC().is_fixed(res) ||   // not fixed, or
                !GetMC().fixed_value(res) ) // fixed to 0
    { // TODO We are in MIP so doing algebra, not DisjunctiveConstr. Why?
      // Well in party1.mod, although this results in more fixed variables,
      // Gurobi 9.5 runs 31s vs 91s.
      auto ae = eq0c.GetArguments();
      auto newvars = GetMC().AddVars_returnIds(2, 0.0, 1.0, var::INTEGER);
      newvars.push_back( res );
      GetMC().AddConstraint( LinConGE(   // b1+b2+resvar >= 1
                                         {{1.0, 1.0, 1.0}, newvars},
                                         1.0 ) );
      auto bNt = GetMC().ComputeBoundsAndType(ae);
      double cmpEps = GetMC().ComparisonEps( bNt.get_result_type() );
      {
        GetMC().AddConstraint(IndicatorConstraintLinLE(
                                newvars[0], 1,
                              { {ae.coefs(), ae.vars()},
                                -ae.constant_term() - cmpEps }));
      }
      ae.negate();
      GetMC().AddConstraint(IndicatorConstraintLinLE(
                              newvars[1], 1,
                            { {ae.coefs(), ae.vars()},
                              -ae.constant_term() - cmpEps }));
    } // else, skip
  }

  /// Reuse the stored ModelConverter
  using Base::GetMC;
};

} // namespace mp

#endif // EQ0_H
