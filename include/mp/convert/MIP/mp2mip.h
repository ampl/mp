#ifndef MP2MIP_H
#define MP2MIP_H

#include "mp/convert/converter_flat.h"

namespace mp {

/// MPToMIPConverter: one of the converters requiring a "minimal" output interface
template <class Impl, class Backend,
          class Model = BasicModel<std::allocator<char> > >
class MPToMIPConverter
    : public BasicMPFlatConverter<Impl, Backend, Model>
{
public:
  using BaseConverter = BasicMPFlatConverter<Impl, Backend, Model>;

  ///////////////////// SPECIALIZED CONSTRAINT CONVERTERS //////////////////
  USE_BASE_CONSTRAINT_CONVERTERS( BaseConverter )           // reuse default ones

  template <int sense, class MinOrMaxConstraint>
  void ConvertMinOrMax(const MinOrMaxConstraint& mc) {
    const auto& args = mc.GetArguments();
    const int nargs = args.size();
    const auto flags = this->AddVars(nargs, 0.0, 1.0, var::Type::INTEGER);   // binary flags
    this->AddConstraint(LinearConstraint(std::vector<double>(nargs, 1.0),    // sum of the flags >= 1
                        flags, 1.0, this->PlusInfinity()));
    const auto resvar = mc.GetResultVar();
    for (int i=0; i<nargs; ++i) {
      this->AddConstraint(LinearConstraint({1.0*sense, -1.0*sense},
                          {args[i], resvar}, this->MinusInfinity(), 0.0));
      this->AddConstraint(IndicatorConstraintLinLE{flags[i], 1,
                          {1.0*sense, -1.0*sense}, {resvar, args[i]}, 0.0});
    }
  }

  void Convert(const MaximumConstraint& mc) {
    ConvertMinOrMax<1>(mc);
  }

  void Convert(const MinimumConstraint& mc) {
    ConvertMinOrMax<-1>(mc);
  }

  void Convert(const LEConstraint& lec) {
    auto& m = this->GetModel();
    if (m.is_fixed(lec.GetResultVar()))
      throw std::logic_error("LEConstraint: result fixed, not implemented");
    assert(!lec.GetContext().IsNone());
    if (lec.GetContext().HasPositive())
      ConvertImplied(lec);
    if (lec.GetContext().HasNegative())
      ConvertReverseImplied(lec);
  }

  void ConvertImplied(const LEConstraint& lec) {
    auto& m = this->GetModel();
    double d=0.0;
    LinearExprUnzipper le;
    const auto& args = lec.GetArguments();
    if (m.is_fixed(args[0]))
      d -= m.fixed_value(args[0]);
    else
      le.AddTerm(args[0], 1.0);
    if (m.is_fixed(args[1]))
      d += m.fixed_value(args[1]);
    else
      le.AddTerm(args[1], -1.0);
    if (0==le.num_terms()) {
      if (d<0.0)
        m.narrow_var_bounds(lec.GetResultVar(), 0.0, 0.0);
    } else {
      this->AddConstraint(IndicatorConstraintLinLE(
                            lec.GetResultVar(), 1, le.c_, le.v_, d));
    }
  }

  void ConvertReverseImplied(const LEConstraint& lec) {
    auto& m = this->GetModel();
    double d=-1.0;
    LinearExprUnzipper le;
    const auto& args = lec.GetArguments();
    if (m.is_fixed(args[0]))
      d += m.fixed_value(args[0]);
    else
      le.AddTerm(args[0], -1.0);
    if (m.is_fixed(args[1]))
      d -= m.fixed_value(args[1]);
    else
      le.AddTerm(args[1], 1.0);
    if (0==le.num_terms()) {
      if (d<0.0)
        m.narrow_var_bounds(lec.GetResultVar(), 1.0, 1.0);
    } else {
      this->AddConstraint(IndicatorConstraintLinLE(
                            lec.GetResultVar(), 0, le.c_, le.v_, d));
    }
  }

  void Convert(const IndicatorConstraintLinLE& indc) {
    auto binvar=indc.get_binary_var();
    if (indc.is_binary_value_1())                  /// If binval==1, complement the variable
      binvar = this->MakeComplementVar(binvar);
    /// Convert indc's linear inequality to 'cmpvar<=0'
    int cmpvar = MP_DISPATCH( Convert2Var(indc.to_lhs_affine_expr()) );
    this->AddConstraint(LinearConstraint(          /// Big-M constraint cmpvar <= ub(cmpvar)*binvar
        {1.0, -this->ub(cmpvar)}, {cmpvar, binvar}, this->MinusInfinity(), 0.0));
  }


};

} // namespace mp

#endif // MP2MIP_H
