#ifndef MUL_H
#define MUL_H

#include "mp/error.h"
#include "mp/flat/redef/redef_base.h"
#include "mp/flat/constr_std.h"

namespace mp {

/// Convert all quadratic terms in the body
template <class ModelConverter, int sens>
class QCConverter_MIP :
    public BasicItemConverter<ModelConverter> {
public:
  /// Base class
  using Base = BasicItemConverter<ModelConverter>;
  /// Constructor
  QCConverter_MIP(ModelConverter& mc) : Base(mc) { }
  /// Converted item type
  using ItemType = QuadConRhs<sens>;

  /// Conversion
  void Convert(const ItemType& qc, int ) {
    LinearizeQPTerms(qc);
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
    // Sort linear body
    lin_terms.sort_terms();
    GetMC().AddConstraint( LinConRhs< sens >{
                             lin_terms, qc.GetRhsOrRange() } );
  }

  LinTerms LinearizeQPTerm(double c, int x, int y) {
    return LinearizeProductWithBinaryVar(c, x, y);
  }

  LinTerms LinearizeProductWithBinaryVar(double c, int x, int y) {
    LinTerms lt;
    MP_ASSERT_ALWAYS(GetMC().is_binary_var(x) ||
                     GetMC().is_binary_var(y),
                     "Can only convert product with a binary variable");
    bool is_x_bin = GetMC().is_binary_var(x);
    auto i_bin = is_x_bin ? x : y;       // index of the (chosen) binary var
    auto i_other = is_x_bin ? y : x;                  // the other var index
    auto x_ifthen = GetMC().template
        AssignResultVar2Args( IfThenConstraint {{      // no context as of now
                                                      i_bin,
                                                      i_other,
                                                      GetMC().MakeFixedVar(0.0)
                              }} );
    lt.add_term(c, x_ifthen);
    return lt;
  }
};


/// Typedef MulCvt<QuadConLE>
template <class MC>
using MulCvtLE_MIP = QCConverter_MIP<MC, -1>;

/// Typedef MulCvt<QuadConEQ>
template <class MC>
using MulCvtEQ_MIP = QCConverter_MIP<MC,  0>;

/// Typedef MulCvt<QuadConGE>
template <class MC>
using MulCvtGE_MIP = QCConverter_MIP<MC,  1>;


} // namespace mp

#endif // MUL_H
