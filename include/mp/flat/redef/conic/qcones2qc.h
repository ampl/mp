#ifndef QCONES2QC_H
#define QCONES2QC_H

#include "mp/flat/redef/redef_base.h"
#include "mp/flat/constr_std.h"

namespace mp {

/// Converts QCone
template <class ModelConverter>
class QConeConverter :
    public BasicFuncConstrCvt<
      QConeConverter<ModelConverter>, ModelConverter> {
public:
  /// Base class
  using Base = BasicFuncConstrCvt<
    QConeConverter<ModelConverter>, ModelConverter>;
  /// Constructor
  QConeConverter(ModelConverter& mc) : Base(mc) { }
  /// Converted item type
  using ItemType = QuadraticConeConstraint;

  /// Check whether the constraint
  /// needs to be converted despite being accepted by ModelAPI.
  bool IfNeedsConversion(const ItemType& , int ) {
    return GetMC().IfConvertSOCP2QC();
  }

  /// Convert to
  /// (c[0]*x[0])^2 >= sum(i)((c[i]*x[i])^2)
  /// with x[0]>=0.
  void Convert(const ItemType& ac, int ) {
    const auto& x = ac.GetArguments();
    auto c = ac.GetParameters();
    for (auto& coef: c)
      coef *= coef;
    c[0] = -c[0];
    if (!GetMC().is_fixed(x[0])) {
      GetMC().NarrowVarBounds(x[0], 0.0, GetMC().Infty());
      auto qc {QuadConLE{ {{}, {c, x, x}}, {0.0} }};
      GetMC().AddConstraint(std::move(qc));
    } else {
      // Reproduce fixed RHS, better for Mosek & COPT. conic/socp_10
      auto rhs = -c[0] * GetMC().fixed_value(x[0]);
      c.erase(c.begin());
      auto x0 = x;
      x0.erase(x0.begin());
      auto qc {QuadConLE{ {{}, {c, x0, x0}}, {rhs} }};
      GetMC().AddConstraint(std::move(qc));
    }
  }

  /// Reuse the stored ModelConverter
  using Base::GetMC;
};


/// Converts RQCone
template <class ModelConverter>
class RQConeConverter :
    public BasicFuncConstrCvt<
      RQConeConverter<ModelConverter>, ModelConverter> {
public:
  /// Base class
  using Base = BasicFuncConstrCvt<
    RQConeConverter<ModelConverter>, ModelConverter>;
  /// Constructor
  RQConeConverter(ModelConverter& mc) : Base(mc) { }
  /// Converted item type
  using ItemType = RotatedQuadraticConeConstraint;

  /// Check whether the constraint
  /// needs to be converted despite being accepted by ModelAPI.
  bool IfNeedsConversion(const ItemType& , int ) {
    return GetMC().IfConvertSOCP2QC();
  }

  /// Convert to 2(c[0]*x[0]*c[1]*x[1]) >= sum(i>=2)((c[i]*x[i])^2).
  void Convert(const ItemType& ac, int ) {
    const auto& x = ac.GetArguments();
    const auto& c = ac.GetParameters();
    std::vector<int> x1{x.begin()+1, x.end()};
    auto x2 = x1;
    x2[0] = x[0];
    std::vector<double> c12{c.begin()+1, c.end()};
    c12[0] *= -2.0*c[0];
    for (auto i=c12.size(); --i; )
      c12[i] *= c12[i];
    // Reproduce linear term, for Mosek & COPT. conic/socp_11.
    if (GetMC().is_fixed(x1[0])) {
      std::vector<double> clin
          = { c12[0] * GetMC().fixed_value(x1[0]) };
      std::vector<int> xlin = { x2[0] };
      x1.erase(x1.begin());
      x2.erase(x2.begin());
      c12.erase(c12.begin());
      GetMC().AddConstraint(
            QuadConLE{ {{clin, xlin}, {c12, x1, x2}}, {0.0} });
    } else
      if (GetMC().is_fixed(x2[0])) {
        std::vector<double> clin
            = { c12[0] * GetMC().fixed_value(x2[0]) };
        std::vector<int> xlin = { x1[0] };
        x1.erase(x1.begin());
        x2.erase(x2.begin());
        c12.erase(c12.begin());
        GetMC().AddConstraint(
              QuadConLE{ {{clin, xlin}, {c12, x1, x2}}, {0.0} });
      } else {
        auto qc = QuadConLE{ {{}, {c12, x1, x2}}, {0.0} };
        GetMC().AddConstraint(std::move(qc));
      }
  }

  /// Reuse the stored ModelConverter
  using Base::GetMC;
};

} // namespace mp

#endif // QCONES2QC_H
