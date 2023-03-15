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
    return 0==GetMC().IfPassSOCPCones() &&
        0!=GetMC().IfPassQuadCon();
  }

  /// Convert to (c[0]*x[0])^2 >= sum(i)((c[i]*x[i])^2).
  void Convert(const ItemType& ac, int ) {
    const auto& x = ac.GetArguments();
    auto c = ac.GetParameters();
    for (auto& coef: c)
      coef *= coef;
    c[0] = -c[0];
    GetMC().NarrowVarBounds(x[0], 0.0, GetMC().Infty());
    auto qc {QuadConLE{ {{}, {c, x, x}}, {0.0} }};
    GetMC().AddConstraint(std::move(qc));
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
    return 0==GetMC().IfPassSOCPCones() &&
        0!=GetMC().IfPassQuadCon();
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
    auto qc = QuadConLE{ {{}, {c12, x1, x2}}, {0.0} };
    GetMC().AddConstraint(std::move(qc));
  }

  /// Reuse the stored ModelConverter
  using Base::GetMC;
};

} // namespace mp

#endif // QCONES2QC_H
