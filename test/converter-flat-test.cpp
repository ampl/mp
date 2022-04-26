
#include "mp/easy-modeler.h"

#include "converter-flat-test.h"

namespace {

using InterfaceTester_MaxConstraint =
             InterfaceTesterWithBackendAcceptingConstraints<mp::MaxConstraint>;

TEST_F(InterfaceTester_MaxConstraint, MaximumConstraintIsPassedToBackend) {
  auto con=GetModel().AddCon(5.0, 5.0);
  const auto args = GetInterface().AddVars(3, -1.0, 11.0);
  con.set_nonlinear_expr(MakeIterated(GetModel(), mp::expr::MAX, args ));
  GetInterface().ConvertModel();
  ASSERT_HAS_CONSTRAINT( GetBackend(), mp::MaxConstraint(3, args) );
}


/////////////////////////////// Quadratics //////////////////////////////////
using InterfaceTester_QuadraticConstraint =
             InterfaceTesterWithBackendAcceptingConstraints<mp::QuadConRange>;

/// EasyModeler syntax
TEST_F(InterfaceTester_QuadraticConstraint, QuadExprIsMultipliedOutAndInlined) {
  auto modeler = mp::MakeEasyModeler(GetModel());
  auto x = modeler.AddVars(3, -1.0, 11.0);
  modeler.AddAlgCon(5.0,
                    (5*x[0]+3) * (6*x[1]+2) + x[2],
                    INFINITY);
  GetInterface().ConvertModel();
  const auto xi = modeler.GetVarIndices(x);
  auto lt = mp::LinTerms{ {10.0, 18.0, 1.0}, {xi[0], xi[1], xi[2]} };
  auto qt = mp::QuadTerms{ {30.0}, {xi[0]}, {xi[1]} };
  ASSERT_HAS_CONSTRAINT( GetBackend(), mp::QuadConRange(
    { std::move(lt), std::move(qt) },
    {-1.0, INFINITY} ) );
}

TEST_F(InterfaceTester_QuadraticConstraint, Pow2_isMultipliedOutAndInlined) {
  auto modeler = mp::MakeEasyModeler(GetModel());
  auto x = modeler.AddVars(3, -1.0, 11.0);
  modeler.AddAlgCon(5.0,
                    ((8*x[0] + 2*x[1] + 3)^2)            // C++ operator precedence
                      + 3.5*x[2],
                    INFINITY);
  GetInterface().ConvertModel();
  const auto xi = modeler.GetVarIndices(x);
  auto lt = mp::LinTerms{ {48, 12, 3.5}, {xi[0], xi[1], xi[2]} };
  auto qt = mp::QuadTerms{ {64, 4, 32}, {xi[0], xi[1], xi[0]}, {xi[0], xi[1], xi[1]} };
  ASSERT_HAS_CONSTRAINT( GetBackend(), mp::QuadConRange(
    { std::move(lt), std::move(qt) },
    {-4.0, INFINITY} ) );
}


TEST_F(InterfaceTester_QuadraticConstraint, QuadConstraintIsPassedToBackend__OldSyntax) {
  auto con=GetModel().AddCon(5.0, 5.0);
  const auto args = GetInterface().AddVars(2, -1.0, 11.0);
  con.set_nonlinear_expr(GetModel().MakeBinary(mp::expr::MUL,
                                               GetModel().MakeVariable( args[0] ),
                                               GetModel().MakeVariable( args[1] ) ) );
  GetInterface().ConvertModel();
  auto qt = mp::QuadTerms{ {1.0}, {args[0]}, {args[1]} };
  ASSERT_HAS_CONSTRAINT( GetBackend(), mp::QuadConRange(
                           { mp::LinTerms{}, std::move(qt) },
    {5.0, 5.0} ) );
}

} // namespace
