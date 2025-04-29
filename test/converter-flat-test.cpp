
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
///
/// Probably this tests more the Problem Flattener,
/// or the whole chain up to ModelAPI
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

TEST_F(InterfaceTester_QuadraticConstraint, QP2Pass0) {
  auto modeler = mp::MakeEasyModeler(GetModel());
  auto vars0 = modeler.AddVars(4, -1.0, 11.0);
  auto x = vars0[0], y = vars0[1], z = vars0[2], t = vars0[3];
  modeler.AddAlgCon(-INFINITY,
                    (((5*x-2)^2) - (4*x-3)*(2*y-1) + ((3*x+2*z+8)^2))
                    * ((5*t-2)^2) / 4             // C++ operator precedence
                        + ((2*x-8)^2) - (y-z)*z + (x-3)*(x-2*z+5),
                    200);
  GetEnv().ParseOptionString("qp2pass=0");
  GetInterface().ConvertModel();
  {
    const auto vars0i = modeler.GetVarIndices(vars0);
    auto x = vars0i[0], y = vars0i[1], z = vars0i[2], t = vars0i[3];
    int u = t+1, v = t+2;            // aux vars
    // Check that we have the left factor of the 4th-degree term
    ASSERT_HAS_CONSTRAINT(
        GetBackend(),
        mp::QuadConEQ(
            { { {32, 6, 32, -1}, {x, y, z, v} },    // v: do we flatten the 2nd factor 1st?
             { {34, 4, -8, 12}, {x, z, x, x}, {x, z, y, z} } },
            {-65} ) );
    // The right factor
    ASSERT_HAS_CONSTRAINT(
        GetBackend(),
        mp::QuadConEQ(
            { { {-20, -1}, {t, u} },
             { {25}, {t}, {t} } },
            {-4} ) );
    // The upper-level constraint
    ASSERT_HAS_CONSTRAINT(
        GetBackend(),
        mp::QuadConRange(
            { { {-30, 6}, {x, z} },
             { {0.25, 5, 1, -2, -1}, {u, x, z, x, y}, {v, x, z, z, z} } },
            {-INFINITY, 151} ) );
  }
}

/// @todo How to join common code?
TEST_F(InterfaceTester_QuadraticConstraint, QP2Pass1) {
  auto modeler = mp::MakeEasyModeler(GetModel());
  auto vars0 = modeler.AddVars(4, -1.0, 11.0);
  auto x = vars0[0], y = vars0[1], z = vars0[2], t = vars0[3];
  modeler.AddAlgCon(-INFINITY,
                    (((5*x-2)^2) - (4*x-3)*(2*y-1) + ((3*x+2*z+8)^2))
                            * ((5*t-2)^2) / 4             // C++ operator precedence
                        + ((2*x-8)^2) - (y-z)*z + (x-3)*(x-2*z+5),
                    200);
  GetEnv().ParseOptionString("qp2pass=1");
  GetInterface().ConvertModel();
  {
    const auto vars0i = modeler.GetVarIndices(vars0);
    auto x = vars0i[0], y = vars0i[1], z = vars0i[2], t = vars0i[3];
    int u = t+1, v = t+2;            // aux vars
    // Check that we have the left factor of the 4th-degree term
    ASSERT_HAS_CONSTRAINT(
        GetBackend(),
        mp::QuadConEQ(
            { { {32, 6, 32, -1}, {x, y, z, v} },    // v: do we flatten the 2nd factor 1st?
             { {34, 4, -8, 12}, {x, z, x, x}, {x, z, y, z} } },
            {-65} ) );
    // The right factor
    ASSERT_HAS_CONSTRAINT(
        GetBackend(),
        mp::QuadConEQ(
            { { {-20, -1}, {t, u} },
             { {25}, {t}, {t} } },
            {-4} ) );
    // The upper-level constraint
    ASSERT_HAS_CONSTRAINT(
        GetBackend(),
        mp::QuadConRange(
            { { {-30, 6}, {x, z} },
             { {0.25, 5, 1, -2, -1}, {u, x, z, x, y}, {v, x, z, z, z} } },
            {-INFINITY, 151} ) );
  }
}


} // namespace
