
#include "mp/easy-modeler.h"

#include "converter-flat-test.h"

namespace {

using InterfaceTester_MaxConstraint =
             InterfaceTesterWithBackendAcceptingConstraints<mp::MaxConstraint>;

TEST_F(InterfaceTester_MaxConstraint, AcceptanceQueriesKeepVirtualAndDirectResults) {
  const auto& keeper = GetInterface().GetFlatCvt().GetConstraintKeeper(
      (const mp::MaxConstraint*)nullptr);
  const mp::BasicConstraintKeeper& base = keeper;
  const auto& api = GetBackend();
  EXPECT_EQ(mp::ConstraintAcceptanceLevel::Recommended,
            keeper.GetModelAPIAcceptance());
  EXPECT_EQ(keeper.GetModelAPIAcceptance(), base.GetModelAPIAcceptance(api));
  EXPECT_EQ(mp::ExpressionAcceptanceLevel::NotAccepted,
            keeper.GetModelAPIAcceptanceEXPR());
  EXPECT_EQ(keeper.GetModelAPIAcceptanceEXPR(), base.GetModelAPIAcceptanceEXPR(api));
}

TEST_F(InterfaceTester_MaxConstraint, MaximumConstraintIsPassedToBackend) {
  auto con=GetModel().AddCon(5.0, 5.0);
  const auto args = GetInterface().AddVars(3, -1.0, 11.0);
  con.set_nonlinear_expr(MakeIterated(GetModel(), mp::expr::MAX, args ));
  GetInterface().ConvertModel();
  ASSERT_HAS_CONSTRAINT( GetBackend(), mp::MaxConstraint(3, args) );
}

/// A repeated solve (e.g., a second AMPLSSolve()) starts the solve
/// iterations again; before, the second run did not solve at all.
TEST_F(InterfaceTester_MaxConstraint, RepeatedSolveRunsAgain) {
  GetInterface().AddVars(1, 0.0, 1.0);
  GetInterface().ConvertModel();
  auto stt = []() { return mp::sol::SOLVED; };
  auto sol = []() { return mp::Solution{}; };
  auto& itf = GetInterface();
  itf.StartSolveIterations();
  EXPECT_TRUE(itf.PrepareSolveIteration(stt, sol));
  EXPECT_FALSE(itf.PrepareSolveIteration(stt, sol));
  itf.StartSolveIterations();                  // the second solve
  EXPECT_TRUE(itf.PrepareSolveIteration(stt, sol));
  EXPECT_FALSE(itf.PrepareSolveIteration(stt, sol));
}


/// Bounds of integer variables are rounded inwards within
/// max(cvt:pre:roundintbnd, 16 ulp(b)): rounding errors of computed
/// bounds are removed, real fractional parts are kept.
TEST_F(InterfaceTester_MaxConstraint, RoundIntBoundUsesUlpTolerance) {
  auto& m = GetInterface().GetFlatCvt().GetModel();
  const int v = m.AddVar__basic(-1e10, 1e10, mp::var::INTEGER);
  const double b = 16800000.7 / 0.1;           // 168000006.99999997
  EXPECT_EQ(168000007.0, m.RoundIntBound(v, b, true));
  EXPECT_EQ(168000007.0, m.RoundIntBound(v, -(-b), false));
  EXPECT_EQ(-168000007.0, m.RoundIntBound(v, -b, false));
  // 1e9 + 5e-5 is far above rounding error: lb 1e9+1, ub 1e9
  EXPECT_EQ(1e9 + 1.0, m.RoundIntBound(v, 1e9 + 5e-5, true));
  EXPECT_EQ(1e9, m.RoundIntBound(v, 1e9 + 5e-5, false));
  EXPECT_EQ(600000000.0, m.RoundIntBound(v, 599999999.5, true));
  EXPECT_EQ(599999999.0, m.RoundIntBound(v, 599999999.5, false));
  // Small values: the absolute tolerance (default 1e-9)
  EXPECT_EQ(3.0, m.RoundIntBound(v, 3.0 + 1e-10, false));
  EXPECT_EQ(3.0, m.RoundIntBound(v, 3.0 - 1e-10, true));
  EXPECT_EQ(4.0, m.RoundIntBound(v, 3.0 + 1e-8, true));
  EXPECT_EQ(2.0, m.RoundIntBound(v, 3.0 - 1e-8, false));
  // An exact integer is returned as is, also where the tolerance (0.5)
  // is below the spacing of doubles: 2^52 + 1 -/+ 0.5 is not
  // representable and would be rounded to the even neighbour
  const double big = 4503599627370497.0;       // 2^52 + 1
  EXPECT_EQ(big, m.RoundIntBound(v, big, true));
  EXPECT_EQ(big, m.RoundIntBound(v, big, false));
  EXPECT_EQ(-big, m.RoundIntBound(v, -big, true));
  EXPECT_EQ(-big, m.RoundIntBound(v, -big, false));
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


/// Backend accepting a^x natively
using InterfaceTester_ExpAConstraint =
             InterfaceTesterWithBackendAcceptingConstraints<mp::ExpAConstraint>;

/// 1^x is the constant 1. As a^x with a = 1, its piecewise-linear
/// approximation would divide by log(a) = 0. 2^x stays a^x.
TEST_F(InterfaceTester_ExpAConstraint, PowOfBase1IsConstant) {
  auto modeler = mp::MakeEasyModeler(GetModel());
  auto x = modeler.AddVars(2, -2.0, 3.0);
  modeler.AddAlgCon(-INFINITY, x[0] + (1.0 ^ x[1]), 5.0);
  modeler.AddAlgCon(-INFINITY, x[0] + (2.0 ^ x[1]), 9.0);
  GetInterface().ConvertModel();
  const auto& expa = GetBackend().GetCons(
      mp::ExpAConstraint(mp::ExpAConstraint::Arguments{ 0 },
                         mp::ExpAConstraint::Parameters{ 2.0 }));
  ASSERT_EQ(1u, expa.size());
  EXPECT_EQ(2.0, expa[0].GetParameters()[0]);
}

/// The same through VisitPow(): a base that is a constant only after
/// conversion (x0 - x0 + 1).
TEST_F(InterfaceTester_ExpAConstraint, PowOfConvertedBase1IsConstant) {
  auto modeler = mp::MakeEasyModeler(GetModel());
  auto x = modeler.AddVars(2, -2.0, 3.0);
  modeler.AddAlgCon(-INFINITY,
                    x[0] + ((x[0] - x[0] + 1.0) ^ (x[1]*x[1])), 5.0);
  modeler.AddAlgCon(-INFINITY,
                    x[0] + ((x[0] - x[0] + 2.0) ^ (x[1]*x[1])), 9.0);
  GetInterface().ConvertModel();
  const auto& expa = GetBackend().GetCons(
      mp::ExpAConstraint(mp::ExpAConstraint::Arguments{ 0 },
                         mp::ExpAConstraint::Parameters{ 2.0 }));
  ASSERT_EQ(1u, expa.size());
  EXPECT_EQ(2.0, expa[0].GetParameters()[0]);
}

} // namespace
