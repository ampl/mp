
#include "mp/easy-modeler.h"

#include "converter-flat-test.h"

namespace {

using InterfaceTester_MaxConstraint =
             InterfaceTesterWithBackendAcceptingConstraints<mp::MaximumConstraint>;

TEST_F(InterfaceTester_MaxConstraint, MaximumConstraintIsPassedToBackend) {
  auto con=GetModel().AddCon(5.0, 5.0);
  const auto args = GetInterface().AddVars(3, -1.0, 11.0);
  con.set_nonlinear_expr(MakeIterated(GetModel(), mp::expr::MAX, args ));
  GetInterface().ConvertModelAndUpdateBackend();
  ASSERT_HAS_CONSTRAINT( GetBackend(), mp::MaximumConstraint(3, args) );
}


using InterfaceTester_QuadraticConstraint =
             InterfaceTesterWithBackendAcceptingConstraints<mp::QuadraticConstraint>;

TEST_F(InterfaceTester_QuadraticConstraint, QuadConstraintIsPassedToBackend) {
  auto con=GetModel().AddCon(5.0, 5.0);
  const auto args = GetInterface().AddVars(2, -1.0, 11.0);
  con.set_nonlinear_expr(GetModel().MakeBinary(mp::expr::MUL,
                                               GetModel().MakeVariable( args[0] ),
                                               GetModel().MakeVariable( args[1]) ) );
  GetInterface().ConvertModelAndUpdateBackend();
  ASSERT_HAS_CONSTRAINT( GetBackend(), mp::QuadraticConstraint(
    {},
    { {1.0, args[0], args[1]} },
    5.0, 5.0 ) );
}

/// This is to illustrate how fun it was to model with ProblemBuilder syntax.
/// You may skip to the next test
TEST_F(InterfaceTester_QuadraticConstraint, QuadExprIsMultipliedOutAndInlinedAndPassedToBackend) {
  auto con=GetModel().AddCon(5.0, GetInterface().Infty());
  const auto args = GetInterface().AddVars(3, -1.0, 11.0);
  auto expr_5x = GetModel().MakeBinary(mp::expr::MUL,
                                       GetModel().MakeNumericConstant( 5.0 ),
                                       GetModel().MakeVariable( args[0] ) );
  auto expr_5xPlus3 = GetModel().MakeBinary(mp::expr::ADD,
                                       expr_5x,
                                       GetModel().MakeNumericConstant( 3.0 ) );
  auto expr_6y = GetModel().MakeBinary(mp::expr::MUL,
                                       GetModel().MakeNumericConstant( 6.0 ),
                                       GetModel().MakeVariable( args[1] ) );
  auto expr_6xPlus2 = GetModel().MakeBinary(mp::expr::ADD,
                                       expr_6y,
                                       GetModel().MakeNumericConstant( 2.0 ) );
  auto expr_5x3Times6y2 = GetModel().MakeBinary(mp::expr::MUL, expr_5xPlus3, expr_6xPlus2 );
  auto expr_5x3Times6y2PlusZ = GetModel().MakeBinary(mp::expr::ADD, expr_5x3Times6y2,
                                                     GetModel().MakeVariable( args[2]) );
  con.set_nonlinear_expr( expr_5x3Times6y2PlusZ );
  GetInterface().ConvertModelAndUpdateBackend();
  ASSERT_HAS_CONSTRAINT( GetBackend(), mp::QuadraticConstraint(
  { {10.0, args[0]}, {18.0, args[1]}, {1.0, args[2]} },
    { {30.0, args[0], args[1]} },
    -1.0, GetInterface().Infty() ) );
}

/// EasyModeler syntax
TEST_F(InterfaceTester_QuadraticConstraint, QuadExprIsMultipliedOutAndInlinedAndPassedToBackend__EasyModel) {
  auto modeler = mp::MakeEasyModeler(GetModel());
  auto x = modeler.AddVars(3, -1.0, 11.0);
  modeler.AddAlgCon(5.0,
                    (5*x[0]+3) * (6*x[1]+2) + x[2],
                    GetInterface().Infty());
  GetInterface().ConvertModelAndUpdateBackend();
  const auto xi = modeler.GetVarIndices(x);
  ASSERT_HAS_CONSTRAINT( GetBackend(), mp::QuadraticConstraint(
  { {10.0, xi[0]}, {18.0, xi[1]}, {1.0, xi[2]} },
    { {30.0, xi[0], xi[1]} },
    -1.0, GetInterface().Infty() ) );
}


} // namespace
