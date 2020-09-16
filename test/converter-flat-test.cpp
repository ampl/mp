
#include "gtest/gtest.h"

#include "converter-flat-test.h"
#include "mp/convert/converter_flat.h"

namespace {

template <class Constraint>
class InterfaceTesterWithBackendAcceptingConstraints : public ::testing::Test {
  using Interface = InterfaceWithBackendAcceptingConstraints<
      mp::BasicMPFlatConverter, Constraint>;
  Interface interface_;
public:
  Interface& GetInterface() { return interface_; }
  typename Interface::ModelType& GetModel() { return interface_.GetModel(); }
  typename Interface::BackendType& GetBackend() { return interface_.GetBackend(); }
};

using InterfaceTester_MaxConstraint =
             InterfaceTesterWithBackendAcceptingConstraints<mp::MaximumConstraint>;

TEST_F(InterfaceTester_MaxConstraint, MaximumConstraintIsPassedToBackend) {
  auto con=GetModel().AddCon(5.0, 5.0);
  const auto args = GetInterface().AddVars(3, -1.0, 11.0);
  con.set_nonlinear_expr(MakeIterated(GetModel(), mp::expr::MAX, args ));
  GetInterface().ConvertModelAndUpdateBackend();
  ASSERT_TRUE( GetBackend().HasConstraint( mp::MaximumConstraint(3, args)) );
}


/// TODO replace actual constraint etc.
using InterfaceTester_QuadraticConstraint =
             InterfaceTesterWithBackendAcceptingConstraints<mp::MaximumConstraint>;

TEST_F(InterfaceTester_QuadraticConstraint, QuadConstraintIsPassedToBackend) {
  auto con=GetModel().AddCon(5.0, 5.0);
  const auto args = GetInterface().AddVars(2, -1.0, 11.0);
  con.set_nonlinear_expr(GetModel().MakeBinary(mp::expr::MUL,
                                               GetModel().MakeVariable( args[0] ),
                         GetModel().MakeVariable( args[1]) ) );
  GetInterface().ConvertModelAndUpdateBackend();
  ASSERT_TRUE( GetBackend().HasConstraint( mp::MaximumConstraint(3, args)) );
}

}
