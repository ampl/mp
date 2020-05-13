
#include "gtest/gtest.h"

#include "converter-flat-test.h"
#include "mp/convert/converter_flat.h"

namespace {

TEST(ConverterFlatConstraintsTest, MaximumConstraintIsPassedToBackend) {
  InterfaceWithBackendAcceptingConstraints<
      mp::BasicMPFlatConverter, mp::MaximumConstraint> interface;
  auto& model = interface.GetModel();
  auto con=model.AddCon(5.0, 5.0);
  const auto args = interface.AddVars(3, -1.0, 11.0);
  con.set_nonlinear_expr(MakeIterated(model, mp::expr::MAX, args ));
  interface.ConvertModelAndUpdateBackend();
  ASSERT_TRUE( interface.GetBackend().HasConstraint(
                 mp::MaximumConstraint(3, args)) );
}

}
