#ifndef CONVERTERFLATTEST_H
#define CONVERTERFLATTEST_H

#include <vector>

#include "gtest/gtest.h"

#include "mp/convert/backend.h"
#include "mp/convert/basic_converters.h"
#include "mp/convert/converter_flat.h"

template <class Constraint>
class TestBackendAcceptingConstraints :
    public mp::BasicBackend<TestBackendAcceptingConstraints<Constraint> > {
  using Base = mp::BasicBackend<TestBackendAcceptingConstraints<Constraint> >;
  /// VARIABLES
  struct Var {
    double lb_, ub_;
    mp::var::Type type_;
  };
  std::vector<Var> vars_;

  std::vector<mp::LinearConstraint> lin_constr_;
public:
  TestBackendAcceptingConstraints() { }
  void AddVariable(typename Base::Variable var) {
    vars_.push_back( {var.lb(), var.ub(), var.type()} );
  }
  int NumVars() const { return (int)vars_.size(); }
  void AddConstraint(const mp::LinearConstraint& lc) {
    lin_constr_.push_back( lc );
  }

  /// ACCEPTING THE CUSTOM CONSTRAINT
public:
  USE_BASE_CONSTRAINT_HANDLERS(Base)
private:
  std::vector<Constraint> constr_;
public:
  ACCEPT_CONSTRAINT(Constraint, mp::Recommended)
  void AddConstraint(const Constraint& con) {
    constr_.push_back(con);
  }
  bool HasConstraint(const Constraint& con) {
    return constr_.end() != std::find(constr_.begin(), constr_.end(), con);
  }
  std::string GetConstraintsPrintout() const {
    std::ostringstream oss;
    int i=0;
    for (const auto& con: constr_) {
      oss << con.GetConstraintName() << ' ' << (i++)
          << ":  ";
//      con.print(oss);
      oss << std::endl;
    }
    return oss.str();
  }
};

template <template <class, class, class> class ConverterTemplate, class Constraint>
using InterfaceWithBackendAcceptingConstraints =
        mp::Interface<ConverterTemplate, TestBackendAcceptingConstraints<Constraint> >;


////////////////////////// SERVICE STUFF ////////////////////////////
inline
mp::IteratedExpr MakeIterated(mp::ExprFactory& ef, mp::expr::Kind kind,
                          const std::vector<int> &args) {
  const auto N = args.size();
  auto builder = ef.BeginIterated(kind, N);
  for (size_t i = 0; i < N; ++i)
    builder.AddArg( ef.MakeVariable(args[i]) );
  return ef.EndIterated(builder);
}

////////////////////////// TEST FIXTURE /////////////////////////////
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

#define ASSERT_HAS_CONSTRAINT( backend, constr ) \
  ASSERT_TRUE( (backend).HasConstraint( constr ) ) \
    << (backend).GetConstraintsPrintout()

} // namespace

#endif // CONVERTERFLATTEST_H
