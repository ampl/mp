#ifndef CONVERTERFLATTEST_H
#define CONVERTERFLATTEST_H

#include <vector>

#include "gtest/gtest.h"

#include "mp/flat/backend.h"
#include "mp/flat/expr_flattener.h"
#include "mp/flat/converter.h"

template <class Constraint>
class TestBackendAcceptingConstraints :
    public mp::Backend< TestBackendAcceptingConstraints<Constraint> > {
  using Base = mp::Backend< TestBackendAcceptingConstraints<Constraint> >;
  /// VARIABLES
  mp::VarArrayDef vars_;

  std::vector<mp::RangeLinCon> lin_constr_;
public:
  TestBackendAcceptingConstraints() { }
  void AddVariables(const mp::VarArrayDef& v) { vars_ = v; }
  int NumVars() const { return (int)vars_.size(); }
  void AddConstraint(const mp::RangeLinCon& lc) {
    lin_constr_.push_back( lc );
  }

  /// ACCEPTING THE CUSTOM CONSTRAINT
public:
  USE_BASE_CONSTRAINT_HANDLERS(Base)
private:
  std::vector<Constraint> constr_;
public:
  ACCEPT_CONSTRAINT(Constraint, mp::Recommended, mp::CG_Default)
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
        mp::ExprFlattenerImpl<mp::ExprFlattener, mp::Problem,
          mp::ConverterImpl<ConverterTemplate,
            TestBackendAcceptingConstraints<Constraint> > >;


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
      mp::FlatConverter, Constraint>;
  using Backend = TestBackendAcceptingConstraints<Constraint>;
  Interface interface_;
public:
  Interface& GetInterface() { return interface_; }
  typename Interface::ModelType& GetModel() { return interface_.GetModel(); }
  Backend& GetBackend()
  { return dynamic_cast<Backend&>( interface_.GetBasicBackend() ); }
};

#define ASSERT_HAS_CONSTRAINT( backend, constr ) \
  ASSERT_TRUE( (backend).HasConstraint( constr ) ) \
    << (backend).GetConstraintsPrintout()

} // namespace

#endif // CONVERTERFLATTEST_H
