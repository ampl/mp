#ifndef CONVERTERFLATTEST_H
#define CONVERTERFLATTEST_H

#include <vector>

#include "gtest/gtest.h"

#include "mp/flat/backend_model_api_base.h"
#include "mp/flat/model_flattener.h"
#include "mp/flat/converter.h"

template <class Constraint>
class TestBackendAcceptingConstraints :
    public mp::BasicBackendFlatModelAPI {
  using Base = mp::BasicBackendFlatModelAPI;
  /// VARIABLES
  mp::VarArrayDef vars_;

  std::vector<mp::LinConEQ> lin_constr_;
public:
  TestBackendAcceptingConstraints() { }
  TestBackendAcceptingConstraints(mp::Env& ) { }

  static constexpr const char* GetTypeName() { return "tester"; }

  void AddVariables(const mp::VarArrayDef& v) { vars_ = v; }
  int NumVars() const { return (int)vars_.size(); }
  void AddConstraint(const mp::LinConEQ& lc) {
    lin_constr_.push_back( lc );
  }

  /// ACCEPTING THE CUSTOM CONSTRAINT
public:
  USE_BASE_CONSTRAINT_HANDLERS(Base)
private:
  std::vector<Constraint> constr_;
public:
  ACCEPT_CONSTRAINT(Constraint, mp::Recommended, mp::CG_Default)
  ACCEPT_CONSTRAINT(mp::LinConEQ, mp::Recommended, mp::CG_Default)
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
      oss << con.GetTypeName() << ' ' << (i++)
          << ":  ";
//      con.print(oss);
      oss << std::endl;
    }
    return oss.str();
  }

public:
};

template <template <class, class, class> class ConverterTemplate, class Constraint>
using InterfaceWithBackendAcceptingConstraints =
        mp::ModelFltImpl<mp::ModelFlattener, mp::Problem,
          mp::FlatCvtImpl<ConverterTemplate,
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
  mp::Env env_;
public:
  InterfaceTesterWithBackendAcceptingConstraints() :
    interface_(env_) { }
  InterfaceTesterWithBackendAcceptingConstraints(mp::Env& e) :
    interface_(e) { }
  Interface& GetInterface() { return interface_; }
  typename Interface::ModelType& GetModel() { return interface_.GetModel(); }
  Backend& GetBackend()
  { return interface_.GetFlatCvt().GetBasicBackend(); }
};

#define ASSERT_HAS_CONSTRAINT( backend, constr ) \
  ASSERT_TRUE( (backend).HasConstraint( constr ) ) \
    << (backend).GetConstraintsPrintout()

} // namespace

#endif // CONVERTERFLATTEST_H
