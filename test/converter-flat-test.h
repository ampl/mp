#ifndef CONVERTERFLATTEST_H
#define CONVERTERFLATTEST_H

#include <vector>

#include "gtest/gtest.h"

#include "mp/flat/model_api_base.h"
#include "mp/flat/problem_flattener.h"
#include "mp/flat/converter.h"
#include "mp/flat/constr_algebraic.h"

using namespace mp;

template <class Constraint>
class TestBackendAcceptingConstraints :
    public mp::BasicFlatModelAPI {
  using Base = mp::BasicFlatModelAPI;
  /// VARIABLES
  mp::VarArrayDef vars_;

public:
  TestBackendAcceptingConstraints() { }
  TestBackendAcceptingConstraints(mp::Env& ) { }

  static constexpr const char* GetTypeName() { return "tester"; }

  void AddVariables(const mp::VarArrayDef& v) { vars_ = v; }
  int NumVars() const { return (int)vars_.size(); }

public:
  USE_BASE_CONSTRAINT_HANDLERS(Base)

#define STORE_CONSTR(Type, accLevel, grp)  \
private: \
  std::vector<Type> con_ ## Type ## _;  \
public:  \
  const std::vector<Type>& GetCons(const Type& ) const  \
  { return con_ ## Type ## _; }  \
  ACCEPT_CONSTRAINT(Type, accLevel, grp)  \
  void AddConstraint(const Type& con) {  \
    con_ ## Type ## _.push_back(con);  \
  }   \
  bool HasConstraint(const Type& con) {  \
    return con_ ## Type ## _.end()  \
      != std::find(con_ ## Type ## _.begin(),  \
           con_ ## Type ## _.end(), con);  \
  }

  /// ACCEPTING THE CUSTOM CONSTRAINT
  STORE_CONSTR(Constraint, Recommended, mp::CG_Default)

  /// ACCEPT LINEAR CONS
  STORE_CONSTR(LinConEQ, Recommended, mp::CG_Default)
  STORE_CONSTR(LinConLE, Recommended, mp::CG_Default)
  STORE_CONSTR(LinConGE, Recommended, mp::CG_Default)
  STORE_CONSTR(LinConRange, Recommended, mp::CG_Default)

  /// ACCEPT Q CONS
  STORE_CONSTR(QuadConEQ, Recommended, mp::CG_Default)
  STORE_CONSTR(QuadConLE, Recommended, mp::CG_Default)
  STORE_CONSTR(QuadConGE, Recommended, mp::CG_Default)
  // STORE_CONSTR(QuadConRange, Recommended, mp::CG_Default)

  mutable ItemNamer in_ {"x"};

  template <class ConType>
  std::string GetConstraintsPrintout(const ConType& con) const {
    std::ostringstream oss;
    int i=0;
    for (const auto& c: GetCons(con)) {
      oss << c.GetTypeName() << ' ' << (i++)
          << ":  ";
      fmt::MemoryWriter wrt;
      WriteModelItem(wrt, c, in_);
      oss << wrt.str() << std::endl;
    }
    oss << "    ====================\n  <== Searched for:  ";
    fmt::MemoryWriter wrt;
    WriteModelItem(wrt, con, in_);
    oss << wrt.str() << std::endl;
    return oss.str();
  }

public:
};

template <template <class, class, class> class ConverterTemplate, class Constraint>
using InterfaceWithBackendAcceptingConstraints =
        mp::ProblemFltImpl<mp::ProblemFlattener, mp::Problem,
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
    interface_(env_) { interface_.InitOptions(); }
  InterfaceTesterWithBackendAcceptingConstraints(mp::Env& e) :
    interface_(e) { }
  Env& GetEnv() { return env_; }
  Interface& GetInterface() { return interface_; }
  typename Interface::ModelType& GetModel() { return interface_.GetModel(); }
  Backend& GetBackend()
  { return interface_.GetFlatCvt().GetModelAPI(); }
};

#define ASSERT_HAS_CONSTRAINT( backend, constr ) \
  ASSERT_TRUE( (backend).HasConstraint( constr ) ) \
    << (backend).GetConstraintsPrintout(constr)

} // namespace

#endif // CONVERTERFLATTEST_H
