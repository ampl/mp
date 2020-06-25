#ifndef CONVERTERFLATTEST_H
#define CONVERTERFLATTEST_H

#include <vector>

#include "mp/backend.h"
#include "mp/convert/basic_converters.h"

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
  TestBackendAcceptingConstraints() : Base("TestBackend") { }
  void AddVariable(typename Base::Variable var) {
    vars_.push_back( {var.lb(), var.ub(), var.type()} );
  }
  int NumVars() const { return (int)vars_.size(); }
  /// The linear constraints
  void AddLinearConstraint(int nnz, const double* c, const int* v,
                           double lb, double ub) {
    lin_constr_.push_back( { {c, c+nnz}, {v, v+nnz}, lb, ub } );
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
  for (int i = 0; i < N; ++i)
    builder.AddArg( ef.MakeVariable(args[i]) );
  return ef.EndIterated(builder);
}


#endif // CONVERTERFLATTEST_H
