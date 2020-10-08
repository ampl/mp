#ifndef EASYMODELER_H
#define EASYMODELER_H

#include <vector>

#include "mp/expr.h"

namespace mp {

/// This class facilitates C++ operator overloading
/// and other convenient stuff for mp expressions
/// and ProblemBuilder functionality
template <class ProblemBuilder>
class EasyModeler {
public:
  using Builder = ProblemBuilder;
  using MP_NumExpr = typename Builder::NumericExpr;

  EasyModeler(ProblemBuilder& bld) : pb_(bld) { }

  /// Wrap mp::NumericExpr to also keep the Builder
  class NumExpr {
  public:
    NumExpr(Builder& bld, MP_NumExpr expr) : bld_(bld), expr_(expr) { }

    /// The underlying mp::NumericExpr object
    MP_NumExpr GetImpl() { return expr_; }
    /// The underlying builder
    Builder& GetBuilder() { return bld_; }

    /// inject operator* etc directly into the out-of-class namespace
#define EXPORT_BINARY_OPERATOR( symbol, NAMEvv, NAMEvc, NAMEcv ) \
    friend NumExpr operator symbol (NumExpr e1, NumExpr e2) { \
      return {e1.GetBuilder(), \
        e1.GetBuilder().MakeBinary(mp::expr::NAMEvv, e1.GetImpl(), e2.GetImpl())}; \
    } \
    friend NumExpr operator symbol (NumExpr e1, double c) { \
      return {e1.GetBuilder(), \
        e1.GetBuilder().MakeBinary(mp::expr::NAMEvc, e1.GetImpl(), \
        e1.GetBuilder().MakeNumericConstant(c))}; \
    } \
    friend NumExpr operator symbol (double c, NumExpr e2) { \
      return {e2.GetBuilder(), \
        e2.GetBuilder().MakeBinary(mp::expr::NAMEcv, \
        e2.GetBuilder().MakeNumericConstant(c), e2.GetImpl())}; \
    }
#define EXPORT_COMMUTATIVE_BINARY_OPERATOR( symbol, NAME ) \
  EXPORT_BINARY_OPERATOR( symbol, NAME, NAME, NAME )

    /// Operator definitions
    EXPORT_COMMUTATIVE_BINARY_OPERATOR( +, ADD )
    EXPORT_COMMUTATIVE_BINARY_OPERATOR( *, MUL )
    EXPORT_BINARY_OPERATOR( ^,   POW,  POW_CONST_EXP,  POW_CONST_BASE )


  private:
    Builder& bld_;
    MP_NumExpr expr_;
  };

  /// Model building: variables
  std::vector<NumExpr> AddVars(int nvars, double lb, double ub, var::Type type = var::CONTINUOUS) {
    std::vector<NumExpr> newVars;
    newVars.reserve(nvars);
    for (int i=0; i<nvars; ++i) {
      auto v = GetBuilder().AddVar(lb, ub, type).index();
      newVars.push_back( NumExpr(GetBuilder(), GetBuilder().MakeVariable(v)) );
    }
    return newVars;
  }

  /// Variable indices
  std::vector<int> GetVarIndices(std::vector<NumExpr>& vars) {
    std::vector<int> result;
    result.reserve(vars.size());
    for (auto& e: vars) {
      auto evar = Cast<Reference>(e.GetImpl());
      assert(expr::VARIABLE == evar.kind());
      result.push_back( evar.index() );
    }
    return result;
  }

  /// Model building: constraints
  typename Builder::MutAlgebraicCon AddAlgCon(double lb, NumExpr expr, double ub) {
    auto con = GetBuilder().AddCon(lb, ub);
    con.set_nonlinear_expr(expr.GetImpl());
    return con;
  }

protected:
  const Builder& GetBuilder() const { return pb_; }
  Builder& GetBuilder() { return pb_; }

private:
  Builder& pb_;
};

template <class ProblemBuilder>
EasyModeler<ProblemBuilder> MakeEasyModeler(ProblemBuilder& bld) {
  return EasyModeler<ProblemBuilder>(bld);
}


} // namespace mp

#endif // EASYMODELER_H
