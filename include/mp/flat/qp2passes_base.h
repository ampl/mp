#ifndef QP2PASSES_BASE_H
#define QP2PASSES_BASE_H

#include <memory>

#include "mp/expr.h"
#include "mp/flat/problem_flattener_base.h"
#include "mp/flat/eexpr.h"

namespace mp {

/// A 2-pass parser of QP expressions
/// (interface)
class BasicQP2Passes {
public:
  /// Construct
  BasicQP2Passes(BasicProblemFlattener& flt)
      : flattener_(flt) { }
  /// Destruct
  virtual ~BasicQP2Passes() { }

  /// Process an expression, converting into an EExpr
  /// @note Currently only SumExpr
  virtual void Process(Expr ) = 0;

  /// Obtain the result - moves out
  virtual EExpr GetResult() = 0;

  /// The Flattener
  const BasicProblemFlattener& GetFlattener() const
  { return flattener_; }
  BasicProblemFlattener& GetFlattener()
  { return flattener_; }

private:
  BasicProblemFlattener& flattener_;
};

/// QP2Passes factory
std::unique_ptr<BasicQP2Passes>
MakeQP2Passes(BasicProblemFlattener& );

}  // namespace mp

#endif // QP2PASSES_BASE_H
