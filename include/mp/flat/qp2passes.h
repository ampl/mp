#ifndef QP2PASSES_H
#define QP2PASSES_H

#include <memory>

#include "mp/expr.h"
#include "mp/flat/problem_flattener_base.h"
#include "mp/flat/eexpr.h"
#include "mp/flat/bucketaccum.h"

namespace mp {

class QP2PassVisitor;

/// A 2-pass parser of QP expressions.
/// See ASL's mqpcheck().
///
/// Pass 1: mark sum terms which are at most QP.
///   Mark linear and QP variables.
/// Pass 2: sum constant term
///   and linear / QP term coefficients (can out-multiply).
/// Postprocess: join the above with over-QP terms.
class QP2Passes {
public:
  /// Construct
  QP2Passes(BasicProblemFlattener& flt, QP2PassVisitor& v)
      : flattener_(flt), visitor_(v) { }

  /// Process an expression, converting into an EExpr
  /// @note Currently only SumExpr
  void Process(Expr );

  /// Obtain the result - moves out
  EExpr GetResult();

protected:
  internal::ExprTypes::SumExpr GetTopExpr() const
  { return top_expr_; }

  void InitPass1();
  void RunPass1();
  bool Pass2SeemsWorth() const;
  void InitPass2Full();
  void RunPass2Full();
  void CollectMarkedTerms();
  void ExtractPass2ResultIntoBuckets();
  /// When normal pass 2 seems not worth
  void RunPass2Buckets();

  /// The Flattener
  const BasicProblemFlattener& GetFlattener() const
  { return flattener_; }
  BasicProblemFlattener& GetFlattener()
  { return flattener_; }

private:
  BasicProblemFlattener& flattener_;
  QP2PassVisitor& visitor_;

  internal::ExprTypes::SumExpr top_expr_;
  std::vector<bool> is_term_qp_;
  unsigned int n_qp_terms_ {};     // N terms of degree <=2
  EExpr result_;
};

/// QP2PassVisitor factory.
/// Flattener should have a member variable of this,
/// to avoid reallocation of timestamps.
QP2PassVisitor*
MakeQP2PassVisitor(BasicProblemFlattener& );

/// Delete
void DeleteQP2PassVisitor(QP2PassVisitor* );

/// Shrink QP2PassVisitor's memory,
/// after flattening the model.
void Shrink(QP2PassVisitor& );

// #define DEBUG_QP2PASSES

}  // namespace mp

#endif // QP2PASSES_H
