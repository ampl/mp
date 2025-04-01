#ifndef QP2PASSES_H
#define QP2PASSES_H

#include "mp/flat/qp2passes_base.h"
#include "mp/expr-visitor.h"

namespace mp {

class QP2Passes;

/// Node visitor result.
/// Polynomial degree of the node.
using QP2PassNodeResult = int;

/// Linear/QP expression visitor
class QP2PassVisitor
    : public
      ExprVisitor<
          QP2PassVisitor,
          QP2PassNodeResult> {
public:
  /// Typedef base class
  using Base = ExprVisitor<
      QP2PassVisitor,
      QP2PassNodeResult>;
  /// Construct
  QP2PassVisitor(QP2Passes& qp2p) : qp2p_(qp2p) { }

  /// Init pass 1
  void InitPass1();

  /// Reuse dispatcher
  using Base::Visit;
  /// Overload the dispatching Visit()
  /// to multiply the constant factor for its subtree
  QP2PassNodeResult Visit(Expr e, double f);

  /// Any unsupported expr - we process this top-level term
  /// via buckets
  QP2PassNodeResult VisitUnsupported(Expr e);

  /// Unary minus
  QP2PassNodeResult VisitMinus(UnaryExpr e);
  /// Add
  QP2PassNodeResult VisitAdd(BinaryExpr e);
  /// Sub
  QP2PassNodeResult VisitSub(BinaryExpr e);
  /// Div
  QP2PassNodeResult VisitDiv(BinaryExpr e);
  /// Visit SumExpr
  QP2PassNodeResult VisitSum(internal::ExprTypes::SumExpr expr);

  /// Constant
  QP2PassNodeResult VisitNumericConstant(NumericConstant );

protected:
  const QP2Passes& GetQP2P() const { return qp2p_; }
  QP2Passes& GetQP2P() { return qp2p_; }

  bool CheckDegree(int degree);
  bool CheckAndMax(int degree_from_node, int& deg_max);

  std::pair<bool, double> IsConst(Expr );

private:
  QP2Passes& qp2p_;

  unsigned int timestamp_ {0};
  std::vector<unsigned int> var_timestamps_;

  int mode_ {};         // 2: root (accepts QP), 1: only linear
  double factor_ {1.0};
};

/// A 2-pass parser of QP expressions
/// (implementation)
class QP2Passes : public BasicQP2Passes {
public:
  /// Construct
  QP2Passes(BasicProblemFlattener& flt)
      : BasicQP2Passes(flt) { }

  /// Process an expression, converting into an EExpr
  /// @note Currently only SumExpr
  void Process(Expr ) override;

  /// Obtain the result - moves out
  EExpr GetResult() override;

protected:
  internal::ExprTypes::SumExpr GetTopExpr() const
  { return top_expr_; }

  void InitPass1();
  void RunPass1();
  bool Pass2SeemsWorth() const;
  void InitPass2Full();
  void RunPass2Full();
  void ExtractPass2ResultIntoBuckets();
  void RunPass2Buckets();

private:
  internal::ExprTypes::SumExpr top_expr_;

  QP2PassVisitor visitor_ {*this};

  std::vector<bool> is_term_qp_;

  BucketAccumulator<EExpr> buckets_ {0};
};

}  // namespace mp

#endif // QP2PASSES_H
