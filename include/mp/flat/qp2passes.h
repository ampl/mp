#ifndef QP2PASSES_H
#define QP2PASSES_H

#include "mp/flat/qp2passes_base.h"
#include "mp/flat/bucketaccum.h"
#include "mp/utils-matrix.h"
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
  /// Init pass 2
  void InitPass2();

  /// Reuse dispatcher
  using Base::Visit;
  /// Overload the dispatching Visit()
  /// to multiply the constant factor for its subtree
  QP2PassNodeResult Visit(Expr e, double f);
  /// Visit at most linear
  /// @param ae: the affine expression to store the subtree
  /// @note resets factor_=1.0 for the subtree
  QP2PassNodeResult VisitAtmostAffine(Expr e, AffineExpr& ae);

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
  /// Mul
  QP2PassNodeResult VisitMul(BinaryExpr e);
  /// PowConstExp
  QP2PassNodeResult VisitPowConstExp(BinaryExpr e);
  /// Pow2
  QP2PassNodeResult VisitPow2(UnaryExpr e);
  /// Pow
  QP2PassNodeResult VisitPow(BinaryExpr e);

  /// Constant
  QP2PassNodeResult VisitNumericConstant(NumericConstant );
  /// Variable
  QP2PassNodeResult VisitVariable(Reference );
  /// Defined variable
  QP2PassNodeResult VisitCommonExpr(Reference );

  /// An estimate on the number of source QP terms
  unsigned long long NumSourceTermsQP() const
  { return n_source_terms_qp_; }

  /// Number of QP vars
  auto NumQPVars() const { return vars_qp_.size(); }

protected:
  const QP2Passes& GetQP2P() const { return qp2p_; }
  QP2Passes& GetQP2P() { return qp2p_; }

  bool CheckDegree(int degree);
  bool CheckAndMax(int degree_from_node, int& deg_max);

  QP2PassNodeResult DoVisitPow2(Expr e);
  /// Add to the top-level or current ae's constant term
  void DoAddConst(long double c);

  /// @todo consider binary terms (cvt:prod)
  bool ProcessAffineFactors(AffineExpr& aeL, AffineExpr& aeR);
  bool EstimateOutmultiplication(AffineExpr& aeL, AffineExpr& aeR);
  void MultiplyOut(AffineExpr& aeL, AffineExpr& aeR);

  std::pair<bool, double> IsConst(Expr );
  /// when the other factor is const
  QP2PassNodeResult VisitFactor(Expr e, double f);

  const AffineExpr* GetPAffineExpr() const { return p_ae_; }
  AffineExpr* GetPAffineExpr() { return p_ae_; }

  /// Note top-level var
  void NoteLinVar(int v);
  void NoteQPVar(int v);
  /// Add top-level lin term
  void AddLinTerm(double c, int v);
  void AddQPTerm(double c, int v1, int v2);

  unsigned int GetTimeStamp() const { return timestamp_; }

private:
  QP2Passes& qp2p_;

  int pass_ {};         // 1 or 2
  int mode_ {};         // 2: root (accepts QP), 1: only linear
  long double factor_ {1.0};

  AffineExpr* p_ae_{};  // pointer to chosen ae being filled

  unsigned int timestamp_ {0};
  std::vector<unsigned int> ts_lin_, ts_qp_;
  SmallVec<int, 64> vars_lin_, vars_qp_;
  unsigned long long n_source_terms_qp_ {};

  /// Pass 2
  long double const_term_ {};   // the top-level constant term
  std::vector<double> coefs_lin_;

  TMatrix<double, 16> coefs_qp_;
  std::vector<int> vperm_qp_;   // inverse of vars_qp_
};


/// A 2-pass parser of QP expressions
/// (implementation).
/// See ASL's mqpcheck().
///
/// Pass 1: mark sum terms which are at most QP.
///   Mark linear and QP variables.
/// Pass 2: sum constant term
///   and linear / QP term coefficients (can out-multiply).
/// Postprocess: join the above with over-QP terms.
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
  /// Obtain the buckets
  const BucketAccumulator<EExpr>& GetBuckets() const
  { return buckets_; }
  /// Obtain the buckets
  BucketAccumulator<EExpr>& GetBuckets()
  { return buckets_; }

  internal::ExprTypes::SumExpr GetTopExpr() const
  { return top_expr_; }

  void InitPass1();
  void RunPass1();
  bool Pass2SeemsWorth() const;
  void InitPass2Full();
  void RunPass2Full();
  void CollectMarkedTerms();
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
