/*
 Expression writer

 Copyright (C) 2014 AMPL Optimization Inc

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization Inc disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Victor Zverovich
 */

#ifndef MP_EXPR_WRITER_H_
#define MP_EXPR_WRITER_H_

#include "mp/basic-expr-visitor.h"
#include "precedence.h"

namespace mp {

// An expression visitor that writes AMPL expressions in a textual form
// to fmt::Writer. It takes into account precedence and associativity
// of operators avoiding unnecessary parentheses except for potentially
// confusing cases such as "!x = y" which is written as "!(x = y) instead.
template <typename ExprTypes>
class ExprWriter :
    public BasicExprVisitor<ExprWriter<ExprTypes>, void, void, ExprTypes> {
 private:
  fmt::Writer &writer_;
  int precedence_;

  MP_DEFINE_EXPR_TYPES(ExprTypes);

  typedef BasicExprVisitor<ExprWriter<ExprTypes>, void, void, ExprTypes> Base;

  static int precedence(Expr e) { return internal::precedence(e.kind()); }

  // Writes an argument list surrounded by parentheses.
  template <typename Iter>
  void WriteArgs(Iter begin, Iter end, const char *sep = ", ",
      int precedence = prec::UNKNOWN);

  template <typename Expr>
  void WriteArgs(Expr e, const char *sep = ", ",
      int precedence = prec::UNKNOWN) {
    WriteArgs(e.begin(), e.end(), sep, precedence);
  }

  // Writes a function or an expression that has a function syntax.
  template <typename Expr>
  void WriteFunc(Expr e) {
    writer_ << str(e.kind());
    WriteArgs(e);
  }

  template <typename ExprType>
  void WriteBinary(ExprType e);

  void WriteCallArg(Expr arg);

  class Parenthesizer {
   private:
    ExprWriter &writer_;
    int saved_precedence_;
    bool write_paren_;

   public:
    Parenthesizer(ExprWriter &w, Expr e, int precedence);
    ~Parenthesizer();
  };

 public:
  explicit ExprWriter(fmt::Writer &w)
  : writer_(w), precedence_(prec::UNKNOWN) {}

  void Visit(NumericExpr e, int precedence = -1) {
    Parenthesizer p(*this, e, precedence);
    Base::Visit(e);
  }

  void Visit(LogicalExpr e, int precedence = -1) {
    Parenthesizer p(*this, e, precedence);
    Base::Visit(e);
  }

  void VisitNumericConstant(NumericConstant c) { writer_ << c.value(); }

  void VisitUnary(UnaryExpr e) {
    writer_ << str(e.kind()) << '(';
    Visit(e.arg(), prec::UNKNOWN);
    writer_ << ')';
  }

  void VisitMinus(UnaryExpr e) {
    writer_ << '-';
    Visit(e.arg());
  }

  void VisitPow2(UnaryExpr e) {
    Visit(e.arg(), prec::EXPONENTIATION + 1);
    writer_ << " ^ 2";
  }

  void VisitBinary(BinaryExpr e) { WriteBinary(e); }
  void VisitBinaryFunc(BinaryExpr e);
  void VisitIf(IfExpr e);
  void VisitVarArg(VarArgExpr e) { WriteFunc(e); }
  void VisitSum(SumExpr e);
  void VisitCount(CountExpr e) { WriteFunc(e); }
  void VisitNumberOf(NumberOfExpr e);
  void VisitPLTerm(PLTerm e);
  void VisitCall(CallExpr e);
  void VisitVariable(Variable v) { writer_ << 'x' << (v.index() + 1); }

  void VisitNot(NotExpr e) {
     writer_ << '!';
     // Use a precedence higher then relational to print expressions
     // as "!(x = y)" instead of "!x = y".
     LogicalExpr arg = e.arg();
     Visit(arg,
           precedence(arg) == prec::RELATIONAL ? prec::RELATIONAL + 1 : -1);
  }

  void VisitBinaryLogical(BinaryLogicalExpr e) { WriteBinary(e); }
  void VisitRelational(RelationalExpr e) { WriteBinary(e); }
  void VisitLogicalCount(LogicalCountExpr e);
  void VisitIteratedLogical(IteratedLogicalExpr e);
  void VisitImplication(ImplicationExpr e);
  void VisitAllDiff(PairwiseExpr e) { WriteFunc(e); }
  void VisitLogicalConstant(LogicalConstant c) { writer_ << c.value(); }
};

template <typename ExprTypes>
ExprWriter<ExprTypes>::Parenthesizer::Parenthesizer(
    ExprWriter<ExprTypes> &w, Expr e, int prec)
: writer_(w), write_paren_(false) {
  saved_precedence_ = w.precedence_;
  if (prec == -1)
    prec = w.precedence_;
  write_paren_ = precedence(e) < prec;
  if (write_paren_)
    w.writer_ << '(';
  w.precedence_ = precedence(e);
}

template <typename ExprTypes>
ExprWriter<ExprTypes>::Parenthesizer::~Parenthesizer() {
  writer_.precedence_ = saved_precedence_;
  if (write_paren_)
    writer_.writer_ << ')';
}

template <typename ExprTypes>
template <typename Iter>
void ExprWriter<ExprTypes>::WriteArgs(
    Iter begin, Iter end, const char *sep, int precedence) {
  writer_ << '(';
  if (begin != end) {
    Visit(*begin, precedence);
    for (++begin; begin != end; ++begin) {
      writer_ << sep;
      Visit(*begin, precedence);
    }
  }
  writer_ << ')';
}

template <typename ExprTypes>
template <typename ExprType>
void ExprWriter<ExprTypes>::WriteBinary(ExprType e) {
  int prec = precedence(e);
  bool right_associative = prec == prec::EXPONENTIATION;
  Visit(e.lhs(), prec + (right_associative ? 1 : 0));
  writer_ << ' ' << str(e.kind()) << ' ';
  Visit(e.rhs(), prec + (right_associative ? 0 : 1));
}

template <typename ExprTypes>
void ExprWriter<ExprTypes>::WriteCallArg(Expr arg) {
  if (NumericExpr e = ExprTypes::template Cast<NumericExpr>(arg)) {
    Visit(e, prec::UNKNOWN);
    return;
  }
  assert(arg.kind() == expr::STRING);
  writer_ << "'";
  const char *s = ExprTypes::template Cast<StringLiteral>(arg).value();
  for ( ; *s; ++s) {
    char c = *s;
    switch (c) {
    case '\n':
      writer_ << '\\' << c;
      break;
    case '\'':
      // Escape quote by doubling.
      writer_ << c;
      // Fall through.
    default:
      writer_ << c;
    }
  }
  writer_ << "'";
}

template <typename ExprTypes>
void ExprWriter<ExprTypes>::VisitBinaryFunc(BinaryExpr e) {
  writer_ << str(e.kind()) << '(';
  Visit(e.lhs(), prec::UNKNOWN);
  writer_ << ", ";
  Visit(e.rhs(), prec::UNKNOWN);
  writer_ << ')';
}

template <typename ExprTypes>
void ExprWriter<ExprTypes>::VisitIf(IfExpr e) {
  writer_ << "if ";
  Visit(e.condition(), prec::UNKNOWN);
  writer_ << " then ";
  NumericExpr false_expr = e.false_expr();
  bool has_else = !IsZero(false_expr);
  Visit(e.true_expr(), prec::CONDITIONAL + (has_else ? 1 : 0));
  if (has_else) {
    writer_ << " else ";
    Visit(false_expr);
  }
}

template <typename ExprTypes>
void ExprWriter<ExprTypes>::VisitSum(SumExpr e) {
  writer_ << "/* sum */ (";
  typename SumExpr::iterator i = e.begin(), end = e.end();
  if (i != end) {
    Visit(*i);
    for (++i; i != end; ++i) {
      writer_ << " + ";
      Visit(*i);
    }
  }
  writer_ << ')';
}

template <typename ExprTypes>
void ExprWriter<ExprTypes>::VisitNumberOf(NumberOfExpr e) {
  writer_ << "numberof ";
  typename NumberOfExpr::iterator i = e.begin();
  Visit(*i++, prec::UNKNOWN);
  writer_ << " in ";
  WriteArgs(i, e.end());
}

template <typename ExprTypes>
void ExprWriter<ExprTypes>::VisitPLTerm(PLTerm e) {
  writer_ << "<<" << e.breakpoint(0);
  for (int i = 1, n = e.num_breakpoints(); i < n; ++i)
    writer_ << ", " << e.breakpoint(i);
  writer_ << "; " << e.slope(0);
  for (int i = 1, n = e.num_slopes(); i < n; ++i)
    writer_ << ", " << e.slope(i);
  writer_ << ">> ";
  NumericExpr arg = e.arg();
  if (Variable var = ExprTypes::template Cast<Variable>(arg))
    writer_ << "x" << (var.index() + 1);
  else
    writer_ << "e" << ((ExprTypes::template Cast<CommonExpr>(arg)).index() + 1);
}

template <typename ExprTypes>
void ExprWriter<ExprTypes>::VisitCall(CallExpr e) {
  writer_ << e.function().name() << '(';
  int num_args = e.num_args();
  if (num_args > 0) {
    typename CallExpr::iterator i = e.begin();
    WriteCallArg(*i);
    for (typename CallExpr::iterator end = e.end(); i != end; ++i) {
      writer_ << ", ";
      WriteCallArg(*i);
    }
  }
  writer_ << ')';
}

template <typename ExprTypes>
void ExprWriter<ExprTypes>::VisitLogicalCount(LogicalCountExpr e) {
  writer_ << str(e.kind()) << ' ';
  Visit(e.lhs());
  writer_ << ' ';
  WriteArgs(e.rhs());
}

template <typename ExprTypes>
void ExprWriter<ExprTypes>::VisitIteratedLogical(IteratedLogicalExpr e) {
  // There is no way to produce an AMPL forall/exists expression because
  // its indexing is not available any more. So we write a count expression
  // instead with a comment about the original expression.
  writer_ << "/* " << str(e.kind()) << " */ ";
  int prec = prec::LOGICAL_AND + 1;
  const char *op = " && ";
  if (e.kind() == expr::EXISTS) {
    prec = prec::LOGICAL_OR + 1;
    op = " || ";
  }
  WriteArgs(e, op, prec);
}

template <typename ExprTypes>
void ExprWriter<ExprTypes>::VisitImplication(ImplicationExpr e) {
  Visit(e.condition());
  writer_ << " ==> ";
  Visit(e.true_expr(), prec::IMPLICATION + 1);
  LogicalExpr false_expr = e.false_expr();
  LogicalConstant c = ExprTypes::template Cast<LogicalConstant>(false_expr);
  if (!c || c.value() != 0) {
    writer_ << " else ";
    Visit(false_expr);
  }
}
}  // namespace mp

// TODO: test

#endif  // MP_EXPR_WRITER_H_
