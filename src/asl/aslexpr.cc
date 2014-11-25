/*
 A C++ interface to AMPL expressions.

 Copyright (C) 2012 AMPL Optimization Inc

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

#include "aslexpr.h"
#include "precedence.h"

#include "asl/aslexpr-visitor.h"

#include <cstdio>
#include <cstring>

using std::size_t;

using mp::asl::Cast;
using mp::asl::Expr;
using mp::asl::NumericConstant;
using mp::asl::NumericExpr;
using mp::asl::LogicalExpr;
namespace prec = mp::prec;
namespace asl = mp::asl;

namespace {

// An expression visitor that writes AMPL expressions in a textual form
// to fmt::Writer. It takes into account precedence and associativity
// of operators avoiding unnecessary parentheses except for potentially
// confusing cases such as "!x = y" which is written as "!(x = y) instead.
class ExprWriter :
    public mp::ExprVisitor<ExprWriter, void, void, asl::ExprTypes> {
 private:
  fmt::Writer &writer_;
  int precedence_;

  typedef mp::ExprVisitor<ExprWriter, void, void, asl::ExprTypes> ExprVisitor;

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
    writer_ << e.opstr();
    WriteArgs(e);
  }

  template <typename Expr>
  void WriteBinary(Expr e);

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
    ExprVisitor::Visit(e);
  }

  void Visit(LogicalExpr e, int precedence = -1) {
    Parenthesizer p(*this, e, precedence);
    ExprVisitor::Visit(e);
  }

  void VisitUnary(asl::UnaryExpr e) {
    writer_ << e.opstr() << '(';
    Visit(e.arg(), prec::UNKNOWN);
    writer_ << ')';
  }

  void VisitUnaryMinus(asl::UnaryExpr e) {
    writer_ << '-';
    Visit(e.arg());
  }

  void VisitPow2(asl::UnaryExpr e) {
    Visit(e.arg(), prec::EXPONENTIATION + 1);
    writer_ << " ^ 2";
  }

  void VisitBinary(asl::BinaryExpr e) { WriteBinary(e); }
  void VisitBinaryFunc(asl::BinaryExpr e);
  void VisitVarArg(asl::VarArgExpr e) { WriteFunc(e); }
  void VisitIf(asl::IfExpr e);
  void VisitSum(asl::SumExpr e);
  void VisitCount(asl::CountExpr e) { WriteFunc(e); }
  void VisitNumberOf(asl::NumberOfExpr e);
  void VisitPiecewiseLinear(asl::PiecewiseLinearExpr e);
  void VisitCall(asl::CallExpr e);
  void VisitNumericConstant(NumericConstant c) { writer_ << c.value(); }
  void VisitVariable(asl::Variable v) { writer_ << 'x' << (v.index() + 1); }

  void VisitNot(asl::NotExpr e) {
     writer_ << '!';
     // Use a precedence higher then relational to print expressions
     // as "!(x = y)" instead of "!x = y".
     LogicalExpr arg = e.arg();
     Visit(arg,
         arg.precedence() == prec::RELATIONAL ? prec::RELATIONAL + 1 : -1);
  }

  void VisitBinaryLogical(asl::BinaryLogicalExpr e) { WriteBinary(e); }
  void VisitRelational(asl::RelationalExpr e) { WriteBinary(e); }
  void VisitLogicalCount(asl::LogicalCountExpr e);
  void VisitIteratedLogical(asl::IteratedLogicalExpr e);
  void VisitImplication(asl::ImplicationExpr e);
  void VisitAllDiff(asl::PairwiseExpr e) { WriteFunc(e); }
  void VisitLogicalConstant(asl::LogicalConstant c) { writer_ << c.value(); }
};

ExprWriter::Parenthesizer::Parenthesizer(ExprWriter &w, Expr e, int precedence)
: writer_(w), write_paren_(false) {
  saved_precedence_ = w.precedence_;
  if (precedence == -1)
    precedence = w.precedence_;
  write_paren_ = e.precedence() < precedence;
  if (write_paren_)
    w.writer_ << '(';
  w.precedence_ = e.precedence();
}

ExprWriter::Parenthesizer::~Parenthesizer() {
  writer_.precedence_ = saved_precedence_;
  if (write_paren_)
    writer_.writer_ << ')';
}

template <typename Iter>
void ExprWriter::WriteArgs(
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

template <typename Expr>
void ExprWriter::WriteBinary(Expr e) {
  int precedence = e.precedence();
  bool right_associative = precedence == prec::EXPONENTIATION;
  Visit(e.lhs(), precedence + (right_associative ? 1 : 0));
  writer_ << ' ' << e.opstr() << ' ';
  Visit(e.rhs(), precedence + (right_associative ? 0 : 1));
}

void ExprWriter::WriteCallArg(Expr arg) {
  if (NumericExpr e = Cast<NumericExpr>(arg)) {
    Visit(e, prec::UNKNOWN);
    return;
  }
  assert(arg.kind() == mp::expr::STRING);
  writer_ << "'";
  const char *s = Cast<asl::StringLiteral>(arg).value();
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

void ExprWriter::VisitBinaryFunc(asl::BinaryExpr e) {
  writer_ << e.opstr() << '(';
  Visit(e.lhs(), prec::UNKNOWN);
  writer_ << ", ";
  Visit(e.rhs(), prec::UNKNOWN);
  writer_ << ')';
}

void ExprWriter::VisitIf(asl::IfExpr e) {
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

void ExprWriter::VisitSum(asl::SumExpr e) {
  writer_ << "/* sum */ (";
  asl::SumExpr::iterator i = e.begin(), end = e.end();
  if (i != end) {
    Visit(*i);
    for (++i; i != end; ++i) {
      writer_ << " + ";
      Visit(*i);
    }
  }
  writer_ << ')';
}

void ExprWriter::VisitNumberOf(asl::NumberOfExpr e) {
  writer_ << "numberof ";
  asl::NumberOfExpr::iterator i = e.begin();
  Visit(*i++, prec::UNKNOWN);
  writer_ << " in ";
  WriteArgs(i, e.end());
}

void ExprWriter::VisitPiecewiseLinear(asl::PiecewiseLinearExpr e) {
  writer_ << "<<" << e.breakpoint(0);
  for (int i = 1, n = e.num_breakpoints(); i < n; ++i)
    writer_ << ", " << e.breakpoint(i);
  writer_ << "; " << e.slope(0);
  for (int i = 1, n = e.num_slopes(); i < n; ++i)
    writer_ << ", " << e.slope(i);
  writer_ << ">> " << "x" << (e.var_index() + 1);
}

void ExprWriter::VisitCall(asl::CallExpr e) {
  writer_ << e.function().name() << '(';
  int num_args = e.num_args();
  if (num_args > 0) {
    WriteCallArg(e[0]);
    for (int i = 1; i < num_args; ++i) {
      writer_ << ", ";
      WriteCallArg(e[i]);
    }
  }
  writer_ << ')';
}

void ExprWriter::VisitLogicalCount(asl::LogicalCountExpr e) {
  writer_ << e.opstr() << ' ';
  Visit(e.lhs());
  writer_ << ' ';
  WriteArgs(e.rhs());
}

void ExprWriter::VisitIteratedLogical(asl::IteratedLogicalExpr e) {
  // There is no way to produce an AMPL forall/exists expression because
  // its indexing is not available any more. So we write a count expression
  // instead with a comment about the original expression.
  writer_ << "/* " << e.opstr() << " */ ";
  int precedence = prec::LOGICAL_AND + 1;
  const char *op = " && ";
  if (e.kind() == mp::expr::EXISTS) {
    precedence = prec::LOGICAL_OR + 1;
    op = " || ";
  }
  WriteArgs(e, op, precedence);
}

void ExprWriter::VisitImplication(asl::ImplicationExpr e) {
  Visit(e.condition());
  writer_ << " ==> ";
  Visit(e.true_expr(), prec::IMPLICATION + 1);
  LogicalExpr false_expr = e.false_expr();
  asl::LogicalConstant c = Cast<asl::LogicalConstant>(false_expr);
  if (!c || c.value() != 0) {
    writer_ << " else ";
    Visit(false_expr);
  }
}

// Compares expressions for equality.
class ExprEqual :
    public mp::ExprVisitor<ExprEqual, bool, bool, asl::ExprTypes> {
 private:
  Expr expr_;

 public:
  explicit ExprEqual(Expr e) : expr_(e) {}

  template <typename T>
  bool VisitNumericConstant(T c) { return Cast<T>(expr_).value() == c.value(); }

  bool VisitVariable(asl::Variable v) {
    return Cast<asl::Variable>(expr_).index() == v.index();
  }

  template <typename E>
  bool VisitUnary(E e) {
    return Equal(Cast<E>(expr_).arg(), e.arg());
  }

  template <typename E>
  bool VisitBinary(E e) {
    E binary = Cast<E>(expr_);
    return Equal(binary.lhs(), e.lhs()) && Equal(binary.rhs(), e.rhs());
  }

  template <typename E>
  bool VisitIf(E e) {
    E if_expr = Cast<E>(expr_);
    return Equal(if_expr.condition(), e.condition()) &&
           Equal(if_expr.true_expr(), e.true_expr()) &&
           Equal(if_expr.false_expr(), e.false_expr());
  }

  bool VisitPiecewiseLinear(asl::PiecewiseLinearExpr e) {
    asl::PiecewiseLinearExpr pl = Cast<asl::PiecewiseLinearExpr>(expr_);
    int num_breakpoints = pl.num_breakpoints();
    if (num_breakpoints != e.num_breakpoints())
      return false;
    for (int i = 0; i < num_breakpoints; ++i) {
      if (pl.slope(i) != e.slope(i) || pl.breakpoint(i) != e.breakpoint(i))
        return false;
    }
    return pl.slope(num_breakpoints) == e.slope(num_breakpoints) &&
           pl.var_index() == e.var_index();
  }

  bool VisitCall(asl::CallExpr e) {
    asl::CallExpr call = Cast<asl::CallExpr>(expr_);
    int num_args = call.num_args();
    if (call.function() != e.function() || num_args != e.num_args())
      return false;
    for (int i = 0; i < num_args; ++i) {
      Expr arg = call[i], other_arg = e[i];
      if (arg.kind() != other_arg.kind())
        return false;
      if (NumericExpr num_arg = Cast<NumericExpr>(arg)) {
        if (!Equal(num_arg, Cast<NumericExpr>(other_arg)))
          return false;
      } else if (std::strcmp(
              Cast<asl::StringLiteral>(arg).value(),
              Cast<asl::StringLiteral>(other_arg).value()) != 0)
        return false;
    }
    return true;
  }

  template <typename E>
  bool VisitVarArg(E e) {
    E vararg = Cast<E>(expr_);
    typename E::iterator i = vararg.begin(), iend = vararg.end();
    typename E::iterator j = e.begin(), jend = e.end();
    for (; i != iend; ++i, ++j) {
      if (j == jend || !Equal(*i, *j))
        return false;
    }
    return j == jend;
  }

  bool VisitSum(asl::SumExpr e) { return VisitVarArg(e); }
  bool VisitCount(asl::CountExpr e) { return VisitVarArg(e); }
  bool VisitNumberOf(asl::NumberOfExpr e) { return VisitVarArg(e); }

  bool VisitLogicalConstant(asl::LogicalConstant c) {
    return VisitNumericConstant(c);
  }

  bool VisitNot(asl::NotExpr e) { return VisitUnary(e); }

  bool VisitBinaryLogical(asl::BinaryLogicalExpr e) { return VisitBinary(e); }
  bool VisitRelational(asl::RelationalExpr e) { return VisitBinary(e); }
  bool VisitLogicalCount(asl::LogicalCountExpr e) { return VisitBinary(e); }

  bool VisitImplication(asl::ImplicationExpr e) { return VisitIf(e); }

  bool VisitIteratedLogical(asl::IteratedLogicalExpr e) {
    return VisitVarArg(e);
  }

  bool VisitAllDiff(asl::PairwiseExpr e) { return VisitVarArg(e); }
};
}  // namespace

#ifdef MP_USE_UNORDERED_MAP

using asl::internal::HashCombine;

namespace {
// Computes a hash value for an expression.
class ExprHasher :
    public mp::ExprVisitor<ExprHasher, size_t, size_t, asl::ExprTypes> {
 private:
  static size_t Hash(Expr e) {
    return HashCombine<int>(0, e.kind());
  }

  template <typename T>
  static size_t Hash(Expr e, const T &value) {
    return HashCombine(Hash(e), value);
  }

 public:
  size_t VisitNumericConstant(NumericConstant c) { return Hash(c, c.value()); }
  size_t VisitVariable(asl::Variable v) { return Hash(v, v.index()); }

  size_t VisitUnary(asl::UnaryExpr e) { return Hash(e, e.arg()); }

  template <typename E>
  size_t VisitBinary(E e) { return HashCombine(Hash(e, e.lhs()), e.rhs()); }

  template <typename E>
  size_t VisitIf(E e) {
    size_t hash = HashCombine(Hash(e), e.condition());
    return HashCombine(HashCombine(hash, e.true_expr()), e.false_expr());
  }

  size_t VisitPiecewiseLinear(asl::PiecewiseLinearExpr e) {
    size_t hash = Hash(e);
    int num_breakpoints = e.num_breakpoints();
    for (int i = 0; i < num_breakpoints; ++i) {
      hash = HashCombine(hash, e.slope(i));
      hash = HashCombine(hash, e.breakpoint(i));
    }
    hash = HashCombine(hash, e.slope(num_breakpoints));
    return HashCombine(hash, e.var_index());
  }

  size_t VisitCall(asl::CallExpr e) {
    // Function name is hashed as a pointer. This works because the function
    // object is the same for all calls to the same function.
    size_t hash = Hash(e, e.function().name());
    for (int i = 0, n = e.num_args(); i < n; ++i)
      hash = HashCombine(hash, e[i]);
    return hash;
  }

  template <typename E>
  size_t VisitVarArg(E e) {
    size_t hash = Hash(e);
    for (typename E::iterator i = e.begin(), end = e.end(); i != end; ++i)
      hash = HashCombine(hash, *i);
    return hash;
  }

  size_t VisitSum(asl::SumExpr e) { return VisitVarArg(e); }
  size_t VisitCount(asl::CountExpr e) { return VisitVarArg(e); }
  size_t VisitNumberOf(asl::NumberOfExpr e) { return VisitVarArg(e); }

  size_t VisitLogicalConstant(asl::LogicalConstant c) {
    return Hash(c, c.value());
  }

  size_t VisitNot(asl::NotExpr e) { return Hash(e, e.arg()); }

  size_t VisitBinaryLogical(asl::BinaryLogicalExpr e) { return VisitBinary(e); }
  size_t VisitRelational(asl::RelationalExpr e) { return VisitBinary(e); }

  size_t VisitLogicalCount(asl::LogicalCountExpr e) {
    NumericExpr rhs = e.rhs();
    return HashCombine(Hash(e, e.lhs()), rhs);
  }

  size_t VisitImplication(asl::ImplicationExpr e) { return VisitIf(e); }

  size_t VisitIteratedLogical(asl::IteratedLogicalExpr e) {
    return VisitVarArg(e);
  }

  size_t VisitAllDiff(asl::PairwiseExpr e) { return VisitVarArg(e); }

  size_t VisitStringLiteral(asl::StringLiteral s) {
    size_t hash = Hash(s);
    for (const char *value = s.value(); *value; ++value)
      hash = HashCombine(hash, *value);
    return hash;
  }
};
}

namespace std {
template <>
struct hash<Expr> {
  std::size_t operator()(Expr e) const;
};
}

size_t std::hash<Expr>::operator()(Expr expr) const {
  ExprHasher hasher;
  NumericExpr n = Cast<NumericExpr>(expr);
  return n ? hasher.Visit(n) :
             hasher.VisitStringLiteral(Cast<asl::StringLiteral>(expr));
}

size_t std::hash<NumericExpr>::operator()(NumericExpr expr) const {
  return ExprHasher().Visit(expr);
}

size_t std::hash<LogicalExpr>::operator()(LogicalExpr expr) const {
  return ExprHasher().Visit(expr);
}
#else
# if defined(AMPL_NO_UNORDERED_MAP_WARNING)
  // Do nothing.
# elif defined(_MSC_VER)
#  pragma message("warning: unordered_map not available, numberof may be slow")
# else
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wpedantic"
#  warning "unordered_map not available, numberof may be slow"
#  pragma clang diagnostic pop
# endif
#endif

namespace mp {
namespace asl {

const de VarArgExpr::END = de();

bool Equal(NumericExpr e1, NumericExpr e2) {
  if (e1.kind() != e2.kind())
    return false;
  return ExprEqual(e1).Visit(e2);
}

bool Equal(LogicalExpr e1, LogicalExpr e2) {
  if (e1.kind() != e2.kind())
    return false;
  return ExprEqual(e1).Visit(e2);
}

#ifdef MP_USE_UNORDERED_MAP
size_t internal::HashNumberOfArgs::operator()(NumberOfExpr e) const {
  size_t hash = 0;
  for (int i = 1, n = e.num_args(); i < n; ++i)
    hash = HashCombine(hash, e[i]);
  return hash;
}
#endif

bool internal::EqualNumberOfArgs::operator()(
    NumberOfExpr lhs, NumberOfExpr rhs) const {
  int num_args = lhs.num_args();
  if (num_args != rhs.num_args())
    return false;
  for (int i = 1; i < num_args; ++i) {
    if (!Equal(lhs[i], rhs[i]))
      return false;
  }
  return true;
}

template <typename LinearExpr>
void WriteExpr(fmt::Writer &w, LinearExpr linear, NumericExpr nonlinear) {
  bool have_terms = false;
  typedef typename LinearExpr::iterator Iterator;
  for (Iterator i = linear.begin(), e = linear.end(); i != e; ++i) {
    double coef = i->coef();
    if (coef != 0) {
      if (have_terms)
        w << " + ";
      else
        have_terms = true;
      if (coef != 1)
        w << coef << " * ";
      w << "x" << (i->var_index() + 1);
    }
  }
  if (!nonlinear || IsZero(nonlinear)) {
    if (!have_terms)
      w << "0";
    return;
  }
  if (have_terms)
    w << " + ";
  ExprWriter(w).Visit(nonlinear);
}

template
void WriteExpr<LinearObjExpr>(
    fmt::Writer &w, LinearObjExpr linear, NumericExpr nonlinear);

template
void WriteExpr<LinearConExpr>(
    fmt::Writer &w, LinearConExpr linear, NumericExpr nonlinear);
}
}
