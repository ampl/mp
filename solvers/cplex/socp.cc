/*
 CPLEX solver with SOCP transformations

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

#include "mp/expr-visitor.h"
#include "mp/nl.h"
#include "mp/problem.h"
#include "asl.h"

#define SKIP_NL2_DEFINES
#include "asl_pfg.h"
#include "asl_pfgh.h"
#undef cde

#include "asl/aslbuilder.h"
#include "asl/aslproblem.h"

namespace asl = mp::asl;
using asl::internal::ASLBuilder;

using mp::Problem;

namespace mp {

// Returns true if e is a positive constant.
inline bool IsPosConstant(NumericExpr e) {
  NumericConstant n = Cast<NumericConstant>(e);
  return n && n.value() > 0;
}

// A detector of special forms of expressions.
template <typename Impl>
class ExprDetector : public ExprVisitor<Impl, bool> {
 public:
  bool VisitUnhandledNumericExpr(NumericExpr) { return false; }
  bool VisitUnhandledLogicalExpr(LogicalExpr) { return false; }
};

// A detector of sums of expressions possibly multipled by positive constants.
template <typename Impl>
class WeightedSumDetector : public ExprDetector<Impl> {
 public:
  bool VisitAdd(BinaryExpr e) {
    return this->Visit(e.lhs()) && this->Visit(e.rhs());
  }

  bool VisitMul(BinaryExpr e) {
    return ((IsPosConstant(e.lhs()) && this->Visit(e.rhs())) ||
            (this->Visit(e.lhs()) && IsPosConstant(e.rhs())));
  }

  bool VisitSum(IteratedExpr e) {
    for (IteratedExpr::iterator i = e.begin(), end = e.end(); i != end; ++i) {
      if (!this->Visit(*i))
        return false;
    }
    return true;
  }
};

// A detector of affine expressions.
// TODO: allow multiplication by any constant, not only positive
class AffineExprDetector : public WeightedSumDetector<AffineExprDetector> {
 public:
  bool VisitNumericConstant(NumericConstant) { return true; }
  bool VisitVariable(Variable) { return true; }
};

// A detector of sums of squares of affine expressions.
class SumOfSquaresDetector : public WeightedSumDetector<SumOfSquaresDetector> {
 private:
  // The number of quadratic terms.
  int num_terms_;

 public:
  SumOfSquaresDetector() : num_terms_(0) {}

  int num_terms() const { return num_terms_; }

  bool VisitPow2(UnaryExpr e) {
    ++num_terms_;
    return AffineExprDetector().Visit(e.arg());
  }
};

// A detector of sums of norms.
class SumOfNormsDetector : public WeightedSumDetector<SumOfNormsDetector> {
 private:
  // The number of quadratic terms in each sqrt expression.
  std::vector<int> num_terms_;

 public:
  const int *num_terms() const { return &num_terms_[0]; }

  bool VisitSqrt(UnaryExpr e) {
    SumOfSquaresDetector detector;
    if (!detector.Visit(e.arg()))
      return false;
    num_terms_.push_back(detector.num_terms());
    return true;
  }
};

template <typename Impl>
class SumConverter : public ExprVisitor<Impl, void> {
 private:
  Problem &problem_;
  double coef_;

 public:
  explicit SumConverter(Problem &p) : problem_(p), coef_(1) {}

  Problem &problem() { return problem_; }
  double coef() const { return coef_; }

  void VisitAdd(BinaryExpr e) {
    this->Visit(e.lhs());
    this->Visit(e.rhs());
  }

  // TODO: handle SUM and MUL
};

// Adds a variable ``x`` and a constraint ``x = <affine-expr>``.
class AffineExprConverter : public SumConverter<AffineExprConverter> {
 private:
  Problem::LinearConBuilder builder_;
  double constant_;
  int var_index_;
  int con_index_;

 public:
  explicit AffineExprConverter(Problem &p)
    : SumConverter<AffineExprConverter>(p),
      builder_(p.AddCon(0, 0)), constant_(0),
      var_index_(0), con_index_(p.num_algebraic_cons() - 1) {
  }

  int var_index() const { return var_index_; }

  void Convert(NumericExpr e);

  void VisitNumericConstant(NumericConstant n) {
    constant_ += coef() * n.value();
  }

  void VisitVariable(Variable v) {
    builder_.AddTerm(v.index(), coef());
  }

  // TODO
};

void AffineExprConverter::Convert(NumericExpr e) {
  // Add a free variable ``x`` to represent the affine expression.
  double inf = std::numeric_limits<double>::infinity();
  var_index_ = problem().AddVar(-inf, inf).index();
  // Build the constraint ``x = expr``.
  builder_.AddTerm(var_index_, -1);
  Visit(e);
  if (constant_ != 0) {
    // Adjust the right-hand side.
    Problem::MutAlgebraicCon con = problem().algebraic_con(con_index_);
    con.set_lb(-constant_);
    con.set_ub(-constant_);
  }
}

class SumOfSquaresConverter : public SumConverter<SumOfSquaresConverter> {
 private:
  Problem::IteratedExprBuilder &sum_;

 public:
  explicit SumOfSquaresConverter(Problem &p, Problem::IteratedExprBuilder &sum)
  : SumConverter<SumOfSquaresConverter>(p), sum_(sum) {}

  void VisitPow2(UnaryExpr e);
};

void SumOfSquaresConverter::VisitPow2(UnaryExpr e) {
  Problem &p = problem();
  AffineExprConverter converter(p);
  converter.Convert(e.arg());
  // Replace the term ``coef * expr ^ 2`` with ``coef * x ^ 2``.
  NumericExpr term = p.MakeUnary(expr::POW2,
                                 p.MakeVariable(converter.var_index()));
  if (coef() != 1)
    term = p.MakeBinary(expr::MUL, p.MakeNumericConstant(coef()), term);
  sum_.AddArg(term);
}

// Converts a sum of norms into the SOCP form.
class SumOfNormsConverter : public SumConverter<SumOfNormsConverter> {
 private:
  Problem::MutObjective obj_;
  const int *num_terms_;
  int term_index_;

 public:
  SumOfNormsConverter(
      Problem &p, Problem::MutObjective obj, const int *num_terms)
    : SumConverter<SumOfNormsConverter>(p), obj_(obj),
      num_terms_(num_terms), term_index_(0) {}

  void VisitSqrt(UnaryExpr e);
};

void SumOfNormsConverter::VisitSqrt(UnaryExpr e) {
  Problem &p = problem();
  double inf = std::numeric_limits<double>::infinity();
  // Add a nonnegative variable ``x`` to represent ``sqrt(expr)``.
  Problem::Variable x = p.AddVar(0, inf);
  obj_.linear_expr().AddTerm(x.index(), coef());
  // Add a constraint ``x ^ 2 >= expr``.
  Problem::IteratedExprBuilder sum =
      p.BeginIterated(expr::SUM, num_terms_[term_index_++] + 1);
  SumOfSquaresConverter(p, sum).Visit(e.arg());
  UnaryExpr term = p.MakeUnary(
        expr::MINUS, p.MakeUnary(expr::POW2, p.MakeVariable(x.index())));
  sum.AddArg(term);
  p.AddCon(-inf, 0, p.EndIterated(sum));
}
}  // namespace mp

namespace {

// Converts a linear expression into the ASL form.
// expr: an expression to convert
// builder: a builder for converted expression
template <typename LinearExprBuilder>
void ConvertLinearExpr(mp::LinearExpr expr, LinearExprBuilder builder) {
  for (mp::LinearExpr::iterator i = expr.begin(), e = expr.end(); i != e; ++i)
    builder.AddTerm(i->var_index(), i->coef());
}

// Converts a nonlinear expression into the ASL form.
class MPToASLExprConverter :
    public mp::ExprVisitor<MPToASLExprConverter,
                           asl::NumericExpr, asl::LogicalExpr> {
 private:
  ASLBuilder &builder_;

 public:
  explicit MPToASLExprConverter(ASLBuilder &b) : builder_(b) {}

  asl::NumericExpr VisitNumericConstant(mp::NumericConstant c) {
    return builder_.MakeNumericConstant(c.value());
  }

  asl::NumericExpr VisitVariable(mp::Reference e) {
    return builder_.MakeVariable(e.index());
  }

  asl::NumericExpr VisitUnary(mp::UnaryExpr e) {
    return builder_.MakeUnary(e.kind(), Visit(e.arg()));
  }

  asl::NumericExpr VisitBinary(mp::BinaryExpr e) {
    return builder_.MakeBinary(e.kind(), Visit(e.lhs()), Visit(e.rhs()));
  }

  asl::NumericExpr VisitSum(mp::IteratedExpr e) {
    ASLBuilder::NumericExprBuilder sum = builder_.BeginSum(e.num_args());
    for (mp::IteratedExpr::iterator i = e.begin(), end = e.end(); i != end; ++i)
      sum.AddArg(Visit(*i));
    return builder_.EndSum(sum);
  }

  // TODO: convert all expresion kinds
};

// Adapts Problem interface for use with .nl reader.
class ProblemBuilder : public Problem {
 public:
  typedef mp::Function Function;
  typedef mp::Expr Expr;
  typedef mp::NumericExpr NumericExpr;
  typedef mp::LogicalExpr LogicalExpr;
  typedef mp::CountExpr CountExpr;
  typedef mp::Reference Reference;

  typedef IteratedExprBuilder NumericExprBuilder;
  typedef IteratedExprBuilder VarArgExprBuilder;
  typedef IteratedExprBuilder NumberOfExprBuilder;

  VarArgExprBuilder BeginVarArg(mp::expr::Kind kind, int num_args) {
    return BeginIterated(kind, num_args);
  }
  NumericExpr EndVarArg(VarArgExprBuilder builder) {
    return EndIterated(builder);
  }

  NumericExprBuilder BeginSum(int num_args) {
    return BeginIterated(mp::expr::SUM, num_args);
  }
  NumericExpr EndSum(NumericExprBuilder builder) {
    return EndIterated(builder);
  }

  struct ColumnSizeHandler {
    void Add(int) {
      // Ignore column sizes as the constraints are stored row-wise.
    }
  };

  // Returns a handler that receives column sizes in Jacobian.
  ColumnSizeHandler GetColumnSizeHandler() {
    return ColumnSizeHandler();
  }

  NumericExpr MakeBinary(
      mp::expr::Kind kind, NumericExpr lhs, NumericExpr rhs) {
    // Translate a POW expression with right-hand side of 2
    // to a POW2 expression.
    if (kind == mp::expr::POW) {
      mp::NumericConstant c = mp::Cast<mp::NumericConstant>(rhs);
      if (c && c.value() == 2)
        return MakeUnary(mp::expr::POW2, lhs);
    }
    return Problem::MakeBinary(kind, lhs, rhs);
  }
};

class SOCPConverter {
 private:
  ProblemBuilder problem_;
  ASLBuilder builder_;
  std::vector<int> col_sizes_;

  // The number of variables in the original problem.
  int num_vars_;

  // The number of algebraic constraints in the original problem.
  int num_algebraic_cons_;

  // Converts a nonlinear expression to ASL format.
  asl::NumericExpr Convert(mp::NumericExpr expr) {
    MPToASLExprConverter converter(builder_);
    return expr ? converter.Visit(expr) : asl::NumericExpr();
  }

  // Convertes an algebraic constraint to ASL format.
  void Convert(Problem::AlgebraicCon con) {
    mp::LinearExpr expr = con.linear_expr();
    ASLBuilder::LinearConBuilder con_builder = builder_.AddCon(
          con.lb(), con.ub(), Convert(con.nonlinear_expr()), expr.num_terms());
    ConvertLinearExpr(expr, con_builder);
  }

  typedef Problem::ObjList ObjList;

 public:
  explicit SOCPConverter(ASL *asl)
    : builder_(asl), num_vars_(0), num_algebraic_cons_(0) {}

  void Run(const char *stub);

  // Converts the problem into ASL format.
  void ConvertToASL();

  // Writes the solution.
  void WriteSol(ASL *asl, const char *msg, double *x, double *y,
                Option_Info *oi) {
    // Adjust the solution.
    asl->i.n_var0 = num_vars_;
    asl->i.n_con0 = num_algebraic_cons_;
    // TODO: reorder constraints if necessary
    write_sol_ASL(asl, msg, x, y, oi);
  }
};

void SOCPConverter::Run(const char *stub) {
  mp::ProblemBuilderToNLAdapter<ProblemBuilder> adapter(problem_);
  ReadNLFile(fmt::format("{}.nl", stub), adapter);
  num_vars_ = problem_.num_vars();
  num_algebraic_cons_ = problem_.num_algebraic_cons();
  if (!problem_.HasComplementarity()) {
    // TODO: check if all the constraints can be converted to SOCP form
    Problem::MutObjective obj = problem_.obj(0);
    mp::NumericExpr obj_expr = obj.nonlinear_expr();
    if (obj.type() == mp::obj::MIN && obj_expr) {
      mp::SumOfNormsDetector detector;
      if (detector.Visit(obj_expr)) {
        obj.set_nonlinear_expr(mp::NumericExpr());
        mp::SumOfNormsConverter(
              problem_, obj, detector.num_terms()).Visit(obj_expr);
        // TODO: convert to SOCP
      }
    }
  }

  mp::ProblemInfo info = mp::ProblemInfo();
  info.num_vars = problem_.num_vars();
  info.num_algebraic_cons = problem_.num_algebraic_cons();
  info.num_objs = problem_.num_objs();

  // Count nonzeros in objectives.
  ObjList objs = problem_.objs();
  for (ObjList::iterator i = objs.begin(), end = objs.end(); i != end; ++i)
    info.num_obj_nonzeros += i->linear_expr().num_terms();

  // Get algebraic constraint information.
  col_sizes_.resize(info.num_vars);
  double inf = std::numeric_limits<double>::infinity();
  Problem::AlgebraicConList cons = problem_.algebraic_cons();
  for (Problem::AlgebraicConList::iterator
       i = cons.begin(), end = cons.end(); i != end; ++i) {
    if (-inf < i->lb() && i->ub() < inf)
      ++info.num_ranges;
    if (i->nonlinear_expr())
      ++info.num_nl_cons;
    mp::LinearExpr expr = i->linear_expr();
    info.num_con_nonzeros += expr.num_terms();
    for (mp::LinearExpr::iterator j = expr.begin(), e = expr.end(); j != e; ++j)
      ++col_sizes_[j->var_index()];
  }

  // TODO: convert all problem info
  builder_.SetInfo(info);
  builder_.set_stub(stub);
}

void SOCPConverter::ConvertToASL() {
  // Convert variables.
  int num_vars = problem_.num_vars();
  for (int i = 0; i < num_vars; ++i) {
    Problem::Variable var = problem_.var(i);
    builder_.AddVar(var.lb(), var.ub(), var.type());
  }

  // Convert objectives.
  ObjList objs = problem_.objs();
  for (ObjList::iterator i = objs.begin(), end = objs.end(); i != end; ++i) {
    mp::LinearExpr expr = i->linear_expr();
    ASLBuilder::LinearObjBuilder obj_builder = builder_.AddObj(
          i->type(), Convert(i->nonlinear_expr()), expr.num_terms());
    ConvertLinearExpr(expr, obj_builder);
  }

  ASLBuilder::ColumnSizeHandler cols = builder_.GetColumnSizeHandler();
  for (int i = 0; i < num_vars - 1; ++i)
    cols.Add(col_sizes_[i]);

  // Convert nonlinear algebraic constraints first as required by ASL.
  Problem::AlgebraicConList algebraic_cons = problem_.algebraic_cons();
  for (Problem::AlgebraicConList::iterator
       i = algebraic_cons.begin(), end = algebraic_cons.end(); i != end; ++i) {
    if (i->nonlinear_expr())
      Convert(*i);
  }

  // Convert linear algebraic constraints.
  for (Problem::AlgebraicConList::iterator
       i = algebraic_cons.begin(), end = algebraic_cons.end(); i != end; ++i) {
    if (!i->nonlinear_expr())
      Convert(*i);
  }

  // Convert logical constraints.
  MPToASLExprConverter converter(builder_);
  Problem::LogicalConList logical_cons = problem_.logical_cons();
  for (Problem::LogicalConList::iterator
       i = logical_cons.begin(), end = logical_cons.end(); i != end; ++i) {
    builder_.AddCon(converter.Visit(i->expr()));
  }

  // TODO: convert common expressions, complementarity conditions and suffixes.

  builder_.EndBuild();
}
}  // namespace

// ASL structure with extra information for SOCP transformations.
struct SOCP_ASL : ASL_fg {
  SOCPConverter *converter;
};

extern "C" {

extern ASLhead ASLhead_ASL;

// Allocate ASL with extra information for SOCP transformations.
ASL *socp_ASL_alloc(int type) {
  // Set Stderr if necessary.
  if (!Stderr)
    Stderr_init_ASL();
  Mach_ASL();
  if (type < 1 || type > 5)
    return 0;
  static int msize[5] = {
    sizeof(SOCP_ASL),
    sizeof(SOCP_ASL),
    sizeof(ASL_fgh),
    sizeof(ASL_pfg),
    sizeof(ASL_pfgh)
  };
  int size = msize[type - 1];
  ASL *asl = static_cast<ASL*>(mymalloc(size));
  memcpy(asl, &edagpars_ASL, sizeof(Edagpars));
  memset(&asl->i, 0, size - sizeof(Edagpars));
  asl->i.ASLtype = type;
  asl->i.n_prob = 1;
  switch (type) {
  case ASL_read_pfg:
    reinterpret_cast<ASL_pfg*>(asl)->P.merge = 1;
    break;
  case ASL_read_pfgh:
    reinterpret_cast<ASL_pfgh*>(asl)->P.merge = 1;
    break;
  }
  ASLhead *head = asl->p.h.next = ASLhead_ASL.next;
  asl->p.h.prev = head->prev;
  head->prev = ASLhead_ASL.next = &asl->p.h;
  return cur_ASL = asl;
}

void socpASL_free(ASL **aslp) {
  delete reinterpret_cast<SOCP_ASL*>(*aslp)->converter;
  ::ASL_free(aslp);
}

void *socp_jac0dim(ASL *asl, const char *stub, ftnlen) {
  SOCP_ASL *socp_asl = reinterpret_cast<SOCP_ASL*>(asl);
  socp_asl->converter = new SOCPConverter(asl);
  socp_asl->converter->Run(stub);
  return 0;
}

int socp_qp_read(ASL *asl, void *, int) {
  reinterpret_cast<SOCP_ASL*>(asl)->converter->ConvertToASL();
  return 0;
}

void socp_write_sol(ASL *asl, const char *msg,
                    double *x, double *y, Option_Info *oi) {
  reinterpret_cast<SOCP_ASL*>(asl)->converter->WriteSol(asl, msg, x, y, oi);
}
}  // extern "C"
