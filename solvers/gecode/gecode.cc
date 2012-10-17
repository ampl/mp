#include "gecode.h"

#include <gecode/int.hh>
#include <gecode/minimodel.hh>
#include <gecode/search.hh>
#include <gecode/gist.hh>

#include <iostream>
#include <memory>
#include <ctime>

#include "solvers/util/util.h"
#include "solvers/util/expr.h"
#include "solvers/getstub.h"
#include "solvers/nlp.h"
#include "solvers/opcode.hd"

using std::cerr;
using std::endl;
using std::vector;

using ampl::Expr;
using ampl::GetOpName;

using Gecode::BoolExpr;
using Gecode::DFS;
using Gecode::IntArgs;
using Gecode::IntVarArgs;
using Gecode::IntVar;
using Gecode::IntVarArray;
using Gecode::IntRelType;
using Gecode::IRT_NQ;
using Gecode::linear;
using Gecode::Space;

namespace {

class GecodeProblem: public Space, public ampl::ExprVisitor<GecodeProblem, IntVar> {
private:
  IntVarArray vars_;
  IntVar obj_;
  IntRelType obj_irt_; // IRT_NQ - no objective,
                       // IRT_LE - minimization, IRT_GQ - maximization
  static const BoolExpr DUMMY_EXPR;

public:
  GecodeProblem(int num_vars) :
      vars_(*this, num_vars), obj_irt_(IRT_NQ) {}

  GecodeProblem(bool share, GecodeProblem &s) :
    Space(share, s), obj_irt_(s.obj_irt_) {
    vars_.update(*this, share, s.vars_);
    if (obj_irt_ != IRT_NQ)
      obj_.update(*this, share, s.obj_);
  }

  Space *copy(bool share);

  IntVarArray &vars() { return vars_; }
  IntVar &obj() { return obj_; }

  void SetObjective(bool max) {
    obj_irt_ = max ? Gecode::IRT_GQ : Gecode::IRT_LE;
  }

  virtual void constrain(const Space& best);

  // Converts the specified logical ASL expression such as 'or', '<=', or
  // 'alldiff' into an equivalent Gecode expression.
  BoolExpr ConvertLogicalExpr(const expr *e);

  IntVar VisitNumber(ampl::Number n) {
    double value = n.value();
    IntVar var(*this, value, value);
    rel(*this, var == value);
    return var;
  }

  IntVar VisitVariable(ampl::Variable v) {
    return vars_[v.index()];
  }

  // Converts the specified arithmetic ASL expression into an equivalent
  // Concert expression.
  IntVar ConvertArithmeticExpr(const expr *e);
};

const BoolExpr GecodeProblem::DUMMY_EXPR((Gecode::BoolVar()));

Space *GecodeProblem::copy(bool share) {
  return new GecodeProblem(share, *this);
}

void GecodeProblem::constrain(const Space& best) {
  if (obj_irt_ != IRT_NQ)
    rel(*this, obj_, obj_irt_, static_cast<const GecodeProblem&>(best).obj_);
}

#define PR if (0) Printf

BoolExpr GecodeProblem::ConvertLogicalExpr(const expr *e) {
  size_t opnum = reinterpret_cast<size_t>(e->op);
  PR("op %d  optype %d  ", opnum, optype[opnum]);

  switch(opnum) {
  case NE:
    PR("!=\n");
    return Visit(Expr(e->L.e)) != Visit(Expr(e->R.e));

  case OPALLDIFF: {
    PR("all different\n");
    expr **ep = e->L.ep, **end = e->R.ep;
    IntVarArgs x(end - ep);
    for (unsigned i = 0; ep != end; ++ep, ++i) {
      x[i] = reinterpret_cast<size_t>((*ep)->op) == OPVARVAL ?
          vars_[(*ep)->a] : Visit(Expr(*ep));
    }
    distinct(*this, vars_);
    return DUMMY_EXPR;
  }

  default:
    throw ampl::IncompleteConstraintExprError(GetOpName(opnum));
  }
}
}

namespace ampl {

Driver::Driver() :
    asl(reinterpret_cast<ASL_fg*>(ASL_alloc(ASL_read_fg))) {
  oinfo_.reset(new Option_Info());
}

Driver::~Driver() {
}

int Driver::run(char **argv) {
  // Get the name of the .nl file and read problem sizes.
  char *stub = getstub(&argv, oinfo_.get());
  if (!stub) {
    usage_noexit_ASL(oinfo_.get(), 1);
    return 1;
  }
  FILE *nl = jac0dim(stub, strlen(stub));

  // Read coefficients, bounds & expression tree from the .nl file.
  Uvx = static_cast<real*>(Malloc(n_var * sizeof(real)));
  Urhsx = static_cast<real*>(Malloc(n_con * sizeof(real)));

  efunc *r_ops_int[N_OPS];
  for (int i = 0; i < N_OPS; ++i)
    r_ops_int[i] = reinterpret_cast<efunc*>(i);
  asl->I.r_ops_ = r_ops_int;
  want_derivs = 0;
  fg_read(nl, ASL_allow_CLP);
  asl->I.r_ops_ = 0;

  // TODO: parse options
  /*if (!parse_options(argv))
   return 1;*/

  // Set up an optimization problem in Gecode.
  int n_var_int = nbv + niv + nlvbi + nlvci + nlvoi;
  int n_var_cont = n_var - n_var_int;
  if (n_var_cont != 0) {
    cerr << "Gecode doesn't support continuous variables" << endl;
    return 1;
  }
  std::auto_ptr<GecodeProblem> problem(new GecodeProblem(n_var));
  IntVarArray &vars = problem->vars();
  for (int j = 0; j < n_var; ++j)
    vars[j] = IntVar(*problem, LUv[j], Uvx[j]);

  // Post branching.
  branch(*problem, vars, Gecode::INT_VAR_SIZE_MIN, Gecode::INT_VAL_MIN);

  if (n_obj > 0) {
    // TODO: convert the objective expr
    problem->SetObjective(objtype[0] == 0);
    /*IloExpr objExpr(env_, objconst0(asl));
    if (0 < nlo)
    objExpr += build_expr (obj_de[0].e);
    for (ograd *og = Ograd[0]; og; og = og->next)
    objExpr += (og -> coef) * vars_[og -> varno];
    IloObjective MinOrMax(env_, objExpr,
    objtype[0] == 0 ? IloObjective::Minimize : IloObjective::Maximize);
    optimizer_->set_obj(MinOrMax);
    IloAdd (mod_, MinOrMax);*/
  }

  // Convert constraints.
  for (int i = 0; i < n_con; ++i) {
    int num_terms = 0;
    for (cgrad *cg = Cgrad[i]; cg; cg = cg->next)
      ++num_terms;
    IntArgs c(num_terms);
    IntVarArgs x(num_terms);
    int index = 0;
    for (cgrad *cg = Cgrad[i]; cg; cg = cg->next) {
      c[index] = cg->coef;
      x[index] = vars[cg->varno];
      ++index;
    }
    double lb = LUrhs[i];
    double ub = Urhsx[i];
    if (lb <= negInfinity) {
      linear(*problem, c, x, Gecode::IRT_LE, ub);
    } else if (ub >= Infinity) {
      linear(*problem, c, x, Gecode::IRT_GQ, lb);
    } else if (lb == ub) {
      linear(*problem, c, x, Gecode::IRT_EQ, lb);
    } else {
      linear(*problem, c, x, Gecode::IRT_GQ, lb);
      linear(*problem, c, x, Gecode::IRT_LE, ub);
    }
    // TODO: nonlinear part
    //if (i < nlc)
    //   conExpr += model->ConvertExpr(con_de[i].e);
    //Con[i] = (LUrhs[i] <= conExpr <= Urhsx[i]);
  }

  // Convert logical constraints.
  for (int i = 0; i < n_lcon; ++i) {
    const expr *asl_expr = lcon_de[i].e;
    BoolExpr gecode_expr(problem->ConvertLogicalExpr(asl_expr));
    if (reinterpret_cast<size_t>(asl_expr->op) != OPALLDIFF)
      rel(*problem, gecode_expr);
  }

  // TODO
  // finish_building_numberof();

  DFS<GecodeProblem> e(problem.get());
  problem.reset();

  // Solve the problem.
  std::auto_ptr<GecodeProblem> solution(e.next());

  // Convert solution status.
  const char *status = 0;
  vector<real> primal;
  if (solution.get()) {
    solve_result_num = 100;
    status = "feasible solution";
    IntVarArray &vars = solution->vars();
    primal.resize(n_var);
    for (int j = 0; j < n_var; ++j)
      primal[j] = vars[j].val();
  } else {
    solve_result_num = 200;
    status = "infeasible problem";
  }

  char message[256];
  Sprintf(message, "%s: %s\n", oinfo_->bsname, status);
  write_sol(message, primal.empty() ? 0 : &primal[0], 0, oinfo_.get());
  return 0;
}
}
