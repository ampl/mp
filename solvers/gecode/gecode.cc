#include "gecode.h"

#include <gecode/int.hh>
#include <gecode/search.hh>
#include <gecode/gist.hh>

#include <iostream>
#include <memory>

#include "getstub.h"
#include "nlp.h"
#include "opcode.hd"

using namespace Gecode;
using namespace std;

class GecodeModel: public Space {
private:
  IntVarArray vars_;

public:
  GecodeModel(int num_vars) :
      vars_(*this, num_vars) {
    // no leading zeros
    /*rel(*this, s, IRT_NQ, 0);
    rel(*this, m, IRT_NQ, 0);
    // all letters distinct
    distinct(*this, vars_);
    // post branching
    branch(*this, vars_, INT_VAR_SIZE_MIN, INT_VAL_MIN);*/
  }

  GecodeModel(bool share, GecodeModel &s) : Space(share, s) {
    vars_.update(*this, share, s.vars_);
  }

  Space *copy(bool share);

  IntVarArray& vars() {
    return vars_;
  }
};

Space *GecodeModel::copy(bool share) {
  return new GecodeModel(share,*this);
}

Driver::Driver() :
    asl(reinterpret_cast<ASL_fg*>(ASL_alloc(ASL_read_fg))) {
  oinfo_.reset(new Option_Info());
}

Driver::~Driver() {
}

int Driver::run(char **argv) {
  // Initialize timers.
  /*IloTimer timer(env_);
  timer.start();

  IloNum Times[5];
  Times[0] = timer.getTime();*/

  // Get name of .nl file and read problem sizes.
  char *stub = getstub(&argv, oinfo_.get());
  if (!stub) {
    usage_noexit_ASL(oinfo_.get(), 1);
    return 1;
  }
  FILE *nl = jac0dim(stub, strlen(stub));

  // Read coefficients, bounds & expression tree from .nl file.
  Uvx = static_cast<real*>(Malloc(n_var * sizeof(real)));
  Urhsx = static_cast<real*>(Malloc(n_con * sizeof(real)));

  efunc *r_ops_int[N_OPS];
  for (int i = 0; i < N_OPS; i++)
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
  std::auto_ptr<GecodeModel> model(new GecodeModel(n_var));
  IntVarArray &vars = model->vars();
  for (int j = 0; j < n_var; ++j)
    vars[j] = IntVar(*model, LUv[j], Uvx[j]);

  // TODO: objective
  /*if (n_obj > 0) {
   IloExpr objExpr(env_, objconst0(asl));
   if (0 < nlo)
   objExpr += build_expr (obj_de[0].e);
   for (ograd *og = Ograd[0]; og; og = og->next)
   objExpr += (og -> coef) * vars_[og -> varno];
   IloObjective MinOrMax(env_, objExpr,
   objtype[0] == 0 ? IloObjective::Minimize : IloObjective::Maximize);
   optimizer_->set_obj(MinOrMax);
   IloAdd (mod_, MinOrMax);
   }*/

  // TODO: constraints
  for (int i = 0; i < n_con; i++) {
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
    // TODO: ranges
    if (LUrhs[i] == Urhsx[i])
      linear(*model, c, x, IRT_EQ, LUrhs[i]);
    // TODO: nonlinear part
    //if (i < nlc)
    //   conExpr += model->ConvertExpr(con_de[i].e);
    //Con[i] = (LUrhs[i] <= conExpr <= Urhsx[i]);
  }

  /*IloConstraintArray LCon(env_,n_lcon);

   for (int i = 0; i < n_lcon; i++)
   LCon[i] = build_constr (lcon_de[i].e);

   if (n_con > 0) mod_.add (Con);
   if (n_lcon > 0) mod_.add (LCon);

   finish_building_numberof ();

   int timing = get_option(TIMING);
   Times[1] = timing ? timer.getTime() : 0;*/

  // Solve the problem.
  // TODO
  DFS<GecodeModel> e(model.get());
  model.reset();

  // TODO
  /*Times[2] = timing ? timer.getTime() : 0;*/
  std::auto_ptr<GecodeModel> solution(e.next());
  // TODO
  /* Times[3] = timing ? timer.getTime() : 0;*/

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
    solve_result_num = 501;
    status = "unknown solution status";
  }

  char message[256];
  int end = Sprintf(message, "%s: %s\n", oinfo_->bsname, status);
  // TODO
  //optimizer_->get_solution(asl, sMsg + sSoFar, primal, dual);
  write_sol(message, primal.empty() ? 0 : &primal[0], 0, oinfo_.get());

  // TODO
  /*if (timing) {
   Times[4] = timer.getTime();
   cerr << endl
   << "Define = " << Times[1] - Times[0] << endl
   << "Setup =  " << Times[2] - Times[1] << endl
   << "Solve =  " << Times[3] - Times[2] << endl
   << "Output = " << Times[4] - Times[3] << endl;
   }*/

  // TODO
  // search and print all solutions
  return 0;
}
