/*-------------------------------------------------------------------------*/
/* AMPL/Concert driver                                       Robert Fourer */
/*                                                                         */
/* Name           : concert.cpp                                            */
/* Title          : AMPL/ILOG Concert driver                               */
/* By             : Robert Fourer                                          */
/* Date           : October 2000                                           */
/*                                                                         */
/* A driver to link AMPL linear integer programs with ILOG Concert 1.0     */
/* October 2000: Linear/Nonlinear version                                  */
/*-------------------------------------------------------------------------*/

#include "concert.h"

#include <algorithm>
#include <iostream>

#include <ilcplex/ilocplex.h>
#include <ilcp/cp.h>

#include "solvers/asl.h"
#include "solvers/nlp.h"
#include "solvers/getstub.h"
#include "solvers/opcode.hd"

using namespace std;

/*----------------------------------------------------------------------

  Process directives accessible from AMPL

----------------------------------------------------------------------*/

int debugexpr;
int ilogopttype;
int timing;
int usenumberof;

namespace {

// for suppressing "String literal to char*" warnings
#define CSTR(s) const_cast<char*>(s)

keyword keywds[] = { /* must be alphabetical */
   KW(CSTR("debugexpr"), I_val, &debugexpr,
      CSTR("print debugging information for expression trees")),
   KW(CSTR("ilogcplex"), IK1_val, &ilogopttype,
      CSTR("use ILOG CPLEX optimizer")),
   KW(CSTR("ilogsolver"), IK0_val, &ilogopttype,
      CSTR("use ILOG Constraint Programming optimizer")),
   KW(CSTR("timing"), I_val, &timing, CSTR("display timings for the run")),
   KW(CSTR("usenumberof"), I_val, &usenumberof,
      CSTR("consolidate 'numberof' expressions"))
};

Option_Info Oinfo = { CSTR("concert"), 0,
   CSTR("concert_options"), keywds, nkeywds, 0, 0,
   0, 0, 0, 0, 0, 20120521, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

// Returns the constant term in the first objective.
real objconst0(ASL_fg *a) {
  expr *e = a->I.obj_de_->e;
  return reinterpret_cast<size_t>(e->op) == OPNUM ?
      reinterpret_cast<expr_n*>(e)->v : 0;
}
}

Driver::Driver() : mod_(env_) {
   version_.resize(strlen(IloConcertVersion::_ILO_NAME) + 100);
   snprintf(&version_[0], version_.size() - 1,
       "%s %d.%d.%d", IloConcertVersion::_ILO_NAME,
       IloConcertVersion::_ILO_MAJOR_VERSION,
       IloConcertVersion::_ILO_MINOR_VERSION,
       IloConcertVersion::_ILO_TECH_VERSION);
   Oinfo.bsname = &version_[0];
   Oinfo.version = &version_[0];

   debugexpr = 0;
   ilogopttype = -1;
   timing = 0;
   usenumberof = 1;
}

Driver::~Driver() {
   env_.end();
}

/*----------------------------------------------------------------------

  Main Program

----------------------------------------------------------------------*/

int Driver::run(int argc, char **argv) {
   /*** Initialize timers ***/

   IloTimer timer(env_);
   timer.start();

   IloNum Times[4];
   Times[0] = timer.getTime();

   /*** Get name of .nl file; read problem sizes ***/

   ASL_fg *asl = reinterpret_cast<ASL_fg*>(ASL_alloc(ASL_read_fg));
   char *stub = getstub(&argv, &Oinfo);
   if (!stub)
     usage_ASL(&Oinfo, 1);
   FILE *nl = jac0dim(stub, strlen(stub));

   /*** Read coefficients & bounds & expression tree from .nl file ***/

   Uvx = static_cast<real*>(Malloc(n_var * sizeof(real)));
   Urhsx = static_cast<real*>(Malloc(n_con * sizeof(real)));

   efunc *r_ops_int[N_OPS];
   for (int i = 0; i < N_OPS; i++)
      r_ops_int[i] = reinterpret_cast<efunc*>(i);
   asl->I.r_ops_ = r_ops_int;
   want_derivs = 0;
   fg_read(nl, ASL_allow_CLP);
   asl->I.r_ops_ = 0;

   int n_var_int = nbv + niv + nlvbi + nlvci + nlvoi;

   /*** Get and process ILOG Concert options ***/

   if (getopts(argv, &Oinfo)) exit(1);

   int n_badvals = 0;

   switch(timing) {
      case 0:
         break;
      case 1:
         break;
      default:
         cerr << "Invalid value " << timing
            << " for directive timing" << endl;
         n_badvals++;
         }

   switch(ilogopttype) {
      case -1:
         nlo + nlc + n_lcon > 0 ? ilogopttype = 0 : ilogopttype = 1;
         break;
      case 0:   // ILOG Solver
         break;
      case 1:   // ILOG CPLEX
         break;
      default:
         cerr << "Invalid value " << ilogopttype
            << " for driver's ilogopttype setting" << endl;
         n_badvals++;
         }

   switch(usenumberof) {
      case 0:
         break;
      case 1:
         break;
      default:
         cerr << "Invalid value " << usenumberof
            << " for directive usenumberof" << endl;
         n_badvals++;
         }

   if (n_badvals)
      exit(1);

   /*-------------------------------------------------------------------

     Set up optimization problem in ILOG Concert

   -------------------------------------------------------------------*/

   vars_ = IloNumVarArray(env_,n_var);

   for (int j = 0; j < n_var - n_var_int; j++)
      vars_[j] = IloNumVar(env_, LUv[j], Uvx[j], ILOFLOAT);
   for (int j = n_var - n_var_int; j < n_var; j++)
      vars_[j] = IloNumVar(env_, LUv[j], Uvx[j], ILOINT);

   IloObjective MinOrMax(env_);

   if (n_obj > 0) {
      IloExpr objExpr(env_, objconst0(asl));
      if (0 < nlo)
         objExpr += build_expr (obj_de[0].e);
      for (ograd *og = Ograd[0]; og; og = og->next)
         objExpr += (og -> coef) * vars_[og -> varno];
      MinOrMax = IloObjective (env_, objExpr,
         objtype[0] == 0 ? IloObjective::Minimize : IloObjective::Maximize);
      IloAdd (mod_, MinOrMax);
   }

   IloRangeArray Con(env_,n_con);

   for (int i = 0; i < n_con; i++) {
      IloExpr conExpr(env_);
      for (cgrad *cg = Cgrad[i]; cg; cg = cg->next)
         conExpr += (cg -> coef) * vars_[cg -> varno];
      if (i < nlc) 
         conExpr += build_expr (con_de[i].e);
      Con[i] = (LUrhs[i] <= conExpr <= Urhsx[i]);
   }

   IloConstraintArray LCon(env_,n_lcon);

   for (int i = 0; i < n_lcon; i++)
      LCon[i] = build_constr (lcon_de[i].e);

   if (n_con > 0) mod_.add (Con);
   if (n_lcon > 0) mod_.add (LCon);

   finish_building_numberof ();

   if (timing) Times[1] = timer.getTime();

   /*-------------------------------------------------------------------

     Solve integer/linear program in CPLEX

   -------------------------------------------------------------------*/

   if (ilogopttype == 1) {
      IloCplex cplex (env_);
      cplex.extract (mod_);

      if (timing) Times[2] = timer.getTime();

      cplex.solve();

      if (timing) Times[3] = timer.getTime();

      IloNum objValue = cplex.getObjValue();

      int sSoFar = 0;
      char sMsg[256];

      sSoFar += Sprintf(sMsg,
         "\n%s: optimal solution found\n", Oinfo.bsname);

      if (nbv + niv > 0) {
         sSoFar += g_fmtop(sMsg+sSoFar,cplex.getNnodes());
         sSoFar += Sprintf(sMsg+sSoFar, " nodes, ");
         sSoFar += g_fmtop(sMsg+sSoFar,cplex.getNiterations());
         sSoFar += Sprintf(sMsg+sSoFar, " iterations, objective ");
         g_fmtop(sMsg+sSoFar, objValue);

         real *Xopt = new real [n_var];
         for (int j = 0; j < n_var; j++) Xopt[j] = cplex.getValue(vars_[j]);
         write_sol(sMsg, Xopt, 0, &Oinfo);
         delete [] Xopt;
      }

      else {
         sSoFar += g_fmtop(sMsg+sSoFar,cplex.getNiterations());
         sSoFar += Sprintf(sMsg+sSoFar, " iterations, objective ");
         g_fmtop(sMsg+sSoFar, objValue);

         real *Xopt = new real [n_var];
         real *Piopt = new real [n_con];
         for (int j = 0; j < n_var; j++) Xopt[j] = cplex.getValue(vars_[j]);
         for (int i = 0; i < n_con; i++) Piopt[i] = cplex.getDual(Con[i]);
         write_sol(sMsg, Xopt, Piopt, &Oinfo);
         delete [] Xopt;
         delete [] Piopt;
      }
   }

   /*-------------------------------------------------------------------

     Solve problem in ILOG Solver

   -------------------------------------------------------------------*/

   else {
      IloSolver solver (env_);
      solver.extract (mod_);

      if (timing) Times[2] = timer.getTime();

      IloBool successful = solver.solve();

      if (timing) Times[3] = timer.getTime();

      int sSoFar = 0;
      char sMsg[256];

      if (successful) {
         sSoFar += Sprintf(sMsg,
            "\n%s: solution found\n", Oinfo.bsname);

         sSoFar += g_fmtop(sMsg+sSoFar,solver.getNumberOfChoicePoints());
         sSoFar += Sprintf(sMsg+sSoFar, " choice points, ");
         sSoFar += g_fmtop(sMsg+sSoFar,solver.getNumberOfFails());
         sSoFar += Sprintf(sMsg+sSoFar, " fails");
         if (n_obj > 0) {
            sSoFar += Sprintf(sMsg+sSoFar, ", objective ");
            g_fmtop(sMsg+sSoFar, solver.getValue(MinOrMax));
         }

         real *Xopt = new real [n_var];

         for (int j = 0; j < n_var; j++) Xopt[j] = solver.getValue(vars_[j]);
         write_sol(sMsg, Xopt, 0, &Oinfo);
         delete [] Xopt;
      }
      else {
         sSoFar += Sprintf(sMsg,
            "\n%s: no solution found!\n", Oinfo.bsname);
         write_sol(sMsg, 0, 0, &Oinfo);
      }
   }

   if (timing) {
      Times[4] = timer.getTime();
      cerr << endl
           << "Define = " << Times[1] - Times[0] << endl
           << "Setup =  " << Times[2] - Times[1] << endl
           << "Solve =  " << Times[3] - Times[2] << endl
           << "Output = " << Times[4] - Times[3] << endl;
   }
   return 0;
}
