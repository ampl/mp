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
#include <stdio.h>
#include <assert.h>
#include "string.h"

#include <ilconcert/ilomodel.h>
#include <ilcplex/ilocplex.h>
#include <ilcp/cp.h>

#include "solvers/asl.h"
#include "solvers/nlp.h"
#include "solvers/getstub.h"
#include "solvers/opcode.hd"

using namespace std;

#define CHR (char*)  // for suppressing "String literal to char*" warnings

#define R_OPS asl->I.r_ops_
#define asl ((ASL_fg*)a)

extern "C" real objconst0(ASL_fg*, int);
#undef objconst
#define objconst(n) objconst0(asl, n)

/*----------------------------------------------------------------------

  Initialize CONCERT environment
  Define array of decision variables

----------------------------------------------------------------------*/

IloEnv env;
IloModel mod(env);
IloNumVarArray Var;

/*----------------------------------------------------------------------

  Process directives accessible from AMPL

----------------------------------------------------------------------*/

int debugexpr = 0;
int ilogopttype = -1;
int timing = 0;
int usenumberof = 1;

namespace {

keyword keywds[] = { /* must be alphabetical */
   KW(CHR"debugexpr", I_val, &debugexpr,
      CHR"print debugging information for expression trees"),
   KW(CHR"ilogcplex", IK1_val, &ilogopttype,
      CHR"use ILOG CPLEX optimizer"),
   KW(CHR"ilogsolver", IK0_val, &ilogopttype,
      CHR"use ILOG Constraint Programming optimizer"),
   KW(CHR"timing", I_val, &timing, CHR"display timings for the run"),
   KW(CHR"usenumberof", I_val, &usenumberof,
      CHR"consolidate 'numberof' expressions"),
   };

Option_Info Oinfo = { CHR"concert", CHR"ILOG CONCERT 1.0",
   CHR"concert_options", keywds, nkeywds, 0, CHR"ILOG CONCERT 1.0",
   0, 0, 0, 0, 0, 20120521, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
}

/*----------------------------------------------------------------------

  Main Program

----------------------------------------------------------------------*/

int concert_main(int argc, char **argv) {

   FILE *nl;
   ASL *a;

   int i, j, n_var_int;
   char *stub;

   cgrad *cg;
   ograd *og;

   efunc *r_ops_int[N_OPS];

   /*** Initialize Concert ***/

   IloTimer timer(env);
   timer.start();

   IloNum Times[4];
   Times[0] = timer.getTime();

   /*** Get name of .nl file; read problem sizes ***/

   if (argc < 2) {
      fprintf(stderr, "Usage: %s stub\n", argv[0]);
      return 1;
      }

   a = ASL_alloc(ASL_read_fg);
   stub = getstub(&argv, &Oinfo);
   nl = jac0dim(stub, (fint)strlen(stub));

   /*** Read coefficients & bounds & expression tree from .nl file ***/

   Uvx = (real *)Malloc(n_var*sizeof(real));
   Urhsx = (real *)Malloc(n_con*sizeof(real));

   for (i = 0; i < N_OPS; i++)
      r_ops_int[i] = (efunc*)(unsigned long)i;
   R_OPS = r_ops_int;
   want_derivs = 0;
   fg_read(nl,ASL_allow_CLP);
   R_OPS = 0;

   n_var_int = nbv + niv + nlvbi + nlvci + nlvoi;

   /*** Set up some useful Concert arrays ***/

   IloNumArray loVarBnd(env,n_var), upVarBnd(env,n_var);
   for (j = 0; j < n_var; j++) {
      loVarBnd[j] = LUv[j];
      upVarBnd[j] = Uvx[j];
   }

   IloNumArray loConBnd(env,n_con), upConBnd(env,n_con);
   for (i = 0; i < n_con; i++) {
      loConBnd[i] = LUrhs[i];
      upConBnd[i] = Urhsx[i];
   }

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

   Var = IloNumVarArray(env,n_var);

   for (j = 0; j < n_var - n_var_int; j++)
      Var[j] = IloNumVar(env, loVarBnd[j], upVarBnd[j], ILOFLOAT);
   for (j = n_var - n_var_int; j < n_var; j++)
      Var[j] = IloNumVar(env, loVarBnd[j], upVarBnd[j], ILOINT);

   IloObjective MinOrMax(env);

   if (n_obj > 0) {
      IloExpr objExpr(env,objconst(0));
      if (0 < nlo)
         objExpr += build_expr (obj_de[0].e);
      for (og = Ograd[0]; og; og = og->next)
         objExpr += (og -> coef) * Var[og -> varno];
      MinOrMax = IloObjective (env, objExpr,
         objtype[0] == 0 ? IloObjective::Minimize : IloObjective::Maximize);
      IloAdd (mod, MinOrMax);
   }

   IloRangeArray Con(env,n_con);

   for (i = 0; i < n_con; i++) {
      IloExpr conExpr(env);
      for (cg = Cgrad[i]; cg; cg = cg->next)
         conExpr += (cg -> coef) * Var[cg -> varno];
      if (i < nlc) 
         conExpr += build_expr (con_de[i].e);
      Con[i] = (loConBnd[i] <= conExpr <= upConBnd[i]);
      }

   IloConstraintArray LCon(env,n_lcon);

   for (i = 0; i < n_lcon; i++) {
      LCon[i] = build_constr (lcon_de[i].e);
      }

   if (n_con > 0) mod.add (Con);
   if (n_lcon > 0) mod.add (LCon);

   finish_building_numberof ();

   if (timing) Times[1] = timer.getTime();

   /*-------------------------------------------------------------------

     Solve integer/linear program in CPLEX

   -------------------------------------------------------------------*/

   try {
      if (ilogopttype == 1) {
         IloCplex cplex (env);
         cplex.extract (mod);

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
            for(j = 0; j < n_var; j++) Xopt[j] = cplex.getValue(Var[j]);
            write_sol(sMsg, Xopt, 0, &Oinfo);
            delete [] Xopt;
         } 

         else {
            sSoFar += g_fmtop(sMsg+sSoFar,cplex.getNiterations());
            sSoFar += Sprintf(sMsg+sSoFar, " iterations, objective ");
            g_fmtop(sMsg+sSoFar, objValue);

            real *Xopt = new real [n_var];
            real *Piopt = new real [n_con];
            for(j = 0; j < n_var; j++) Xopt[j] = cplex.getValue(Var[j]);
            for(i = 0; i < n_con; i++) Piopt[i] = cplex.getDual(Con[i]);
            write_sol(sMsg, Xopt, Piopt, &Oinfo);
            delete [] Xopt;
            delete [] Piopt;
         }
      }

   /*-------------------------------------------------------------------

     Solve problem in ILOG Solver

   -------------------------------------------------------------------*/

      else {
         IloSolver solver (env);
         solver.extract (mod);

         if (timing) Times[2] = timer.getTime();

         IloBool successful = solver.solve();   // solve(IloSplit(env,Var))

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

            for(j = 0; j < n_var; j++) Xopt[j] = solver.getValue(Var[j]);
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
   }
   catch (const IloException& ex) {
      cerr << "Error: " << ex << endl;
      exit (-1);
   }
   env.end();
   return 0;
}
