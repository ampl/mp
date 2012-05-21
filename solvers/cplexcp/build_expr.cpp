/*-------------------------------------------------------------------------*/
/* AMPL/Concert expression builder                           Robert Fourer */
/*                                                                         */
/* Name           : build_expr.cpp                                         */
/* Title          : AMPL/ILOG Concert expression builder                   */
/* By             : Robert Fourer                                          */
/* Date           : October 2000                                           */
/*                                                                         */
/* Expression-building component of the AMPL/ILOG Concert driver           */
/* October 2000: Linear/Nonlinear version                                  */
/*-------------------------------------------------------------------------*/

#include <iostream>
#include <stdio.h>
#include <assert.h>
#include "string.h"

#include <ilconcert/ilomodel.h>
#include <ilcp/cp.h>

#include "asl.h"
#include "nlp.h"
#include "getstub.h"
#include "opnames.hd"

extern IloEnv env;
extern IloModel mod;
extern IloNumVarArray Var;

extern int usenumberof;
extern int debugexpr;
#define PR if(debugexpr)Printf

extern IloConstraint build_constr (expr*);
extern IloNumVar build_numberof (expr*);

/*----------------------------------------------------------------------

  Walk expression tree and construct function,
  returning a Concert IloExpr

----------------------------------------------------------------------*/

IloExpr build_expr (expr *e)
{
   expr **ep;
   expr_if *eif;
   de *d;

   IloExpr sumExpr, targetExpr;
   IloConstraint ifCond;
   IloNumVar minmaxVar, ifVar;
   IloNumVarArray minmaxArray;

   IloNumArray loSubBnd, upSubBnd;
   IloNumVar resultVar;

   plterm *p;
   int npce, i, j;
   real *pce;
   IloNumArray bkps, slps;

   size_t opnum = reinterpret_cast<size_t>(e->op);
   PR ("op %d  optype %2d  ", opnum, optype[opnum]);

   switch(opnum) {

      case PLUS_opno:
         PR ("+\n");
         return build_expr (e->L.e) + build_expr (e->R.e);

      case MINUS_opno:
         PR ("-\n");
         return build_expr (e->L.e) - build_expr (e->R.e);

      case MULT_opno:
         PR ("*\n");
         return build_expr (e->L.e) * build_expr (e->R.e);

      case DIV_opno:
         PR ("/\n");
         return build_expr (e->L.e) / build_expr (e->R.e);

      case REM_opno:
         Printf ("remainder -- not implemented\n");
         exit(1);

      case LESS_opno:
         PR ("less\n");
         return IloMax (build_expr (e->L.e) - build_expr (e->R.e), 0.0);

      case MINLIST_opno:
         PR ("min\n");

         minmaxArray = IloNumVarArray(env);
         for (d = ((expr_va*)e)->L.d; d->e; d++) {
            minmaxVar = IloNumVar (env, -IloInfinity, IloInfinity);
            mod.add (minmaxVar == build_expr (d->e));
            minmaxArray.add (minmaxVar);
         }
         return IloMin(minmaxArray);

      case MAXLIST_opno:
         PR ("max\n");

         minmaxArray = IloNumVarArray(env);
         for (d = ((expr_va*)e)->L.d; d->e; d++) {
            minmaxVar = IloNumVar (env, -IloInfinity, IloInfinity);
            mod.add (minmaxVar == build_expr (d->e));
            minmaxArray.add (minmaxVar);
         }
         return IloMax(minmaxArray);

      case FLOOR_opno:
         Printf ("floor -- not implemented\n");
         exit(1);

      case CEIL_opno:
         Printf ("ceil -- not implemented\n");
         exit(1);

      case ABS_opno:
         PR ("abs\n");
         return IloAbs (build_expr (e->L.e));

      case UMINUS_opno:
         PR ("unary -\n");
         return - build_expr (e->L.e);

      case IFnl_opno:
         PR ("if\n");

         eif = (expr_if*)e;
         ifCond = build_constr (eif->e);
         ifVar = IloNumVar (env, -IloInfinity, IloInfinity);
         mod.add (IloIfThen (env, ifCond,  ifVar == build_expr (eif->T)));
         mod.add (IloIfThen (env, !ifCond, ifVar == build_expr (eif->F)));
         return ifVar;

      case tanh_opno:
         Printf ("tanh -- not implemented\n");
         exit(1);

      case tan_opno:
         PR ("tan\n");
         return IloTan (build_expr (e->L.e));

      case sqrt_opno:
         PR ("sqrt\n");
         return IloPower (IloExprBase(build_expr (e->L.e)), 0.5);

      case sinh_opno:
         Printf ("sinh -- not implemented\n");
         exit(1);

      case sin_opno:
         PR ("sin\n");
         return IloSin (build_expr (e->L.e));

      case log10_opno:
         PR ("log10\n");
         return IloLog (build_expr (e->L.e)) / IloLog(10);

      case log_opno:
         PR ("log\n");
         return IloLog (build_expr (e->L.e));

      case exp_opno:
         PR ("exp\n");
         return IloExponent (build_expr (e->L.e));

      case cosh_opno:
         Printf ("cosh -- not implemented\n");
         exit(1);

      case cos_opno:
         PR ("cos\n");
         return IloCos (build_expr (e->L.e));

      case atanh_opno:
         Printf ("atanh -- not implemented\n");
         exit(1);

      case atan2_opno:
         Printf ("atan2 -- not implemented\n");
         exit(1);

      case atan_opno:
         PR ("atan\n");
         return IloArcTan (build_expr (e->L.e));

      case asinh_opno:
         Printf ("asin -- not implemented\n");
         exit(1);

      case asin_opno:
         PR ("asin\n");
         return IloArcSin (build_expr (e->L.e));

      case acosh_opno:
         Printf ("acos -- not implemented\n");
         exit(1);

      case acos_opno:
         PR ("acos\n");
         return IloArcCos (build_expr (e->L.e));

      case SUMLIST_opno:
         PR ("summation\n");

         sumExpr = IloExpr(env);
         for (ep = e->L.ep; ep < e->R.ep; ep++)
            sumExpr += build_expr (*ep);
         return sumExpr;

      case intDIV_opno:
         Printf ("int division -- not implemented\n");
         exit(1);

      case precision_opno:
         Printf ("precision -- not implemented\n");
         exit(1);

      case round_opno:
         Printf ("round -- not implemented\n");
         exit(1);

      case trunc_opno:
         Printf ("trunc -- not implemented\n");
         exit(1);

      case POWBAS_opno:
         PR ("1pow %e\n", e->R.en->v);
	 return IloPower (IloExprBase(build_expr (e->L.e)), e->R.en->v);

      case POW2_opno:
         PR ("^2\n");
         return IloSquare (build_expr (e->L.e));

      case POWEXP_opno:
         PR ("cpow %e\n", e->L.en->v);
	 return IloPower (e->L.en->v, IloExprBase(build_expr (e->R.e)));

      case POW_opno:
         PR ("^\n");
         return IloPower (IloExprBase(build_expr (e->L.e)), 
                          IloExprBase(build_expr (e->R.e)));

      case FUNCALL_opno:
         Printf ("function call -- not implemented\n");
         exit(1);

      case NUM_opno:
         PR ("%e\n", ((expr_n*)e)->v);
         return IloExpr (env, ((expr_n*)e)->v);

      case PLTERM_opno:
         p = e->L.p;
         npce = p->n;
         pce = p->bs;
         j = ((expr_v *)e->R.e)->a;

         PR ("pl ");
         for (i = 0; i < npce-1; i++)
            PR ("slp %f bkp %f ", pce[2*i], pce[2*i+1]);
         PR ("slp %f ", pce[2*(npce-1)]);
         PR ("X[%d]\n", j+1);

         bkps = IloNumArray(env);
         slps = IloNumArray(env);
         for (i = 0; i < npce-1; i++) {
            slps.add (pce[2*i]);
            bkps.add (pce[2*i+1]);
         }
         slps.add (pce[2*(npce-1)]);
         return IloPiecewiseLinear (Var[j],bkps,slps,0,0);

      case IFSYM_opno:
         Printf ("if sym -- not implemented\n");
         exit(1);

      case HOL_opno:
         Printf ("string argument -- not implemented\n");
         exit(1);

      case VARVAL_opno:
         PR ("X[%d]\n", e->a + 1);
         return Var[e->a];

      /*----------------------------------------------------------------
        Logic extensions
      ----------------------------------------------------------------*/

      case COUNT_opno:
         PR ("count\n");

         sumExpr = IloExpr(env);
         for (ep = e->L.ep; ep < e->R.ep; ep++)
            sumExpr += build_constr (*ep);
         return sumExpr;

      case NUMBEROF_opno:
         PR ("number of\n");

         ep = e->L.ep;
         if (reinterpret_cast<size_t>((*ep)->op) != NUM_opno || !usenumberof) {
            sumExpr = IloExpr(env);
            targetExpr = build_expr (*ep);
            for (ep++; ep < e->R.ep; ep++)
               sumExpr += (build_expr (*ep) == targetExpr);
            return sumExpr;
         }
         else
            return build_numberof (e);

      case /* VARSUBVAR_opno */ 99: {
         PR ("vars in subscript of var\n");

         IloIntVar selectVar = IloIntVar (env,loSubBnd[e->a],upSubBnd[e->a]); 
         mod.add (selectVar == build_expr (e->L.e)); 

         return Var[selectVar];
      }

      case ATMOST_opno:
         Printf ("invalid atmost in expression\n");
         exit(1);

      case ATLEAST_opno:
         Printf ("invalid atleast in expression\n");
         exit(1);

      case EXACTLY_opno:
         Printf ("invalid exactly in expression\n");
         exit(1);

      case OR_opno:
         Printf ("invalid logical OR in expression\n");
         exit(1);

      case AND_opno:
         Printf ("invalid logical AND in expression\n");
         exit(1);

      case LT_opno:
         Printf ("invalid < in expression\n");
         exit(1);

      case LE_opno:
         Printf ("invalid <= in expression\n");
         exit(1);

      case EQ_opno:
         Printf ("invalid = in expression\n");
         exit(1);

      case GE_opno:
         Printf ("invalid >= in expression\n");
         exit(1);

      case GT_opno:
         Printf ("invalid > in expression\n");
         exit(1);

      case NE_opno:
         Printf ("invalid != in expression\n");
         exit(1);

      case ALLDIFF_opno:
         Printf ("invalid alldiff in expression\n");
         exit(1);

      default:
         Printf ("other -- not implemented\n");
         exit(1);
         return IloExpr();
   }
}
