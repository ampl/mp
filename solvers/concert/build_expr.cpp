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

#include "concert.h"

#include <iostream>
#include <stdio.h>
#include <assert.h>
#include "string.h"

#include <ilconcert/ilomodel.h>
#include <ilcp/cp.h>

#include "solvers/asl.h"
#include "solvers/nlp.h"
#include "solvers/getstub.h"
#include "solvers/opcode.hd"
#include "opnames.hd"

#define PR if(debugexpr)Printf

IloConstraint build_constr (expr*);
IloNumVar build_numberof (expr*);

// Builds an array of expressions from the argument list of e.
IloNumExprArray build_minmax_array(expr *e)
{
   IloNumExprArray array(env);
   for (de *d = reinterpret_cast<expr_va*>(e)->L.d; d->e; ++d)
      array.add (build_expr (d->e));
   return array;
}

/*----------------------------------------------------------------------

  Walk expression tree and construct function,
  returning a Concert IloExpr

----------------------------------------------------------------------*/

IloExpr build_expr (expr *e)
{
   expr **ep;

   IloExpr sumExpr, targetExpr;

   IloNumArray loSubBnd, upSubBnd;
   IloNumVar resultVar;

   plterm *p;
   int npce, i, j;
   real *pce;
   IloNumArray bkps, slps;

   size_t opnum = reinterpret_cast<size_t>(e->op);
   PR ("op %d  optype %2d  ", opnum, optype[opnum]);

   switch(opnum) {
      case OPPLUS:
         PR ("+\n");
         return build_expr (e->L.e) + build_expr (e->R.e);

      case OPMINUS:
         PR ("-\n");
         return build_expr (e->L.e) - build_expr (e->R.e);

      case OPMULT:
         PR ("*\n");
         return build_expr (e->L.e) * build_expr (e->R.e);

      case OPDIV:
         PR ("/\n");
         return build_expr (e->L.e) / build_expr (e->R.e);

      case OPREM: {
         PR ("remainder\n");
         // a mod b = a - trunc(a / b) * b
         IloNumExpr lhs = build_expr (e->L.e), rhs = build_expr (e->R.e);
         return lhs - IloTrunc(lhs / rhs) * rhs;
      }

      case OPPOW:
         PR ("^\n");
         return IloPower (IloExprBase(build_expr (e->L.e)),
                          IloExprBase(build_expr (e->R.e)));

      case OPLESS:
         PR ("less\n");
         return IloMax (build_expr (e->L.e) - build_expr (e->R.e), 0.0);

      case MINLIST:
         PR ("min\n");
         return IloMin(build_minmax_array(e));

      case MAXLIST:
         PR ("max\n");
         return IloMax(build_minmax_array(e));

      case FLOOR:
         PR ("floor\n");
         return IloFloor(build_expr(e->L.e));

      case CEIL:
         PR ("ceil\n");
         return IloCeil(build_expr(e->L.e));

      case ABS:
         PR ("abs\n");
         return IloAbs (build_expr (e->L.e));

      case OPUMINUS:
         PR ("unary -\n");
         return - build_expr (e->L.e);

      case OPIFnl: {
         PR ("if\n");
         expr_if *eif = reinterpret_cast<expr_if*>(e);
         IloConstraint ifCond = build_constr (eif->e);
         IloNumVar ifVar = IloNumVar (env, -IloInfinity, IloInfinity);
         mod.add (IloIfThen (env, ifCond,  ifVar == build_expr (eif->T)));
         mod.add (IloIfThen (env, !ifCond, ifVar == build_expr (eif->F)));
         return ifVar;
      }

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
         PR ("int division\n");
         return IloTrunc (build_expr (e->L.e) / build_expr (e->R.e));

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

      case FUNCALL_opno:
         Printf ("function call -- not implemented\n");
         exit(1);

      case NUM_opno: {
         double n = reinterpret_cast<expr_n*>(e)->v;
         PR ("%e\n", n);
         return IloExpr (env, n);
      }

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

      case OPOR:
         throw Error("invalid logical OR in expression");

      case OPAND:
         throw Error("invalid logical AND in expression");

      case LT:
         throw Error("invalid < in expression");

      case LE:
         throw Error("invalid <= in expression");

      case EQ:
         throw Error("invalid = in expression");

      case GE:
         throw Error("invalid >= in expression");

      case GT:
         throw Error("invalid > in expression");

      case NE:
         throw Error("invalid != in expression");

      case OPNOT:
         throw Error("invalid logical NOT in expression");

      case ALLDIFF_opno:
         Printf ("invalid alldiff in expression\n");
         exit(1);

      default:
         Printf ("other -- not implemented\n");
         exit(1);
         return IloExpr();
   }
}
