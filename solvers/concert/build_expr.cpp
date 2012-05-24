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
#include <cmath>

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
         IloNumExpr lhs(build_expr (e->L.e)), rhs(build_expr (e->R.e));
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
         IloConstraint ifCond(build_constr (eif->e));
         IloNumVar ifVar(env, -IloInfinity, IloInfinity);
         mod.add (IloIfThen (env, ifCond,  ifVar == build_expr (eif->T)));
         mod.add (IloIfThen (env, !ifCond, ifVar == build_expr (eif->F)));
         return ifVar;
      }

      case OP_tanh: {
         PR ("tanh\n");
         IloNumExpr exp(IloExponent(2 * build_expr(e->L.e)));
         return (exp - 1) / (exp + 1);
      }

      case OP_tan:
         PR ("tan\n");
         return IloTan (build_expr (e->L.e));

      case OP_sqrt:
         PR ("sqrt\n");
         return IloPower (IloExprBase(build_expr (e->L.e)), 0.5);

      case OP_sinh: {
         PR ("sinh\n");
         IloNumExpr arg(build_expr(e->L.e));
         return (IloExponent(arg) - IloExponent(-arg)) / 2;
      }

      case OP_sin:
         PR ("sin\n");
         return IloSin (build_expr (e->L.e));

      case OP_log10:
         PR ("log10\n");
         return IloLog (build_expr (e->L.e)) / IloLog(10);

      case OP_log:
         PR ("log\n");
         return IloLog (build_expr (e->L.e));

      case OP_exp:
         PR ("exp\n");
         return IloExponent (build_expr (e->L.e));

      case OP_cosh: {
         PR ("cosh\n");
         IloNumExpr arg(build_expr(e->L.e));
         return (IloExponent(arg) + IloExponent(-arg)) / 2;
      }

      case OP_cos:
         PR ("cos\n");
         return IloCos (build_expr (e->L.e));

      case OP_atanh: {
         PR ("atanh\n");
         IloNumExpr arg(build_expr(e->L.e));
         return (IloLog(1 + arg) - IloLog(1 - arg)) / 2;
      }

      case OP_atan2: {
         PR ("atan2\n");
         IloNumExpr y(build_expr(e->L.e)), x(build_expr(e->R.e));
         IloNumExpr atan(IloArcTan(y / x));
         IloNumVar result(env, -IloInfinity, IloInfinity);
         mod.add(IloIfThen(env, x >= 0, result == atan));
         mod.add(IloIfThen(env, x <= 0 && y >= 0, result == atan + M_PI));
         mod.add(IloIfThen(env, x <= 0 && y <= 0, result == atan - M_PI));
         return result;
      }

      case OP_atan:
         PR ("atan\n");
         return IloArcTan (build_expr (e->L.e));

      case OP_asinh: {
         PR ("asinh\n");
         IloNumExpr arg(build_expr(e->L.e));
         return IloLog(arg + IloPower(IloPower(arg, 2) + 1, 0.5));
      }

      case OP_asin:
         PR ("asin\n");
         return IloArcSin (build_expr (e->L.e));

      case OP_acosh: {
         PR ("acosh\n");
         IloNumExpr arg(build_expr(e->L.e));
         return IloLog(arg + IloPower(arg + 1, 0.5) * IloPower(arg - 1, 0.5));
      }

      case OP_acos:
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
