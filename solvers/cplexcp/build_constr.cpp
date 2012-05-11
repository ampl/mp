/*-------------------------------------------------------------------------*/
/* AMPL/Concert constraint builder                           Robert Fourer */
/*                                                                         */
/* Name           : build_expr.cpp                                         */
/* Title          : AMPL/ILOG Concert constraint builder                   */
/* By             : Robert Fourer                                          */
/* Date           : October 2000                                           */
/*                                                                         */
/* Constraint-building component of the AMPL/ILOG Concert driver           */
/* October 2000: Linear/Nonlinear version                                  */
/*-------------------------------------------------------------------------*/

#include <iostream.h>
#include <stdio.h>
#include <assert.h>
#include "string.h"

#include <ilconcert/ilomodel.h>
#include <ilcplex/ilocplex.h>
#include <ilsolver/ilosolver.h>

#include "asl.h"
#include "nlp.h"
#include "getstub.h"
#include "opnames.hd"

extern IloEnv env;
extern IloModel mod;
extern IloNumVarArray Var;

extern int debugexpr;
#define PR if(debugexpr)Printf

extern IloExpr build_expr (expr*);

/*----------------------------------------------------------------------

  Walk expression tree and construct constraint,
  returning a Concert IloConstraint

----------------------------------------------------------------------*/

IloConstraint build_constr (expr *e)
{
   efunc *op;
   expr **ep;
   expr_if *eif;

   IloInt opnum;
   IloOr disjunction;
   IloAnd conjunction;
   IloConstraint ifCond;
   IloNumVar dummy, alldiffVar;
   IloNumVarArray alldiffArray;

   opnum = (int) e->op;
   PR ("op %d  optype %d  ", opnum, optype[opnum]);

   switch(opnum) {

      case PLUS_opno:
         Printf ("incomplete constraint expression using +\n");
         exit(1);

      case MINUS_opno:
         Printf ("incomplete constraint expression using -\n");
         exit(1);

      case MULT_opno:
         Printf ("incomplete constraint expression using *\n");
         exit(1);

      case DIV_opno:
         Printf ("incomplete constraint expression using /\n");
         exit(1);

      case REM_opno:
         Printf ("incomplete constraint expression using remainder\n");
         exit(1);

      case POW_opno:
         Printf ("incomplete constraint expression using ^\n");
         exit(1);

      case LESS_opno:
         Printf ("incomplete constraint expression using less\n");
         exit(1);

      case MINLIST_opno:
         Printf ("incomplete constraint expression using min\n");
         exit(1);

      case MAXLIST_opno:
         Printf ("incomplete constraint expression using max\n");
         exit(1);

      case FLOOR_opno:
         Printf ("incomplete constraint expression using floor\n");
         exit(1);

      case CEIL_opno:
         Printf ("incomplete constraint expression using ceil\n");
         exit(1);

      case ABS_opno:
         Printf ("incomplete constraint expression using abs\n");
         exit(1);

      case UMINUS_opno:
         Printf ("incomplete constraint expression using unary -\n");
         exit(1);

      case IFnl_opno:
         Printf ("incomplete constraint expression using if nl\n");
         exit(1);

      case tanh_opno:
         Printf ("incomplete constraint expression using tanh\n");
         exit(1);

      case tan_opno:
         Printf ("incomplete constraint expression using tan\n");
         exit(1);

      case sqrt_opno:
         Printf ("incomplete constraint expression using sqrt\n");
         exit(1);

      case sinh_opno:
         Printf ("incomplete constraint expression using sinh\n");
         exit(1);

      case sin_opno:
         Printf ("incomplete constraint expression using sin\n");
         exit(1);

      case log10_opno:
         Printf ("incomplete constraint expression using log10\n");
         exit(1);

      case log_opno:
         Printf ("incomplete constraint expression using log\n");
         exit(1);

      case exp_opno:
         Printf ("incomplete constraint expression using exp\n");
         exit(1);

      case cosh_opno:
         Printf ("incomplete constraint expression using cosh\n");
         exit(1);

      case cos_opno:
         Printf ("incomplete constraint expression using cos\n");
         exit(1);

      case atanh_opno:
         Printf ("incomplete constraint expression using atanh\n");
         exit(1);

      case atan2_opno:
         Printf ("incomplete constraint expression using atan2\n");
         exit(1);

      case atan_opno:
         Printf ("incomplete constraint expression using atan\n");
         exit(1);

      case asinh_opno:
         Printf ("incomplete constraint expression using asin\n");
         exit(1);

      case asin_opno:
         Printf ("incomplete constraint expression using asin\n");
         exit(1);

      case acosh_opno:
         Printf ("incomplete constraint expression using acos\n");
         exit(1);

      case acos_opno:
         Printf ("incomplete constraint expression using acos\n");
         exit(1);

      case SUMLIST_opno:
         Printf ("incomplete constraint expression using summation\n");
         exit(1);

      case intDIV_opno:
         Printf ("incomplete constraint expression using int division\n");
         exit(1);

      case precision_opno:
         Printf ("incomplete constraint expression using precision\n");
         exit(1);

      case round_opno:
         Printf ("incomplete constraint expression using round\n");
         exit(1);

      case trunc_opno:
         Printf ("incomplete constraint expression using trunc\n");
         exit(1);

      case POWBAS_opno:
         Printf ("incomplete constraint expression using 1pow %e\n",
            e->R.en->v);
         exit(1);

      case POW2_opno:
         Printf ("incomplete constraint expression using ^2\n");
         exit(1);

      case POWEXP_opno:
         Printf ("incomplete constraint expression using cpow %e\n", 
            e->L.en->v);
         exit(1);

      case FUNCALL_opno:
         Printf ("incomplete constraint expression using function call\n");
         exit(1);

      case NUM_opno:
         PR ("%e\n", ((expr_n*)e)->v);
         if (((expr_n*)e)->v == 1) {
            dummy = IloNumVar (env, 1, 1);
            return  dummy == 1;
         }
         else if (((expr_n*)e)->v == 0) {
            dummy = IloNumVar (env, 1, 1);
            return  dummy == 0;
         }
         else {
            Printf ("unexpected use of %e as logical constant\n", 
               ((expr_n*)e)->v);
            exit(1);
         }

      case PLTERM_opno:
         Printf ("incomplete constraint expression using pl term\n");
         exit(1);

      case IFSYM_opno:
         Printf ("incomplete constraint expression using if sym\n");
         exit(1);

      case HOL_opno:
         Printf ("incomplete constraint expression using string argument\n");
         exit(1);

      case VARVAL_opno:
         Printf ("incomplete constraint expression using X[%d]\n", e->a + 1);
         exit(1);

      /*----------------------------------------------------------------
        Logic extensions
      ----------------------------------------------------------------*/

      case COUNT_opno:
         Printf ("incomplete constraint expression using count\n");
         exit(1);

      case NUMBEROF_opno:
         Printf ("incomplete constraint expression using numberof\n");
         exit(1);

      case ATMOST_opno:
         PR ("atmost\n");
         return  build_expr (e->L.e) >=  build_expr (e->R.e);

      case NOTATMOST_opno:
         PR ("not atmost\n");
         return  ! (build_expr (e->L.e) >=  build_expr (e->R.e));

      case ATLEAST_opno:
         PR ("atleast\n");
         return  build_expr (e->L.e) <=  build_expr (e->R.e);

      case NOTATLEAST_opno:
         PR ("not atleast\n");
         return  ! (build_expr (e->L.e) <=  build_expr (e->R.e));

      case EXACTLY_opno:
         PR ("exactly\n");
         return  build_expr (e->L.e) == build_expr (e->R.e);

      case NOTEXACTLY_opno:
         PR ("not exactly\n");
         return  ! (build_expr (e->L.e) == build_expr (e->R.e));

      case OR_opno:
         PR ("logical OR\n");
         return  build_constr (e->L.e) ||  build_constr (e->R.e);

      case ORLIST_opno:
         PR ("logical EXISTS\n");
         disjunction = IloOr(env);
         for (ep = e->L.ep; ep < e->R.ep; *ep++)
            disjunction.add (build_constr (*ep));
         return disjunction;

      case AND_opno:
         PR ("logical AND\n");
         return build_constr (e->L.e) && build_constr (e->R.e);

      case ANDLIST_opno:
         PR ("logical FORALL\n");
         conjunction = IloAnd(env);
         for (ep = e->L.ep; ep < e->R.ep; *ep++)
            conjunction.add (build_constr (*ep));
         return conjunction;

      case NOT_opno:
         PR ("logical NOT\n");
         return ! build_constr (e->L.e);

      case IFF_opno:
         PR ("iff\n");
         return  build_constr (e->L.e) == build_constr (e->R.e);

      case IMPELSE_opno:
         PR ("implies else\n");
         eif = (expr_if*)e;
         ifCond = build_constr (eif->e);
         return IloIfThen (env,  ifCond, build_constr (eif->T))
             && IloIfThen (env, !ifCond, build_constr (eif->F));

      case LT_opno:
         PR ("<\n");
         return  build_expr (e->L.e) <  build_expr (e->R.e);

      case LE_opno:
         PR ("<=\n");
         return  build_expr (e->L.e) <=  build_expr (e->R.e);

      case EQ_opno:
         PR ("=\n");
         return  build_expr (e->L.e) ==  build_expr (e->R.e);

      case GE_opno:
         PR (">=\n");
         return  build_expr (e->L.e) >=  build_expr (e->R.e);

      case GT_opno:
         PR (">\n");
         return  build_expr (e->L.e) >  build_expr (e->R.e);

      case NE_opno:
         PR ("!=\n");
         return  build_expr (e->L.e) !=  build_expr (e->R.e);

      case ALLDIFF_opno:
         PR ("all different\n");
         alldiffArray = IloNumVarArray(env);
         for (ep = e->L.ep; ep < e->R.ep; *ep++) {
            if ((int) (*ep)->op == VARVAL_opno)
               alldiffArray.add (Var[(*ep)->a]); 
            else {
               alldiffVar = IloIntVar (env, -IloInfinity, IloInfinity);
               mod.add (alldiffVar == build_expr (*ep));
               alldiffArray.add (alldiffVar); }
            }
         return IloAllDiff (env, alldiffArray);

      default:
         Printf ("other\n");
         exit(1);
   }
}
