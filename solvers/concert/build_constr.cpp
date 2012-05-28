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

#include "concert.h"

#include <ilconcert/ilomodel.h>

#include "util.h"
#include "solvers/asl.h"
#include "solvers/nlp.h"
#include "solvers/getstub.h"
#include "solvers/opcode.hd"

#define PR if (debugexpr) Printf

using std::size_t;

/*----------------------------------------------------------------------

  Walk expression tree and construct constraint,
  returning a Concert IloConstraint

----------------------------------------------------------------------*/

IloConstraint build_constr (const expr *e)
{
   expr **ep;
   expr_if *eif;

   IloOr disjunction;
   IloAnd conjunction;
   IloConstraint ifCond;
   IloNumVar dummy, alldiffVar;

   size_t opnum = reinterpret_cast<size_t>(e->op);
   PR ("op %d  optype %d  ", opnum, optype[opnum]);

   switch(opnum) {
      case OPNUM:
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

      case OPATMOST:
         PR ("atmost\n");
         return  build_expr (e->L.e) >=  build_expr (e->R.e);

      case OPNOTATMOST:
         PR ("not atmost\n");
         return  ! (build_expr (e->L.e) >=  build_expr (e->R.e));

      case OPATLEAST:
         PR ("atleast\n");
         return  build_expr (e->L.e) <=  build_expr (e->R.e);

      case OPNOTATLEAST:
         PR ("not atleast\n");
         return  ! (build_expr (e->L.e) <=  build_expr (e->R.e));

      case OPEXACTLY:
         PR ("exactly\n");
         return  build_expr (e->L.e) == build_expr (e->R.e);

      case OPNOTEXACTLY:
         PR ("not exactly\n");
         return  ! (build_expr (e->L.e) == build_expr (e->R.e));

      case OPOR:
         PR ("logical OR\n");
         return  build_constr (e->L.e) ||  build_constr (e->R.e);

      case ORLIST:
         PR ("logical EXISTS\n");
         disjunction = IloOr(env);
         for (ep = e->L.ep; ep < e->R.ep; ep++)
            disjunction.add (build_constr (*ep));
         return disjunction;

      case OPAND:
         PR ("logical AND\n");
         return build_constr (e->L.e) && build_constr (e->R.e);

      case ANDLIST:
         PR ("logical FORALL\n");
         conjunction = IloAnd(env);
         for (ep = e->L.ep; ep < e->R.ep; ep++)
            conjunction.add (build_constr (*ep));
         return conjunction;

      case OPNOT:
         PR ("logical NOT\n");
         return ! build_constr (e->L.e);

      case OP_IFF:
         PR ("iff\n");
         return  build_constr (e->L.e) == build_constr (e->R.e);

      case OPIMPELSE:
         PR ("implies else\n");
         eif = (expr_if*)e;
         ifCond = build_constr (eif->e);
         return IloIfThen (env,  ifCond, build_constr (eif->T))
             && IloIfThen (env, !ifCond, build_constr (eif->F));

      case LT:
         PR ("<\n");
         return  build_expr (e->L.e) <  build_expr (e->R.e);

      case LE:
         PR ("<=\n");
         return  build_expr (e->L.e) <=  build_expr (e->R.e);

      case EQ:
         PR ("=\n");
         return  build_expr (e->L.e) ==  build_expr (e->R.e);

      case GE:
         PR (">=\n");
         return  build_expr (e->L.e) >=  build_expr (e->R.e);

      case GT:
         PR (">\n");
         return  build_expr (e->L.e) >  build_expr (e->R.e);

      case NE:
         PR ("!=\n");
         return  build_expr (e->L.e) !=  build_expr (e->R.e);

      case OPALLDIFF: {
         PR ("all different\n");
         IloIntVarArray alldiffArray(env);
         for (ep = e->L.ep; ep < e->R.ep; ep++) {
            if (reinterpret_cast<size_t>((*ep)->op) == OPVARVAL)
               alldiffArray.add (Var[(*ep)->a]);
            else {
               alldiffVar = IloIntVar (env, IloIntMin, IloIntMax);
               mod.add (alldiffVar == build_expr (*ep));
               alldiffArray.add (alldiffVar);
            }
         }
         return IloAllDiff (env, alldiffArray);
      }

      default:
         throw IncompleteConstraintExprError(get_opname(opnum));
   }
}
