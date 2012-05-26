/*-------------------------------------------------------------------------*/
/* AMPL/Concert "number of" expression builder               Robert Fourer */
/*                                                                         */
/* Name           : build_numberof.cpp                                     */
/* Title          : AMPL/ILOG Concert "number of" expression builder       */
/* By             : Robert Fourer                                          */
/* Date           : April 2001                                             */
/*                                                                         */
/* A component of the AMPL/ILOG Concert driver:                            */
/* Processes individual AMPL number-of constraints                         */
/* into ILOG Solver IloDistribute calls.                                   */
/*-------------------------------------------------------------------------*/

#include "concert.h"

#include <cassert>
#include <cstddef>
#include <sstream>

#include <ilconcert/ilomodel.h>

#include "solvers/asl.h"
#include "solvers/nlp.h"
#include "solvers/opcode.hd"
#include "opnames.hd"

using std::size_t;

namespace {

class numberof {
public:
   IloIntVarArray cards;
   IloIntArray values;
   IloIntVarArray vars;
   expr *numberofexpr;
   numberof *next;
};

numberof *numberofstart = 0;

}

/*----------------------------------------------------------------------

  Given a node for a number-of operator
  that has a constant as its first operand,
  add it to the driver's data structure that collects these operators.

----------------------------------------------------------------------*/

IloNumVar build_numberof (expr *e)
{
   assert(reinterpret_cast<size_t>(e->op) == OPNUMBEROF &&
          reinterpret_cast<size_t>((*e->L.ep)->op) == OPNUM);
   expr **ep, **enp;
   int elen, i;

   numberof *np;
   IloIntVar cardVar, listVar;
   // Did we previously see a number-of operator
   // having the same expression-list?

   ep = e->L.ep;
   for (ep++, elen = 0; ep < e->R.ep; ep++, elen++);

   for (np = numberofstart; np; np = np->next) {
      if (np->vars.getSize() != elen)
         continue;

      ep = e->L.ep;
      enp = np->numberofexpr->L.ep;
      for (ep++, enp++; ep < e->R.ep; ep++, enp++)
         if (!same_expr(*ep,*enp))
            break;

      if (ep < e->R.ep)
         continue;
      else
         break;
   }

   // New expression-list:
   // Build a new numberof structure.

   if (!np) {
      np = new numberof;
      np->next = numberofstart;
      numberofstart = np;
      np->numberofexpr = e;

      np->cards = IloIntVarArray(env);
      np->values = IloIntArray(env);
      np->vars = IloIntVarArray(env);

      ep = e->L.ep;
      (np->values).add (((expr_n*)(*ep))->v);

      for (ep++; ep < e->R.ep; ep++) {
         listVar = IloIntVar(env, IloIntMin, IloIntMax);
         (np->vars).add (listVar);
         mod.add (listVar == build_expr (*ep));
      }

      cardVar = IloIntVar(env, IloIntMin, IloIntMax);
      (np->cards).add (cardVar);
      return cardVar;
   }

   // Previously seen expression-list:
   // Add to its numberof structure.

   else {
      for (i=0; i < (np->values).getSize(); i++)
         if (np->values[i] == ((expr_n*)(*e->L.ep))->v)
            return np->cards[i];

      (np->values).add (((expr_n*)(*e->L.ep))->v);

      cardVar = IloIntVar(env, IloIntMin, IloIntMax);
      (np->cards).add (cardVar);
      return cardVar;
   }
}

bool same_expr (expr *e1, expr *e2)
{
   size_t opnum = reinterpret_cast<size_t>(e1->op);
   if (opnum != reinterpret_cast<size_t>(e2->op))
      return false;

   int type = optype[opnum];
   switch (type) {
      case OPTYPE_UNARY:
         return same_expr (e1->L.e, e2->L.e);

      case OPTYPE_BINARY:
         return same_expr (e1->L.e, e2->L.e) && 
                same_expr (e1->R.e, e2->R.e);

      case OPTYPE_VARARG: {
         de *d1 = reinterpret_cast<expr_va*>(e1)->L.d;
         de *d2 = reinterpret_cast<expr_va*>(e2)->L.d;
         for (; d1->e && d2->e; d1++, d2++)
            if (!same_expr (d1->e, d2->e))
               return false;
         return !d1->e && !d2->e;
      }

      case OPTYPE_PLTERM: {
         plterm *p1 = e1->L.p, *p2 = e2->L.p;
         if (p1->n != p2->n)
            return false;
         real *pce1 = p1->bs, *pce2 = p2->bs;
         for (int i = 0, n = p1->n * 2 - 1; i < n; i++) {
            if (pce1[i] != pce2[i])
               return false;
         }
         return same_expr (e1->R.e, e2->R.e);
      }

      case OPTYPE_IF: {
         expr_if *eif1 = reinterpret_cast<expr_if*>(e1);
         expr_if *eif2 = reinterpret_cast<expr_if*>(e2);
         return same_expr (eif1->e, eif2->e) &&
                same_expr (eif1->T, eif2->T) &&
                same_expr (eif1->F, eif2->F);
      }

      case OPTYPE_SUM:
      case OPTYPE_COUNT: {
         expr **ep1 = e1->L.ep;
         expr **ep2 = e2->L.ep;
         for (; ep1 < e1->R.ep && ep2 < e2->R.ep; ep1++, ep2++)
            if (!same_expr (*ep1, *ep2))
               return false;
         return !(ep1 < e1->R.ep) && !(ep2 < e2->R.ep);
      }

      case OPTYPE_FUNCALL:
      case OPTYPE_STRING:
         throw UnsupportedExprError(get_opname(opnum));

      case OPTYPE_NUMBER:
         return reinterpret_cast<expr_n*>(e1)->v ==
                reinterpret_cast<expr_n*>(e2)->v;

      case OPTYPE_VARIABLE:
         return e1->a == e2->a;

      default: {
         std::ostringstream oss;
         oss << "unknown operator type " << type << " in build_numberof";
         throw Error(oss.str());
      }
   }
}

void finish_building_numberof ()
{
   for (numberof *np = numberofstart; np; np = np->next)
      mod.add (IloDistribute (env, np->cards, np->values, np->vars));
}
