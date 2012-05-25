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

#include <ilconcert/ilomodel.h>

#include "solvers/asl.h"
#include "solvers/nlp.h"
#include "solvers/opcode.hd"
#include "opnames.hd"

using std::size_t;

bool same_expr (expr*,expr*);

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
   size_t opnum1 = reinterpret_cast<size_t>(e1->op);
   size_t opnum2 = reinterpret_cast<size_t>(e2->op);

   if (opnum1 != opnum2)
      return false;

   switch(optype[opnum1]) {

      case 1:
         return same_expr (e1->L.e, e2->L.e);

      case 2:
         return same_expr (e1->L.e, e2->L.e) && 
                same_expr (e1->R.e, e2->R.e);

      case 3: {
         de *d1 = ((expr_va*)e1)->L.d;
         de *d2 = ((expr_va*)e2)->L.d;
         for (; d1->e && d2->e; d1++, d2++)
            if (!same_expr (d1->e, d2->e))
               return false;
         return !d1->e && !d2->e;
      }

      case 4:
         Printf ("pl terms not implemented in build_numberof\n");
         exit(1);

      case 5: {
         expr_if *eif1 = (expr_if*)e1;
         expr_if *eif2 = (expr_if*)e2;
         return same_expr (eif1->e, eif2->e) &&
                same_expr (eif1->T, eif2->T) &&
                same_expr (eif1->F, eif2->F);
      }

      case 6:
      case 11: {
         expr **ep1 = e1->L.ep;
         expr **ep2 = e2->L.ep;
         for (; ep1 < e1->R.ep && ep2 < e2->R.ep; ep1++, ep2++)
            if (!same_expr (*ep1, *ep2))
               return false;
         return !(ep1 < e1->R.ep) && !(ep2 < e2->R.ep);
      }

      case 7:
         Printf ("function calls not implemented in build_numberof\n");
         exit(1);

      case 8:
         Printf ("string arguments not implemented in build_numberof\n");
         exit(1);

      case 9:
         return ((expr_n*)e1)->v == ((expr_n*)e2)->v;

      case 10:
         return e1->a == e2->a;

      default:
         Printf ("unknown operator type %d in build_numberof\n",
            optype[opnum1]);
         exit(1);
         return false;
   }
}

void finish_building_numberof () 
{
   numberof *np;

   for (np = numberofstart; np; np = np->next)
      mod.add (IloDistribute (env, np->cards, np->values, np->vars));
}
