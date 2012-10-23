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

#include "ilogcp.h"

#include <cstddef>
#include <cmath>
#include <algorithm>

#include "solvers/util/util.h"
#include "nlp.h"
#include "opcode.hd"

#define PR if (get_option(DEBUGEXPR)) Printf

#ifndef M_PI
# define M_PI 3.14159265358979323846
#endif

using std::size_t;
using std::vector;

namespace {

class SameExpr {
 private:
  const expr *e;
  int elen;

 public:
  SameExpr(const expr *e);

  // Returns true if the stored expression is the same as the argument's
  // expression.
  bool operator()(const ampl::NumberOf& nof) const;
};

SameExpr::SameExpr(const expr *e) : e(e), elen(0)
{
  for (expr **ep = e->L.ep + 1; ep < e->R.ep; ep++, elen++) ;
}

bool SameExpr::operator()(const ampl::NumberOf& nof) const
{
  if (nof.num_vars() != elen)
    return false;

  for (expr **ep = e->L.ep + 1, **enp = nof.numberofexpr()->L.ep + 1;
      ep != e->R.ep; ep++, enp++) {
    if (!ampl::Equal(*ep, *enp))
      return false;
  }
  return true;
}
}

namespace ampl {

IloIntVar NumberOf::add(real value, IloEnv env) {
  for (int i = 0, n = values_.getSize(); i < n; i++)
     if (values_[i] == value)
        return cards_[i];
  values_.add (value);
  IloIntVar cardVar(env, IloIntMin, IloIntMax);
  cards_.add (cardVar);
  return cardVar;
}

IloNumVar Driver::build_numberof (const expr *e)
{
   assert(reinterpret_cast<size_t>(e->op) == OPNUMBEROF &&
          reinterpret_cast<size_t>((*e->L.ep)->op) == OPNUM);

   // Did we previously see a number-of operator
   // having the same expression-list?

   vector<NumberOf>::reverse_iterator np =
     find_if(numberofs_.rbegin(), numberofs_.rend(), SameExpr(e));

   // New expression-list:
   // Build a new numberof structure.

   if (np == numberofs_.rend()) {
      expr **ep = e->L.ep;
      IloIntArray values(env_);
      values.add (reinterpret_cast<expr_n*>(*ep)->v);

      IloIntVarArray vars(env_);
      for (ep++; ep < e->R.ep; ep++) {
        IloIntVar var(env_, IloIntMin, IloIntMax);
        vars.add(var);
        mod_.add(var == Visit(Expr(*ep)));
      }

      IloIntVar cardVar(env_, IloIntMin, IloIntMax);
      IloIntVarArray cards(env_);
      cards.add (cardVar);
      numberofs_.push_back(NumberOf(cards, values, vars, e));
      return cardVar;
   }

   // Previously seen expression-list:
   // Add to its numberof structure.
   return np->add(reinterpret_cast<expr_n*>(*e->L.ep)->v, env_);
}

void Driver::finish_building_numberof()
{
   for (vector<NumberOf>::const_iterator
        i = numberofs_.begin(), end = numberofs_.end(); i != end; ++i) {
      mod_.add (i->to_distribute(env_));
   }
   numberofs_.clear();
}

IloExpr Driver::VisitAtan2(BinaryExpr e) {
  IloNumExpr y(Visit(e.lhs())), x(Visit(e.rhs()));
  IloNumExpr atan(IloArcTan(y / x));
  IloNumVar result(env_, -IloInfinity, IloInfinity);
  mod_.add(IloIfThen(env_, x >= 0, result == atan));
  mod_.add(IloIfThen(env_, x <= 0 && y >= 0, result == atan + M_PI));
  mod_.add(IloIfThen(env_, x <= 0 && y <= 0, result == atan - M_PI));
  return result;
}
}
