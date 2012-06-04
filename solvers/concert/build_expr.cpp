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

#include <cstddef>
#include <cmath>
#include <algorithm>

#include "util.h"
#include "solvers/asl.h"
#include "solvers/nlp.h"
#include "solvers/opcode.hd"

#define PR if (get_option(DEBUGEXPR)) Printf

#ifndef M_PI
# define M_PI 3.14159265358979323846
#endif

using std::size_t;
using std::vector;

namespace {

bool has_zero_rhs(const expr *e)
{
   expr_n *rhs = e->R.en;
   return reinterpret_cast<size_t>(rhs->op) == OPNUM && rhs->v == 0;
}

class SameExpr {
 private:
  const expr *e;
  int elen;

 public:
  SameExpr(const expr *e);

  // Returns true if the stored expression is the same as the argument's
  // expression.
  bool operator()(const NumberOf& nof) const;
};

SameExpr::SameExpr(const expr *e) : e(e), elen(0)
{
   for (expr **ep = e->L.ep + 1; ep < e->R.ep; ep++, elen++);
}

bool SameExpr::operator()(const NumberOf& nof) const
{
   if (nof.num_vars() != elen)
      return false;

   for (expr **ep = e->L.ep + 1, **enp = nof.numberofexpr()->L.ep + 1;
        ep != e->R.ep; ep++, enp++) {
      if (!same_expr(*ep, *enp))
         return false;
   }
   return true;
}
}

IloIntVar NumberOf::add(real value, IloEnv env) {
  for (int i = 0, n = values_.getSize(); i < n; i++)
     if (values_[i] == value)
        return cards_[i];
  values_.add (value);
  IloIntVar cardVar(env, IloIntMin, IloIntMax);
  cards_.add (cardVar);
  return cardVar;
}

IloNumExprArray Driver::build_minmax_array(const expr *e)
{
   IloNumExprArray array(env_);
   for (de *d = reinterpret_cast<const expr_va*>(e)->L.d; d->e; ++d)
      array.add (build_expr (d->e));
   return array;
}

IloExpr Driver::build_expr (const expr *e)
{
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
         const expr_if *eif = reinterpret_cast<const expr_if*>(e);
         IloConstraint ifCond(build_constr (eif->e));
         IloNumVar ifVar(env_, -IloInfinity, IloInfinity);
         mod_.add (IloIfThen (env_, ifCond,  ifVar == build_expr (eif->T)));
         mod_.add (IloIfThen (env_, !ifCond, ifVar == build_expr (eif->F)));
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
         return (IloExponent(arg) - IloExponent(-arg)) * 0.5;
      }

      case OP_sin:
         PR ("sin\n");
         return IloSin (build_expr (e->L.e));

      case OP_log10:
         PR ("log10\n");
         return IloLog10 (build_expr (e->L.e));

      case OP_log:
         PR ("log\n");
         return IloLog (build_expr (e->L.e));

      case OP_exp:
         PR ("exp\n");
         return IloExponent (build_expr (e->L.e));

      case OP_cosh: {
         PR ("cosh\n");
         IloNumExpr arg(build_expr(e->L.e));
         return (IloExponent(arg) + IloExponent(-arg)) * 0.5;
      }

      case OP_cos:
         PR ("cos\n");
         return IloCos (build_expr (e->L.e));

      case OP_atanh: {
         PR ("atanh\n");
         IloNumExpr arg(build_expr(e->L.e));
         return (IloLog(1 + arg) - IloLog(1 - arg)) * 0.5;
      }

      case OP_atan2: {
         PR ("atan2\n");
         IloNumExpr y(build_expr(e->L.e)), x(build_expr(e->R.e));
         IloNumExpr atan(IloArcTan(y / x));
         IloNumVar result(env_, -IloInfinity, IloInfinity);
         mod_.add(IloIfThen(env_, x >= 0, result == atan));
         mod_.add(IloIfThen(env_, x <= 0 && y >= 0, result == atan + M_PI));
         mod_.add(IloIfThen(env_, x <= 0 && y <= 0, result == atan - M_PI));
         return result;
      }

      case OP_atan:
         PR ("atan\n");
         return IloArcTan (build_expr (e->L.e));

      case OP_asinh: {
         PR ("asinh\n");
         IloNumExpr arg(build_expr(e->L.e));
         return IloLog(arg + IloPower(IloSquare(arg) + 1, 0.5));
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

      case OPSUMLIST: {
         PR ("summation\n");
         IloExpr sumExpr(env_);
         for (expr **ep = e->L.ep, **end = e->R.ep; ep != end; ep++)
            sumExpr += build_expr (*ep);
         return sumExpr;
      }

      case OPintDIV:
         PR ("int division\n");
         return IloTrunc (build_expr (e->L.e) / build_expr (e->R.e));

      case OPround:
         PR ("round\n");
         if (!has_zero_rhs(e))
            throw UnsupportedExprError("round with nonzero second parameter");
         // Note that IloOplRound rounds half up.
         return IloOplRound(build_expr(e->L.e));

      case OPtrunc:
         PR ("trunc\n");
         if (!has_zero_rhs(e))
            throw UnsupportedExprError("trunc with nonzero second parameter");
         return IloTrunc(build_expr(e->L.e));

      case OP1POW:
         PR ("1pow %e\n", e->R.en->v);
         return IloPower (build_expr (e->L.e), e->R.en->v);

      case OP2POW:
         PR ("^2\n");
         return IloSquare (build_expr (e->L.e));

      case OPCPOW:
         PR ("cpow %e\n", e->L.en->v);
         return IloPower (e->L.en->v, build_expr (e->R.e));

      case OPNUM: {
         real n = reinterpret_cast<const expr_n*>(e)->v;
         PR ("%e\n", n);
         return IloExpr (env_, n);
      }

      case OPPLTERM: {
         plterm *p = e->L.p;
         int npce = p->n - 1;
         real *pce = p->bs;
         int j = reinterpret_cast<expr_v*>(e->R.e)->a;

         PR ("pl ");
         for (int i = 0; i < npce; i++)
            PR ("slp %f bkp %f ", pce[2*i], pce[2*i+1]);
         PR ("slp %f ", pce[2 * npce]);
         PR ("X[%d]\n", j+1);

         IloNumArray bkps(env_), slps(env_);
         for (int i = 0; i < npce; i++) {
            slps.add (pce[2*i]);
            bkps.add (pce[2*i+1]);
         }
         slps.add (pce[2 * npce]);
         return IloPiecewiseLinear (vars_[j], bkps, slps, 0, 0);
      }

      case OPVARVAL:
         PR ("X[%d]\n", e->a + 1);
         return vars_[e->a];

      /*----------------------------------------------------------------
        Logic extensions
      ----------------------------------------------------------------*/

      case OPCOUNT: {
         PR ("count\n");
         IloExpr sumExpr(env_);
         for (expr **ep = e->L.ep, **end = e->R.ep; ep != end; ep++)
            sumExpr += build_constr (*ep);
         return sumExpr;
      }

      case OPNUMBEROF: {
         PR ("number of\n");
         expr **ep = e->L.ep;
         if (reinterpret_cast<size_t>((*ep)->op) == OPNUM &&
             get_option(USENUMBEROF)) {
            return build_numberof(e);
         }
         IloExpr sumExpr(env_);
         IloExpr targetExpr(build_expr (*ep++));
         for (expr **end = e->R.ep; ep != end; ep++)
            sumExpr += (build_expr (*ep) == targetExpr);
         return sumExpr;
      }

      case OPVARSUBVAR: {
         PR ("vars in subscript of var\n");
         // AMPL should provide bounds for selectVar, use arbitrary for now.
         IloIntVar selectVar = IloIntVar (env_, 0, vars_.getSize() - 1);
         mod_.add (selectVar == build_expr (e->L.e));
         return vars_[selectVar];
      }

      default:
         throw UnsupportedExprError(get_opname(opnum));
   }
}

IloConstraint Driver::build_constr (const expr *e)
{
   size_t opnum = reinterpret_cast<size_t>(e->op);
   PR ("op %d  optype %d  ", opnum, optype[opnum]);

   switch(opnum) {
      case OPNUM: {
         real value = reinterpret_cast<const expr_n*>(e)->v;
         PR ("%e\n", value);
         IloNumVar dummy(env_, 1, 1);
         if (value == 1)
            return dummy == 1;
         if (value == 0)
            return dummy == 0;
         std::ostringstream oss;
         oss << "unexpected use of " << value << " as logical constant";
         throw Error(oss.str());
      }

      case LT:
         PR ("<\n");
         return !(build_expr (e->L.e) >= build_expr (e->R.e));

      case LE:
         PR ("<=\n");
         return build_expr (e->L.e) <= build_expr (e->R.e);

      case EQ:
         PR ("=\n");
         return build_expr (e->L.e) == build_expr (e->R.e);

      case GE:
         PR (">=\n");
         return build_expr (e->L.e) >= build_expr (e->R.e);

      case GT:
         PR (">\n");
         return !(build_expr (e->L.e) <= build_expr (e->R.e));

      case NE:
         PR ("!=\n");
         return build_expr (e->L.e) != build_expr (e->R.e);

      case OPATMOST:
         PR ("atmost\n");
         return build_expr (e->L.e) >= build_expr (e->R.e);

      case OPNOTATMOST:
         PR ("not atmost\n");
         return !(build_expr (e->L.e) >= build_expr (e->R.e));

      case OPATLEAST:
         PR ("atleast\n");
         return build_expr (e->L.e) <= build_expr (e->R.e);

      case OPNOTATLEAST:
         PR ("not atleast\n");
         return !(build_expr (e->L.e) <= build_expr (e->R.e));

      case OPEXACTLY:
         PR ("exactly\n");
         return build_expr (e->L.e) == build_expr (e->R.e);

      case OPNOTEXACTLY:
         PR ("not exactly\n");
         return build_expr (e->L.e) != build_expr (e->R.e);

      case OPOR:
         PR ("logical OR\n");
         return IloIfThen (env_, !build_constr (e->L.e), build_constr (e->R.e));

      case ORLIST: {
         PR ("logical EXISTS\n");
         IloOr disjunction(env_);
         for (expr **ep = e->L.ep; ep < e->R.ep; ep++)
            disjunction.add (build_constr (*ep));
         return disjunction;
      }

      case OPAND:
         PR ("logical AND\n");
         return build_constr (e->L.e) && build_constr (e->R.e);

      case ANDLIST: {
         PR ("logical FORALL\n");
         IloAnd conjunction(env_);
         for (expr **ep = e->L.ep, **end = e->R.ep; ep < end; ep++)
            conjunction.add (build_constr (*ep));
         return conjunction;
      }

      case OPNOT:
         PR ("logical NOT\n");
         return !build_constr (e->L.e);

      case OP_IFF:
         PR ("iff\n");
         return build_constr (e->L.e) == build_constr (e->R.e);

      case OPIMPELSE: {
         PR ("implies else\n");
         const expr_if *eif = reinterpret_cast<const expr_if*>(e);
         IloConstraint ifCond(build_constr (eif->e));
         return IloIfThen (env_,  ifCond, build_constr (eif->T))
             && IloIfThen (env_, !ifCond, build_constr (eif->F));
      }

      case OPALLDIFF: {
         PR ("all different\n");
         IloIntVarArray alldiffArray(env_);
         for (expr **ep = e->L.ep, **end = e->R.ep; ep != end; ep++) {
            if (reinterpret_cast<size_t>((*ep)->op) == OPVARVAL) {
               alldiffArray.add (vars_[(*ep)->a]);
            } else {
               IloIntVar alldiffVar(env_, IloIntMin, IloIntMax);
               mod_.add (alldiffVar == build_expr (*ep));
               alldiffArray.add (alldiffVar);
            }
         }
         return IloAllDiff (env_, alldiffArray);
      }

      default:
         throw IncompleteConstraintExprError(get_opname(opnum));
   }
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
         IloIntVar listVar(env_, IloIntMin, IloIntMax);
         vars.add (listVar);
         mod_.add (listVar == build_expr (*ep));
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
