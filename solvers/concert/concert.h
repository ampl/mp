#ifndef AMPL_SOLVERS_CONCERT_H
#define AMPL_SOLVERS_CONCERT_H

struct expr;

class IloConstraint;
class IloExpr;
class IloEnv;
class IloModel;
class IloNumVar;
class IloNumVarArray;

extern IloEnv env;
extern IloModel mod;
extern IloNumVarArray Var;

extern int usenumberof;
extern int debugexpr;

// Variable subscripted by a variable - not implemented in AMPL yet.
enum { OPVARSUBVAR = 99 };

IloExpr build_expr(expr *e);
IloConstraint build_constr (expr*);
void finish_building_numberof ();

int concert_main(int argc, char **argv);

#endif // AMPL_SOLVERS_CONCERT_H
