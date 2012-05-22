#ifndef AMPL_SOLVERS_CONCERT_H
#define AMPL_SOLVERS_CONCERT_H

struct expr;

class IloExpr;
class IloEnv;
class IloModel;
class IloNumVarArray;

extern IloEnv env;
extern IloModel mod;
extern IloNumVarArray Var;

extern int usenumberof;
extern int debugexpr;

IloExpr build_expr(expr *e);
int concert_main(int argc, char **argv);

#endif // AMPL_SOLVERS_CONCERT_H
