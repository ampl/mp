#ifndef AMPL_SOLVERS_CONCERT_H
#define AMPL_SOLVERS_CONCERT_H

#include <stdexcept>

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
enum {OPVARSUBVAR = 99};

IloExpr build_expr(expr *e);
IloConstraint build_constr (expr*);
IloNumVar build_numberof (expr*);

const char *get_opname(int opcode);

int concert_main(int argc, char **argv);

// Thrown when an ASL expression cannot be converted to a Concert one.
class UnsupportedExprError : public std::runtime_error {
 public:
  UnsupportedExprError(const char *expr) :
    std::runtime_error(std::string("unsupported expression: ") + expr) {}
};


#endif // AMPL_SOLVERS_CONCERT_H
